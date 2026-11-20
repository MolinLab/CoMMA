# SPDX-License-Identifier: LGPL-2.1-only
# Exact target-block execution of the existing TASSEL-rule native algorithm.
gbs_parallel_sources <- function() c("R/utilities.R", "R/gbs_imputation_validation.R", "R/gbs_native_imputation.R",
  "R/gbs_imputation_io.R", "R/gbs_tassel_native.R", "R/gbs_tassel_search.R", "R/gbs_parallel.R")

gbs_parallel_budget <- function(resources, shared_gib, worker_gib, tasks) {
  values <- unlist(resources[c("global_cpu_budget", "reserved_cpus", "workers", "worker_threads", "memory_budget_gib", "reported_ram_gb")])
  if (length(values) != 6L || any(!is.finite(values)) || resources$worker_threads != 1L ||
      resources$global_cpu_budget < 1L || resources$reserved_cpus < 0L || resources$workers < 1L ||
      any(values[1:4] != floor(values[1:4])) || resources$memory_budget_gib <= 0 ||
      resources$memory_budget_gib > resources$reported_ram_gb * 1e9/1024^3 ||
      !is.finite(shared_gib) || shared_gib < 0 || !is.finite(worker_gib) || worker_gib <= 0 || tasks < 1L) stop("Invalid CPU/RAM budget.")
  scheduler <- suppressWarnings(as.numeric(Sys.getenv(c("SLURM_CPUS_PER_TASK","NSLOTS","PBS_NP"))))
  scheduler <- scheduler[is.finite(scheduler) & scheduler >= 1]
  effective_cpu <- min(c(resources$global_cpu_budget,scheduler))
  cpu <- effective_cpu - resources$reserved_cpus
  memory <- floor((resources$memory_budget_gib - shared_gib)/worker_gib)
  workers <- min(resources$workers, cpu, memory, tasks)
  if (workers < 1L) stop("Allocation cannot support one worker with the estimated shared and worker memory.")
  list(workers = as.integer(workers), threads = 1L, global_cpu_budget = resources$global_cpu_budget,effective_cpu_limit=effective_cpu,
    estimated_shared_gib = shared_gib, estimated_per_worker_gib = worker_gib,
    estimated_total_gib = shared_gib + workers*worker_gib, memory_budget_gib = resources$memory_budget_gib)
}

gbs_parallel_atomic <- function(object, path) {
  temporary <- paste0(path, ".", Sys.getpid(), ".tmp")
  saveRDS(object, temporary)
  if (file.exists(path) || !file.rename(temporary, path)) stop("Checkpoint cannot be committed: ", path)
  invisible(path)
}

gbs_parallel_task <- function(task, context, worker, directory, fingerprint) {
  path <- file.path(directory, sprintf("block_%06d.rds", task$id))
  if (file.exists(path)) {
    old <- readRDS(path)
    if (!identical(old$fingerprint, fingerprint) || !identical(old$task, task)) stop("Stale task checkpoint.")
    return(old$value)
  }
  value <- worker(task, context)
  gbs_parallel_atomic(list(fingerprint = fingerprint, task = task, value = value), path)
  value
}

gbs_parallel_map <- function(tasks, worker, context, allocation, checkpoint, fingerprint) {
  directory <- gbs_project_path(checkpoint, must_exist = FALSE)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  # Workers never receive genotype truth ledgers or phenotype tables.
  execute <- function(task) gbs_parallel_task(task, context, worker, directory, fingerprint)
  if (allocation$workers == 1L) result <- lapply(tasks, execute) else if (.Platform$OS.type != "windows") {
    result <- parallel::mclapply(tasks, execute, mc.cores = allocation$workers,
      mc.preschedule = TRUE, mc.set.seed = FALSE, mc.allow.recursive = FALSE)
  } else {
    executable <- Sys.getenv("GWAS_R_EXECUTABLE")
    if (!nzchar(executable)) stop("Use the authorized project R launcher for Windows worker tests.")
    # The project authorizes R.exe rather than Rscript.exe. R.exe needs --args
    # before the worker connection parameters; its first -e starts the worker.
    cluster <- parallel::makePSOCKcluster(allocation$workers, rscript = executable,
      rscript_args = c("--vanilla", "--slave", "-e", shQuote("parallel:::.workRSOCK()"), "--args"),
      setup_timeout = 30, outfile = file.path(directory, "workers.log"))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterCall(cluster, function(root, paths, libraries) {
      setwd(root); .libPaths(libraries)
      Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
      data.table::setDTthreads(1L)
      for (path in paths) sys.source(path, envir = .GlobalEnv)
      NULL
    }, getwd(), gbs_parallel_sources(), .libPaths())
    # Export the immutable context once, rather than once per marker block.
    parallel::clusterExport(cluster, c("context", "worker", "directory", "fingerprint"), envir = environment())
    result <- parallel::parLapply(cluster, tasks, function(task)
      gbs_parallel_task(task, context, worker, directory, fingerprint))
  }
  if (length(result) != length(tasks) || any(vapply(result, function(x) is.null(x) || inherits(x, "try-error"), logical(1L)))) stop("A parallel task failed; incomplete results cannot be scored.")
  result
}

gbs_parallel_ld_block <- function(task, context) {
  map <- context$map; p <- nrow(map); targets <- task$indices
  allowed <- lapply(targets, function(j) {
    ok <- seq_len(p) != j
    if (context$distance >= 0) ok <- ok & abs(map$position_bp-map$position_bp[[j]]) <= context$distance
    if (context$policy == "same_chromosome") ok <- ok & map$chromosome == map$chromosome[[j]]
    which(ok)
  })
  pool <- sort(unique(unlist(allowed, use.names = FALSE)))
  neighbours <- r2_values <- vector("list", length(targets)); ties <- logical(length(targets))
  if (length(pool)) {
    a <- crossprod(context$major[, targets, drop=FALSE], context$major[, pool, drop=FALSE])
    b <- crossprod(context$minor[, targets, drop=FALSE], context$major[, pool, drop=FALSE])
    c <- crossprod(context$major[, targets, drop=FALSE], context$minor[, pool, drop=FALSE])
    d <- crossprod(context$minor[, targets, drop=FALSE], context$minor[, pool, drop=FALSE])
    r2 <- gbs_tassel_r2_counts(a,b,c,d)
    for (i in seq_along(targets)) {
      values <- r2[i,match(allowed[[i]],pool)]; valid <- is.finite(values)
      candidates <- allowed[[i]][valid]; values <- values[valid]
      rank <- order(-values,map$marker_id[candidates])
      if (length(rank)>context$max_ld) ties[[i]] <- values[rank[[context$max_ld]]] == values[rank[[context$max_ld+1L]]]
      rank <- head(rank,context$max_ld); neighbours[[i]] <- candidates[rank]; r2_values[[i]] <- values[rank]
    }
  }
  list(indices=targets,neighbours=neighbours,r2=r2_values,ties=ties)
}

fit_gbs_tassel_parallel <- function(reference, markers, settings, resources, target_ids, checkpoint) {
  gbs_validate_calls(reference)
  map <- markers[match(colnames(reference),markers$marker_id),,drop=FALSE]; rownames(map)<-NULL
  if (anyNA(map$marker_id) || anyDuplicated(map$marker_id) || anyNA(map$position_bp) || any(map$position_bp<=0) ||
      !settings$chromosome_policy %in% c("tassel_position_only","same_chromosome") ||
      !is.finite(settings$max_distance_bp) || (settings$max_distance_bp<0 && settings$max_distance_bp != -1) ||
      settings$block_size<1 || settings$block_size != floor(settings$block_size) ||
      settings$max_ld_sites<2 || settings$max_ld_sites>200 || settings$max_ld_sites != floor(settings$max_ld_sites)) stop("Invalid parallel LD contract.")
  targets <- match(sort(unique(target_ids)),map$marker_id)
  if (!length(targets)||anyNA(targets)) stop("Invalid target marker IDs.")
  byte <- gbs_tassel_byte_codes(map$alleles); rownames(byte)<-map$marker_id
  statistics <- gbs_reference_marker_statistics(reference)
  alt_minor <- statistics$alt_frequency<.5 | (statistics$alt_frequency==.5 & byte[,1]<byte[,3])
  alt_minor[!is.finite(statistics$alt_frequency)]<-TRUE
  major <- reference<2; minor <- reference>0
  flip<-which(!alt_minor)
  if(length(flip)) {major[,flip]<-reference[,flip,drop=FALSE]>0;minor[,flip]<-reference[,flip,drop=FALSE]<2}
  major[is.na(major)]<-FALSE;minor[is.na(minor)]<-FALSE
  groups <- split(targets,ceiling(seq_along(targets)/settings$block_size))
  tasks <- lapply(seq_along(groups),function(i)list(id=as.integer(i),indices=unname(groups[[i]])))
  shared <- as.double(length(reference))*8*12/1024^3+.5
  per_worker <- as.double(nrow(reference))*ncol(reference)*8*4/1024^3+.25
  allocation <- gbs_parallel_budget(resources,shared,per_worker,length(tasks))
  reference_hash<-digest::digest(reference,algo="sha256"); map_hash<-digest::digest(map,algo="sha256")
  fingerprint<-digest::digest(list(reference_hash,map_hash,settings,targets,
    digest::digest(file="R/gbs_parallel.R",algo="sha256")),algo="sha256")
  context<-list(map=map,major=major,minor=minor,distance=settings$max_distance_bp,
    policy=settings$chromosome_policy,max_ld=settings$max_ld_sites)
  started<-proc.time()
  cat("LD:",nrow(reference),"lines,",ncol(reference),"predictors,",length(targets),"target SNPs;",allocation$workers,"workers.\n");flush.console()
  results<-gbs_parallel_map(tasks,gbs_parallel_ld_block,context,allocation,checkpoint,fingerprint)
  neighbours<-r2<-vector("list",ncol(reference));ties<-logical(ncol(reference))
  names(neighbours)<-names(r2)<-map$marker_id
  for(result in results) {neighbours[result$indices]<-result$neighbours;r2[result$indices]<-result$r2;ties[result$indices]<-result$ties}
  structure(list(reference=reference,markers=map,marker_ids=map$marker_id,marker_statistics=statistics,
    neighbours=neighbours,r2=r2,byte_codes=byte,LD_cutoff_tie=ties,max_ld_sites=as.integer(settings$max_ld_sites),
    settings=list(max_distance_bp=settings$max_distance_bp,chromosome_policy=settings$chromosome_policy,
      maximum_memory_gib=resources$memory_budget_gib,block_size=settings$block_size),
    provenance=list(algorithm="native_R_TASSEL_5.2.96_parallel_target_blocks_v1",fitted_reference_sha256=reference_hash,
      marker_metadata_sha256=map_hash,fitted_target_ids=map$marker_id[targets],phenotype_used=FALSE,query_used_for_fit=FALSE,
      estimated_fit_memory_gib=allocation$estimated_total_gib,allocation=allocation,
      fingerprint=fingerprint,elapsed_fit_seconds=unname((proc.time()-started)[["elapsed"]]))),class="gbs_tassel_native")
}

gbs_parallel_predict_block <- function(task, context) {
  target<-context$targets[task$indices,,drop=FALSE]; model<-context$model; query<-context$query
  prediction<-matrix(NA_real_,nrow(target),nrow(context$grid),dimnames=list(NULL,context$grid$id))
  for(id in unique(target$marker_id)) {
    j<-match(id,model$marker_ids);positions<-which(target$marker_id==id)
    for(ld in unique(context$grid$ld)) {
      predictors<-head(model$neighbours[[j]],ld)
      predictors<-predictors[order(model$markers$chromosome[predictors],model$markers$position_bp[predictors],model$marker_ids[predictors])]
      donor<-model$reference[,predictors,drop=FALSE]
      coverage<-if(length(predictors))rowMeans(!is.na(donor))else rep(0,nrow(donor))
      for(i in positions) {
        frequency<-model$marker_statistics$alt_frequency[[j]]
        columns<-which(context$grid$ld==ld)
        if(is.finite(frequency)&&frequency %in% c(0,1)) {prediction[i,columns]<-2*frequency;next}
        if(!length(predictors))next
        q<-query[target$sample_id[[i]],predictors]
        eligible<-which(coverage*mean(!is.na(q))*length(predictors)>10 & !is.na(model$reference[,j]) &
          rownames(model$reference)!=target$sample_id[[i]])
        if(!length(eligible))next
        distance<-gbs_tassel_distance(q,donor[eligible,,drop=FALSE])$distance
        good<-is.finite(distance);eligible<-eligible[good];distance<-distance[good]
        if(!length(eligible))next
        calls<-model$reference[eligible,j]
        rank<-head(order(distance,model$byte_codes[j,calls+1L],rownames(model$reference)[eligible]),max(context$grid$k[columns]))
        for(column in columns) {
          selected<-head(rank,context$grid$k[[column]])
          prediction[i,column]<-gbs_tassel_vote(calls[selected],distance[selected],ld,model$byte_codes[j,])
        }
      }
    }
  }
  list(targets=target,prediction=prediction)
}

predict_gbs_targets_parallel <- function(model,query,target_cells,grid,resources,checkpoint) {
  gbs_validate_calls(query)
  if(!all(model$marker_ids %in% colnames(query))||!nrow(target_cells)||
    anyDuplicated(target_cells[c("sample_id","marker_id")])||anyNA(target_cells[c("sample_id","marker_id")])||
    any(!target_cells$marker_id %in% model$provenance$fitted_target_ids)||any(grid$ld>model$max_ld_sites))stop("Invalid parallel prediction targets/grid.")
  query<-query[,model$marker_ids,drop=FALSE]
  shared<-intersect(rownames(model$reference),rownames(query))
  if(length(shared)&&!identical(unname(model$reference[shared,,drop=FALSE]),unname(query[shared,,drop=FALSE])))stop("Shared IDs have different raw calls.")
  target<-target_cells[c("sample_id","marker_id")];target<-target[order(target$marker_id,target$sample_id),,drop=FALSE];rownames(target)<-NULL
  index<-cbind(match(target$sample_id,rownames(query)),match(target$marker_id,colnames(query)))
  if(anyNA(index)||any(!is.na(query[index])))stop("Targets must be missing observed-query positions.")
  groups<-split(seq_len(nrow(target)),ceiling(match(target$marker_id,unique(target$marker_id))/16L))
  tasks<-lapply(seq_along(groups),function(i)list(id=as.integer(i),indices=unname(groups[[i]])))
  shared_memory<-as.double(length(model$reference)+length(query))*8*6/1024^3+.5
  worker_memory<-as.double(nrow(model$reference))*model$max_ld_sites*8*12/1024^3+.25
  allocation<-gbs_parallel_budget(resources,shared_memory,worker_memory,length(tasks))
  fingerprint<-digest::digest(list(model$provenance$fingerprint,digest::digest(query,algo="sha256"),target,grid,
    digest::digest(file="R/gbs_parallel.R",algo="sha256")),algo="sha256")
  cat("Masked predictions:",nrow(target),"cells,",nrow(grid),"settings;",allocation$workers,"workers.\n");flush.console()
  results<-gbs_parallel_map(tasks,gbs_parallel_predict_block,list(model=model,query=query,targets=target,grid=grid),allocation,checkpoint,fingerprint)
  list(targets=do.call(rbind,lapply(results,`[[`,"targets")),prediction=do.call(rbind,lapply(results,`[[`,"prediction")),
    allocation=allocation,truth_values_received=FALSE,fingerprint=fingerprint)
}
