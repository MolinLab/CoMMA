# Common coordinate-defined Shapley games for saved breeding predictors.
# Source-map coordinates define the universe; model fitting and outcomes are
# deliberately absent. Existing interpretation implementations stay unchanged.
breeding_coordinate_map<-function(markers,width_bp=10000000){
 required<-c("marker_id","chromosome","position_bp")
 if(!is.data.frame(markers)||!all(required%in%names(markers))||!nrow(markers)||
    anyNA(markers[,required])||anyDuplicated(markers$marker_id)||
    any(!nzchar(markers$marker_id))||any(!nzchar(as.character(markers$chromosome)))||
    any(!is.finite(markers$position_bp))||any(markers$position_bp<1)||
    any(markers$position_bp!=floor(markers$position_bp))||length(width_bp)!=1L||
    !is.finite(width_bp)||width_bp<1||width_bp!=floor(width_bp))stop("Invalid coordinate-only marker map.")
 start<-floor((markers$position_bp-1)/width_bp)*width_bp+1
 data.frame(marker_id=as.character(markers$marker_id),chromosome=as.character(markers$chromosome),
            position_bp=as.numeric(markers$position_bp),region_id=paste0(markers$chromosome,":",sprintf("%.0f",start),"-",sprintf("%.0f",start+width_bp-1)),
            region_start=start,region_end=start+width_bp-1,stringsAsFactors=FALSE)
}

breeding_coordinate_universe<-function(source_map,width_bp=10000000){
 m<-breeding_coordinate_map(source_map,width_bp)
 u<-unique(m[,c("region_id","chromosome","region_start","region_end")])
 u<-u[order(u$chromosome,u$region_start),];rownames(u)<-NULL
 u$source_marker_count<-as.integer(table(factor(m$region_id,levels=u$region_id)))
 u
}

breeding_coordinate_groups<-function(marker_ids,fitted_map,source_map,width_bp=10000000){
 if(anyDuplicated(marker_ids)||!identical(as.character(fitted_map$marker_id),marker_ids))stop("Fitted coordinate map does not preserve model marker identity.")
 source<-breeding_coordinate_map(source_map,width_bp);m<-breeding_coordinate_map(fitted_map,width_bp)
 idx<-match(marker_ids,source$marker_id)
 if(anyNA(idx)||!identical(m,source[idx,,drop=FALSE])&&
    !isTRUE(all.equal(m,source[idx,,drop=FALSE],check.attributes=FALSE)))stop("Fitted marker coordinates differ from the fixed source map.")
 universe<-breeding_coordinate_universe(source_map,width_bp)
 groups<-setNames(lapply(universe$region_id,function(id)which(m$region_id==id)),universe$region_id)
 counts<-lengths(groups);groups<-groups[counts>0]
 universe$fitted_marker_count<-as.integer(counts);universe$active_in_fitted_model<-counts>0
 if(!length(groups)||anyDuplicated(unlist(groups))||!setequal(unlist(groups),seq_along(marker_ids)))stop("Coordinate partition does not cover every fitted marker exactly once.")
 list(groups=groups,coverage=universe,marker_map=m)
}

breeding_exact_background_groups<-function(model,training,query,background,groups,tolerance=1e-8){
 if(!identical(rownames(training),model$training_ids)||!identical(colnames(training),model$marker_ids)||
    !identical(colnames(query),model$marker_ids)||!identical(colnames(background),model$marker_ids)||
    !nrow(query)||!nrow(background)||anyDuplicated(rownames(background))||
    !all(rownames(background)%in%rownames(training))||
    !identical(background,training[rownames(background),,drop=FALSE])||
    any(!is.finite(query))||anyDuplicated(unlist(groups))||
    !setequal(unlist(groups),seq_len(ncol(query))))stop("Invalid exact background-game identities.")
 linear<-breeding_explanation_linear(model)
 marker<-sweep(sweep(query,2,colMeans(background),"-"),2,linear$beta,"*")
 phi<-matrix(unlist(lapply(groups,function(columns)rowSums(marker[,columns,drop=FALSE])),use.names=FALSE),nrow=nrow(query),dimnames=list(rownames(query),names(groups)))
 baseline<-rep(mean(breeding_predict_learner(model,background)),nrow(query))
 prediction<-breeding_predict_learner(model,query)
 error<-max(abs(rowSums(phi)+baseline-prediction))
 if(!is.finite(error)||error>tolerance)stop("Exact background-game decomposition differs from saved prediction.")
 summary<-data.frame(group_id=names(groups),n_markers=as.integer(lengths(groups)),mean_absolute_shapley=colMeans(abs(phi)),mean_signed_shapley=colMeans(phi),mean_local_MCSE=0,split_half_first_mean_abs=colMeans(abs(phi)),split_half_second_mean_abs=colMeans(abs(phi)),is_other_regions=FALSE)
 summary<-summary[order(-summary$mean_absolute_shapley,summary$group_id),]
 list(contributions=phi,MCSE=phi*0,baseline=baseline,prediction=prediction,summary=summary,
      maximum_additivity_error=error,maximum_path_error=error,paired_permutations=0L,paths=0L,
      background_ids=rownames(background),query_ids=rownames(query),
      estimated_working_gib=8*length(marker)*4/1024^3,
      method="exact_linear_group_Shapley_fixed_empirical_background",
      uncertainty_scope="no_permutation_noise_exact_linear_game_background_and_model_fixed",refitted=FALSE)
}

breeding_coordinate_shapley<-function(model,training,query,background,groups,pairs=16L,seed=1L,
                                     checkpoint,tolerance=1e-8,maximum_working_gib=2,
                                     maximum_pair_evidence_gib=.25){
 pair_gib<-8*nrow(query)*length(groups)*pairs*3/1024^3
 if(!is.finite(pair_gib)||pair_gib>maximum_pair_evidence_gib)stop("Retained permutation evidence exceeds the memory estimate limit.")
 if(model$learner%in%c("gblup_reml","elastic_net"))return(breeding_exact_background_groups(model,training,query,background,groups,tolerance))
 if (identical(model$learner, "random_forest")) {
  # ranger stores the fitting call, which can embed live R objects. Response
  # prediction uses the fitted forest and prediction settings, not that call.
  # Keep the persisted model intact and hash only the prediction copy here.
  model$fit$call <- NULL
 }
 result<-breeding_group_shapley(model,training,query,background,groups,pairs,seed,checkpoint,tolerance,maximum_working_gib)
 draws<-lapply(seq_len(pairs),function(i){
  path<-file.path(checkpoint,sprintf("pair_%04d.rds",i))
  if(!file.exists(paste0(path,".sha256"))||readLines(paste0(path,".sha256"))!=digest::digest(path,file=TRUE,algo="sha256"))stop("Permutation evidence checksum differs.")
  d<-readRDS(path)
  if(!identical(dimnames(d$phi),dimnames(result$contributions))||any(!is.finite(d$phi))||
     max(abs(d$baseline-result$baseline),abs(d$prediction-result$prediction),
         abs(rowSums(d$phi)+d$baseline-d$prediction))>tolerance)stop("Permutation evidence identity/endpoint differs.")
  d$phi
 })
 result$pair_contributions<-array(unlist(draws,use.names=FALSE),
  dim=c(nrow(query),length(groups),pairs),dimnames=c(dimnames(result$contributions),list(paste0("pair",seq_len(pairs)))))
 result$estimated_pair_evidence_gib<-pair_gib
 breeding_validate_coordinate_result(result,groups,tolerance)
 result
}

breeding_validate_coordinate_result<-function(result,groups,tolerance=1e-8){
 phi<-result$contributions;nq<-length(result$query_ids)
 if(!identical(dimnames(phi),list(result$query_ids,names(groups)))||
    length(result$baseline)!=nq||length(result$prediction)!=nq||
    any(!is.finite(result$baseline))||any(!is.finite(result$prediction))||
    !identical(dimnames(result$MCSE),dimnames(phi))||any(!is.finite(phi))||
    any(!is.finite(result$MCSE))||any(result$MCSE<0)||
    max(abs(rowSums(phi)+result$baseline-result$prediction))>tolerance)stop("Invalid coordinate-Shapley decomposition.")
 s<-result$summary;idx<-match(names(groups),s$group_id)
 if(anyNA(idx)||nrow(s)!=length(groups)||anyDuplicated(s$group_id)||
    !identical(as.integer(s$n_markers[idx]),as.integer(lengths(groups)))||
    max(abs(s$mean_absolute_shapley[idx]-colMeans(abs(phi))),
        abs(s$mean_signed_shapley[idx]-colMeans(phi)),
        abs(s$mean_local_MCSE[idx]-colMeans(result$MCSE)))>1e-12)stop("Global coordinate-Shapley summaries differ from local contributions.")
 if(result$paired_permutations>0){
  draws<-result$pair_contributions;p<-result$paired_permutations
  if(!identical(dim(draws),c(nrow(phi),ncol(phi),p))||any(!is.finite(draws)))stop("Missing or invalid permutation-pair evidence.")
  mean_phi<-apply(draws,c(1,2),mean);dimnames(mean_phi)<-dimnames(phi)
  mcse<-apply(draws,c(1,2),stats::sd)/sqrt(p)
  half<-floor(p/2);first<-apply(draws[,,seq_len(half),drop=FALSE],c(1,2),mean)
  second<-apply(draws[,,seq.int(half+1L,p),drop=FALSE],c(1,2),mean)
  if(max(abs(mean_phi-phi),abs(mcse-result$MCSE),
         abs(colMeans(abs(first))-s$split_half_first_mean_abs[idx]),
         abs(colMeans(abs(second))-s$split_half_second_mean_abs[idx]))>1e-12)stop("Independently recomputed pair means/MCSE/split halves differ.")
 }else if(any(result$MCSE!=0)||
           max(abs(s$split_half_first_mean_abs[idx]-colMeans(abs(phi))),
               abs(s$split_half_second_mean_abs[idx]-colMeans(abs(phi))))>1e-12)stop("Exact linear decomposition has spurious sampling uncertainty.")
 invisible(TRUE)
}

breeding_coordinate_precision<-function(result,top=10L){
 s<-result$summary;k<-min(top,nrow(s));selected<-head(s,k)
 ratio<-selected$mean_local_MCSE/pmax(selected$mean_absolute_shapley,.Machine$double.eps)
 first<-s$group_id[order(-s$split_half_first_mean_abs,s$group_id)]
 second<-s$group_id[order(-s$split_half_second_mean_abs,s$group_id)]
 rho<-if(stats::sd(s$split_half_first_mean_abs)>0&&stats::sd(s$split_half_second_mean_abs)>0)stats::cor(s$split_half_first_mean_abs,s$split_half_second_mean_abs,method="spearman") else NA_real_
 data.frame(top=k,top_overlap=length(intersect(head(first,k),head(second,k))),
            rank_spearman=rho,top_mean_local_MCSE_to_mean_abs=mean(ratio),
            top_max_local_MCSE_to_mean_abs=max(ratio),
            constant_region_importance=all(s$mean_absolute_shapley==s$mean_absolute_shapley[1]))
}
