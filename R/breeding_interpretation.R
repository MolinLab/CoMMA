# Saved-model explanations only. No outcomes, fitting, tuning or GWAS are inputs.
breeding_interpretation_require<-function(){
 for(p in c("glmnet","ranger","xgboost"))if(!requireNamespace(p,quietly=TRUE))stop("Saved-model explanation requires installed package: ",p)
 invisible(TRUE)
}
breeding_explanation_ids<-function(ids,metadata,n,seed){
 if(anyDuplicated(ids)||anyNA(ids)||anyNA(match(ids,metadata$sample_id))||anyDuplicated(metadata$sample_id)||n<1L)stop("Invalid explanation sampling identities.")
 ids<-sort(ids);m<-metadata[match(ids,metadata$sample_id),]
 groups<-interaction(m$source_group,m$genetic_stratum,drop=TRUE,lex.order=TRUE)
 with_seed_preserved(seed,{
  queues<-lapply(split(ids,groups),function(x)sample(x,length(x)));out<-character()
  while(length(out)<min(n,length(ids)))for(i in seq_along(queues))if(length(queues[[i]])&&length(out)<min(n,length(ids))){out<-c(out,queues[[i]][1]);queues[[i]]<-queues[[i]][-1]}
  out
 })
}
breeding_explanation_regions<-function(markers,width_bp){
 if(anyDuplicated(markers$marker_id)||anyNA(markers$chromosome)||any(!is.finite(markers$position_bp))||width_bp<=0)stop("Invalid physical region map.")
 bin<-floor((markers$position_bp-1)/width_bp);x<-data.frame(marker_id=markers$marker_id,chromosome=as.character(markers$chromosome),position_bp=markers$position_bp,region_id=paste0(markers$chromosome,":",bin*width_bp+1,"-",(bin+1)*width_bp),region_start=bin*width_bp+1,region_end=(bin+1)*width_bp)
 x
}
breeding_explanation_linear<-function(model){
 if(model$learner=="gblup_reml"){
  beta<-model$fit$marker_effects;intercept<-model$fit$intercept
 }else if(model$learner=="elastic_net"){
  cf<-as.matrix(stats::coef(model$fit$fit,s=model$fit$lambda));intercept<-as.numeric(cf[1,1]);beta<-setNames(as.numeric(cf[-1,1]),rownames(cf)[-1])
 }else stop("Exact additive contributions require a linear fitted learner.")
 if(!identical(names(beta),model$marker_ids)||any(!is.finite(beta))||!is.finite(intercept))stop("Linear saved-model marker coefficients differ.")
 list(beta=beta,intercept=intercept)
}
breeding_explanation_prepare_model<-function(model,query,training,tolerance=1e-8){
 if(!identical(colnames(query),model$marker_ids)||!identical(colnames(training),model$marker_ids)||!identical(rownames(training),model$training_ids)||any(!is.finite(query))||any(!is.finite(training)))stop("Explanation matrices differ from fitted sample/marker identities.")
 context<-list(model=model,query=query,training=training)
 context$baseline<-breeding_predict_learner(model,query)
 if(model$learner%in%c("gblup_reml","elastic_net")){
  context$linear<-breeding_explanation_linear(model);v<-as.numeric(context$linear$intercept+query%*%context$linear$beta)
  if(max(abs(v-context$baseline))>tolerance)stop("Linear fast prediction differs from saved learner.")
 }else if(model$learner=="gaussian_rkhs"){
  if(!identical(model$fit$z_train,training))stop("RKHS training reference differs.")
  context$distances<-breeder_squared_distances(query,training)
 }else if(model$learner=="xgboost")context$booster<-xgboost::xgb.load.raw(model$fit$serialized_booster)
 context
}
breeding_explanation_changed_prediction<-function(context,columns,donor_rows){
 q<-context$query;t<-context$training;model<-context$model
 if(length(donor_rows)!=nrow(q)||anyNA(donor_rows)||any(donor_rows<1|donor_rows>nrow(t))||anyDuplicated(columns)||any(!columns%in%seq_len(ncol(q))))stop("Invalid regional donor or marker indices.")
 replacement<-t[donor_rows,columns,drop=FALSE]
 if(model$learner%in%c("gblup_reml","elastic_net"))return(as.numeric(context$baseline+(replacement-q[,columns,drop=FALSE])%*%context$linear$beta[columns]))
 if(model$learner=="gaussian_rkhs"){
  original<-breeder_squared_distances(q[,columns,drop=FALSE],t[,columns,drop=FALSE]);changed<-breeder_squared_distances(replacement,t[,columns,drop=FALSE])
  distances<-context$distances-original+changed
  if(any(distances< -1e-7))stop("Updated RKHS distances became invalid.")
  kernel<-exp(-model$fit$gamma*pmax(distances,0)/model$fit$distance_scale)
  return(as.numeric(model$fit$intercept+kernel%*%model$fit$alpha))
 }
 q[,columns]<-replacement
 if(model$learner=="xgboost")return(as.numeric(stats::predict(context$booster,newdata=xgboost::xgb.DMatrix(q))))
 breeding_predict_learner(model,q)
}
breeding_regional_sensitivity<-function(model,training,query,regions,donors,checkpoint,tolerance=1e-8){
 context<-breeding_explanation_prepare_model(model,query,training,tolerance)
 if(!identical(regions$marker_id,colnames(query))||nrow(donors)!=nrow(query))stop("Region/donor identities do not align.")
 dir.create(checkpoint,recursive=TRUE,showWarnings=FALSE)
 request<-list(model_sha256=digest::digest(model,algo="sha256"),training_sha256=digest::digest(training,algo="sha256"),query_sha256=digest::digest(query,algo="sha256"),regions=regions,donors=donors,method="training_donor_grouped_prediction_sensitivity",tolerance=tolerance,implementation_sha256=digest::digest(list(body(breeding_explanation_prepare_model),body(breeding_explanation_changed_prediction),body(breeding_regional_sensitivity)),algo="sha256"))
 path<-file.path(checkpoint,"request.rds")
 if(file.exists(path)){if(!identical(readRDS(path),request))stop("Regional sensitivity resume inputs changed.")}else gbs_parallel_atomic(request,path)
 groups<-split(seq_len(nrow(regions)),regions$region_id);rows<-vector("list",length(groups))
 for(i in seq_along(groups)){
  path<-file.path(checkpoint,sprintf("region_%05d.rds",i));columns<-groups[[i]]
  if(file.exists(path)){
   if(!file.exists(paste0(path,".sha256"))||readLines(paste0(path,".sha256"))!=digest::digest(path,file=TRUE,algo="sha256"))stop("Regional checkpoint hash differs.")
   rows[[i]]<-readRDS(path)
  }else {
   changes<-vapply(seq_len(ncol(donors)),function(r)breeding_explanation_changed_prediction(context,columns,donors[,r])-context$baseline,numeric(nrow(query)))
   rows[[i]]<-data.frame(learner=model$learner,region_id=names(groups)[i],chromosome=regions$chromosome[columns[1]],region_start=regions$region_start[columns[1]],region_end=regions$region_end[columns[1]],n_markers=length(columns),n_query=nrow(query),donor_replicates=ncol(donors),mean_absolute_prediction_change=mean(abs(changes)),root_mean_squared_prediction_change=sqrt(mean(changes^2)),mean_signed_prediction_change=mean(changes),stringsAsFactors=FALSE)
   gbs_parallel_atomic(rows[[i]],path);writeLines(digest::digest(path,file=TRUE,algo="sha256"),paste0(path,".sha256"))
  }
  if(rows[[i]]$learner!=model$learner||rows[[i]]$region_id!=names(groups)[i]||rows[[i]]$n_markers!=length(columns)||rows[[i]]$n_query!=nrow(query)||rows[[i]]$donor_replicates!=ncol(donors)||!is.finite(rows[[i]]$mean_absolute_prediction_change))stop("Regional checkpoint identity or sensitivity differs.")
 }
 result<-do.call(rbind,rows);result[order(-result$mean_absolute_prediction_change,result$region_id),]
}
breeding_marker_contributions<-function(model,z,training_means,chunk_rows=64L,maximum_gib=.25,tolerance=1e-5,maximum_output_gib=1){
 if(!nrow(z)||!identical(colnames(z),model$marker_ids)||!identical(names(training_means),model$marker_ids)||any(!is.finite(z))||any(!is.finite(training_means))||chunk_rows<1L)stop("Invalid attribution identities or chunk size.")
 allocation_gib<-8*min(chunk_rows,nrow(z))*(ncol(z)+1)*4/1024^3
 if(allocation_gib>maximum_gib)stop("Attribution chunk exceeds the declared memory estimate.")
 retained_gib<-8*nrow(z)*(ncol(z)+1)*3/1024^3
 if(!is.finite(maximum_output_gib)||maximum_output_gib<=0||retained_gib>maximum_output_gib)stop("Retained attribution matrices exceed the declared output-memory estimate.")
 linear<-if(model$learner%in%c("gblup_reml","elastic_net"))breeding_explanation_linear(model) else NULL
 if(is.null(linear)&&model$learner!="xgboost")stop("Marker contributions implemented for additive models and XGBoost only.")
 booster<-if(is.null(linear))xgboost::xgb.load.raw(model$fit$serialized_booster) else NULL
 sums<-abs_sums<-setNames(numeric(ncol(z)),colnames(z));values<-vector("list",ceiling(nrow(z)/chunk_rows));bias<-numeric(nrow(z));maximum_error<-0
 for(i in seq_along(values)){
  rows<-seq.int((i-1L)*chunk_rows+1L,min(nrow(z),i*chunk_rows));q<-z[rows,,drop=FALSE]
  if(is.null(linear)){
   predicted<-as.numeric(stats::predict(booster,newdata=xgboost::xgb.DMatrix(q)))
   # Use the existing project TreeSHAP primitive and its additivity check.
   native<-xgboost_shap_contributions(booster,q,predicted,tolerance);contribution<-native$marker_contributions
   maximum_error<-max(maximum_error,native$max_abs_additivity_error)
   b<-predicted-rowSums(contribution)
  }else{
   contribution<-sweep(sweep(q,2,training_means,"-"),2,linear$beta,"*")
   b<-rep(linear$intercept+sum(training_means*linear$beta),length(rows));predicted<-as.numeric(linear$intercept+q%*%linear$beta)
  }
  error<-max(abs(rowSums(contribution)+b-predicted));maximum_error<-max(maximum_error,error)
  if(!is.finite(error)||error>tolerance)stop("Saved-model attribution additivity failed.")
  sums<-sums+colSums(contribution);abs_sums<-abs_sums+colSums(abs(contribution));values[[i]]<-contribution;bias[rows]<-b
 }
 list(summary=data.frame(marker_id=colnames(z),mean_absolute_contribution=as.numeric(abs_sums/nrow(z)),mean_signed_contribution=as.numeric(sums/nrow(z))),contributions=do.call(rbind,values),bias=bias,maximum_additivity_error=maximum_error,estimated_chunk_gib=allocation_gib,estimated_retained_gib=retained_gib,full_interaction_array_gib=8*nrow(z)*(ncol(z)+1)^2/1024^3,method=if(is.null(linear))"exact_XGBoost_TreeSHAP_default_tree_path_expectations" else "exact_additive_contribution_relative_to_training_marker_means")
}
