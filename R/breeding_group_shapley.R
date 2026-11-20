# Permutation-sampled group Shapley values for a fixed empirical training
# background. No fitting, response values, GWAS results or conditional imputer.
breeding_group_partition<-function(marker_ids,region_map,focus_regions){
 if(!identical(region_map$marker_id,marker_ids)||anyDuplicated(marker_ids)||anyNA(region_map$region_id)||anyDuplicated(focus_regions)||!all(focus_regions%in%region_map$region_id))stop("Invalid complete marker/group partition.")
 groups<-setNames(lapply(sort(focus_regions),function(id)which(region_map$region_id==id)),sort(focus_regions))
 remaining<-which(!region_map$region_id%in%focus_regions)
 if(length(remaining))groups[["OTHER_REGIONS"]]<-remaining
 if(!length(groups)||anyDuplicated(unlist(groups))||!setequal(unlist(groups),seq_along(marker_ids)))stop("Every fitted marker must belong to exactly one Shapley group.")
 groups
}
breeding_group_path<-function(model,training,query,background,groups,permutation,tolerance=1e-8,reference_context=NULL,target_prediction=NULL){
 if(!identical(sort(permutation),seq_along(groups))||!identical(colnames(query),colnames(background))||!identical(colnames(query),model$marker_ids)||nrow(background)<1L||nrow(query)<1L)stop("Invalid Shapley path inputs.")
 nq<-nrow(query);nb<-nrow(background);qi<-rep(seq_len(nq),each=nb);bi<-rep(seq_len(nb),times=nq)
 current<-background[bi,,drop=FALSE];rownames(current)<-paste0("coalition_",seq_len(nrow(current)))
 context<-if(is.null(reference_context))breeding_explanation_prepare_model(model,background,training,tolerance) else reference_context
 previous<-context$baseline[bi]
 if(model$learner=="gaussian_rkhs")context$distances<-context$distances[bi,,drop=FALSE]
 baseline<-colMeans(matrix(previous,nrow=nb));phi<-matrix(0,nq,length(groups),dimnames=list(rownames(query),names(groups)))
 for(g in permutation){
  columns<-groups[[g]];replacement<-query[qi,columns,drop=FALSE]
  if(model$learner=="gaussian_rkhs"){
   old<-breeder_squared_distances(current[,columns,drop=FALSE],training[,columns,drop=FALSE]);new<-breeder_squared_distances(replacement,training[,columns,drop=FALSE])
   context$distances<-context$distances-old+new
   if(any(context$distances< -1e-7))stop("Invalid incremental coalition distance.")
   value<-as.numeric(model$fit$intercept+exp(-model$fit$gamma*pmax(context$distances,0)/model$fit$distance_scale)%*%model$fit$alpha)
  }else if(model$learner%in%c("gblup_reml","elastic_net"))value<-as.numeric(previous+(replacement-current[,columns,drop=FALSE])%*%context$linear$beta[columns]) else {
   current[,columns]<-replacement
   value<-if(model$learner=="xgboost")as.numeric(stats::predict(context$booster,newdata=xgboost::xgb.DMatrix(current))) else breeding_predict_learner(model,current)
  }
  phi[,g]<-colMeans(matrix(value-previous,nrow=nb));current[,columns]<-replacement;previous<-value
 }
 target<-if(is.null(target_prediction))breeding_predict_learner(model,query) else target_prediction;error<-max(abs(colMeans(matrix(previous,nrow=nb))-target),abs(rowSums(phi)+baseline-target))
 if(!is.finite(error)||error>tolerance)stop("Shapley path fails endpoint prediction/additivity checks.")
 list(phi=phi,baseline=baseline,prediction=target,maximum_error=error)
}
breeding_group_shapley<-function(model,training,query,background,groups,pairs=8L,seed=1L,checkpoint,tolerance=1e-8,maximum_working_gib=2){
 if(pairs<2L||pairs!=as.integer(pairs)||!identical(rownames(training),model$training_ids)||!identical(colnames(training),model$marker_ids)||!identical(colnames(query),model$marker_ids)||!identical(colnames(background),model$marker_ids)||!all(rownames(background)%in%rownames(training))||anyDuplicated(rownames(background))||!identical(background,training[rownames(background),,drop=FALSE])||any(!is.finite(query))||anyDuplicated(unlist(groups))||!setequal(unlist(groups),seq_len(ncol(query))))stop("Invalid saved-model Shapley reference, query or partition.")
 estimate<-8*nrow(query)*nrow(background)*(ncol(query)*6+nrow(training)*4)/1024^3+8*nrow(query)*length(groups)*pairs/1024^3
 if(!is.finite(estimate)||estimate>maximum_working_gib)stop("Coalition working-memory estimate exceeds budget.")
 dir.create(checkpoint,recursive=TRUE,showWarnings=FALSE)
 permutations<-with_seed_preserved(seed,lapply(seq_len(pairs),function(i)sample.int(length(groups))))
 request<-list(model_hash=digest::digest(model,algo="sha256"),training_hash=digest::digest(training,algo="sha256"),query=query,background=background,groups=groups,permutations=permutations,tolerance=tolerance,algorithm_hash=digest::digest(list(body(breeding_group_path),body(breeding_group_shapley)),algo="sha256"))
 path<-file.path(checkpoint,"request.rds");if(file.exists(path)){if(!identical(readRDS(path),request))stop("Changed group-Shapley resume request.")}else gbs_parallel_atomic(request,path)
 draws<-vector("list",pairs)
 reference_context<-breeding_explanation_prepare_model(model,background,training,tolerance);target_prediction<-breeding_predict_learner(model,query)
 for(i in seq_len(pairs)){
  path<-file.path(checkpoint,sprintf("pair_%04d.rds",i))
  if(file.exists(path)){
   if(!file.exists(paste0(path,".sha256"))||readLines(paste0(path,".sha256"))!=digest::digest(path,file=TRUE,algo="sha256"))stop("Group-Shapley checkpoint hash differs.")
   draws[[i]]<-readRDS(path)
  }else{
   forward<-breeding_group_path(model,training,query,background,groups,permutations[[i]],tolerance,reference_context,target_prediction)
   reverse<-breeding_group_path(model,training,query,background,groups,rev(permutations[[i]]),tolerance,reference_context,target_prediction)
   draws[[i]]<-list(phi=(forward$phi+reverse$phi)/2,baseline=forward$baseline,prediction=forward$prediction,maximum_error=max(forward$maximum_error,reverse$maximum_error))
   gbs_parallel_atomic(draws[[i]],path);writeLines(digest::digest(path,file=TRUE,algo="sha256"),paste0(path,".sha256"))
  }
  if(!identical(dimnames(draws[[i]]$phi),list(rownames(query),names(groups))))stop("Group-Shapley checkpoint identities differ.")
 }
 mean_phi<-Reduce(`+`,lapply(draws,`[[`,"phi"))/pairs
 mcse<-sqrt(Reduce(`+`,lapply(draws,function(d)(d$phi-mean_phi)^2))/(pairs-1)/pairs)
 baseline<-draws[[1]]$baseline;prediction<-draws[[1]]$prediction
 additivity<-max(abs(rowSums(mean_phi)+baseline-prediction))
 if(additivity>tolerance)stop("Averaged Shapley decomposition does not reconstruct predictions.")
 half<-floor(pairs/2);first<-Reduce(`+`,lapply(draws[seq_len(half)],`[[`,"phi"))/half;second<-Reduce(`+`,lapply(draws[seq.int(half+1L,pairs)],`[[`,"phi"))/(pairs-half)
 summary<-data.frame(group_id=names(groups),n_markers=lengths(groups),mean_absolute_shapley=colMeans(abs(mean_phi)),mean_signed_shapley=colMeans(mean_phi),mean_local_MCSE=colMeans(mcse),split_half_first_mean_abs=colMeans(abs(first)),split_half_second_mean_abs=colMeans(abs(second)),is_other_regions=names(groups)=="OTHER_REGIONS")
 summary<-summary[order(-summary$mean_absolute_shapley,summary$group_id),]
 list(contributions=mean_phi,MCSE=mcse,baseline=baseline,prediction=prediction,summary=summary,maximum_additivity_error=additivity,maximum_path_error=max(vapply(draws,`[[`,numeric(1L),"maximum_error")),paired_permutations=pairs,paths=2L*pairs,background_ids=rownames(background),query_ids=rownames(query),estimated_working_gib=estimate,method="permutation_sampled_group_Shapley_fixed_empirical_background",uncertainty_scope="Monte_Carlo_over_permutation_pairs_conditional_on_fixed_background_and_groups",refitted=FALSE)
}
breeding_plot_group_shapley<-function(result,display_groups,learner,output){
 if(!length(display_groups)||anyDuplicated(display_groups)||!all(display_groups%in%colnames(result$contributions)))stop("Invalid preselected Shapley display groups.")
 dir.create(output,recursive=TRUE,showWarnings=FALSE)
 grDevices::png(file.path(output,"region_importance.png"),width=2100,height=1450,res=150)
 tryCatch({par(mar=c(5,14,4,2));x<-result$summary[match(rev(display_groups),result$summary$group_id),];barplot(x$mean_absolute_shapley,horiz=TRUE,names.arg=x$group_id,las=1,cex.names=.55,col="#287B8E",xlab="Mean absolute region Shapley contribution",main=paste(learner,"| training-selected regions"))},finally=dev.off())
 i<-1L;v<-result$contributions[i,display_groups];remaining<-setdiff(colnames(result$contributions),display_groups)
 if(length(remaining))v<-c(v,"Remaining groups (sum)"=sum(result$contributions[i,remaining]))
 grDevices::png(file.path(output,"local_contributions.png"),width=2100,height=1450,res=150)
 tryCatch({par(mar=c(5,14,5,2));v<-rev(v);barplot(v,horiz=TRUE,names.arg=names(v),las=1,cex.names=.55,col=ifelse(v>=0,"#B2182B","#2166AC"),xlab="Signed contribution to this line's prediction (trait units)",main=paste(learner,"|",result$query_ids[i]));abline(v=0,col="grey50");mtext(sprintf("Fixed background prediction %.4f + all contributions = prediction %.4f",result$baseline[i],result$prediction[i]),side=3,line=.5,cex=.8)},finally=dev.off())
 invisible(c("region_importance.png","local_contributions.png"))
}
