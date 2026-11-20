# Explanation-only native TreeSHAP validation. No contributions, predictions or
# native bias values are adjusted to force additivity.
breeding_native_shap_audit<-function(native,predicted,absolute_floor=1e-5,relative_tolerance=2e-6,absolute_cap=1e-4){
 policy<-c(absolute_floor,relative_tolerance,absolute_cap)
 if(length(policy)!=3L||any(!is.finite(policy))||!is.matrix(native)||!nrow(native)||nrow(native)!=length(predicted)||ncol(native)<2L||any(!is.finite(native))||any(!is.finite(predicted))||absolute_floor<=0||relative_tolerance<0||absolute_cap<absolute_floor)stop("Invalid native TreeSHAP validation inputs.")
 reconstructed<-rowSums(native);residual<-reconstructed-predicted
 allowed<-pmin(absolute_cap,absolute_floor+relative_tolerance*pmax(1,abs(predicted),abs(reconstructed)))
 data.frame(predicted=predicted,native_sum=reconstructed,native_bias=native[,ncol(native)],signed_residual=residual,absolute_residual=abs(residual),allowed_residual=allowed,legacy_absolute_pass=abs(residual)<=absolute_floor,scale_aware_pass=abs(residual)<=allowed)
}
breeding_xgb_precision_contributions<-function(model,z,training_means,chunk_rows=64L,maximum_gib=.25,tolerance=1e-5,maximum_output_gib=1,audit_directory,relative_tolerance=2e-6,absolute_cap=1e-4){
 limits<-c(chunk_rows,maximum_gib,maximum_output_gib)
 if(model$learner!="xgboost"||!is.matrix(z)||is.null(rownames(z))||anyDuplicated(rownames(z))||anyNA(rownames(z))||!identical(colnames(z),model$marker_ids)||!identical(names(training_means),model$marker_ids)||!nrow(z)||any(!is.finite(z))||any(!is.finite(training_means))||length(limits)!=3L||any(!is.finite(limits))||any(limits<=0)||chunk_rows!=as.integer(chunk_rows))stop("Saved XGBoost attribution schema or memory limits differ.")
 chunk_gib<-8*min(chunk_rows,nrow(z))*(ncol(z)+1)*4/1024^3;retained_gib<-8*nrow(z)*(ncol(z)+1)*3/1024^3
 if(chunk_gib>maximum_gib||retained_gib>maximum_output_gib)stop("Native attribution memory estimate exceeds budget.")
 dir.create(audit_directory,recursive=TRUE,showWarnings=FALSE)
 role<-if(all(rownames(z)%in%model$training_ids))"training" else if(!any(rownames(z)%in%model$training_ids))"heldout_genotype" else stop("Training/query attribution rows were mixed.")
 booster<-xgboost::xgb.load.raw(model$fit$serialized_booster)
 sums<-abs_sums<-setNames(numeric(ncol(z)),colnames(z));values<-list();audits<-list();bias<-numeric(nrow(z))
 for(i in seq_len(ceiling(nrow(z)/chunk_rows))){
  rows<-seq.int((i-1L)*chunk_rows+1L,min(nrow(z),i*chunk_rows));q<-z[rows,,drop=FALSE];dm<-xgboost::xgb.DMatrix(q)
  predicted<-as.numeric(stats::predict(booster,newdata=dm))
  native<-as.matrix(stats::predict(booster,newdata=dm,predcontrib=TRUE,approxcontrib=FALSE))
  if(!identical(dim(native),c(nrow(q),ncol(q)+1L)))stop("Native TreeSHAP dimensions differ.")
  if(!is.null(colnames(native))&&!identical(colnames(native)[seq_len(ncol(q))],colnames(q)))stop("Native TreeSHAP marker order differs.")
  audit<-breeding_native_shap_audit(native,predicted,tolerance,relative_tolerance,absolute_cap);audit$sample_id<-rownames(q);audit$chunk<-i;audits[[i]]<-audit
  # Persist the measured discrepancy even if the tight, bounded guard rejects it.
  write.table(do.call(rbind,audits),file.path(audit_directory,paste0(role,"_native_precision.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
  if(!all(audit$scale_aware_pass))stop(sprintf("Native TreeSHAP residual exceeds bounded precision policy: max error %.9g; max permitted %.9g; %d/%d rows failed (chunk %d).",max(audit$absolute_residual),max(audit$allowed_residual),sum(!audit$scale_aware_pass),nrow(audit),i))
  contribution<-native[,seq_len(ncol(z)),drop=FALSE];colnames(contribution)<-colnames(z);rownames(contribution)<-rownames(q)
  values[[i]]<-contribution;bias[rows]<-audit$native_bias;sums<-sums+colSums(contribution);abs_sums<-abs_sums+colSums(abs(contribution))
 }
 if(diff(range(bias))>1e-8)stop("Unexpected row-varying native bias for this fixed regression model.")
 audit<-do.call(rbind,audits)
 list(summary=data.frame(marker_id=colnames(z),mean_absolute_contribution=as.numeric(abs_sums/nrow(z)),mean_signed_contribution=as.numeric(sums/nrow(z))),contributions=do.call(rbind,values),bias=bias,maximum_additivity_error=max(audit$absolute_residual),estimated_chunk_gib=chunk_gib,estimated_retained_gib=retained_gib,full_interaction_array_gib=8*nrow(z)*(ncol(z)+1)^2/1024^3,method="exact_native_XGBoost_TreeSHAP_with_original_bias_and_bounded_FP32_validation",precision_policy=list(absolute_floor=tolerance,relative_tolerance=relative_tolerance,absolute_cap=absolute_cap,native_values_adjusted=FALSE,legacy_failed_rows=sum(!audit$legacy_absolute_pass),all_scale_aware_checks_passed=all(audit$scale_aware_pass)))
}
