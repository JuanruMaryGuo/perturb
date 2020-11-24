library(Biostrings)
library(clusterProfiler)
library(devtools)
library(SAVER)
library(MUSIC)
library(Seurat)

setwd("/Users/guojuanru/Desktop/Yale/MUSIC/mock")

crop_seq_list_mock<-Input_preprocess_10X_modified("/Users/guojuanru/Desktop/Yale/MUSIC/mock/")

save(crop_seq_list_mock,file = "crop_seq_list_mock.Rdata")
expression_profile = crop_seq_list_mock$expression_profile
perturb_information = crop_seq_list_mock$perturb_information

dim(expression_profile)

expression_profile[1:3,1:3]

# cell quality control
crop_seq_qc<-Cell_qc_modified(crop_seq_list_mock$expression_profile,crop_seq_list_mock$perturb_information,plot=T)
save(crop_seq_qc,file = "crop_seq_qc.Rdata")

# data imputation, it may take a little long time without parallel computation.
crop_seq_imputation<-Data_imputation(crop_seq_qc$expression_profile,crop_seq_qc$perturb_information,cpu_num=5)# cell filtering, it may take a little long time without parallel computation.
save(crop_seq_imputation,file = "crop_seq_imputation.Rdata")

crop_seq_filtered<-Cell_filtering_modified(crop_seq_imputation$expression_profile,crop_seq_imputation$perturb_information,cpu_num=6,cell_num_threshold = 15,plot = TRUE)# obtain highly dispersion differentially expressed genes.
save(crop_seq_filtered,file = "crop_seq_filtered.Rdata")

crop_seq_vargene<-Get_high_varGenes_modified(crop_seq_filtered$expression_profile,crop_seq_filtered$perturb_information,plot=T)# get topics. 
save(crop_seq_vargene,file = "crop_seq_vargene.Rdata")

topic_model_list<-Get_topics_modified(crop_seq_vargene$expression_profile,crop_seq_vargene$perturb_information,topic_number=c(4:6))
save(topic_model_list,file = "topic_model_list.Rdata")

# select the optimal topic number.  
optimalModel<-Select_topic_number(topic_model_list$models,plot=T)

#If you just calculated one topic number, you can skip this step, just run the following:
optimalModel<-topic_model_list$models[[1]]
save(optimalModel,file = "optimalModel.Rdata")

# calculate topic distribution for each cell.
distri_diff<-Diff_topic_distri_modified(optimalModel,topic_model_list$perturb_information,plot=T)
write.csv(distri_diff, file = "distri_diff.csv")

# calculate the overall perturbation effect ranking list without "offTarget_Info".
rank_overall_result<-Rank_overall_modified(distri_diff)
write.csv(rank_overall_result, file = "rank_overall_result.csv")

# calculate the topic-specific ranking list.
rank_topic_specific_result<-Rank_specific_modified(distri_diff)
write.csv(rank_topic_specific_result, file = "rank_topic_specific_result.csv")

# calculate the perturbation correlation.
perturb_cor<-Correlation_perturbation(distri_diff,plot=T)
write.csv(perturb_cor, file = "perturb_cor.csv")

setwd("/Users/guojuanru/Desktop/Yale/MUSIC/sars2_result")

crop_seq_list_sars<-Input_preprocess_10X_modified("/Users/guojuanru/Desktop/Yale/MUSIC/sars2")


save(crop_seq_list_sars,file = "crop_seq_list_sars.Rdata")
expression_profile = crop_seq_list_sars$expression_profile
perturb_information = crop_seq_list_sars$perturb_information

dim(expression_profile)

expression_profile[1:3,1:3]

# cell quality control
crop_seq_qc<-Cell_qc_modified(crop_seq_list_sars$expression_profile,crop_seq_list_sars$perturb_information,plot=T)
save(crop_seq_qc,file = "crop_seq_qc.Rdata")

# data imputation, it may take a little long time without parallel computation.
crop_seq_imputation<-Data_imputation(crop_seq_qc$expression_profile,crop_seq_qc$perturb_information,cpu_num=4)# cell filtering, it may take a little long time without parallel computation.
save(crop_seq_imputation,file = "crop_seq_imputation.Rdata")

crop_seq_filtered<-Cell_filtering_modified(crop_seq_imputation$expression_profile,crop_seq_imputation$perturb_information,cpu_num=4,cell_num_threshold = 15,plot = TRUE)# obtain highly dispersion differentially expressed genes.
save(crop_seq_filtered,file = "crop_seq_filtered.Rdata")

crop_seq_vargene<-Get_high_varGenes_modified(crop_seq_filtered$expression_profile,crop_seq_filtered$perturb_information,plot=T)# get topics. 
save(crop_seq_vargene,file = "crop_seq_vargene.Rdata")

topic_model_list<-Get_topics_modified(crop_seq_vargene$expression_profile,crop_seq_vargene$perturb_information,topic_number=c(4:6))
save(topic_model_list,file = "topic_model_list.Rdata")

# select the optimal topic number.  
optimalModel<-Select_topic_number(topic_model_list$models,plot=T)

#If you just calculated one topic number, you can skip this step, just run the following:
optimalModel<-topic_model_list$models[[1]]
save(optimalModel,file = "optimalModel.Rdata")

# calculate topic distribution for each cell.
distri_diff<-Diff_topic_distri_modified(optimalModel,topic_model_list$perturb_information,plot=T)
write.csv(distri_diff, file = "distri_diff.csv")

# calculate the overall perturbation effect ranking list without "offTarget_Info".
rank_overall_result<-Rank_overall_modified(distri_diff)
write.csv(rank_overall_result, file = "rank_overall_result.csv")

# calculate the topic-specific ranking list.
rank_topic_specific_result<-Rank_specific_modified(distri_diff)
write.csv(rank_topic_specific_result, file = "rank_topic_specific_result.csv")

# calculate the perturbation correlation.
perturb_cor<-Correlation_perturbation(distri_diff,plot=T)
write.csv(perturb_cor, file = "perturb_cor.csv")

