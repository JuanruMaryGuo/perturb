Input_preprocess_10X_modified<-function(directory){
  require(hash)
  require(stringr)
  perturb_seq <- Read10X(directory)$`Gene Expression`
  perturb_seq <- as.matrix(perturb_seq)
  cbc_gbc <- read.table(paste(directory, "cbc_gbc_dict.tsv",sep = "/"), stringsAsFactors = FALSE)
  cbc_gbc <- unique(cbc_gbc)
  data_preprocess <- function(perturb_seq, cbc_gbc) {
    cell_KO_hash = hash()
    for (i in 1:nrow(cbc_gbc)) {
      cbc = cbc_gbc[i, 1]
      gbc = cbc_gbc[i, 2]
      if (has.key(cbc, cell_KO_hash)) {
        cell_KO_hash[cbc] = paste(cell_KO_hash[[cbc]],gbc, sep = ",")
      }
      else {
        cell_KO_hash[cbc] = gbc
      }
    }
    perturb_information = c()
    j = 1
    k = 1
    nogbc_col = c()
    for (i in 1:ncol(perturb_seq)) {
      if (!is.null(cell_KO_hash[[colnames(perturb_seq)[i]]])) {
        perturb_information[k] = cell_KO_hash[[colnames(perturb_seq)[i]]]
        k = k + 1
      }
      else {
        nogbc_col[j] <- i
        j = j + 1
      }
    }
    perturb_seq <- perturb_seq[, -nogbc_col]
    names(perturb_information) = colnames(perturb_seq)
    for (i in 1:length(perturb_information)) {
      sample_info_arr <- unlist(strsplit(perturb_information[i],","))
      if (length(sample_info_arr) > 1) {
        sortedMultiKO<-c()
        sample_info_arr <- sort(sample_info_arr)
        if(sample_info_arr[1]!="Non-Targeting"){
          sortedMultiKO=sample_info_arr[1]
          for (j in 2:length(sample_info_arr)) {
            if(sample_info_arr[j]!="Non-Targeting"){
              sortedMultiKO = paste(sortedMultiKO, sample_info_arr[j],sep = ",")
            }
          }
        }
        else{
          sortedMultiKO=sample_info_arr[2]
          if(length(sample_info_arr)>=3){
            for (j in 3:length(sample_info_arr)) {
              if(sample_info_arr[j]!="Non-Targeting"){
                sortedMultiKO = paste(sortedMultiKO, sample_info_arr[j],sep = ",")
              }
            }
          }
        }
        perturb_information[i] = sortedMultiKO
      }
    }
    perturb_seq <- perturb_seq[!str_detect(row.names(perturb_seq),"^MRP"), ]
    perturb_seq <- perturb_seq[!str_detect(row.names(perturb_seq),"^RP"), ]
    return(list("perturb_data" = perturb_seq, "perturb_information" = perturb_information))
  }
  perturb_seq_list <- data_preprocess(perturb_seq, cbc_gbc)
  perturb_list = list("expression_profile" = perturb_seq_list$perturb_data, "perturb_information" = perturb_seq_list$perturb_information)
  return(perturb_list)
}

Cell_qc_modified<-function(expression_profile,perturb_information,species="Hs",gene_low=500,gene_high=10000,mito_high=0.1,umi_low=1000,umi_high=Inf,plot=FALSE,plot_path="./quality_control.pdf"){
  require(stringr)
  nUMI<-apply(expression_profile,2,sum)
  nGene<-apply(expression_profile,2,function(x){length(x[x>0])})
  expression_profile<-apply(expression_profile,2,function(x){x/(sum(x)/10000)})
  expression_profile<-log(expression_profile+1)    
  filter_retain<-rep("filter",ncol(expression_profile))
  for(i in 1:length(filter_retain)){
    if(nGene[i]>gene_low & nGene[i]<gene_high & nUMI[i]>umi_low & nUMI[i]<umi_high){
      filter_retain[i]<-"retain"
    }
  }
  if(plot){
    pdf(file=plot_path)
    par(mfrow=c(1,3))
    hist(nGene,breaks = 20,freq=FALSE,xlab = "Gene numbers",ylab = "Density",main = "Gene numbers distribution")
    lines(nGene,col="red",lwd=1)
    hist(nUMI,breaks = 20,freq=FALSE,xlab = "UMI numbers",ylab = "Density",main = "UMI numbers distribution")
    lines(density(nUMI),col="red",lwd=1)
    dev.off()
  }
  SQ_filter<-as.matrix(expression_profile[,which(filter_retain=="retain")])
  SQ_data_qc<-SQ_filter#have adopted total_expr=10000 and log
  perturb_information_qc<-perturb_information[colnames(SQ_data_qc)]
  perturb_information_abandon<-perturb_information[setdiff(names(perturb_information),names(perturb_information_qc))]
  return(list("expression_profile"=SQ_data_qc,"perturb_information"=perturb_information_qc,"perturb_information_abandon"=perturb_information_abandon))
}

Cell_filtering_modified<-function(expression_profile,perturb_information,cpu_num=4,cell_num_threshold=15,umi=0.01,pvalue=0.05,vargene_min_num=5,filtered_rate=0.9,plot=FALSE,plot_path="./invalid_rate.pdf"){
  library(parallel)
  library(stringr)
  options(warn=-1)
  get_varGene<-function(ex,s){
    a=ks.test(as.numeric(ex[1:s]),as.numeric(ex[(s+1):length(ex)]))
    return(as.numeric(a$p.value))
  }
  perturb_information_delete_name<-function(perturb_information,delete_name){
    require(stringr)
    if(delete_name=="*"){
      perturb_information<-perturb_information[perturb_information!="*"]
      return(perturb_information)
    }
    delete_name_f<-paste(",",delete_name,sep="")
    delete_name_b<-paste(delete_name,",",sep="")
    perturb_information[perturb_information==delete_name]<-"wait_delete"
    for(i in 1:length(perturb_information)){
      if(str_detect(perturb_information[i],delete_name_b)){
        perturb_information[i]<-sub(delete_name_b,"",perturb_information[i])
      }
      if(str_detect(perturb_information[i],delete_name_f)){
        perturb_information[i]<-sub(delete_name_f,"",perturb_information[i])
      }
    }
    perturb_information<-perturb_information[perturb_information!="wait_delete"]
    return(perturb_information)
  }
  perturb_information_ko<-perturb_information[perturb_information!="Non-Targeting"]
  perturb_information_ko_split<-c()
  for(i in 1:length(perturb_information_ko)){
    ko_split<-unlist(str_split(perturb_information_ko[i],","))
    names(ko_split)<-rep(names(perturb_information_ko)[i],length(ko_split))
    perturb_information_ko_split<-c(perturb_information_ko_split,ko_split)
  }
  #filter KO genes who express little in ctrl samples, because it makes no sense.
  #calculate percent of zero value of control sample
  ko_names<-unique(perturb_information_ko_split)
  Zero_ra<-c()
  if(length(ko_names[!(ko_names %in% row.names(expression_profile))])>0){
    print(paste("Warning! ",ko_names[!(ko_names %in% row.names(expression_profile))],"can't be found in the expression profile, the names of this knockout or knockdown maybe not official gene name, please check and use official gene name instead. Then run this function again. If it is already official gene name, then just go on!"))
  }
  for(i in 1:length(ko_names)){
    if(ko_names[i] %in% row.names(expression_profile)){
      zero_ratio<-length(expression_profile[ko_names[i],][which(expression_profile[ko_names[i],]==0)])/ncol(expression_profile)
      Zero_ra[i]<-zero_ratio
      if(zero_ratio==1){
        print(paste(ko_names[i],"doesn't express and will be filtered.",sep=" "))
        perturb_information_ko<-perturb_information_delete_name(perturb_information_ko,ko_names[i])
      }
    }else{
      Zero_ra[i]<-NA
      #print(paste(ko_names[i],"is missing and will be filtered.",sep=" "))
      #perturb_information_ko<-perturb_information_delete_name(perturb_information_ko,ko_names[i])
    }
  }
  names(Zero_ra)<-ko_names
  expression_profile_ko<-expression_profile[,names(perturb_information_ko)]
  perturb_information_ctrl<-perturb_information[perturb_information=="Non-Targeting"]
  expression_profile_ctrl<-expression_profile[,names(perturb_information_ctrl)]
  
  expression_profile_ko<-expression_profile[,names(perturb_information_ko)]
  cellNum_eachKo<-table(perturb_information_ko)
  deci_m_choose=c()
  filter_record<-matrix(rep(NA,length(cellNum_eachKo)*3),length(cellNum_eachKo))
  filter_record[,1]<-cellNum_eachKo
  colnames(filter_record)<-c("original_num","valid_num","invalid_ratio")
  row.names(filter_record)<-names(cellNum_eachKo)
  for(i in 1:length(cellNum_eachKo)){
    print(paste("filtering for perturbation:",names(cellNum_eachKo)[i],sep=" "))
    if(cellNum_eachKo[i]>=cell_num_threshold){
      vargene<-c()
      ko_barcode<-names(perturb_information_ko[perturb_information_ko==names(cellNum_eachKo)[i]])
      expression_profile_ko_each<-expression_profile_ko[,ko_barcode]
      expr_c<-cbind(expression_profile_ko_each,expression_profile_ctrl)
      gene_highUMI<-apply(expr_c,1,mean)
      names_highUMI<-names(gene_highUMI[gene_highUMI>umi])
      expr_c<-expr_c[names_highUMI,]
      
      s1=length(ko_barcode)
      vargene<-apply(expr_c,1,get_varGene,s=s1)
      vargene<-vargene[vargene<pvalue]
      vargene_name<-names(vargene)
      if(length(vargene)<vargene_min_num){
        print(paste("The number of variable genes of",names(cellNum_eachKo)[i],"is less than ",vargene_min_num,",this perturbation will be filtered directory."))
        filter_record[i,2]<-0
        filter_record[i,3]<-(filter_record[i,1]-filter_record[i,2])/filter_record[i,1]
        next
      }
      expression_profile_ctrl_var<-t(expression_profile_ctrl[vargene_name,])
      expression_profile_ko_var<-t(expression_profile_ko_each[vargene_name,])
      a=cosin_dis_diffMatrix(expression_profile_ko_var,expression_profile_ko_var,cpu_num = cpu_num)
      a_median<-apply(a,1,median)
      b=cosin_dis_diffMatrix(expression_profile_ko_var,expression_profile_ctrl_var,cpu_num=cpu_num)
      b_median<-apply(b,2,median)
      perturb_like<-a_median-b_median
      perturb_like<-perturb_like[perturb_like>0]
      filter_record[i,2]<-length(perturb_like)
      filter_record[i,3]<-(filter_record[i,1]-filter_record[i,2])/filter_record[i,1]
      deci_m_choose<-c(deci_m_choose,names(perturb_like))
    }
    else{
      print("The number of cells with this perturbation is less than 30, this perturbation will be filtered directory.")
    }
  }
  perturb_information_filter<-c(perturb_information[deci_m_choose],perturb_information_ctrl)
  cellNum_eachKo_filter<-as.matrix(table(perturb_information_filter))
  filter_record2<-na.omit(filter_record)
  ko_save<-row.names(filter_record2[filter_record2[,3]<filtered_rate,])
  cellNum_eachKo_filter<-cellNum_eachKo_filter[c(ko_save,"Non-Targeting"),]
  cellNum_eachKo_filter<-cellNum_eachKo_filter[cellNum_eachKo_filter>=cell_num_threshold]
  perturb_information_filter<-perturb_information_filter[perturb_information_filter %in% names(cellNum_eachKo_filter)]
  expression_profile_filter<-expression_profile[,names(perturb_information_filter)]
  perturb_information_filter_abandon<-perturb_information[setdiff(names(perturb_information),names(perturb_information_filter))]
  if(plot){
    pdf(plot_path)
    forPlot<-sort(filter_record2[,3])
    barplot(forPlot,names.arg=names(forPlot),xlab="Perturbation",ylab="Invalid_rate",las=2,,ylim=c(0,1))
    dev.off()
  }
  return(list("expression_profile"=expression_profile_filter,"perturb_information"=perturb_information_filter,"perturb_information_abandon"=perturb_information_filter_abandon,"filter_record"=filter_record,"zero_rate"=Zero_ra))
}

Get_high_varGenes_modified<-function(expression_profile,perturb_information,x.low.cutoff=0.01,y.cutoff=0,num.bin=30,plot=FALSE,plot_path="./get_high_var_genes.pdf"){
  logVarDivMean=function(x) return(log(var(x)/mean(x)))
  expMean=function(x) return(log(mean(x)+1))
  data_norm<-function(xy,num.bin){
    seperate<-seq(0,max(xy[,1]),length.out=num.bin+1)
    for(i in 2:length(seperate)){
      xy[which(xy[,1]>seperate[i-1] & xy[,1]<=seperate[i]),3]<-(xy[which(xy[,1]>seperate[i-1] & xy[,1]<=seperate[i]),2]-mean(xy[which(xy[,1]>seperate[i-1] & xy[,1]<=seperate[i]),2]))/sd(xy[which(xy[,1]>seperate[i-1] & xy[,1]<=seperate[i]),2])
    }
    return(xy)
  }
  label=c()
  j=1
  for(i in 1:length(perturb_information)){
    if(perturb_information[i]=="Non-Targeting"){
      label[j]="Non-Targeting"
      j=j+1
    }
    else{
      label[j]="ko"
      j=j+1
    }
  }
  nCountsEndo_ko<-expression_profile[,which(label=="ko")]
  nCountsEndo_ctrl<-expression_profile[,which(label=="Non-Targeting")]
  data.x_ko<-apply(expression_profile,1,expMean)
  data.y_ko<-apply(nCountsEndo_ko,1,logVarDivMean)
  data.y_ko[data.y_ko==-Inf]<-NaN
  mid<-data.y_ko
  mid[is.nan(mid)]<-0
  data.y_ko[is.nan(data.y_ko)]<-min(mid)-1
  data.y_ctrl<-apply(nCountsEndo_ctrl,1,logVarDivMean)
  data.y_ctrl[data.y_ctrl==-Inf]<-NaN
  mid<-data.y_ctrl
  mid[is.nan(mid)]<-0
  data.y_ctrl[is.nan(data.y_ctrl)]<-min(mid)-1
  data.y_diff<-abs(data.y_ko-data.y_ctrl)
  
  diff_xy<-cbind(data.x_ko,data.y_diff,0,1)
  row.names(diff_xy)<-row.names(nCountsEndo_ko)
  colnames(diff_xy)<-c("data.x","data.y","data.norm.y","vargene")
  diff_xy<-data_norm(diff_xy,num.bin)
  diff_xy[is.na(diff_xy[,3]),3]<-0
  diff_xy[which(diff_xy[,1]>x.low.cutoff  & diff_xy[,3]>y.cutoff),4]<-2
  y_choose<-diff_xy[which(diff_xy[,1]>x.low.cutoff),]
  y_choose<-y_choose[order(-y_choose[,3]),]
  if(plot){
    pdf(file=plot_path)
    par(mfrow=c(1,2))
    plot(diff_xy[,1],diff_xy[,3],type="p",xlab="Average expression",pch=16,ylab="Dispersion difference",col=diff_xy[,4])
    plot(1:nrow(y_choose),y_choose[,3],type="b",xlab="Gene number",ylab="Dispersion difference",col=y_choose[,4])
    grid(nx=NA,ny=50,lwd=1,lty=2,col="blue")
    dev.off()
  }
  expression_profile_varGene<-expression_profile[which(diff_xy[,4]==2),]
  return(list("expression_profile"=expression_profile_varGene,"perturb_information"=perturb_information))
}

Get_topics_modified<-function(expression_profile,perturb_information,topic_number=c(4:6),seed_num=2018,burnin=0,thin=500,iter=500){
  require(slam)
  require(topicmodels)
  require(ggplot2)
  print("Adjusting data for topic model inputting ...")
  #adjust data for topic model
  label=c()
  j=1
  for(i in 1:length(perturb_information)){
    if(perturb_information[i]=="Non-Targeting"){
      label[j]="Non-Targeting"
      j=j+1
    }
    else{
      label[j]="ko"
      j=j+1
    }
  }
  gene_counts_sum<-apply(expression_profile,1,sum)
  expression_profile<-expression_profile[which(gene_counts_sum>0),]
  nCountsEndo_ctrl<-expression_profile[,which(label=="Non-Targeting")]
  for(i in 1:nrow(expression_profile)){
    if(mean(nCountsEndo_ctrl[i,])>0){
      expression_profile[i,]=(expression_profile[i,]-mean(nCountsEndo_ctrl[i,]))/mean(nCountsEndo_ctrl[i,])
    }
  }
  expression_profile<-round((expression_profile+abs(min(expression_profile)))*10)
  #
  topic_model_list<-list()
  control=list(seed=seed_num,burnin=burnin,thin=thin,iter=iter)
  dtm<-as.simple_triplet_matrix(t(expression_profile))
  i=1
  print("It may take a few hours. Please wait patiently.")
  for(k in topic_number){
    print(paste("now the calculating topic number is",k,sep=" "))
    topic_model=LDA(dtm,k=k,method="Gibbs",control=control)
    topic_model_list[[i]]=topic_model
    i<-i+1
  }
  return(list("models"=topic_model_list,"perturb_information"=perturb_information))
}

Diff_topic_distri_modified<-function(model,perturb_information,plot=FALSE,plot_path="./distribution_of_topic.pdf"){
  require(reshape2)
  require(dplyr)
  require(entropy)
  options(warn = -1)
  pmatrix<-model@gamma
  row.names(pmatrix)<-model@documents
  topicNum<-ncol(pmatrix)
  topicName<-paste('Topic',1:topicNum,sep='')
  colnames(pmatrix)<-paste('Topic',1:topicNum,sep='')
  p.matrix<-data.frame(pmatrix,samples=rownames(pmatrix),knockout=perturb_information)
  p.matrix <- melt(p.matrix,id=c('samples','knockout'))
  ko_name<-as.character(unique(p.matrix$knockout))
  t_D<-as.data.frame(matrix(rep(0,length(ko_name)*ncol(pmatrix)*4),ncol = 4))
  colnames(t_D)<-c("knockout","topic","t_D","pvalue")
  p.step1=p.matrix %>% group_by(knockout,variable) %>% summarise(number=sum(value))
  total_number=sum(p.step1$number)
  p.step2=p.step1 %>% group_by(knockout) %>% summarise(cellNum=sum(number))
  p.step1=merge(p.step1,p.step2,by='knockout')
  p.step3<-(p.step1$number)/(p.step1$cellNum)
  p.step4<-data.frame(p.step1,ratio=p.step3)
  p.step4$ctrlNum<-p.step4[which(p.step4$knockout=="Non-Targeting"),"cellNum"]
  p.step4$ctrl_ratio<-p.step4[which(p.step4$knockout=="Non-Targeting"),"ratio"]
  p.step4$diff_index<-(p.step4$ratio-p.step4$ctrl_ratio)
  cell_num_min<-round(min(p.step4$cellNum)*0.9)
  k=1
  for(i in topicName){
    p.matrix.topic<-p.matrix[p.matrix$variable==i,]
    ctrl_topic<-p.matrix.topic[p.matrix.topic$knockout=="Non-Targeting",4]
    ctrl_topic_z<-(ctrl_topic-mean(ctrl_topic))/sqrt(var(ctrl_topic))
    for(j in ko_name){
      ko_topic<-p.matrix.topic[p.matrix.topic$knockout==j,4]
      ko_topic_z<-(ko_topic-mean(ctrl_topic))/sqrt(var(ctrl_topic))
      t_D[k,1]<-j
      t_D[k,2]<-i
      test_s<-matrix(rep(0,2*1000),1000)
      for(t in 1:1000){
        ko_topic_s<-sample(ko_topic,cell_num_min)
        test<-t.test(ko_topic_s,ctrl_topic)
        test_s[t,1]<-test$statistic
      }
      t_D[k,3]<-sort(test_s[,1])[round(length(test_s[,1])/2)]
      k<-k+1
    }
  }
  t_D<-t_D[order(t_D$knockout),]
  p.step4$t_D<-t_D$t_D
  p.step4$t_D_ctrl<-p.step4[p.step4$knockout=="Non-Targeting","t_D"]
  p.step4$t_D_diff<-p.step4$t_D-p.step4$t_D_ctrl
  t_D_sum<-p.step4 %>% group_by(knockout) %>% summarise(t_D_sum=sum(abs(t_D_diff)))
  t_D_sum<-as.data.frame(t_D_sum)
  p.step5=merge(p.step4,t_D_sum,by='knockout')
  t_D_sum2<-p.step4 %>% group_by(knockout) %>% summarise(t_D_sum2=sum(abs(t_D)))
  t_D_sum2<-as.data.frame(t_D_sum2)
  p.step5=merge(p.step5,t_D_sum2,by='knockout')
  p.step5$perturb_percent<-abs(p.step5$t_D)/p.step5$t_D_sum2
  if(plot){
    topic_perturb<-function(distri_diff,plot_path){
      require(gplots)
      require(reshap2)
      perturb_topic<-distri_diff[,c("knockout","variable","t_D_diff")]
      perturb_topic_matrix<-dcast(perturb_topic,knockout ~ variable)
      row.names(perturb_topic_matrix)<-perturb_topic_matrix$knockout
      perturb_topic_matrix<-as.matrix(perturb_topic_matrix[,-1])
      pdf(plot_path)
      heatmap.2(perturb_topic_matrix,col=bluered,dendrogram="both",adjCol=c(NA,1),cexCol = 0.5,cexRow = 0.5,srtCol = 45, srtRow=-45,key=TRUE,trace="none", breaks=seq.int(from = min(perturb_topic_matrix), to = max(perturb_topic_matrix), length.out = 100),hclustfun=hclust)
      dev.off()
    }
    topic_perturb(p.step5,plot_path = plot_path)
  }
  return(p.step5)
}

Rank_overall_modified<-function(distri_diff,offTarget_hash=hash(),output=FALSE,file_path="./rank_overall.txt"){
  require(dplyr)
  require(entropy)
  require(hash)
  KO_offTarget_hash<-hash()
  rank_overall<-distri_diff[,c("knockout","t_D_sum")]
  rank_overall<-unique(rank_overall)
  rank_overall<-rank_overall[order(rank_overall$t_D_sum,decreasing=T),]
  rank_overall$ranking=1:nrow(rank_overall)
  row.names(rank_overall)=1:nrow(rank_overall)
  rank_overall$off_target<-"none"
  for(i in 1:nrow(rank_overall)){
    ko_gene_arr<-unlist(split(rank_overall$knockout,","))
    off_target_arr<-c()
    k=1
    for(j in ko_gene_arr){
      if(has.key(j,offTarget_hash)){
        off_target_arr[k]=offTarget_hash[[j]]
        k=k+1
      }
    }
    off_target=off_target_arr[1]
    if(length(off_target_arr)>1){
      for(k in 2:length(off_target_arr)){
        off_target<-paste(off_target,off_target_arr[k],sep=",")
      }
    }
    if(!is.null(off_target)){
      KO_offTarget_hash[rank_overall$off_target[i]]=off_target
      rank_overall$off_target[i]=off_target
    }
  }
  rankOverall_result<-rank_overall[,c("knockout","ranking","t_D_sum","off_target")]
  rankOverall_result<-rankOverall_result[rankOverall_result$knockout!="Non-Targeting",]
  colnames(rankOverall_result)<-c("perturbation","ranking","Score","off_target")
  if(output){
    write.table(rankOverall_result,file_path,col.names=T,row.names=F,quote=F,sep="\t")
  }
  return(rankOverall_result)
}

Rank_specific_modified<-function(distri_diff,output=FALSE,file_path="./rank_specific.txt"){
  distri_diff<-distri_diff[order(distri_diff$variable,-abs(distri_diff$t_D)),c("variable","knockout","t_D","perturb_percent")]
  t_D_max<-distri_diff %>% group_by(variable) %>% summarise(t_D_max=max(abs(t_D)))
  t_D_min<-distri_diff %>% group_by(variable) %>% summarise(t_D_min=min(abs(t_D)))
  perturb_percent_max<-distri_diff %>% group_by(variable) %>% summarise(perturb_percent_max=max(perturb_percent))
  perturb_percent_min<-distri_diff %>% group_by(variable) %>% summarise(perturb_percent_min=min(perturb_percent))
  rank_topic_specific<-merge(distri_diff,merge(t_D_min,merge(t_D_max,merge(perturb_percent_max,perturb_percent_min))))
  rank_topic_specific$t_D_standard<-(abs(rank_topic_specific$t_D)-rank_topic_specific$t_D_min)/(rank_topic_specific$t_D_max-rank_topic_specific$t_D_min)
  rank_topic_specific$perturb_percent_standard<-(rank_topic_specific$perturb_percent-rank_topic_specific$perturb_percent_min)/(rank_topic_specific$perturb_percent_max-rank_topic_specific$perturb_percent_min)
  rank_topic_specific$recommand_score<-rank_topic_specific$t_D_standard+rank_topic_specific$perturb_percent_standard
  rank_topic_specific<-rank_topic_specific[order(rank_topic_specific$knockout),]
  rank_topic_specific$recommand_score_ctrl<-rank_topic_specific[rank_topic_specific$knockout=="Non-Targeting","recommand_score"]
  rank_topic_specific<-rank_topic_specific[rank_topic_specific$recommand_score>rank_topic_specific$recommand_score_ctrl,]
  rank_topic_specific<-rank_topic_specific[order(rank_topic_specific$variable,-rank_topic_specific$recommand_score),]
  topic_times<-table(as.character(rank_topic_specific$variable))
  ranking<-c()
  for(i in topic_times){
    ranking<-c(ranking,1:i)
  }
  rank_topic_specific$ranking<-ranking
  rank_topic_specific<-rank_topic_specific[,c("variable","knockout","ranking")]
  row.names(rank_topic_specific)<-1:nrow(rank_topic_specific)
  colnames(rank_topic_specific)<-c("topic","perturbation","ranking")
  if(output){
    write.table(rank_topic_specific,file_path,col.names=T,row.names=F,quote=F,sep="\t")
  }
  return(rank_topic_specific)
}

Correlation_perturbation<-function(distri_diff,cutoff=0.9,gene="all",plot=FALSE,plot_path="./correlation_network.pdf",output=FALSE,file_path="./correlation_perturbation.txt"){
  require(gplots)
  require(reshape2)
  distri_diff<-distri_diff[which(distri_diff$knockout!="Non-Targeting"),]
  correlation<-function(matrix,method="pearson"){
    sample_num<-nrow(matrix)
    cor_matrix<-matrix(rep(1,sample_num^2),sample_num)
    colnames(cor_matrix)<-row.names(matrix)
    row.names(cor_matrix)<-row.names(matrix)
    for(i in 1:(sample_num-1)){
      for(j in (i+1):sample_num){
        cor_matrix[i,j]<-cor(matrix[i,],matrix[j,],method = method)
        cor_matrix[j,i]<-cor_matrix[i,j]
      }
    }
    return(cor_matrix)
  }
  perturb_topic<-distri_diff[,c("knockout","variable","t_D_diff")]
  perturb_topic_matrix<-dcast(perturb_topic,knockout ~ variable)
  row.names(perturb_topic_matrix)<-perturb_topic_matrix$knockout
  perturb_topic_matrix<-as.matrix(perturb_topic_matrix[,-1])
  perturb_topic_cor<-correlation(perturb_topic_matrix)
  perturb_cor<-melt(perturb_topic_cor)
  colnames(perturb_cor)<-c("Perturbation_1","Perturbation_2","Correlation")
  perturb_cor<-perturb_cor[perturb_cor$Perturbation_1!=perturb_cor$Perturbation_2,]
  perturb_cor<-perturb_cor[order(perturb_cor$Correlation),]
  perturb_cor<-perturb_cor[seq(1,nrow(perturb_cor),by=2),]
  perturb_cor<-perturb_cor[order(-abs(perturb_cor$Correlation)),]
  perturb_cor<-perturb_cor[order(perturb_cor$Perturbation_1),]
  if(gene!="all"){
    for(i in 1:length(gene)){
      if(!(gene[i] %in% distri_diff$knockout)){
        print(paste("Warning! Can't find",gene[i],",please check it again!"))
        stop()
      }
    }
    perturb_cor<-perturb_cor[(perturb_cor$Perturbation_1 %in% gene) | (perturb_cor$Perturbation_2 %in% gene),]
  }
  if(output){
    write.table(perturb_cor,file_path,col.names=T,row.names=F,quote=F,sep="\t")
  }
  if(plot){
    require(igraph)
    pdf(plot_path)
    #cutoff<-quantile(abs(perturb_cor$Correlation),quanti)
    t<-perturb_cor[abs(perturb_cor$Correlation)>cutoff,]
    t$direction<-sign(t$Correlation)
    color_e<-c()
    for(i in 1:nrow(t)){
      if(t$direction[i]==-1){
        color_e[i]="blue"
      }else{
        color_e[i]="red"
      }
    }
    opar <- par(no.readonly = TRUE)
    par(mar = c(0,0,0,0))
    g <- graph.data.frame(t, directed = FALSE) 
    E(g)$color=color_e
    plot(g, layout = layout.fruchterman.reingold)
    par(opar)
    dev.off()
  }
  return(perturb_cor)
}
