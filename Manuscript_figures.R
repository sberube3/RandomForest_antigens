##################
####Libraries#####
##################

library(tidyverse)
library(randomForest)
library(ggplot2)
library (cowplot)
library(pROC)
library(multiROC)
#########################
#Read in and Clean Data#
########################

#######################
####ICEMR Site Data###
######################

site_data<- read.csv("ICEMR_SurveillanceData_20211215.csv")

site_data<- site_data%>%
  mutate(Location= )

summary_rdt_counts<- site_data%>%
  group_by(hcid_ymd,wk_start_date)%>%
  summarize(avg_num_pos= mean(c(rdtp_u5+rdtp_a5),na.rm=T))


setwd("~/Desktop/Hopkins/Hopkins/Research_Protein_Arrays/Tamaki_arrays")

bar_ranks<- read.csv("Malaria_array_ranks_bar.csv")
meta_data<- read.csv("Raw_covariates.csv")

meta_data<- as.data.frame(t(meta_data))
colnames(meta_data)<- meta_data[1,]
meta_data<- meta_data[-1,]


bar_ranks<- bar_ranks[c(1:488),]



bar_ranks<- bar_ranks%>%
  select(V1:V1038)%>%
  mutate(Location=meta_data$Location)%>%
  filter(Location=="Macha "|Location=="Choma District " |Location=="Nchelenge"|Location=="Honde Valley")%>%
  mutate(Location=ifelse(Location=="Macha "|Location=="Choma District ", "Macha",Location))

colnames(meta_data)[37]<- "Unknown"
meta_data<- meta_data%>%
  filter(Location=="Macha "|Location=="Choma District " |Location=="Nchelenge"|Location=="Honde Valley")%>%
  mutate(Location=ifelse(Location=="Macha "|Location=="Choma District ", "Macha",Location))




index_neg_ctrl_empty_M<- c(574,863,866,1152,1155)
index_neg_ctrl_blank_M<- c(289,578,867,1156)
index_neg_ctrl_blank1_M<- c(285,349,695,1144,1147)
index_neg_ctrl_noDNA_M<- c(17,18,19,20,21,22,306,307,308,309,310,311,
                           595,596,597,598,599,600,884,885,886,887,888,889)
index_neg_ctrl_ttbs_M<- c(13,14,15,16,29,30,31,32,302,303,304,305,
                          318,319,320,321,591,592,593,594,607,608,609,610,
                          880,881,882,883,896,897,898,899)

nc_index_M<- c(index_neg_ctrl_empty_M,index_neg_ctrl_blank_M,
               index_neg_ctrl_blank1_M,index_neg_ctrl_noDNA_M,index_neg_ctrl_ttbs_M)


index_Igg003<-c(11,300,589,878)
index_Igg03<- c(9,298,587,876)
index_Igg3<- c(7,296,585,874)
index_MIgg003<-c(5,294,583,872)
index_MIgg03<- c(3,292,581,870)
index_MIgg3<- c(1,290,579,868)
index_Igg001<- c(12,301,590,879)
index_Igg01<- c(10,299,588,877)
index_Igg1<- c(8,297,586,875)
index_MIgg001<- c(6,295,584,873)
index_MIgg01<- c(4,293,582,871)
index_MIgg1<- c(2,291,580,869)

pc_index_M<- c(index_Igg003,index_Igg03,index_Igg3,index_MIgg003,
               index_MIgg03,index_MIgg3,index_Igg001,index_Igg01,
               index_Igg1,index_MIgg001,index_MIgg01,index_MIgg1)


array_clean<-function(x,data){
  Array_x_whole<-data
  
  Array_x<- Array_x_whole[60:9308,]
  
  colnames(Array_x)<-c("Index","Array_Row","Array_Column","Spot_Row","Spot_Column",
                       "Name","ID","X","Y","Diameter","F_Pixels","B_Pixels","Footprint",
                       "Flags","Ch1_Median","Ch1_Mean","Ch1_SD","Ch1_B_Median","Ch1_B_Mean",
                       "Ch1_B_SD","Ch1_%_>_B_+_1_SD","Ch1_%_>_B+2_SD", 
                       "Ch1_F_%_Sat.", "Ch1_Median_-_B", "Ch1_Mean_-B",
                       "Ch1_SignalNoiseRatio")
  
  a<-(8*x)-7
  b<-(8*x)-6
  c<-(8*x)-5
  d<-(8*x)-4
  e<-(8*x)-3
  f<-(8*x)-2
  g<-(8*x)-1
  h<-(8*x)
  
  
  Array_a<- Array_x%>%
    filter(Array_Row==1)
  
  Array_b<- Array_x%>%
    filter(Array_Row==2)
  
  Array_c<- Array_x%>%
    filter(Array_Row==3)
  
  Array_d<- Array_x%>%
    filter(Array_Row==4)
  
  Array_e<- Array_x%>%
    filter(Array_Row==5)
  
  Array_f<- Array_x%>%
    filter(Array_Row==6)
  
  Array_g<- Array_x%>%
    filter(Array_Row==7)
  
  Array_h<- Array_x%>%
    filter(Array_Row==8)
  
  
  return(list(Array_a,Array_b,Array_c,Array_d,Array_e,Array_f,Array_g,Array_h))
  
}



data_1<-read.csv("./tamaki_csv/1_280_100.csv",header=TRUE,stringsAsFactors = F)



Array_1<-array_clean(x=1,data=data_1)

Sample_1<- Array_1[[1]]

array_proteins<- as.data.frame(cbind(Sample_1$Name[-c(pc_index_M,nc_index_M)],Sample_1$ID[-c(pc_index_M,nc_index_M)],rep(1:1038)))

colnames(array_proteins)<-c("Name", "ID", "Index")

antigen_id_table<- read.csv("Antigen_Profile.csv")


############################################################################
#Run base model, do nested CV, select top antigens, run top antigen models#
###########################################################################
bar_ranks$Location<- as.factor(bar_ranks$Location)
set.seed(33)
m1<- randomForest(Location~.,data=bar_ranks,importance=T,classwt=c(1/3,1/3,1/3))


m_select<- rfcv(bar_ranks[,-1039],as.factor(bar_ranks$Location), cv.fold=10)


importance_df<- as.data.frame(cbind(importance(m1,type=2),array_proteins))

importance_df<- importance_df%>%
  mutate(Protein=rownames(importance_df))%>%
  mutate(Index=paste0("V",Index))


importance_of_vars<- importance_df[order(importance_df$MeanDecreaseGini,decreasing = T),]

importance_of_vars_ID<- full_join(importance_of_vars,antigen_id_table,by="ID")

importance_of_vars_30<- importance_of_vars_ID[1:30,c(1,2,12)]
importance_of_vars_65<- importance_of_vars_ID[1:65,c(1,2,12)]
importance_of_vars_130<- importance_of_vars_ID[1:130,c(1,2,12)]
importance_of_vars_260<- importance_of_vars_ID[1:260,c(1,2,12)]


bar_ranks_top_vars_20<- bar_ranks[,c(which(colnames(bar_ranks)%in%importance_of_vars$Index[1:20]),1039)]
bar_ranks_top_vars_30<- bar_ranks[,c(which(colnames(bar_ranks)%in%importance_of_vars$Index[1:30]),1039)]
bar_ranks_top_vars_65<- bar_ranks[,c(which(colnames(bar_ranks)%in%importance_of_vars$Index[1:65]),1039)]
bar_ranks_top_vars_130<- bar_ranks[,c(which(colnames(bar_ranks)%in%importance_of_vars$Index[1:130]),1039)]
bar_ranks_top_vars_260<- bar_ranks[,c(which(colnames(bar_ranks)%in%importance_of_vars$Index[1:260]),1039)]

set.seed(33)
bar_ranks_top_vars_20$Location<- as.factor(bar_ranks_top_vars_20$Location)
m20_top<- randomForest(Location~.,data=bar_ranks_top_vars_20,importance=T, classwt= c(1/3,1/3,1/3))

bar_ranks_top_vars_30$Location<- as.factor(bar_ranks_top_vars_30$Location)
m30_top<- randomForest(Location~.,data=bar_ranks_top_vars_30,importance=T, classwt= c(1/3,1/3,1/3))

bar_ranks_top_vars_65$Location<- as.factor(bar_ranks_top_vars_65$Location)
m65_top<- randomForest(Location~.,data=bar_ranks_top_vars_65,importance=T, classwt= c(1/3,1/3,1/3))

bar_ranks_top_vars_130$Location<- as.factor(bar_ranks_top_vars_130$Location)
m130_top<- randomForest(Location~.,data=bar_ranks_top_vars_130,importance=T, classwt= c(1/3,1/3,1/3))

bar_ranks_top_vars_260$Location<- as.factor(bar_ranks_top_vars_260$Location)
m260_top<- randomForest(Location~.,data=bar_ranks_top_vars_260,importance=T, classwt= c(1/3,1/3,1/3))

###########################################
#Repeat with a random subset of antigens #
##########################################
set.seed(33)
bar_ranks_random_vars_20<- bar_ranks[,c(sample(1:1038,size=20),1039)]
bar_ranks_random_vars_30<- bar_ranks[,c(sample(1:1038,size=30),1039)]
bar_ranks_random_vars_65<- bar_ranks[,c(sample(1:1038,size=65),1039)]
bar_ranks_random_vars_130<- bar_ranks[,c(sample(1:1038,size=130),1039)]
bar_ranks_random_vars_260<- bar_ranks[,c(sample(1:1038,size=260),1039)]

set.seed(33)
bar_ranks_random_vars_20$Location<- as.factor(bar_ranks_random_vars_20$Location)
m20_random<- randomForest(Location~.,data=bar_ranks_random_vars_20,importance=T, classwt= c(1/3,1/3,1/3))

bar_ranks_random_vars_30$Location<- as.factor(bar_ranks_random_vars_30$Location)
m30_random<- randomForest(Location~.,data=bar_ranks_random_vars_30,importance=T, classwt= c(1/3,1/3,1/3))

bar_ranks_random_vars_65$Location<- as.factor(bar_ranks_random_vars_65$Location)
m65_random<- randomForest(Location~.,data=bar_ranks_random_vars_65,importance=T, classwt= c(1/3,1/3,1/3))

bar_ranks_random_vars_130$Location<- as.factor(bar_ranks_random_vars_130$Location)
m130_random<- randomForest(Location~.,data=bar_ranks_random_vars_130,importance=T, classwt= c(1/3,1/3,1/3))

bar_ranks_random_vars_260$Location<- as.factor(bar_ranks_random_vars_260$Location)
m260_random<- randomForest(Location~.,data=bar_ranks_random_vars_260,importance=T, classwt= c(1/3,1/3,1/3))



###########
#Figure 1#
##########

global_oob_figure1<- c(6.89,6.47,6.47,7.10,5.85,21.92,16.08,11.06,9.81,10.44)
mutasa_oob_figure1<-c(3.68,2.63,3.16,3.68,2.63,17.89,7.37,5.26,5.26,4.21)
macha_oob_figure1<- c(9.58,9.58,9.58,10.78,8.38,24.55,20.36,17.37,15.57,15.57)
nchelenge_oob_figure1<-c(8.20,8.20,7.38,7.38,7.38,24.59,23.77,11.48,9.02,13.11)

figure_1_df<- as.data.frame(cbind(c(global_oob_figure1,mutasa_oob_figure1,macha_oob_figure1,nchelenge_oob_figure1),
                                  c(rep(c("Top", "Random"),each=5,times=4)),
                                  c(rep(c("Global","Mutasa", "Choma", "Nchelenge"),each=10)),
                                  c(rep(c(20,30,65,130,260),times=8))))

colnames(figure_1_df)<- c("OOB Error", "Subset","Region", "Antigens")
figure_1_df$`OOB Error`<- as.numeric(figure_1_df$`OOB Error`)
figure_1_df$Antigens<- as.numeric(figure_1_df$Antigens)


figure_1<- ggplot(figure_1_df, aes(x=Antigens,y=`OOB Error`, group=Subset))+
  geom_line(aes(color=Subset))+
  geom_point(aes(shape=Subset,color=Subset))+
  theme_classic()+scale_color_brewer(palette = "Set1")+labs(x="Number of Antigens", y="Out of Bag Error Rate (%)")+
  facet_wrap(~Region)+theme(legend.key.size = unit(1,"cm"),legend.text=element_text(size=10))


###########################################
#Split test/train and re-run random forests
###########################################
set.seed(33)
training_set<- c(sample(which(bar_ranks$Location=="Macha"),size=111),
                 sample(which(bar_ranks$Location=="Nchelenge"),size=81),
                 sample(which(bar_ranks$Location=="Honde Valley"),size=126))

m20_top_train<- randomForest(Location~.,data=bar_ranks_top_vars_20[training_set,],importance=T, classwt= c(1/3,1/3,1/3))

m30_top_train<- randomForest(Location~.,data=bar_ranks_top_vars_30[training_set,],importance=T, classwt= c(1/3,1/3,1/3))

m65_top_train<- randomForest(Location~.,data=bar_ranks_top_vars_65[training_set,],importance=T, classwt= c(1/3,1/3,1/3))

m130_top_train<- randomForest(Location~.,data=bar_ranks_top_vars_130[training_set,],importance=T, classwt= c(1/3,1/3,1/3))

m260_top_train<- randomForest(Location~.,data=bar_ranks_top_vars_260[training_set,],importance=T, classwt= c(1/3,1/3,1/3))


pred_m20_train<- predict(m20_top_train,bar_ranks_top_vars_20[-training_set,],type="class")
conf_m20<- table(pred_m20_train,bar_ranks_top_vars_20[-training_set,21])
honde_error<- sum(conf_m20[1,2],conf_m20[1,3])/(sum(conf_m20[1,]))
macha_error<- sum(conf_m20[2,1],conf_m20[2,3])/(sum(conf_m20[2,]))
nchelenge_error<- sum(conf_m20[3,1],conf_m20[3,2])/(sum(conf_m20[3,]))
global_err<- mean(c(honde_error,macha_error,nchelenge_error))


pred_m30_train<- predict(m30_top_train,bar_ranks_top_vars_30[-training_set,],type="class")
conf_m30<- table(pred_m30_train,bar_ranks_top_vars_30[-training_set,31])
honde_error<- sum(conf_m30[1,2],conf_m30[1,3])/(sum(conf_m30[1,]))
macha_error<- sum(conf_m30[2,1],conf_m30[2,3])/(sum(conf_m30[2,]))
nchelenge_error<- sum(conf_m30[3,1],conf_m30[3,2])/(sum(conf_m30[3,]))
global_err<- mean(c(honde_error,macha_error,nchelenge_error))

pred_m65_train<- predict(m65_top_train,bar_ranks_top_vars_65[-training_set,],type="class")
conf_m65<- table(pred_m65_train,bar_ranks_top_vars_65[-training_set,66])
honde_error<- sum(conf_m65[1,2],conf_m65[1,3])/(sum(conf_m65[1,]))
macha_error<- sum(conf_m65[2,1],conf_m65[2,3])/(sum(conf_m65[2,]))
nchelenge_error<- sum(conf_m65[3,1],conf_m65[3,2])/(sum(conf_m65[3,]))
global_err<- mean(c(honde_error,macha_error,nchelenge_error))

pred_m130_train<- predict(m130_top_train,bar_ranks_top_vars_130[-training_set,],type="class")
conf_m130<- table(pred_m130_train,bar_ranks_top_vars_130[-training_set,131])
honde_error<- sum(conf_m130[1,2],conf_m130[1,3])/(sum(conf_m130[1,]))
macha_error<- sum(conf_m130[2,1],conf_m130[2,3])/(sum(conf_m130[2,]))
nchelenge_error<- sum(conf_m130[3,1],conf_m130[3,2])/(sum(conf_m130[3,]))
global_err<- mean(c(honde_error,macha_error,nchelenge_error))


pred_m260_train<- predict(m260_top_train,bar_ranks_top_vars_260[-training_set,],type="class")
conf_m260<- table(pred_m260_train,bar_ranks_top_vars_260[-training_set,261])
honde_error<- sum(conf_m260[1,2],conf_m260[1,3])/(sum(conf_m260[1,]))
macha_error<- sum(conf_m260[2,1],conf_m260[2,3])/(sum(conf_m260[2,]))
nchelenge_error<- sum(conf_m260[3,1],conf_m260[3,2])/(sum(conf_m260[3,]))
global_err<- mean(c(honde_error,macha_error,nchelenge_error))
##################################################################
m20_random_train<- randomForest(Location~.,data=bar_ranks_random_vars_20[training_set,],importance=T, classwt= c(1/3,1/3,1/3))

m30_random_train<- randomForest(Location~.,data=bar_ranks_random_vars_30[training_set,],importance=T, classwt= c(1/3,1/3,1/3))

m65_random_train<- randomForest(Location~.,data=bar_ranks_random_vars_65[training_set,],importance=T, classwt= c(1/3,1/3,1/3))

m130_random_train<- randomForest(Location~.,data=bar_ranks_random_vars_130[training_set,],importance=T, classwt= c(1/3,1/3,1/3))

m260_random_train<- randomForest(Location~.,data=bar_ranks_random_vars_260[training_set,],importance=T, classwt= c(1/3,1/3,1/3))


pred_m20_train<- predict(m20_random_train,bar_ranks_random_vars_20[-training_set,],type="class")
conf_m20<- table(pred_m20_train,bar_ranks_random_vars_20[-training_set,21])
honde_error<- sum(conf_m20[1,2],conf_m20[1,3])/(sum(conf_m20[1,]))
macha_error<- sum(conf_m20[2,1],conf_m20[2,3])/(sum(conf_m20[2,]))
nchelenge_error<- sum(conf_m20[3,1],conf_m20[3,2])/(sum(conf_m20[3,]))
global_err<- mean(c(honde_error,macha_error,nchelenge_error))


pred_m30_train<- predict(m30_random_train,bar_ranks_random_vars_30[-training_set,],type="class")
conf_m30<- table(pred_m30_train,bar_ranks_random_vars_30[-training_set,31])
honde_error<- sum(conf_m30[1,2],conf_m30[1,3])/(sum(conf_m30[1,]))
macha_error<- sum(conf_m30[2,1],conf_m30[2,3])/(sum(conf_m30[2,]))
nchelenge_error<- sum(conf_m30[3,1],conf_m30[3,2])/(sum(conf_m30[3,]))
global_err<- mean(c(honde_error,macha_error,nchelenge_error))

pred_m65_train<- predict(m65_random_train,bar_ranks_random_vars_65[-training_set,],type="class")
conf_m65<- table(pred_m65_train,bar_ranks_random_vars_65[-training_set,66])
honde_error<- sum(conf_m65[1,2],conf_m65[1,3])/(sum(conf_m65[1,]))
macha_error<- sum(conf_m65[2,1],conf_m65[2,3])/(sum(conf_m65[2,]))
nchelenge_error<- sum(conf_m65[3,1],conf_m65[3,2])/(sum(conf_m65[3,]))
global_err<- mean(c(honde_error,macha_error,nchelenge_error))

pred_m130_train<- predict(m130_random_train,bar_ranks_random_vars_130[-training_set,],type="class")
conf_m130<- table(pred_m130_train,bar_ranks_random_vars_130[-training_set,131])
honde_error<- sum(conf_m130[1,2],conf_m130[1,3])/(sum(conf_m130[1,]))
macha_error<- sum(conf_m130[2,1],conf_m130[2,3])/(sum(conf_m130[2,]))
nchelenge_error<- sum(conf_m130[3,1],conf_m130[3,2])/(sum(conf_m130[3,]))
global_err<- mean(c(honde_error,macha_error,nchelenge_error))


pred_m260_train<- predict(m260_random_train,bar_ranks_random_vars_260[-training_set,],type="class")
conf_m260<- table(pred_m260_train,bar_ranks_random_vars_260[-training_set,261])
honde_error<- sum(conf_m260[1,2],conf_m260[1,3])/(sum(conf_m260[1,]))
macha_error<- sum(conf_m260[2,1],conf_m260[2,3])/(sum(conf_m260[2,]))
nchelenge_error<- sum(conf_m260[3,1],conf_m260[3,2])/(sum(conf_m260[3,]))
global_err<- mean(c(honde_error,macha_error,nchelenge_error))


############
#Figure S1#
###########

global_oob_figureS1<- c(8.37,6.86,6.29,7.61,8.24,21.92,17.73,9.39,10.71,9.87)
mutasa_oob_figureS1<-c(7.58,8.70,8.96,10.14,10.29,24.64,15.71,12.68,12.86,10.29)
macha_oob_figureS1<- c(12.28,9.26,7.41,7.69,9.43,25.00,13.04,8.16,11.76,9.80)
nchelenge_oob_figureS1<-c(5.26,2.63,2.50,5.00,5.00,17.50,24.44,7.32,7.50,9.52)

figure_S1_df<- as.data.frame(cbind(c(global_oob_figureS1,mutasa_oob_figureS1,macha_oob_figureS1,nchelenge_oob_figureS1),
                                  c(rep(c("Top", "Random"),each=5,times=4)),
                                  c(rep(c("Global","Mutasa", "Choma", "Nchelenge"),each=10)),
                                  c(rep(c(20,30,65,130,260),times=8))))

colnames(figure_S1_df)<- c("OOB Error", "Subset","Region", "Antigens")
figure_S1_df$`OOB Error`<- as.numeric(figure_S1_df$`OOB Error`)
figure_S1_df$Antigens<- as.numeric(figure_S1_df$Antigens)


figure_S1<- ggplot(figure_S1_df, aes(x=Antigens,y=`OOB Error`, group=Subset))+
  geom_line(aes(color=Subset))+
  geom_point(aes(shape=Subset,color=Subset))+
  theme_classic()+scale_color_brewer(palette = "Set1")+labs(x="Number of Antigens", y="Test Error (%)")+
  facet_wrap(~Region)+theme(legend.key.size = unit(1,"cm"),legend.text=element_text(size=10))


###############################################
######AUC for 30, 65, 130, 260 antigens########
###############################################


pred_m20_all<- predict(m20_top,type="prob")


pred_matrix_m20_all<- as.data.frame(cbind(as.character(bar_ranks_top_vars_20[,21]),pred_m20_all))


pred_matrix_m20_all<- pred_matrix_m20_all%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

roc_m20_all<-multi_roc(pred_matrix_m20_all)
roc_m20_all_proc<- multiclass.roc(bar_ranks_top_vars_20[,21]~pred_m20_all)
ci_roc_m20_all<- roc_auc_with_ci(pred_matrix_m20_all,type="basic",R=1000)

euclid_m20_macha<-sqrt((0-(1-roc_m20_all$Specificity$RF$Macha))^2 + (1-roc_m20_all$Sensitivity$RF$Macha)^2)
macha_max_sens<- roc_m20_all$Sensitivity$RF$Macha[which(euclid_m20_macha==min(euclid_m20_macha))]
macha_max_spec<- roc_m20_all$Specificity$RF$Macha[which(euclid_m20_macha==min(euclid_m20_macha))]

euclid_m20_nchelenge<-sqrt((0-(1-roc_m20_all$Specificity$RF$Nchelenge))^2 + (1-roc_m20_all$Sensitivity$RF$Nchelenge)^2)
nchelenge_max_sens<- roc_m20_all$Sensitivity$RF$Nchelenge[which(euclid_m20_nchelenge==min(euclid_m20_nchelenge))]
nchelenge_max_spec<- roc_m20_all$Specificity$RF$Nchelenge[which(euclid_m20_nchelenge==min(euclid_m20_nchelenge))]

euclid_m20_mutasa<-sqrt((0-(1-roc_m20_all$Specificity$RF$Mutasa))^2 + (1-roc_m20_all$Sensitivity$RF$Mutasa)^2)
mutasa_max_sens<- roc_m20_all$Sensitivity$RF$Mutasa[which(euclid_m20_mutasa==min(euclid_m20_mutasa))]
mutasa_max_spec<- roc_m20_all$Specificity$RF$Mutasa[which(euclid_m20_mutasa==min(euclid_m20_mutasa))]

euclid_m20_macro<-sqrt((0-(1-roc_m20_all$Specificity$RF$macro))^2 + (1-roc_m20_all$Sensitivity$RF$macro)^2)
macro_max_sens<- roc_m20_all$Sensitivity$RF$macro[which(euclid_m20_macro==min(euclid_m20_macro))]
macro_max_spec<- roc_m20_all$Specificity$RF$macro[which(euclid_m20_macro==min(euclid_m20_macro))]


#######################################

pred_m30_all<- predict(m30_top,type="prob")


pred_matrix_m30_all<- as.data.frame(cbind(as.character(bar_ranks_top_vars_30[,31]),pred_m30_all))


pred_matrix_m30_all<- pred_matrix_m30_all%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

roc_m30_all<-multi_roc(pred_matrix_m30_all)
roc_m30_all_proc<- multiclass.roc(bar_ranks_top_vars_30[,31]~pred_m30_all)
ci_roc_m30_all<- roc_auc_with_ci(pred_matrix_m30_all,type="basic",R=1000)

###################################################

pred_m65_all<- predict(m65_top,type="prob")


pred_matrix_m65_all<- as.data.frame(cbind(as.character(bar_ranks_top_vars_65[,66]),pred_m65_all))


pred_matrix_m65_all<- pred_matrix_m65_all%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

roc_m65_all<-multi_roc(pred_matrix_m65_all)
roc_m65_all_proc<- multiclass.roc(bar_ranks_top_vars_65[,66]~pred_m65_all)
ci_roc_m65_all<- roc_auc_with_ci(pred_matrix_m65_all,type="basic",R=1000)

###################################################

pred_m130_all<- predict(m130_top,type="prob")


pred_matrix_m130_all<- as.data.frame(cbind(as.character(bar_ranks_top_vars_130[,131]),pred_m130_all))


pred_matrix_m130_all<- pred_matrix_m130_all%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

roc_m130_all<-multi_roc(pred_matrix_m130_all)
roc_m130_all_proc<- multiclass.roc(bar_ranks_top_vars_130[,131]~pred_m130_all)
ci_roc_m130_all<- roc_auc_with_ci(pred_matrix_m130_all,type="basic",R=1000)


#####################################################
pred_m260_all<- predict(m260_top,type="prob")


pred_matrix_m260_all<- as.data.frame(cbind(as.character(bar_ranks_top_vars_260[,261]),pred_m260_all))


pred_matrix_m260_all<- pred_matrix_m260_all%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

roc_m260_all<-multi_roc(pred_matrix_m260_all)
roc_m260_all_proc<- multiclass.roc(bar_ranks_top_vars_260[,261]~pred_m260_all)
ci_roc_m260_all<- roc_auc_with_ci(pred_matrix_m260_all,type="basic",R=1000)


#####################################################

#############################################
####Repeat AUC for test/train split#########
############################################

pred_m20_train<- predict(m20_top_train,bar_ranks_top_vars_20[-training_set,-21],type="prob")

pred_matrix_m20_all<- as.data.frame(cbind(as.character(bar_ranks_top_vars_20[-training_set,21]),pred_m20_train))


pred_matrix_m20_all<- pred_matrix_m20_all%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

roc_m20_all<-multi_roc(pred_matrix_m20_all)
roc_m20_all_proc<- multiclass.roc(bar_ranks_top_vars_20[,21]~pred_m20_all)
ci_roc_m20_all<- roc_auc_with_ci(pred_matrix_m20_all,type="basic",R=1000)

#######################################

pred_m30_train<- predict(m30_top_train,bar_ranks_top_vars_30[-training_set,-31],type="prob")


pred_matrix_m30_all<- as.data.frame(cbind(as.character(bar_ranks_top_vars_30[-training_set,31]),pred_m30_train))


pred_matrix_m30_all<- pred_matrix_m30_all%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

roc_m30_all<-multi_roc(pred_matrix_m30_all)
roc_m30_all_proc<- multiclass.roc(bar_ranks_top_vars_30[,31]~pred_m30_all)
ci_roc_m30_all<- roc_auc_with_ci(pred_matrix_m30_all,type="basic",R=1000)

###################################################

pred_m65_train<- predict(m65_top_train,bar_ranks_top_vars_65[-training_set,-66],type="prob")


pred_matrix_m65_all<- as.data.frame(cbind(as.character(bar_ranks_top_vars_65[-training_set,66]),pred_m65_train))


pred_matrix_m65_all<- pred_matrix_m65_all%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

roc_m65_all<-multi_roc(pred_matrix_m65_all)
roc_m65_all_proc<- multiclass.roc(bar_ranks_top_vars_65[,66]~pred_m65_all)
ci_roc_m65_all<- roc_auc_with_ci(pred_matrix_m65_all,type="basic",R=1000)

###################################################

pred_m130_train<- predict(m130_top_train,bar_ranks_top_vars_130[-training_set,-131],type="prob")


pred_matrix_m130_all<- as.data.frame(cbind(as.character(bar_ranks_top_vars_130[-training_set,131]),pred_m130_train))


pred_matrix_m130_all<- pred_matrix_m130_all%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

roc_m130_all<-multi_roc(pred_matrix_m130_all)
roc_m130_all_proc<- multiclass.roc(bar_ranks_top_vars_130[,131]~pred_m130_all)
ci_roc_m130_all<- roc_auc_with_ci(pred_matrix_m130_all,type="basic",R=1000)


#####################################################
pred_m260_train<- predict(m260_top_train,bar_ranks_top_vars_260[-training_set,-261],type="prob")


pred_matrix_m260_all<- as.data.frame(cbind(as.character(bar_ranks_top_vars_260[-training_set,261]),pred_m260_train))


pred_matrix_m260_all<- pred_matrix_m260_all%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

roc_m260_all<-multi_roc(pred_matrix_m260_all)
roc_m260_all_proc<- multiclass.roc(bar_ranks_top_vars_260[,261]~pred_m260_all)
ci_roc_m260_all<- roc_auc_with_ci(pred_matrix_m260_all,type="basic",R=1000)


#####################################################

###############################################
########Split data by age/RDT status##########
##############################################

#age 5 cutoff

child_5<- which(as.numeric(meta_data$Age)<=5)
adult_5<- which(as.numeric(meta_data$Age)>5)
child_15<- which(as.numeric(meta_data$Age)<=15)
adult_15<- which(as.numeric(meta_data$Age)>15)
rdt_pos<- which(meta_data$`RDT Result`==1)
rdt_neg<- which(meta_data$`RDT Result`==0)

set.seed(33)

m20_top_child_5<- randomForest(Location~.,data=bar_ranks_top_vars_20[child_5,],importance=T, classwt= c(1/3,1/3,1/3))
m20_top_adult_5<- randomForest(Location~.,data=bar_ranks_top_vars_20[adult_5,],importance=T, classwt= c(1/3,1/3,1/3))
m20_top_child_15<- randomForest(Location~.,data=bar_ranks_top_vars_20[child_15,],importance=T, classwt= c(1/3,1/3,1/3))
m20_top_adult_15<- randomForest(Location~.,data=bar_ranks_top_vars_20[adult_15,],importance=T, classwt= c(1/3,1/3,1/3))
m20_top_rdt_pos<- randomForest(Location~.,data=bar_ranks_top_vars_20[rdt_pos,],importance=T, classwt= c(1/3,1/3,1/3))
m20_top_rdt_neg<- randomForest(Location~.,data=bar_ranks_top_vars_20[rdt_neg,],importance=T, classwt= c(1/3,1/3,1/3))



m30_top_child_5<- randomForest(Location~.,data=bar_ranks_top_vars_30[child_5,],importance=T, classwt= c(1/3,1/3,1/3))
m30_top_adult_5<- randomForest(Location~.,data=bar_ranks_top_vars_30[adult_5,],importance=T, classwt= c(1/3,1/3,1/3))
m30_top_child_15<- randomForest(Location~.,data=bar_ranks_top_vars_30[child_15,],importance=T, classwt= c(1/3,1/3,1/3))
m30_top_adult_15<- randomForest(Location~.,data=bar_ranks_top_vars_30[adult_15,],importance=T, classwt= c(1/3,1/3,1/3))
m30_top_rdt_pos<- randomForest(Location~.,data=bar_ranks_top_vars_30[rdt_pos,],importance=T, classwt= c(1/3,1/3,1/3))
m30_top_rdt_neg<- randomForest(Location~.,data=bar_ranks_top_vars_30[rdt_neg,],importance=T, classwt= c(1/3,1/3,1/3))

m65_top_child_5<- randomForest(Location~.,data=bar_ranks_top_vars_65[child_5,],importance=T, classwt= c(1/3,1/3,1/3))
m65_top_adult_5<- randomForest(Location~.,data=bar_ranks_top_vars_65[adult_5,],importance=T, classwt= c(1/3,1/3,1/3))
m65_top_child_15<- randomForest(Location~.,data=bar_ranks_top_vars_65[child_15,],importance=T, classwt= c(1/3,1/3,1/3))
m65_top_adult_15<- randomForest(Location~.,data=bar_ranks_top_vars_65[adult_15,],importance=T, classwt= c(1/3,1/3,1/3))
m65_top_rdt_pos<- randomForest(Location~.,data=bar_ranks_top_vars_65[rdt_pos,],importance=T, classwt= c(1/3,1/3,1/3))
m65_top_rdt_neg<- randomForest(Location~.,data=bar_ranks_top_vars_65[rdt_neg,],importance=T, classwt= c(1/3,1/3,1/3))

m130_top_child_5<- randomForest(Location~.,data=bar_ranks_top_vars_130[child_5,],importance=T, classwt= c(1/3,1/3,1/3))
m130_top_adult_5<- randomForest(Location~.,data=bar_ranks_top_vars_130[adult_5,],importance=T, classwt= c(1/3,1/3,1/3))
m130_top_child_15<- randomForest(Location~.,data=bar_ranks_top_vars_130[child_15,],importance=T, classwt= c(1/3,1/3,1/3))
m130_top_adult_15<- randomForest(Location~.,data=bar_ranks_top_vars_130[adult_15,],importance=T, classwt= c(1/3,1/3,1/3))
m130_top_rdt_pos<- randomForest(Location~.,data=bar_ranks_top_vars_130[rdt_pos,],importance=T, classwt= c(1/3,1/3,1/3))
m130_top_rdt_neg<- randomForest(Location~.,data=bar_ranks_top_vars_130[rdt_neg,],importance=T, classwt= c(1/3,1/3,1/3))

m260_top_child_5<- randomForest(Location~.,data=bar_ranks_top_vars_260[child_5,],importance=T, classwt= c(1/3,1/3,1/3))
m260_top_adult_5<- randomForest(Location~.,data=bar_ranks_top_vars_260[adult_5,],importance=T, classwt= c(1/3,1/3,1/3))
m260_top_child_15<- randomForest(Location~.,data=bar_ranks_top_vars_260[child_15,],importance=T, classwt= c(1/3,1/3,1/3))
m260_top_adult_15<- randomForest(Location~.,data=bar_ranks_top_vars_260[adult_15,],importance=T, classwt= c(1/3,1/3,1/3))
m260_top_rdt_pos<- randomForest(Location~.,data=bar_ranks_top_vars_260[rdt_pos,],importance=T, classwt= c(1/3,1/3,1/3))
m260_top_rdt_neg<- randomForest(Location~.,data=bar_ranks_top_vars_260[rdt_neg,],importance=T, classwt= c(1/3,1/3,1/3))


############
#Figure 2#
###########

global_oob_figure2_20<- c(10.26,8.04,8.37,9.16,6.31)
mutasa_oob_figure2_20<-c(2.56,3.97,1.92,4.65,2.70)
macha_oob_figure2_20<- c(31.25,9.46,14.29,12.04,6.06)
nchelenge_oob_figure2_20<-c(8.70,12.12,16.28,10.13,15.19)

global_oob_figure2_30<- c(6.41,7.79,5.42,7.33,5.05)
mutasa_oob_figure2_30<-c(2.56,3.97,0.96,4.65,2.70)
macha_oob_figure2_30<- c(25.00,8.78,10.71,8.33,3.79)
nchelenge_oob_figure2_30<-c(0.00,12.12,9.30,8.86,12.66)

global_oob_figure2_65<- c(7.69,5.78,6.90,5.86,4.29)
mutasa_oob_figure2_65<-c(2.56,3.97,0.96,4.65,2.70)
macha_oob_figure2_65<- c(25.00,6.08,14.29,6.48,3.79)
nchelenge_oob_figure2_65<-c(4.35,8.08,11.63,6.33,8.86)

global_oob_figure2_130<- c(6.41,6.78,6.40,6.23,4.29)
mutasa_oob_figure2_130<-c(0.00,3.97,1.92,4.65,1.08)
macha_oob_figure2_130<- c(25.00,7.43,14.29,8.33,5.30)
nchelenge_oob_figure2_130<-c(4.35,10.10,6.98,5.06,10.13)

global_oob_figure2_260<- c(6.41,7.04,5.91,8.79,4.80)
mutasa_oob_figure2_260<-c(0.00,3.97,1.92,4.65,1.62)
macha_oob_figure2_260<- c(25.00,9.46,14.29,11.11,5.30)
nchelenge_oob_figure2_260<-c(4.35,8.08,4.65,10.13,11.39)


figure_2_df<- as.data.frame(cbind(c(global_oob_figure2_20,mutasa_oob_figure2_20,macha_oob_figure2_20,nchelenge_oob_figure2_20,
                                     global_oob_figure2_30,mutasa_oob_figure2_30,macha_oob_figure2_30,nchelenge_oob_figure2_30,
                                     global_oob_figure2_65,mutasa_oob_figure2_65,macha_oob_figure2_65,nchelenge_oob_figure2_65,
                                     global_oob_figure2_130,mutasa_oob_figure2_130,macha_oob_figure2_130,nchelenge_oob_figure2_130,
                                     global_oob_figure2_260,mutasa_oob_figure2_260,macha_oob_figure2_260,nchelenge_oob_figure2_260),
                                   c(rep(c("Children<5","Adults>5","Children<15","Adults>15", "RDT Neg"),times=20)),
                                   c(rep(c("Global","Mutasa", "Choma", "Nchelenge"),each=5,times=5)),
                                   c(rep(c(20,30,65,130,260),each=20))))

colnames(figure_2_df)<- c("OOB Error", "Subset","Region", "Antigens")
figure_2_df$`OOB Error`<- as.numeric(figure_2_df$`OOB Error`)
figure_2_df$Antigens<- as.numeric(figure_2_df$Antigens)


figure_2<- ggplot(figure_2_df, aes(x=Antigens,y=`OOB Error`, group=Subset))+
  geom_line(aes(color=Subset))+
  geom_point(aes(shape=Subset,color=Subset))+
  theme_classic()+scale_color_brewer(palette = "Set1")+labs(x="Number of Antigens", y="Out of Bag Error (%)")+
  facet_wrap(~Region)+theme(legend.key.size = unit(1,"cm"),legend.text=element_text(size=10))

#####################################################################

###############################################
######AUC for 30, 65, 130, 260 antigens########
###############################################


pred_m20_child_5<- predict(m20_top_child_5,type="prob")
pred_m20_adult_5<- predict(m20_top_adult_5,type="prob")
pred_m20_child_15<- predict(m20_top_child_15,type="prob")
pred_m20_adult_15<- predict(m20_top_adult_15,type="prob")
pred_m20_rdt_neg<- predict(m20_top_rdt_neg,type="prob")

pred_matrix_m20_child_5<- as.data.frame(cbind(as.character(bar_ranks_top_vars_20[child_5,21]),pred_m20_child_5))
pred_matrix_m20_adult_5<- as.data.frame(cbind(as.character(bar_ranks_top_vars_20[adult_5,21]),pred_m20_adult_5))
pred_matrix_m20_child_15<- as.data.frame(cbind(as.character(bar_ranks_top_vars_20[child_15,21]),pred_m20_child_15))
pred_matrix_m20_adult_15<- as.data.frame(cbind(as.character(bar_ranks_top_vars_20[adult_15,21]),pred_m20_adult_15))
pred_matrix_m20_rdt_neg<- as.data.frame(cbind(as.character(bar_ranks_top_vars_20[rdt_neg,21]),pred_m20_rdt_neg))


pred_matrix_m20_child_5<- pred_matrix_m20_child_5%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

pred_matrix_m20_adult_5<- pred_matrix_m20_adult_5%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m20_child_15<- pred_matrix_m20_child_15%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m20_adult_15<- pred_matrix_m20_adult_15%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m20_rdt_neg<- pred_matrix_m20_rdt_neg%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


ci_roc_m20_child_5<- roc_auc_with_ci(pred_matrix_m20_child_5,type="basic",R=1000)
ci_roc_m20_adult_5<- roc_auc_with_ci(pred_matrix_m20_adult_5,type="basic",R=1000)
ci_roc_m20_child_15<- roc_auc_with_ci(pred_matrix_m20_child_15,type="basic",R=1000)
ci_roc_m20_adult_15<- roc_auc_with_ci(pred_matrix_m20_adult_15,type="basic",R=1000)
ci_roc_m20_rdt_neg<- roc_auc_with_ci(pred_matrix_m20_rdt_neg,type="basic",R=1000)


#######################################

pred_m30_child_5<- predict(m30_top_child_5,type="prob")
pred_m30_adult_5<- predict(m30_top_adult_5,type="prob")
pred_m30_child_15<- predict(m30_top_child_15,type="prob")
pred_m30_adult_15<- predict(m30_top_adult_15,type="prob")
pred_m30_rdt_neg<- predict(m30_top_rdt_neg,type="prob")

pred_matrix_m30_child_5<- as.data.frame(cbind(as.character(bar_ranks_top_vars_30[child_5,31]),pred_m30_child_5))
pred_matrix_m30_adult_5<- as.data.frame(cbind(as.character(bar_ranks_top_vars_30[adult_5,31]),pred_m30_adult_5))
pred_matrix_m30_child_15<- as.data.frame(cbind(as.character(bar_ranks_top_vars_30[child_15,31]),pred_m30_child_15))
pred_matrix_m30_adult_15<- as.data.frame(cbind(as.character(bar_ranks_top_vars_30[adult_15,31]),pred_m30_adult_15))
pred_matrix_m30_rdt_neg<- as.data.frame(cbind(as.character(bar_ranks_top_vars_30[rdt_neg,31]),pred_m30_rdt_neg))


pred_matrix_m30_child_5<- pred_matrix_m30_child_5%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

pred_matrix_m30_adult_5<- pred_matrix_m30_adult_5%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m30_child_15<- pred_matrix_m30_child_15%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m30_adult_15<- pred_matrix_m30_adult_15%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m30_rdt_neg<- pred_matrix_m30_rdt_neg%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


ci_roc_m30_child_5<- roc_auc_with_ci(pred_matrix_m30_child_5,type="basic",R=1000)
ci_roc_m30_adult_5<- roc_auc_with_ci(pred_matrix_m30_adult_5,type="basic",R=1000)
ci_roc_m30_child_15<- roc_auc_with_ci(pred_matrix_m30_child_15,type="basic",R=1000)
ci_roc_m30_adult_15<- roc_auc_with_ci(pred_matrix_m30_adult_15,type="basic",R=1000)
ci_roc_m30_rdt_neg<- roc_auc_with_ci(pred_matrix_m30_rdt_neg,type="basic",R=1000)


###################################################

pred_m65_child_5<- predict(m65_top_child_5,type="prob")
pred_m65_adult_5<- predict(m65_top_adult_5,type="prob")
pred_m65_child_15<- predict(m65_top_child_15,type="prob")
pred_m65_adult_15<- predict(m65_top_adult_15,type="prob")
pred_m65_rdt_neg<- predict(m65_top_rdt_neg,type="prob")

pred_matrix_m65_child_5<- as.data.frame(cbind(as.character(bar_ranks_top_vars_65[child_5,66]),pred_m65_child_5))
pred_matrix_m65_adult_5<- as.data.frame(cbind(as.character(bar_ranks_top_vars_65[adult_5,66]),pred_m65_adult_5))
pred_matrix_m65_child_15<- as.data.frame(cbind(as.character(bar_ranks_top_vars_65[child_15,66]),pred_m65_child_15))
pred_matrix_m65_adult_15<- as.data.frame(cbind(as.character(bar_ranks_top_vars_65[adult_15,66]),pred_m65_adult_15))
pred_matrix_m65_rdt_neg<- as.data.frame(cbind(as.character(bar_ranks_top_vars_65[rdt_neg,66]),pred_m65_rdt_neg))


pred_matrix_m65_child_5<- pred_matrix_m65_child_5%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

pred_matrix_m65_adult_5<- pred_matrix_m65_adult_5%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m65_child_15<- pred_matrix_m65_child_15%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m65_adult_15<- pred_matrix_m65_adult_15%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m65_rdt_neg<- pred_matrix_m65_rdt_neg%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


ci_roc_m65_child_5<- roc_auc_with_ci(pred_matrix_m65_child_5,type="basic",R=1000)
ci_roc_m65_adult_5<- roc_auc_with_ci(pred_matrix_m65_adult_5,type="basic",R=1000)
ci_roc_m65_child_15<- roc_auc_with_ci(pred_matrix_m65_child_15,type="basic",R=1000)
ci_roc_m65_adult_15<- roc_auc_with_ci(pred_matrix_m65_adult_15,type="basic",R=1000)
ci_roc_m65_rdt_neg<- roc_auc_with_ci(pred_matrix_m65_rdt_neg,type="basic",R=1000)


###################################################
pred_m130_child_5<- predict(m130_top_child_5,type="prob")
pred_m130_adult_5<- predict(m130_top_adult_5,type="prob")
pred_m130_child_15<- predict(m130_top_child_15,type="prob")
pred_m130_adult_15<- predict(m130_top_adult_15,type="prob")
pred_m130_rdt_neg<- predict(m130_top_rdt_neg,type="prob")

pred_matrix_m130_child_5<- as.data.frame(cbind(as.character(bar_ranks_top_vars_130[child_5,131]),pred_m130_child_5))
pred_matrix_m130_adult_5<- as.data.frame(cbind(as.character(bar_ranks_top_vars_130[adult_5,131]),pred_m130_adult_5))
pred_matrix_m130_child_15<- as.data.frame(cbind(as.character(bar_ranks_top_vars_130[child_15,131]),pred_m130_child_15))
pred_matrix_m130_adult_15<- as.data.frame(cbind(as.character(bar_ranks_top_vars_130[adult_15,131]),pred_m130_adult_15))
pred_matrix_m130_rdt_neg<- as.data.frame(cbind(as.character(bar_ranks_top_vars_130[rdt_neg,131]),pred_m130_rdt_neg))


pred_matrix_m130_child_5<- pred_matrix_m130_child_5%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

pred_matrix_m130_adult_5<- pred_matrix_m130_adult_5%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m130_child_15<- pred_matrix_m130_child_15%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m130_adult_15<- pred_matrix_m130_adult_15%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m130_rdt_neg<- pred_matrix_m130_rdt_neg%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


ci_roc_m130_child_5<- roc_auc_with_ci(pred_matrix_m130_child_5,type="basic",R=1000)
ci_roc_m130_adult_5<- roc_auc_with_ci(pred_matrix_m130_adult_5,type="basic",R=1000)
ci_roc_m130_child_15<- roc_auc_with_ci(pred_matrix_m130_child_15,type="basic",R=1000)
ci_roc_m130_adult_15<- roc_auc_with_ci(pred_matrix_m130_adult_15,type="basic",R=1000)
ci_roc_m130_rdt_neg<- roc_auc_with_ci(pred_matrix_m130_rdt_neg,type="basic",R=1000)



#####################################################
pred_m260_child_5<- predict(m260_top_child_5,type="prob")
pred_m260_adult_5<- predict(m260_top_adult_5,type="prob")
pred_m260_child_15<- predict(m260_top_child_15,type="prob")
pred_m260_adult_15<- predict(m260_top_adult_15,type="prob")
pred_m260_rdt_neg<- predict(m260_top_rdt_neg,type="prob")

pred_matrix_m260_child_5<- as.data.frame(cbind(as.character(bar_ranks_top_vars_260[child_5,261]),pred_m260_child_5))
pred_matrix_m260_adult_5<- as.data.frame(cbind(as.character(bar_ranks_top_vars_260[adult_5,261]),pred_m260_adult_5))
pred_matrix_m260_child_15<- as.data.frame(cbind(as.character(bar_ranks_top_vars_260[child_15,261]),pred_m260_child_15))
pred_matrix_m260_adult_15<- as.data.frame(cbind(as.character(bar_ranks_top_vars_260[adult_15,261]),pred_m260_adult_15))
pred_matrix_m260_rdt_neg<- as.data.frame(cbind(as.character(bar_ranks_top_vars_260[rdt_neg,261]),pred_m260_rdt_neg))


pred_matrix_m260_child_5<- pred_matrix_m260_child_5%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))

pred_matrix_m260_adult_5<- pred_matrix_m260_adult_5%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m260_child_15<- pred_matrix_m260_child_15%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m260_adult_15<- pred_matrix_m260_adult_15%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


pred_matrix_m260_rdt_neg<- pred_matrix_m260_rdt_neg%>%
  mutate(Macha_true= ifelse(V1=="Macha",1,0))%>%
  mutate(Nchelenge_true= ifelse(V1=="Nchelenge",1,0))%>%
  mutate(Mutasa_true= ifelse(V1=="Honde Valley",1,0))%>%
  mutate(Macha_pred_RF=Macha)%>%
  mutate(Nchelenge_pred_RF=Nchelenge)%>%
  mutate(Mutasa_pred_RF=`Honde Valley`)%>%
  select(c(Macha_true,Nchelenge_true,Mutasa_true,Macha_pred_RF,Nchelenge_pred_RF,Mutasa_pred_RF))


ci_roc_m260_child_5<- roc_auc_with_ci(pred_matrix_m260_child_5,type="basic",R=1000)
ci_roc_m260_adult_5<- roc_auc_with_ci(pred_matrix_m260_adult_5,type="basic",R=1000)
ci_roc_m260_child_15<- roc_auc_with_ci(pred_matrix_m260_child_15,type="basic",R=1000)
ci_roc_m260_adult_15<- roc_auc_with_ci(pred_matrix_m260_adult_15,type="basic",R=1000)
ci_roc_m260_rdt_neg<- roc_auc_with_ci(pred_matrix_m260_rdt_neg,type="basic",R=1000)


#####################################################
###################Antigen top 30 table#############
####################################################

tamaki_top_30<- c("PF3D7_0220000",
                  "PF3D7_1002100",
                  "PF3D7_0800200",
                  "PF3D7_0532100",
                  "PF3D7_1410400",
                  "PF3D7_0202500",
                  "PF3D7_0930300",
                  "PF3D7_1401400",
                  "PF3D7_1035900",
                  "PF3D7_0223300",
                  "PF3D7_1300300",
                  "PF3D7_0207000",
                  "PF3D7_1007700",
                  "PF3D7_0402400",
                  "PF3D7_0530100",
                  "PF3D7_0220000",
                  "PF3D7_0800300",
                  "PF3D7_1335300",
                  "PF3D7_1036400",
                  "PF3D7_0422100",
                  "PF3D7_0420700",
                  "PF3D7_0903500",
                  "PF3D7_0801000",
                  "PF3D7_1002000",
                  "PF3D7_0207700",
                  "PF3D7_0904900",
                  "PF3D7_0206800",
                  "PF3D7_1036000",
                  "PF3D7_0620400",
                  "PF3D7_0808600"
)





gene_IDs_in_tamaki_top_30<- which(importance_of_vars_ID$Gene.ID%in%tamaki_top_30)

tamaki_top_30_antigens<- which(colnames(bar_ranks)%in%importance_of_vars_ID$Index.x[gene_IDs_in_tamaki_top_30[1:44]])


greenhouse_antigens_days_since<- c("PF3D7_1002000",
                                   "PF3D7_0402400",
                                   "PF3D7_1106300",
                                   "PF3D7_0711700",
                                   "PF3D7_0800300",
                                   "PF3D7_0501100",
                                   "PF3D7_0423700",
                                   "PF3D7_1020800",
                                   "PF3D7_0731600",
                                   "PF3D7_1002100",
                                   "PF3D7_0223300",
                                   "PF3D7_0532100",
                                   "PF3D7_0702300",
                                   "PF3D7_0800200",
                                   "PF3D7_0801000",
                                   "PF3D7_0936300",
                                   "PF3D7_1033200",
                                   "PF3D7_1129100",
                                   "PF3D7_1300300",
                                   "PF3D7_1353100",
                                   "PF3D7_0304600",
                                   "PF3D7_0414700",
                                   "PF3D7_0420700",
                                   "PF3D7_0532400",
                                   "PF3D7_0620400",
                                   "PF3D7_0808600",
                                   "PF3D7_1024800",
                                   "PF3D7_1133400",
                                   "PF3D7_1438100",
                                   "PF3D7_1477500")
greenhouse_antigens_incidence<-c("PF3D7_1002000",
                                 "PF3D7_1020800",
                                 "PF3D7_0731600",
                                 "PF3D7_0532100",
                                 "PF3D7_0801000",
                                 "PF3D7_0304600",
                                 "PF3D7_0206800",
                                 "PF3D7_0800300",
                                 "PF3D7_1106300",
                                 "PF3D7_0423700",
                                 "PF3D7_0207000",
                                 "PF3D7_0402400",
                                 "PF3D7_0501100",
                                 "PF3D7_0532400",
                                 "PF3D7_0702300",
                                 "PF3D7_0711700",
                                 "PF3D7_0936300",
                                 "PF3D7_1002100",
                                 "PF3D7_1036400",
                                 "PF3D7_1129100",
                                 "PF3D7_0223300",
                                 "PF3D7_0620400",
                                 "PF3D7_0800200",
                                 "PF3D7_0930300",
                                 "PF3D7_1024800",
                                 "PF3D7_1100800",
                                 "PF3D7_1133400",
                                 "PF3D7_1300300",
                                 "PF3D7_1401400",
                                 "PF3D7_1410400")



gene_IDs_in_greenhouse_days_since<- which(importance_of_vars_ID$Gene.ID%in%greenhouse_antigens_days_since)

greenhouse_days_since_top_30_antigens<- which(colnames(bar_ranks)%in%importance_of_vars_ID$Index.x[gene_IDs_in_greenhouse_days_since[1:33]])

gene_IDs_in_greenhouse_incidence<- which(importance_of_vars_ID$Gene.ID%in%greenhouse_antigens_incidence)

greenhouse_incidence_top_30_antigens<- which(colnames(bar_ranks)%in%importance_of_vars_ID$Index.x[gene_IDs_in_greenhouse_incidence[1:40]])


colnames(bar_ranks[,tamaki_top_30_antigens])
colnames(bar_ranks_top_vars_30)

#########################################
########Repeat with pf only#############
########################################
library(gdata)



pf_indices<- importance_of_vars$Index[-which(startsWith(importance_of_vars$ID,"PVX"))]

pf_bar_ranks<- bar_ranks[,c(pf_indices,1039)]

set.seed(33)
m1<- randomForest(Location~.,data=bar_ranks,importance=T,classwt=c(1/3,1/3,1/3))




