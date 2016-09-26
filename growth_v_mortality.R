#User provided input
path = "/Users/chrisgraves/Documents/Yeast_data/Microscopy_data/" #Path of data
expt_ID = c('CB_141007')#Experiment ID for analysis
max_t = 20 #Maximum number of time points recorded in all experiments
library(FNN,lib.loc = "/Library/Frameworks/R.framework/Versions/3.2/Resources/library/")
library(ggplot2)
library(reshape2)
library(dplyr)


#Load data from disk
for (j in 1:length(expt_ID)) {
  
  load(paste(path,expt_ID[j],"/d.Rfile",sep=""))
  load(paste(path,expt_ID[j],"/gr.Rfile",sep=""))
  load(paste(path,expt_ID[j],"/md.Rfile",sep=""))
  
  #Create data frame containing raw microcolony growth data
  #Microcolony areas
  areas = as.data.frame(cbind(d$areas,matrix(data=NA,nrow = dim(d$areas)[1],ncol = max_t-dim(d$areas)[2])))
  
  
  #NYU computed parameters
  dist = md$mindist #distance to the nearest microcolony
  growth = gr$colparam #growth rate calculations
  old_coords = cbind(md$cxm,md$cym)
  cenfield = md$cenfield
  
  
  
  #Experiment date
  date = unlist(strsplit(expt_ID[j],""))[4:9]
  date = paste(date[1],date[2],date[3],date[4],date[5],date[6],sep="")
  date = matrix(data = date,nrow=length(dist),ncol=1)
  
  #Extract strain IDs
  condition.variable.names<-c("Treatment",'Strain')
  number.condition.variables<-length(condition.variable.names)
  
  wells<-vector()
  for(i in 1:length(d$well.list)){
    colonies<-as.vector(unlist(d$well.list[i]))
    wells[colonies]<-names(d$well.list[i])
  }
  
  conditions<-as.vector(d$condition.names[wells])
  conditionsMat<-matrix(NA,length(conditions),number.condition.variables)
  
  if(number.condition.variables==1){
    for(i in 1:length(conditions)){
      conditionsMat[i,]<-conditions[i]
    }
  }
  else{
    for(i in 1:length(conditions)){
      spltCond<-strsplit(conditions[i], "-", fixed=T)
      spltCondVector<-unlist(spltCond)
      conditionsMat[i,]<-spltCondVector
    }
  }
  
  colnames(conditionsMat)<-condition.variable.names
  
  if (j==1) {
    df = cbind(as.data.frame(conditionsMat),growth,as.data.frame(date),dist,cenfield,old_coords,areas)
    colnames(df) = c("treatment","strain",'rate',"r",'foldx',"lag","rate_first","rate_last","date","dist",'cenfield','x','y','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','t11','t12','t13','t14','t15','t16','t17','t18','t19','t20')
    head}
  else{
    temp = cbind(as.data.frame(conditionsMat),growth,as.data.frame(date),dist,cenfield,old_coords,areas)
    colnames(temp) = c("treatment","strain",'rate',"r",'foldx',"lag","rate_first","rate_last","date","dist",'cenfield','x','y','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10','t11','t12','t13','t14','t15','t16','t17','t18','t19','t20')
    df = rbind(df,temp)
    
  }
}

head(df)
max(df$rate_last,na.rm=T)


#Determine whether microcolony grew after heat shock

areas_after = df[,(dim(df)[2]-14):dim(df)[2]]#Subset containing measurements after heat shock (after t=5)

good_data = subset(df,rowSums(areas_after!=0)>=8)

areas_after = good_data[,(dim(good_data)[2]-14):dim(good_data)[2]]


#function that determines whether a colony grew based on whether it's area increased by a 
#factor greater than the threshold argument
did_grow = function(areas_after,threshold) {
  area_0 = areas_after[1]
  area_end = areas_after[max(which(areas_after!=0))]
  if ((area_end/area_0) >= threshold) {
    return(1)
  }else {
    return(0)
  }
}



good_data$did_grow = apply(areas_after,1,did_grow,3) #add a column specifying whether colony grew
good_data$bins = cut(good_data$rate,breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),labels=c('0-0.1','0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5','0.5-0.6','0.6-0.7','0.7-0.8')) #bin growth rates 


calc_mort = function(did_grow) {
  num_obs = length(did_grow)
  num_grew = sum(did_grow)
  
  return(c(num_grew/num_obs,num_obs,num_grew))
}


P = subset(good_data,good_data$treatment=='P')

mort = tapply(P$did_grow,INDEX=P$bins,FUN=calc_mort)



mort_mat = matrix(nrow=8,ncol=3)


for (j in 1:dim(mort)) {
  if (length(unlist(mort[j])) == 3){
    mort_mat[j,] = unlist(mort[j])
  }
}

mort_mat = cbind(names(mort),mort_mat)
colnames(mort_mat) = c('growth_rate','p_alive','num_obs','num_grew')
mort_mat




#Re-calculate colony growth rates to inlcude slow-growing cells that are filtered out by NYU pipeline

pre_HS = df[,(dim(df)[2]-(max_t-1)):((dim(df)[2]-(max_t-1))+4)]
trimmed_df = subset(df,rowSums(pre_HS!=0)==5) #only include observations complete data before HS
pre_HS = trimmed_df[,(dim(trimmed_df)[2]-(max_t-1)):((dim(trimmed_df)[2]-(max_t-1))+4)]

x=c(0:4)
var_x = var(x)
mean_x = mean(x)
ls_5 = function(row) {
  row = log(row[1:5])
  b = cov(x,row)/var_x
  a = mean(row)-b*mean_x
  y_pred = a+b*x
  SS_tot = sum((row-mean(row))^2)
  SS_res = sum((row-y_pred)^2)
  r_sq = 1-(SS_res/SS_tot)
  var = var(row)
  foldx = exp(row[5])/exp(row[1])
  return(c(b,r_sq,var,foldx))
  
}

new_rates = apply(pre_HS,1,ls_5)

new_rates=t(new_rates)
colnames(new_rates) = c('new_rate','new_r_sq','new_var','new_foldx')
trimmed_df=cbind(trimmed_df,new_rates)
head(trimmed_df)
trimmed_df$new_rate[trimmed_df$new_rate<0] =0 #Change negative values of rate to zero

new_vs_old = ggplot(trimmed_df,aes(x=rate,y=new_rate))+
  geom_point()+
  geom_smooth(method=lm)

new_vs_old



#Filter out low quality estimates and re-run mortality analysis
good_data = subset(trimmed_df,trimmed_df$new_r_sq>0.9|abs(trimmed_df$new_rate)<0.1) #filter low r_sq
good_data$new_rate[good_data$new_rate<0]=0
areas_after = good_data[,(dim(good_data)[2]-18):(dim(good_data)[2]-4)]
good_data = subset(good_data,rowSums(areas_after!=0)>=8)
areas_after = good_data[,(dim(good_data)[2]-18):(dim(good_data)[2]-4)]

good_data$did_grow = apply(areas_after,1,did_grow,8) #add a column specifying whether colony grew
good_data$bins = cut(good_data$new_rate,breaks=c(-0.0001,0.15,0.3,0.45,0.6,0.75),labels=c('0-0.15','0.15-0.3','0.3-0.45','0.45-0.6','0.6-0.75')) #bin growth rates 


P = subset(good_data,good_data$treatment=="P")
C = subset(good_data,good_data$treatment=="C")
H = subset(good_data,good_data$treatment=="H")

head(P)

mort_P = tapply(P$did_grow,INDEX=P$bins,FUN=calc_mort)
mort_C = tapply(C$did_grow,INDEX=C$bins,FUN=calc_mort)
mort_H = tapply(H$did_grow,INDEX=H$bins,FUN=calc_mort)

mort_P

mort_mat = matrix(nrow=3*length(levels(good_data$bins)),ncol=4)
num_pts = dim(mort_P)



for (j in 1:num_pts) {
  
  
  mort_mat[j,2:4] = c(unlist(mort_P[j]))
  mort_mat[j+num_pts,2:4]=c(unlist(mort_C[j]))
  mort_mat[j+(2*num_pts),2:4] = c(unlist(mort_H[j]))
  
}

mort_mat = as.data.frame(mort_mat)
mort_mat[,1] = c(replicate(num_pts,'P'),replicate(num_pts,'C'),replicate(num_pts,'H'))
colnames(mort_mat) = c('treatment','p_alive','num_obs','num_grew')


mort_mat$SE = sqrt((mort_mat$p_alive*(1 - mort_mat$p_alive))/mort_mat$num_obs)
mort_mat$upper = mort_mat$p_alive + 2*mort_mat$SE
mort_mat$lower = mort_mat$p_alive - 2*mort_mat$SE
mort_mat$bin_center = c(seq(0.075,0.8,by=0.15),seq(0.075,0.8,by=0.15),seq(0.075,0.8,by=0.15))

mort_mat



model_P = glm(data=P,did_grow~new_rate,family='binomial')
model_C = glm(data=C,did_grow~new_rate,family='binomial')
model_H = glm(data=H,did_grow~new_rate,family='binomial')

x = seq(0.075,0.8,by=0.15)
mort_mat$pred = c(predict(model_P,list(new_rate=x),type="response"),predict(model_C,list(new_rate=x),type="response"),predict(model_H,list(new_rate=x),type="response"))

summary(model_P)


mort_plot = ggplot(data=mort_mat,aes(x=bin_center,y=p_alive,color = treatment))+
  geom_errorbar(aes(ymin=lower,ymax=upper,width=0.1))+
  geom_point()+
  geom_line(aes(x=bin_center,y=pred,color=treatment))


mort_plot

log_plot = ggplot(P,aes(x=new_rate,color=as.factor(did_grow),fill=as.factor(did_grow)))+
  geom_histogram(alpha=0.6)+
  scale_y_continuous(labels = percent_format())+
  labs(x='Growth Rate',y='Survival')+
  theme_bw()+
  theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"))+
  theme(legend.text=element_text(size=16,face="bold"),
        legend.title=element_text(size=18,face="bold"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=2))

log_plot

P_plot = ggplot(P,aes(x=new_rate,fill=as.factor(did_grow),color = as.factor(did_grow)))+
  geom_density(alpha = 0.5,aes(y=..scaled..))+
  labs(x='Growth rate',y='Frequency')+
  theme_bw()+
  theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"))+
  theme(legend.text=element_text(size=16,face="bold"),
        legend.title=element_text(size=18,face="bold"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=2))+ 
  scale_fill_manual(values=c("#0000FF", "#00FF33"), 
                    name="",
                    breaks=c("0","1"),
                    labels=c("Died", "Survived"))+
  scale_color_manual(values=c("#0000FF", "#00FF33"), 
                     name="",
                     breaks=c("0","1"),
                     labels=c("Died", "Survived"))
P_plot

H_plot = ggplot(H,aes(x=new_rate,fill=as.factor(did_grow),color = as.factor(did_grow)))+
  geom_density(alpha = 0.5,aes(y=..scaled..))+
  labs(x='Growth rate',y='Frequency')+
  theme_bw()+
  theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"))+
  theme(legend.text=element_text(size=16,face="bold"),
        legend.title=element_text(size=18,face="bold"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=2))+ 
  scale_fill_manual(values=c("#99FF00", "#009900"), 
                    name="",
                    breaks=c("0","1"),
                    labels=c("Died", "Survived"))+
  scale_color_manual(values=c("#99FF00", "#009900"), 
                     name="",
                     breaks=c("0","1"),
                     labels=c("Died", "Survived"))
H_plot

C_plot = ggplot(C,aes(x=new_rate,fill=as.factor(did_grow),color = as.factor(did_grow)))+
  geom_density(alpha = 0.5,aes(y=..scaled..))+
  labs(x='Growth rate',y='Frequency')+
  theme_bw()+
  theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=18,face="bold"))+
  theme(legend.text=element_text(size=16,face="bold"),
        legend.title=element_text(size=18,face="bold"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=2))+ 
  scale_fill_manual(values=c("#00CCFF", "#0000FF"), 
                    name="",
                    breaks=c("0","1"),
                    labels=c("Died", "Survived"))+
  scale_color_manual(values=c("#00CCFF", "#0000FF"), 
                     name="",
                     breaks=c("0","1"),
                     labels=c("Died", "Survived"))
C_plot