R
source("read.admb.R")
library(tidyverse)
library(ggplot2)
library(GGally)
library(ggridges,quietly=TRUE)
setwd("~/Google Drive/Atka tag model/2018 tag model/Main/runs")
mytheme <- theme(panel.grid.major.x = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank() )
mytheme <- mytheme + theme(text=element_text(size=18)) + theme(axis.title.x=element_text(size=24) ,axis.title.y=element_text(size=24))
mytheme <- mytheme + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),panel.grid.minor.y = element_line(colour="grey60", linetype="dashed"), panel.grid.major.y = element_blank() )
mytheme <- mytheme + theme( panel.background = element_rect(fill="white"), panel.border = element_rect(colour="black", fill=NA, size=1))

#----Immigration (recruitment) factor estimated ------------------------
mc4 <- read.table("results/both_mcmc.rep",header=T,row.names=NULL)
mc1 <- read.table("results/both_mcmc.rep",header=T,row.names=NULL)
mc2 <- read.table("results/nearshore_mcmc.rep",header=T,row.names=NULL)
mc3 <- read.table("results/Seamounts_mcmc.rep",header=T,row.names=NULL)
mc4$Biomass <- mc2$Biomass + mc3$Biomass
mc <- rbind(mc1,mc2,mc3,mc4)
#mc$trial <- c(rep('543 Combined',4999),rep('543 Nearshore',4999),rep('543 Seamounts',4999),rep('541 Seguam',4999))
mc$Area <- c(rep('543 Combined',4999),rep('543 Nearshore',4999),rep('543 Seamounts',4999),rep('Nearshore+Seamounts',4999))
ggplot(mc,aes(x=Biomass,fill = Area)) + geom_density(alpha=0.5) + mytheme + xlim(c(0,1000))  + xlab("Biomass (kt)") + ggtitle("Immigration estimated")

#----Immigration (recruitment) factor fixed at 1.0----------------------
mc <- read.table("results/both_no_RF_mcmc.rep",header=T,row.names=NULL)
mc <- rbind(mc,read.table("results/nearshore_no_RF_mcmc.rep",header=T,row.names=NULL))
mc <- rbind(mc,read.table("results/seamounts_no_RF_mcmc.rep",header=T,row.names=NULL))
#mc$trial <- c(rep('543 Combined',4999),rep('543 Nearshore',4999),rep('543 Seamounts',4999),rep('541 Seguam',4999))
mc$Area <- c(rep('543 Combined',4999),rep('543 Nearshore',4999),rep('543 Seamounts',4999))
ggplot(mc,aes(x=Biomass,fill = Area)) + geom_density(alpha=0.5) + mytheme + xlim(c(0,1000))  + xlab("Biomass (kt)") + ggtitle("No immigration estimated") +  scale_y_continuous(breaks=NULL)
#theme(axis.title.y=element_blank())
#----Seguam sensitivities ---------------------------------
mcp <- data.frame(read.table("results/Seguam_M_mcmc.rep",header=T,row.names=NULL))[,c(1,3:6,16)] %>% sample_n(2000) 
ggpairs(mcp,aes(color="blue",alpha=.3))
names(mc)

#mc <- read.table("res_seguam_sensitivity/Seguam_no_RF_mcmc.rep",header=T,row.names=NULL)
#mc <- rbind(mc,read.table("res_seguam_sensitivity/Seguam_mcmc.rep",header=T,row.names=NULL))
#mc <- rbind(mc,read.table("res_seguam_sensitivity/Seguam_M_mcmc.rep",header=T,row.names=NULL))
#mc <- rbind(mc,read.table("res_seguam_sensitivity/Seguam_no_RF_M_mcmc.rep",header=T,row.names=NULL))
#mc <- read.table("results/Seguam_no_RF_mcmc.rep",header=T,row.names=NULL)
#mc <- rbind(mc,read.table("results/Seguam_mcmc.rep",header=T,row.names=NULL))
#mc <- rbind(mc,read.table("results/Seguam_M_mcmc.rep",header=T,row.names=NULL))
#mc <- rbind(mc,read.table("results/Seguam_no_RF_M_mcmc.rep",header=T,row.names=NULL))
##mc$trial <- c(rep('543 Combined',4999),rep('543 Nearshore',4999),rep('543 Seamounts',4999),rep('541 Seguam',4999))
##mc$Area <- c(rep('543 Combined',4999),rep('543 Nearshore',4999),rep('543 Seamounts',4999))
##mc$Model <- c(rep('No immig.',4999), rep('Immig. est.',4999), rep('Immig. est., CV M 0.25',4999), rep('No immig., CV M 0.25',4999) )
#mc$Model <- c(rep('C',4999),
#	rep('I',4999),
#	rep('IM',4999),
#	rep('M',4999)
#	)
mc <- read.table("results/Seguam_mcmc.rep",header=T,row.names=NULL)
#mc <- read.table("results/Seguam1_mcmc.rep",header=T,row.names=NULL)
#mc1 <- read.table("results/Seguam1b_mcmc.rep",header=T,row.names=NULL)
mc1 <- read.table("results/Seguam1_mcmc.rep",header=T,row.names=NULL)
mc34 <- read.table("results/Seguam34_mcmc.rep",header=T,row.names=NULL)
#mc2 <- mc1
#mc2$Biomass <- mc1$Biomass / mc34$Biomass
#ggplot(mc2,aes(x=Biomass)) + geom_density(alpha=0.5,fill="salmon") + mytheme + xlim(c(0,3))  + xlab("Area 1 biomass / areas 3&4 biomass") + ggtitle("Seguam area")+ geom_vline(xintercept = 1.)
#ggplot(mc2,aes(x=Biomass)) + stat_ecdf(size=2,color="salmon") + mytheme + xlim(c(0,2))  + xlab("Area 1 biomass / areas 3&4 biomass") + ylab("Probability")+ ggtitle("Seguam area") + geom_vline(xintercept = 1.)
#mc <- rbind(mc,mc2,mc3,mc1)
mc <- rbind(mc,mc2,mc34)
mc$Model <- c(rep('All strata',4999),
	rep('Stratum 1',4999),
	rep('Strata 3 & 4',4999))
mc$Model <- c(rep('Seguam all strata',4999),
	rep('Seguam strata 3 & 4',4999),
	rep('All - (3+4) ',4999),
	rep('Seguam strata 1 ',4999))
#mc$Area <- c(rep('Seguam no immigration',4999))
mcp <- sample_n(mc,2000)[,c(1,3:6,17)]
ggpairs(mcp,aes(color=Model,alpha=.2))
ggplot(mc,aes(x=Biomass,fill = Model)) + geom_density(alpha=0.5) + mytheme + xlim(c(0,1000))  + xlab("Biomass (kt)") + ggtitle("Seguam area")
ggplot(mc,aes(x=Biomass,color = Model)) + stat_ecdf(size=2) + mytheme + xlim(c(0,1000))  + xlab("Biomass (kt)") + ggtitle("Seguam area")
stat_ecdf()
mapping = NULL, data = NULL, geom = "step",
  position = "identity", ..., n = NULL, pad = TRUE, na.rm = FALSE,
  show.legend = NA, inherit.aes = TRUE)

#----updated datafiles------------------------
#----updated datafiles------------------------
mc <- read.table("results/both2_mcmc.rep",header=T,row.names=NULL)
mc <- rbind(mc,read.table("results/nearshore2_mcmc.rep",header=T,row.names=NULL))
mc <- rbind(mc,read.table("results/seamounts2_mcmc.rep",header=T,row.names=NULL))
mc <- rbind(mc,read.table("results/Seguam2_mcmc.rep",header=T,row.names=NULL))
mc$trial <- c(rep('543 Combined',4999),rep('543 Nearshore',4999),rep('543 Seamounts',4999),rep('541 Seguam',4999))

ggplot(mc,aes(x=Biomass,fill = trial)) + geom_density(alpha=0.5) + mytheme + xlim(c(0,1000))  + xlab("Biomass (kt)") + ggtitle("Up updated datafiles")


ggplot(mc,aes(x=Biomass,fill = trial)) + geom_density(alpha=0.5) + mytheme + xlim(c(0,1000))  + xlab("Biomass (kt)")
ggplot(mc,aes(x=N,fill = trial)) + geom_density() + mytheme + xlim(c(0,9500)) 
ggplot(mc,aes(x=N,fill = trial)) + geom_density() + mytheme + xlim(c(0,5500)) 

mc %>% filter(trial=="Seguam") %>%      sample_n(500) %>% select(1,3,4,5,6,7,8) %>% ggpairs() + mytheme
mc %>% filter(trial=="Both") %>%      sample_n(500) %>% select(1,4,5,6,7) %>% ggpairs() + mytheme
mc %>% filter(trial=="Seamounts") %>%      sample_n(500) %>% select(1,4,5,6,7) %>% ggpairs() + mytheme
mc %>% filter(trial=="Nearshore") %>%      sample_n(500) %>% select(1,4,5,6,7) %>% ggpairs() + mytheme
mc %>% filter(trial=="Nearshore") %>% sample_n(500) %>% select(1,3,4,5,6,7,8) %>% ggpairs() + mytheme
mc %>% filter(trial=="Seamounts") %>% sample_n(500) %>% select(1,3,4,5,6,7,8) %>% ggpairs() + mytheme
ggpairs(mc)
stat = "identity",scale = 1, alpha = .3,
p1  <- ggplot(mc,aes(x=age,y=as.factor(Year),height = sel)) + geom_density_ridges(stat = "identity",scale = 1, alpha = .3,
  	     fill=fill,color="black") + .THEME +
         xlim(c(1,lage))+ ylab("Year") + xlab("Age (years)") + scale_y_discrete(limits=rev(levels(as.factor(sdf$Year))))
names(mc)



getwd()
setwd("R")
system("../src/run.bat pet_all")
mc=read.table("../src/arc/pet_all_s2_mcout.rep",skip=1,header=F)
mc=read.table("../src/arc/pet_all_mcout.rep",skip=1,header=F)
colnames(mc)=c("B","N","p_loss","Rep_Rate_F","Ref_Rate_SF","M","ER","Survival","ObjFun","ObjF_Tags","ObjF_Repr","ObjF_M","x","x","x")
plot(density(mc$B),xlim=c(0,800),xlab="Biomass")
mc=read.table("../src/arc/seg_all_mcout.rep",skip=1,header=F)
mc=read.table("../src/arc/tan_all_mcout.rep",skip=1,header=F)
mc=read.table("../src/arc/pet_011_mcout.rep",skip=1,header=F)
mc=read.table("../src/arc/seg_011_mcout.rep",skip=1,header=F)
mc=read.table("../src/arc/tan_011_mcout.rep",skip=1,header=F)
.MyPlot("../src/arc/pet_all_mcout.rep")
.MyLines("../src/arc/pet_011_mcout.rep",col="green")
.MyLines("../src/arc/pet_all_s3_mcout.rep",col="red")
b
.MyPlot <- function( repObj = "../src/arc/tan_011_mcout.rep",xlim=c(0,700),main="Biomass"){
  mc=read.table(repObj,skip=1,header=F)
  colnames(mc)=c("B","N","p_loss","Rep_Rate_F","Ref_Rate_SF","M","ER","Survival","ObjFun","ObjF_Tags","ObjF_Repr","ObjF_M","x","x","x")
  plot(density(mc$B),xlim=xlim,main=main,xlab="Biomass")
}

.MyLines <- function( repObj = "../src/arc/tan_011_mcout.rep",col=1,lty=1,lwd=1 ){
  mc=read.table(repObj,skip=1,header=F)
  colnames(mc)=c("B","N","p_loss","Rep_Rate_F","Ref_Rate_SF","M","ER","Survival","ObjFun","ObjF_Tags","ObjF_Repr","ObjF_M","x","x","x")
  lines(density(mc$B),lty=lty,lwd=lwd,col=col)
}

quantile(mc$B,.8)
lines(density(mc$B),lty=1,lwd=2,col="red")
density(mc$B)
names(mc)
dim(mc)

p_all=read.rep("../src/arc/pet_all.rep")
p_alls1=read.rep("../src/arc/pet_all_s1.rep")
p_alls2=read.rep("../src/arc/pet_all_s2.rep")
p_alls3=read.rep("../src/arc/pet_all_s3.rep")
p_alls4=read.rep("../src/arc/pet_all_s4.rep")
plot(p_all$Obs_Tags,ylim=c(0,170),xlab="Event",ylab="Tags returned",cex=1.5)
points(p_all$Pred_Tags,pch=19,cex=1.2)
points(p_alls1$Pred_Tags,pch=19,cex=1.2,col="red")
points(p_alls2$Pred_Tags,pch=19,cex=1.2,col="blue")
points(p_alls3$Pred_Tags,pch=2,cex=1.2,col="blue")
points(p_alls4$Pred_Tags,pch=18,cex=1.2,col="red")
# For 2011
p_011=read.rep("../src/arc/pet_011.rep")
plot(p_011$Obs_Tags,ylim=c(0,170),xlab="Event",ylab="Tags returned",cex=1.5)
points(p_011$Pred_Tags,pch=19,cex=1.2)

# Seguam
p_seg=read.rep("../src/arc/seg_all.rep")
p_segs1=read.rep("../src/arc/seg_all_s1.rep")
p_segs2=read.rep("../src/arc/seg_all_s2.rep")
p_segs3=read.rep("../src/arc/seg_all_s3.rep")
p_seg011=read.rep("../src/arc/seg_011.rep")
p_seg011s1=read.rep("../src/arc/seg_011_s1.rep")

plot(p_seg$Obs_Tags,ylim=c(0,170),xlab="Event",ylab="Tags returned",cex=1.5)
points(p_seg$Pred_Tags,pch=19,cex=1.2)

plot(p_segs1$Obs_Tags,ylim=c(0,170),xlab="Event",ylab="Tags returned",cex=1.5)
points(p_segs1$Pred_Tags,pch=19,cex=1.2,col="red")
points(p_segs2$Pred_Tags,pch=18,cex=1.2)
points(p_segs3$Pred_Tags,pch=5,cex=1.8)

plot(p_seg011$Obs_Tags,ylim=c(0,180),xlab="Event",ylab="Tags returned",cex=1.5)
points(p_seg011$Pred_Tags,pch=19,cex=1.2)

plot(p_seg011s1$Obs_Tags,ylim=c(0,180),xlab="Event",ylab="Tags returned",cex=1.5)
points(p_seg011s1$Pred_Tags,pch=19,cex=1.2)

# Tanaga 
p_tan=read.rep("../src/arc/tan_all.rep")



p_all=read.rep("../src/arc/pet_all.rep")
r_2011=read.table("../src/mcout1.rep",header=T)
r_2011_od=read.table("../src/mcout.rep",header=T)
r_all=read.table("../src/mcout_all.rep",header=T)
r_all_od= r_all
#=read.table("../src/mcout.rep",header=T)
pairs(r_all[,1:7],pch=19,cex=.6)

plot(density(r_all$N),ylim=c(0,.013),xlim=c(0,1000),xlab="Population abundance",main="Seguam",col="green")
lines(density(r_2011$N),col="red",lwd=1)
lines(density(r_2011_od$N),col="red",lwd=2,lty=2)
lines(density(r_all_od$N),col="green",lwd=2,lty=2)
windows()
plot(density(r_all$M),xlim=c(0.0,0.4),xlab="M ",main="Seguam",col="green")
lines(density(r_2011$M),col="red",lwd=1)
lines(density(r_2011_od$M),col="red",lwd=2,lty=2)
lines(density(r_all_od$M),col="green",lwd=2,lty=2)


