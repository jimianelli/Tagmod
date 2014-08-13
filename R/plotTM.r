getwd()
setwd("R")
system("../src/run.bat pet_all")
source("read.admb.R")
?read.table
mc=read.table("../src/arc/pet_all_s2_mcout.rep",skip=1,header=F)
mc=read.table("../src/arc/pet_all_mcout.rep",skip=1,header=F)
names(mc)
mc=read.table("../src/arc/seg_all_mcout.rep",skip=1,header=F)
mc=read.table("../src/arc/pet_all_mcout.rep",skip=1,header=F)
colnames(mc)=c("Biom","Nmale","NFemale","p_loss","Rep_Rate_F","Ref_Rate_SF","M","ER","Survival","ObjFun","ObjF_Tags","ObjF_Repr","ObjF_M","x","x")
plot(density(mc$ER))
plot(density(mc$B),xlim=c(0,800),xlab="Biomass")
mcnames=1:10
mcnames
mcnames[1]="../src/arc/seg_all_mcout.rep"
mcnames[2]="../src/arc/tan_all_mcout.rep"
mcnames[3]="../src/arc/pet_all_mcout.rep"
mcnames[2]="../src/arc/pet_all_s1_mcout.rep"
mcnames[3]="../src/arc/pet_all_s2_mcout.rep"
mcnames[4]="../src/arc/pet_all_s3_mcout.rep"
mcnames[5]="../src/arc/pet_all_s4_mcout.rep"
mcnames[6]="../src/arc/pet_all_s5_mcout.rep"
mcnames[7]="../src/arc/pet_011_mcout.rep"
mcnames[8]="../src/arc/tan_all_mcout.rep"
mcnames[9]="../src/arc/tan_011_mcout.rep"
mcnames[10]="../src/arc/seg_all_mcout.rep"
mcnames[11]="../src/arc/seg_011_mcout.rep"
mcnames[12]="../src/arc/seg_all_s1_mcout.rep"
mcnames[13]="../src/arc/seg_all_s2_mcout.rep"
mcnames[14]="../src/arc/seg_all_s3_mcout.rep"
mcnames[15]="../src/arc/seg_all_s4_mcout.rep"
mcnames[16]="../src/arc/seg_all_s5_mcout.rep"

.MyPlot(mcnames[1],xlim=c(0,550))
.MyLines(mcnames[2] ,col="green")
.MyLines(mcnames[3] ,col="red")
.MyLines(mcnames[4] ,col="green")
.MyLines(mcnames[5] ,col="red")
.MyLines(mcnames[6] ,col="red",lty=2,lwd=3)
.MyLines(mcnames[7] ,col="red",lty=2,lwd=1)

.MyLines(mcnames[8] ,col="red",lty=2,lwd=3)
.MyLines(mcnames[9] ,col="red",lty=2,lwd=1)

.MyPlot2(mcnames[1])
.MyPlot2(mcnames[2])
.MyPlot2(mcnames[3])
.MyPlot(mcnames[10])
.MyLines(mcnames[11])
.MyLines(mcnames[12] ,col="green")
.MyLines(mcnames[13] ,col="red")
.MyLines(mcnames[14] ,col="green")
.MyLines(mcnames[15] ,col="red")
.MyLines(mcnames[16] ,col="red",lty=2,lwd=3)




.MyPlot2 <- function( repObj = "../src/arc/tan_011_mcout.rep",xlim=c(0,.2),main="Exploitation rate"){
  mc=read.table(repObj,skip=1,header=F)
  colnames(mc)=c("B","Nmales","NFemale","p_loss","Rep_Rate_F","Ref_Rate_SF","M","ER","Survival","ObjFun","ObjF_Tags","ObjF_Repr","ObjF_M","x","x","x")
  plot(density(mc$ER),xlim=xlim,main=main,xlab="Exploitation rate")
}

.MyPlot <- function( repObj = "../src/arc/tan_011_mcout.rep",xlim=c(0,700),main="Biomass"){
  mc=read.table(repObj,skip=1,header=F)
  colnames(mc)=c("B","Nmales","NFemale","p_loss","Rep_Rate_F","Ref_Rate_SF","M","ER","Survival","ObjFun","ObjF_Tags","ObjF_Repr","ObjF_M","x","x","x")
  plot(density(mc$B),xlim=xlim,main=main,xlab="Biomass")
}

.MyLines <- function( repObj = "../src/arc/tan_011_mcout.rep",col=1,lty=1,lwd=1 ){
  mc=read.table(repObj,skip=1,header=F)
  lines(density(mc$B),lty=lty,lwd=lwd,col=col)
  colnames(mc)=c("B","Nmales","NFemale","p_loss","Rep_Rate_F","Ref_Rate_SF","M","ER","Survival","ObjFun","ObjF_Tags","ObjF_Repr","ObjF_M","x","x","x")
  # print(density(mc$B)) print(density(mc$B)$x) print(density(mc$B)$y) names(density(mc$B))
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


