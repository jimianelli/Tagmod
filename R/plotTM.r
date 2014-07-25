getwd()
setwd("R")
system("../src/run.bat pet_all")
source("read.admb.R")
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


