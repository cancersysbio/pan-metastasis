#######
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


library(ggplot2)
library(gmodels)
library(spatstat)

#setwd("dir/data")

dat <- read.delim("paired_primary_met_data.txt", header=T)
mono <- dat[which(dat$clonality=="monoclonal" & dat$met_site != "LN"),]
met_status <- rep("0",nrow(mono))

met_status[which(mono$time_to_relapse_mons<3)] = "synchronous"
met_status[which(mono$time_to_relapse_mons>=3)] = "metachronous"
mono <- data.frame(mono,met_status=met_status)

######################################################
pdf("hist_metastatic_timing_3cancers.pdf", width=12, height=2.5)
par(mfrow = c(1, 3))

bc <- mono[which(mono$cancer_type=="Breast"),]
crc <- mono[which(mono$cancer_type=="Colorectal"),]
lc <- mono[which(mono$cancer_type=="Lung"),]

#####breast
mtime <- data.frame(cancer_type=bc$subtype, sample=bc$met_name, ratio=bc$LmLp_ratio,met_status=bc$met_status, time_to_relapse=bc$time_to_relapse_mons)
Rt_fraction <- 1-mtime$ratio*0.13
bc_norm_mtime <- Rt_fraction[which(Rt_fraction>0)]
summary(bc_norm_mtime)
tmet <- Rt_fraction*4.6
mtime1 <- data.frame(mtime, tmet=tmet)
print(nrow(mtime1))
mtime1 <- mtime1[which(mtime1$tmet*12 > -mtime1$time_to_relapse),]
print(nrow(mtime1))
print(summary(mtime1$tmet))


p1 = ggplot(mtime1, aes(tmet, fill = met_status)) + 
   geom_histogram(alpha = 1, aes(y = ..count..)) + xlim(-7,5) +
   geom_vline(xintercept=median(mtime1$tmet), color="black",linetype="dashed", size=1)+
   ggtitle("Breast cancer")+
   xlab("ts") +
   theme(axis.text.x = element_text(size=14, margin=margin(5,0,0,0), hjust=0.8, colour="black"),
        axis.text.y = element_text(size=14), axis.title.x =element_text(size=14), axis.title.y =element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
		    panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))		

######colon
mtime <- data.frame(cancer_type=crc$cancer_type, sample=crc$met_name,ratio=crc$LmLp_ratio,met_status=crc$met_status, time_to_relapse=crc$time_to_relapse_mons)
Rt_fraction <- 1-mtime$ratio*0.13
crc_norm_mtime <- Rt_fraction[which(Rt_fraction>0)]
summary(crc_norm_mtime)
tmet <- Rt_fraction*5.2
mtime2 <- data.frame(mtime, tmet=tmet)
print(nrow(mtime2))
mtime2 <- mtime2[which(mtime2$tmet*12 > -mtime2$time_to_relapse),]
print(nrow(mtime2))
print(summary(mtime2$tmet))


p2 = ggplot(mtime2, aes(tmet, fill = met_status)) + 
   geom_histogram(alpha = 1, aes(y = ..count..)) +  xlim(-5,5) +
   geom_vline(xintercept=median(mtime2$tmet), color="black",linetype="dashed", size=1)+
   ggtitle("Colorectal cancer")+
   xlab("ts") +
   theme(axis.text.x = element_text(size=14, margin=margin(5,0,0,0), hjust=0.8, colour="black"),
        axis.text.y = element_text(size=14), axis.title.x =element_text(size=14), axis.title.y =element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
		    panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))		



######lung
mtime <- data.frame(cancer_type=lc$cancer_type, sample=lc$met_name, ratio=lc$LmLp_ratio,met_status=lc$met_status, time_to_relapse=lc$time_to_relapse_mons)
Rt_fraction <- 1-mtime$ratio*0.13
lc_norm_mtime <- Rt_fraction[which(Rt_fraction>0)]
summary(lc_norm_mtime)
tmet <- Rt_fraction*4.3
mtime3 <- data.frame(mtime, tmet=tmet)
print(nrow(mtime3))
mtime3 <- mtime3[which(mtime3$tmet*12 > -mtime3$time_to_relapse),]
print(nrow(mtime3))
print(summary(mtime3$tmet))


p3 = ggplot(mtime3, aes(tmet, fill = met_status)) + 
   geom_histogram(alpha = 1, aes(y = ..count..)) +  xlim(-5,5) +
   geom_vline(xintercept=median(mtime3$tmet), color="black",linetype="dashed", size=1)+
   ggtitle("Lung cancer")+
   xlab("ts") +
   theme(axis.text.x = element_text(size=14, margin=margin(5,0,0,0), hjust=0.8, colour="black"),
        axis.text.y = element_text(size=14), axis.title.x =element_text(size=14), axis.title.y =element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
		    panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))		


multiplot(p2,p3,p1,cols=3)

dev.off()


#######################################
norm_time <- data.frame(cancer_type=c(rep("Breast",length(bc_norm_mtime)),rep("Colorectal",length(crc_norm_mtime)),rep("Lung",length(lc_norm_mtime))),tmet=c(bc_norm_mtime,crc_norm_mtime,lc_norm_mtime))

pdf("normalized_met_timing.pdf", width=5, height=3)

ggplot(norm_time, aes(x=cancer_type, y=tmet, fill=cancer_type)) + 
  geom_jitter(shape=16, size=2, aes(color = cancer_type), position=position_jitter(0.2)) +
  geom_boxplot(width=0.4, outlier.shape = NA, alpha=0) + ylim(0,1) +
  scale_x_discrete(limits=c("Colorectal","Lung","Breast")) +
  #stat_summary(fun.y=mean, geom="point", size=2, color="red") +
  #stat_summary(fun.y=median, geom="point", size=2, color="blue") +
  #scale_fill_manual(values=c("#df65b0","#41ae76","#fe9929"),guide=FALSE) +
  xlab("") + ylab("normalized ts") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14, margin=margin(5,0,0,0), angle=30, hjust=0.8, colour="black"),
        axis.text.y = element_text(size=16), axis.title.x =element_text(size=16), axis.title.y =element_text(size=16), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))   

dev.off()


mtimex <- rbind(mtime1,mtime2)
mtime_all <- rbind(mtimex,mtime3)
syn <- mtime_all[which(mtime_all$met_status=="synchronous"),]
meta <- mtime_all[which(mtime_all$met_status=="metachronous"),]
mtime_all <- rbind(syn,meta)
wilcox.test(syn$tmet, meta$tmet)


pdf("synchronous_vs_metachromous_met_timing.pdf", width=4, height=3)

ggplot(mtime_all, aes(x=met_status, y=tmet, fill=met_status)) + 
  geom_jitter(shape=16, size=2, aes(color = met_status), position=position_jitter(0.2)) +
  geom_boxplot(width=0.5, outlier.shape = NA, alpha=0) + ylim(-7,5) +
  scale_x_discrete(limits=c("synchronous","metachronous")) +
  #stat_summary(fun.y=mean, geom="point", size=2, color="red") +
  #stat_summary(fun.y=median, geom="point", size=2, color="blue") +
  #scale_fill_manual(values=c("#df65b0","#41ae76","#fe9929"),guide=FALSE) +
  xlab("") + ylab("ts") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14, margin=margin(5,0,0,0), angle=30, hjust=0.8, colour="black"),
        axis.text.y = element_text(size=16), axis.title.y =element_text(size=16), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))   

dev.off()


#########breast cancer subtypes##########
mtime1$subtype = mtime1$cancer_type

mtime_bc <- mtime1[which(mtime1$subtype=="ER+/HER2-" | mtime1$subtype=="ER-/HER2+" | mtime1$subtype=="ER+/HER2+" | mtime1$subtype=="TN"),]
pdf("BC_met_timing_by_subtypes_1.pdf", width=5, height=3)

ggplot(mtime_bc, aes(x=subtype, y=tmet, fill=subtype)) + 
  geom_jitter(shape=16, size=3, aes(color = subtype), position=position_jitter(0.2)) +
  geom_boxplot(width=0.5, outlier.shape = NA, alpha=0) + ylim(-7,5) +
  scale_x_discrete(limits=c("ER+/HER2+","ER-/HER2+","ER+/HER2-","TN")) +
  #stat_summary(fun.y=mean, geom="point", size=2, color="red") +
  #stat_summary(fun.y=median, geom="point", size=2, color="blue") +
  #scale_fill_manual(values=c("#df65b0","#41ae76","#fe9929"),guide=FALSE) +
  xlab("") + ylab("ts") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14, margin=margin(5,0,0,0), angle=30, hjust=0.8, colour="black"),
        axis.text.y = element_text(size=16), axis.title.y =element_text(size=16), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
    	panel.background = element_blank(), axis.line = element_line(colour = "black"))   

dev.off()


mtime_bc <- mtime1[which(mtime1$subtype=="ER+/HER2-" | mtime1$subtype=="ER-/HER2+" | mtime1$subtype=="ER+/HER2+" | mtime1$subtype=="TN"),]
her2_neg <- mtime_bc[which(mtime1$subtype=="ER+/HER2-" | mtime1$subtype=="TN"),]
her2_pos <- mtime_bc[which(mtime1$subtype=="ER-/HER2+" | mtime1$subtype=="ER+/HER2+"),]
mtime_bc <- data.frame(subtype=c(rep("HER2-",nrow(her2_neg)),rep("HER2+",nrow(her2_pos))), tmet=c(her2_neg$tmet,her2_pos$tmet))


pdf("BC_met_timing_by_subtypes_2.pdf", width=3, height=3)

ggplot(mtime_bc, aes(x=subtype, y=tmet, fill=subtype)) + 
  geom_jitter(shape=16, size=3, aes(color = subtype), position=position_jitter(0.2)) +
  geom_boxplot(width=0.5, outlier.shape = NA, alpha=0) + ylim(-7,5) +
  scale_x_discrete(limits=c("HER2-","HER2+")) +
  #stat_summary(fun.y=mean, geom="point", size=2, color="red") +
  #stat_summary(fun.y=median, geom="point", size=2, color="blue") +
  #scale_fill_manual(values=c("#df65b0","#41ae76","#fe9929"),guide=FALSE) +
  xlab("") + ylab("ts") +
  theme_bw() +
  theme(axis.text.x = element_text(size=14, margin=margin(5,0,0,0), angle=30, hjust=0.8, colour="black"),
        axis.text.y = element_text(size=16), axis.title.y =element_text(size=16), legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
    	panel.background = element_blank(), axis.line = element_line(colour = "black"))   

dev.off()


################
mtime_all <- read.delim("met_timing_3cancers_final.txt",header=T)

pdf("correlation_ts_vs_diagnosis_span.pdf", width=5, height=2.5)

p1<-ggplot(mtime_all, aes(x=time_to_relapse_years, y=tmet)) + 
  geom_point(aes(colour = factor(cancer_type)), size = 2) +
  scale_color_manual(values=c("#8c6bb1", "#2b8cbe", "#d95f02")) + 
  geom_smooth(method=lm, color='black') + xlim(-1,25)+ylim(-7,5) + theme_bw() +
  xlab("Time span between Dx of P and M") +  ylab("ts") +
  theme(axis.text.x = element_text(size=12, hjust=0.5, colour="black"),
        axis.text.y = element_text(size=12), axis.title.y =element_text(size=12), axis.title.x =element_text(size=12), plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), axis.line = element_line(colour = "black"))	

 print(p1)

 dev.off()


reg <- lm(mtime_all$tmet ~ mtime_all$time_to_relapse_years)
cor.test(mtime_all$tmet, mtime_all$time_to_relapse_years, method="spearman")


