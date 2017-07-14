
library(ggplot2)
library(dplyr)

dat = tbl_df(read.table("2d_F1_scores.tab",header=TRUE))

datAvgF1 = dat %>%  group_by(Tag, Method) %>% summarise(avgF1=mean(F1)) 

ggplot(datAvgF1,aes(Tag,avgF1))+geom_jitter(aes(color=Method,size=3),width=0,height=0.03)+ylim(-0.03,1.03)+ylab("Average F1 Score")+scale_x_discrete("Label Overlap Range (fraction)", labels = c("1_0.5-0.55" = "[0.50,0,55)","2_0.55-0.60" = "[0.55,0.60)","3_0.60-0.65" = "[0.60,0.65)","4_0.65-0.70" = "[0.65,0.70)","5_0.70-0.75" = "[0.70,0.75)","6_0.75-0.80" = "[0.75,0.80)", "7_0.80-0.85" = "[0.80,0.85)", "8_0.85-0.9" = "[0.85,0.90)", "8_0.90-0.95" = "[0.90,0.95)", "8_0.95-1" = "[0.95,1)"))+theme(axis.text = element_text(size=20),axis.text.x = element_text(angle = 60, hjust = 1,size=20),legend.title= element_text(size=25),legend.text= element_text(size=20), title=element_text(size=35),axis.title=element_text(size=25))+ggtitle("F1 Score at \nDifferent Overlaps ")
ggsave("F1.png",width=10,height=7)

dat = tbl_df(read.table("2e_and_2f_FP_and_TP.tab",sep="\t",header=TRUE))

glimpse(dat)

#ggplot(dat,aes(Tag,TP))+geom_dotplot(binaxis="y",dotsize=0.5,aes(fill=Method))



ggplot(dat,aes(Tag,TP))+geom_jitter(aes(color=Method,size=3),width=0,height=0.03)+ylim(-0.03,1.03)+ylab("True Positive Rate")+scale_x_discrete("Label Overlap Range (fraction)", labels = c("1_0.5-0.55" = "[0.50,0,55)","2_0.55-0.60" = "[0.55,0.60)","3_0.60-0.65" = "[0.60,0.65)","4_0.65-0.70" = "[0.65,0.70)","5_0.70-0.75" = "[0.70,0.75)","6_0.75-0.80" = "[0.75,0.80)", "7_0.80-0.85" = "[0.80,0.85)", "8_0.85-0.9" = "[0.85,0.90)", "8_0.90-0.95" = "[0.90,0.95)", "8_0.95-1" = "[0.95,1)"))+ggtitle("True Positive Rates \nat Different Overlaps ")+theme(axis.title=element_text(size=25),axis.text= element_text(size=20), axis.text.x = element_text(angle = 60, hjust = 1,size=20),legend.title= element_text(size=25),legend.text= element_text(size=20), title=element_text(size=35))

ggsave("TP.png",width=10,height=7)
ggplot(dat,aes(Tag,FP))+geom_jitter(aes(color=Method,size=3),width=0,height=0.03)+ylim(-0.03,1.03)+ylab("False Positive Rate")+scale_x_discrete("Label Overlap Range (fraction)", labels = c("1_0.5-0.55" = "[0.50,0,55)","2_0.55-0.60" = "[0.55,0.60)","3_0.60-0.65" = "[0.60,0.65)","4_0.65-0.70" = "[0.65,0.70)","5_0.70-0.75" = "[0.70,0.75)","6_0.75-0.80" = "[0.75,0.80)", "7_0.80-0.85" = "[0.80,0.85)", "8_0.85-0.9" = "[0.85,0.90)", "8_0.90-0.95" = "[0.90,0.95)", "8_0.95-1" = "[0.95,1)"))+ggtitle("False Positive Rates \nat Different Overlaps ")+theme(axis.title=element_text(size=25),axis.text= element_text(size=20), axis.text.x = element_text(angle = 60, hjust = 1,size=20),legend.title= element_text(size=25),legend.text= element_text(size=20), title=element_text(size=35))


ggsave("FP.png",width=10,height=7)

ggplot(dat,aes(Tag,FN))+geom_jitter(aes(color=Method,size=3),width=0,height=0.03)+ylim(-0.03,1.03)+xlab("Overlap Range")+ylab("False Negative Rate")+scale_x_discrete("Label Overlap Range (fraction)", labels = c("1_0.5-0.55" = "[0.50,0,55)","2_0.55-0.60" = "[0.55,0.60)","3_0.60-0.65" = "[0.60,0.65)","4_0.65-0.70" = "[0.65,0.70)","5_0.70-0.75" = "[0.70,0.75)","6_0.75-0.80" = "[0.75,0.80)", "7_0.80-0.85" = "[0.80,0.85)", "8_0.85-0.9" = "[0.85,0.90)", "8_0.90-0.95" = "[0.90,0.95)", "8_0.95-1" = "[0.95,1)"))+ggtitle("False Negative Rates at Different Overlaps ")+theme(axis.title=element_text(size=13,face="bold"),axis.text= element_text(size=12,face="bold"), axis.text.x = element_text(angle = 60, hjust = 1,size=12,face="bold"),legend.title= element_text(size=13,face="bold"),legend.text= element_text(size=12,face="bold"), title=element_text(size=17,face="bold"))
ggsave("FN.png")
