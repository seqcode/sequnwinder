library(dplyr)
library(ggplot2)

early = "#4C2C7B"
active = "#3852A3"
inactive = "#6BBD45"
shared = "#F17E20"
late = "#E92528"

dat = tbl_df(read.table("3d_Oct4_ChIPSeq_tagcounts.mat",header=TRUE))

vals = c("early"=early,"late"=late,"shared"=shared,"inactive"=inactive,"active"=active)

tmpE = dat %>% select(EB,eClass) %>% mutate(class=eClass)
tmpd = dat %>% select(EB,dClass) %>% mutate(class=dClass)
tmp = bind_rows(tmpE,tmpd)

#ggplot(oct_df,aes(eClass,log(EB),fill=dClass,color=eClass))+geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=vals)+scale_color_manual(values=vals)
ggplot(tmp,aes(class,log(EB),fill=class))+geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("early"=early,"late"=late,"shared"=shared,"inactive"=inactive,"active"=active))
ggsave("3d_Oct4_ChIPSeq_tagcounts_boxplot.png")


dat = tbl_df(read.table("3c_Onecut2_ChIPSeq_tagcounts.mat",header=TRUE))

vals = c("early"=early,"late"=late,"shared"=shared,"inactive"=inactive,"active"=active)

tmpE = dat %>% select(EB48h,eClass) %>% mutate(class=eClass)
tmpd = dat %>% select(EB48h,dClass) %>% mutate(class=dClass)
tmp = bind_rows(tmpE,tmpd)

#ggplot(oct_df,aes(eClass,log(EB),fill=dClass,color=eClass))+geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=vals)+scale_color_manual(values=vals)
ggplot(tmp,aes(class,log(EB48h),fill=class))+geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("early"=early,"late"=late,"shared"=shared,"inactive"=inactive,"active"=active))
ggsave("3c_Onecut2_ChIPSeq_tagcounts_boxplot.png")


nil_df = tbl_df(read.table("3c_d_motifHitcounts.mat",header=TRUE))
tmpE = nil_df %>% select(Oct4MotifScore,esClass) %>% mutate(class=esClass)
tmpd = nil_df %>% select(Oct4MotifScore,dynClass) %>% mutate(class=dynClass)
tmp = bind_rows(tmpE,tmpd)
ggplot(tmp,aes(class,Oct4MotifScore,fill=class))+geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=c("early"=early,"late"=late,"shared"=shared,"inactive"=inactive,"active"=active))
ggsave("3d_Oct4_motif_logodds_boxplot.png")


tmpE = nil_df %>% select(OnecutMotifScore,esClass) %>% mutate(class=esClass)
tmpd = nil_df %>% select(OnecutMotifScore,dynClass) %>% mutate(class=dynClass)
tmp = bind_rows(tmpE,tmpd)

ggplot(tmp,aes(class,OnecutMotifScore,fill=class))+geom_boxplot(outlier.shape=NA)+scale_fill_manual(values=vals)
ggsave("3c_Onecut2_motif_logodds_boxplot.png")
