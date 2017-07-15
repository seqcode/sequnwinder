library(dplyr)
library(ggplot2)

early = "#4C2C7B"
active = "#3852A3"
inactive = "#6BBD45"
shared = "#F17E20"
late = "#E92528"

model_scores_for_plotting = read.table("model_scores_for_plotting.tab",header=TRUE)

dat = tbl_df(model_scores_for_plotting)

tmp = dat %>% filter(Fac == "Onecut") %>% filter(Model=="SeqUnwinder")
ggplot(tmp,aes(Class,Score,fill=Class))+geom_bar(stat="identity")+ylim(c(-1,0.8))+scale_fill_manual(values=c("shared"=shared,"inactive"=inactive,"active"=active,"late"=late,"early"=early))
ggsave("OnecutSeqBar.png")
tmp = dat %>% filter(Fac == "Onecut") %>% filter(Model=="MCC")
ggplot(tmp,aes(Class,Score,fill=Class))+geom_bar(stat="identity")+ylim(c(-1,0.8))+scale_fill_manual(values=c("shared"=shared,"inactive"=inactive,"active"=active,"late"=late,"early"=early))
ggsave("OnecutMCCBar.png")

tmp = dat %>% filter(Fac == "Zfp281") %>% filter(Model=="SeqUnwinder")
ggplot(tmp,aes(Class,Score,fill=Class))+geom_bar(stat="identity")+ylim(c(-1,0.8))+scale_fill_manual(values=c("shared"=shared,"inactive"=inactive,"active"=active,"late"=late,"early"=early))
ggsave("Zfp281SeqBar.png")
tmp = dat %>% filter(Fac == "Zfp281") %>% filter(Model=="MCC")
ggplot(tmp,aes(Class,Score,fill=Class))+geom_bar(stat="identity")+ylim(c(-1,0.8))+scale_fill_manual(values=c("shared"=shared,"inactive"=inactive,"active"=active,"late"=late,"early"=early))
ggsave("Zfp281MCCBar.png")

tmp = dat %>% filter(Fac == "Oct4") %>% filter(Model=="SeqUnwinder")
ggplot(tmp,aes(Class,Score,fill=Class))+geom_bar(stat="identity")+ylim(c(-1,0.8))+scale_fill_manual(values=c("shared"=shared,"inactive"=inactive,"active"=active,"late"=late,"early"=early))
ggsave("Oct4SeqBar.png")
tmp = dat %>% filter(Fac == "Oct4") %>% filter(Model=="MCC")
ggplot(tmp,aes(Class,Score,fill=Class))+geom_bar(stat="identity")+ylim(c(-1,0.8))+scale_fill_manual(values=c("shared"=shared,"inactive"=inactive,"active"=active,"late"=late,"early"=early))
ggsave("Oct4MCCBar.png")

tmp = dat %>% filter(Fac == "Lhx3") %>% filter(Model=="SeqUnwinder")
ggplot(tmp,aes(Class,Score,fill=Class))+geom_bar(stat="identity")+ylim(c(-1,0.8))+scale_fill_manual(values=c("shared"=shared,"inactive"=inactive,"active"=active,"late"=late,"early"=early))
ggsave("Lhx3SeqBar.png")
tmp = dat %>% filter(Fac == "Lhx3") %>% filter(Model=="MCC")
ggplot(tmp,aes(Class,Score,fill=Class))+geom_bar(stat="identity")+ylim(c(-1,0.8))+scale_fill_manual(values=c("shared"=shared,"inactive"=inactive,"active"=active,"late"=late,"early"=early))
ggsave("Lhx3MCCBar.png")

