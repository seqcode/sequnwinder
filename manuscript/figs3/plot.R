library(dplyr)
library(ggplot2)


dat = tbl_df(read.table("cognateMotifHitInfoForAllFacs.tab",header=TRUE))

GM12878 = "#3852A4"
K562= "#4C297C"
H1hESC ="#69BD45"
shared="#231F20"
distal="red"
proximal = "green"

ggplot(dat,aes(fac,cognateScore,fill=cClass)) + geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.85))+scale_fill_manual(values = c("GM12878"=GM12878,"K562"=K562,"H1hESC"=H1hESC,"Shared"=shared))
ggsave("cognate_motif_hitinfo_boxplot.png",width=10.7,height=6.3)
