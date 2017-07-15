library(ggplot2)
library(dplyr)

GM12878 = "#3852A4"
K562= "#4C297C"
H1hESC ="#69BD45"
shared="#231F20"

All_TFs_distances <- read.table("All_TFs_distances.tab", header=TRUE)
dat = tbl_df(All_TFs_distances)
ggplot(dat,aes(log(distance),color=type))+geom_density()+facet_wrap(~factor,dir="v",scales = "free_x" )+scale_x_continuous(breaks = c(5,10,15),labels = c("150","22k","3mil"))+scale_color_manual(values = c("shared"=shared,"K562"=K562,"GM12878"=GM12878,"H1hESC"=H1hESC))
ggsave("s2_distance_from_TSS.png")
