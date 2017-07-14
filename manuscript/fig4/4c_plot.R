
library(dplyr)
library(ggplot2)

GM12878_overlaps <- read.delim("4c_GM12878_overlaps.tab", quote="")
K562_overlaps <- read.delim("4c_K562_overlaps.tab", quote="")
H1hESC_overlaps <- read.delim("4c_H1hESC_overlaps.tab", quote="")

gdata = tbl_df(GM12878_overlaps)
kdata = tbl_df(K562_overlaps)
hdata = tbl_df(H1hESC_overlaps)

GM12878 = "#3852A4"
K562= "#4C297C"
H1hESC ="#69BD45"
shared="#231F20"
distal="red"
proximal = "green"

tmpG = gdata %>% filter(tClass == "Distal") %>% mutate(Motif=ifelse(CognateMotifScore < Cognate0.01,"Absent","Present"),CellType="GM12878") %>% select(Peak,tClass,cClass,Fac,numCofacs,Motif,CellType)
tmpK = kdata %>% filter(tClass == "Distal") %>% mutate(Motif=ifelse(CognateMotifScore < Cognate0.01,"Absent","Present"),CellType="K562") %>% select(Peak,tClass,cClass,Fac,numCofacs,Motif,CellType)
tmpH = hdata %>% filter(tClass == "Distal") %>% mutate(Motif=ifelse(CognateMotifScore < Cognate0.01,"Absent","Present"),CellType="H1hESC") %>% select(Peak,tClass,cClass,Fac,numCofacs,Motif,CellType)

ggplot(tmpK,aes(cClass,numCofacs,fill=Motif))+geom_boxplot()
ggsave("4c_k562_cofacdegree.png")
ggplot(tmpG,aes(cClass,numCofacs,fill=Motif))+geom_boxplot()
ggsave("4c_GM12878_cofacdegree.png")
ggplot(tmpH,aes(cClass,numCofacs,fill=Motif))+geom_boxplot()
ggsave("4c_H1hESC_cofacdegree.png")
