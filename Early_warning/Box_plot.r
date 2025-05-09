##### Box plots illustrating the distribution of Kendall's tau values.
# The following code can be use to reproduce the box plot figures presented in the main text.
# Since the code structure is identical across models and indicators,
# we demonstrate it here using the shallow lake model as an example.
# This example can serve as a guide for generating figures for other models as well.


library(reshape2)
library(ggplot2)

data=read.csv(file.choose()) #This will open a file dialog box to select the file you want to open in R.
#choose the proper csv file containing 150 Kendall tau values for each species response correlation in columns
# for instance, to plot Kendall tau values for AR1 in shallow-lake model we choose  "lake150_M_ktau_AR1_k0.0.csv" (for macrophyte)
data1=data.frame(data[,2:6])
colnames(data1)=c("0.0","0.225","0.45","0.675","0.9")
# choose the csv file "lake150_A_ktau_AR1_k0.0.csv" for species 2 (Algae population)
data=read.csv(file.choose())
data2=data.frame(data[,2:6])
colnames(data2)=c("0.0","0.225","0.45","0.675","0.9")

d1=melt(data1)
d2=melt(data2)
d3=cbind("sp"=rep("M",750),d1)
d4=cbind("sp"=rep("A",750),d2)
d5=rbind(d3,d4)

ggplot(d5,aes(x=variable, y=value,fill=sp))+
  geom_boxplot(position = position_dodge(0.8), width=0.65,alpha=0.1,outlier.shape = NA)+
  geom_point(aes(x=variable,y=value,color=sp),position = position_jitterdodge(0.25),alpha=0.3)+
  theme(legend.position = "top",legend.title = element_blank() ,panel.background = element_rect(fill="white"),panel.grid.major.y = element_line(size=0.25,color="grey"),
        panel.grid.minor.y = element_line(size=0.25,color="grey"),panel.grid.minor.x = element_blank(),panel.grid.major.x = element_blank(),
        panel.border = element_rect(fill = "NA",color = "black"),axis.text =element_text(color="black",size = 13),
        plot.title = element_text(size=17,vjust = -8),axis.title.x = element_text(size=17),axis.title.y = element_text(size=17))+
  xlab("Species response correlation")+ylab("Kendall Tau")+ggtitle("AR1")


# For the indicator SD follow the same procedure and choose the files for kendall tau values of SD for each species.
