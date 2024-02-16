#Figures of E_coli Loop analyses
source("scripts/functions/functions.R")

##FFL
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")
plast_medgrowthloops <- read.csv("scripts/data/plast_medium_E_coli_FFL.csv", sep = ",")
plast_ph1 <- read.csv("scripts/data/ph1_plast_E_coli_FFL.csv", sep = ",")
plast_ph2 <- read.csv("scripts/data/ph2_plast_E_coli_FFL.csv", sep = ",")##
plast_mg_c <- read.csv("scripts/data/mg_c_growth_E_coli_FFL.csv", sep = ",")
plast_aero <- read.csv("scripts/data/plast_aero_E_coli_FFL.csv", sep = ",")
plast_stringent <- read.csv("scripts/data/plast_stringent_E_coli_FFL.csv", sep = ",")
plast_tempr1 <- read.csv("scripts/data/plast_temptr1_E_coli_FFL.csv", sep = ",")
plast_tempr2 <- read.csv("scripts/data/plast_temptr2_E_coli_FFL.csv", sep = ",")
plast_juice <- read.csv("scripts/data/plast_juice_E_coli_FFL.csv", sep = ",")
plast_stress <- read.csv("scripts/data/plast_stress_E_coli_FFL.csv", sep = ",")
plast_ox <- read.csv("scripts/data/plast_ox_E_coli_FFL.csv", sep = ",")

df1 <- as.data.frame(rbind(colSums(rbind(non_plast,all_plast)[,3:12])*100/(nrow(rbind(non_plast,all_plast))),
            colSums(non_plast[,3:12])*100/nrow(non_plast),
            colSums(all_plast[,3:12])*100/nrow(all_plast),
            colSums(plast_medgrowthloops[,3:12])*100/nrow(plast_medgrowthloops),
            colSums(plast_ph1[,3:12])*100/nrow(plast_ph1),
            colSums(plast_ph2[,3:12])*100/nrow(plast_ph2),
            colSums(plast_aero[,3:12])*100/nrow(plast_aero),
            colSums(plast_mg_c[,3:12])*100/nrow(plast_mg_c),
            colSums(plast_stringent[,3:12])*100/nrow(plast_stringent),
            colSums(plast_tempr1[,3:12])*100/nrow(plast_tempr1),
            colSums(plast_tempr1[,3:12])*100/nrow(plast_tempr1),
            colSums(plast_juice[,3:12])*100/nrow(plast_juice),
            colSums(plast_stress[,3:12])*100/nrow(plast_stress),
            colSums(plast_ox[,3:12])*100/nrow(plast_ox)
            ))
rownames(df1) <- c("All genes","Non-plastic\ngenes", "Plastic\ngenes","Feugas\n2016", "Tucker\n2002", "Maurer\n2004","Ng\n2018","Caglar\n2017","Durfee\n2008","White-Ziegler\n2007","Kim\n2020","Bergholz\n2009","Bhatia\n2022","Wang\n2009")
df1[,3:10] <- df1[,c(order(as.character(colnames(df1[,3:10])), method = c("radix"))+2)]
colnames(df1) <- c("FFL","No_FFL",colnames(df1[,c(order(as.character(colnames(df1[,3:10])), method = c("radix"))+2)]))



##DIAMONDs
non_plast <- read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")
plast_medgrowthloops <- read.csv("scripts/data/plast_medium_E_coli_diamond.csv", sep = ",")
plast_ph1 <- read.csv("scripts/data/ph1_plast_E_coli_diamond.csv", sep = ",")
plast_ph2 <- read.csv("scripts/data/ph2_plast_E_coli_diamond.csv", sep = ",")##
plast_mg_c <- read.csv("scripts/data/mg_c_growth_E_coli_diamond.csv", sep = ",")
plast_aero <- read.csv("scripts/data/plast_aero_E_coli_diamond.csv", sep = ",")
plast_stringent <- read.csv("scripts/data/plast_stringent_E_coli_diamond.csv", sep = ",")
plast_tempr1 <- read.csv("scripts/data/plast_temptr1_E_coli_diamond.csv", sep = ",")
plast_tempr2 <- read.csv("scripts/data/plast_temptr2_E_coli_diamond.csv", sep = ",")
plast_juice <- read.csv("scripts/data/plast_juice_E_coli_diamond.csv", sep = ",")
plast_stress <- read.csv("scripts/data/plast_stress_E_coli_diamond.csv", sep = ",")
plast_ox <- read.csv("scripts/data/plast_ox_E_coli_diamond.csv", sep = ",")

df2 <- as.data.frame(rbind(colSums(rbind(non_plast,all_plast)[,3:14])*100/(nrow(rbind(non_plast,all_plast))),
                           colSums(non_plast[,3:14])*100/nrow(non_plast),
                           colSums(all_plast[,3:14])*100/nrow(all_plast),
                           colSums(plast_medgrowthloops[,3:14])*100/nrow(plast_medgrowthloops),
                           colSums(plast_ph1[,3:14])*100/nrow(plast_ph1),
                           colSums(plast_ph2[,3:14])*100/nrow(plast_ph2),
                           colSums(plast_aero[,3:14])*100/nrow(plast_aero),
                           colSums(plast_mg_c[,3:14])*100/nrow(plast_mg_c),
                           colSums(plast_stringent[,3:14])*100/nrow(plast_stringent),
                           colSums(plast_tempr1[,3:14])*100/nrow(plast_tempr1),
                           colSums(plast_tempr2[,3:14])*100/nrow(plast_tempr2),
                           colSums(plast_juice[,3:14])*100/nrow(plast_juice),
                           colSums(plast_stress[,3:14])*100/nrow(plast_stress),
                           colSums(plast_ox[,3:14])*100/nrow(plast_ox)
))
rownames(df2) <- c("All genes","Non-plastic\ngenes", "Plastic\ngenes","Feugas\n2016", "Tucker\n2002", "Maurer\n2004","Ng\n2018","Caglar\n2017","Durfee\n2008","White-Ziegler\n2007","Kim\n2020","Bergholz\n2009","Bhatia\n2022","Wang\n2009")
df2[,3:12] <- df2[,c(order(as.character(colnames(df2[,3:12])), method = c("radix"))+2)]
colnames(df2) <- c("FFL","No_FFL",colnames(df2[,c(order(as.character(colnames(df2[,3:12])), method = c("radix"))+2)]))


pdf("figures/Suppfig_3_type_FFL.pdf", width=15, height=10)
par(mar = c(8,4, 2,1), mfrow=c(2,1))
barplot(t(df1[,2:10]), col=c("grey","blue4","blue","dodgerblue3","skyblue3","lightskyblue","deepskyblue","cadetblue2","cyan3"), ylab="%",
        legend.text = c("No FFL",colnames(df1[,3:10])), args.legend = list(ncol=4, x = "topright", inset = c(0.3, 1.22)),
        main = "Feedforward motifs")

barplot(t(df2[,2:12]), col=c("grey","#E65100","#FF6C00","#FB8C00","#FFB300","#FFCC80","#FDD835","#FFEE58","#FFF59D","#EEFF91","#FFFDE7"), ylab="%",
        legend.text = c("No Diamond",colnames(df2[,3:12])), args.legend = list(ncol=4, x = "topright", inset = c(0.27, 1.22)),
        main = "Diamond  motifs")
dev.off()

