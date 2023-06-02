#Figures of E_coli Loop analyses
source("scripts/functions/functions.R")

##FFL
non_plast <- read.csv("scripts/data/nonplast_FFL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_FFL.csv", sep = ",")
np_medgrowthloops <- read.csv("scripts/data/np_medium_growth_FFL.csv", sep = ",")
plast_medgrowthloops <- read.csv("scripts/data/plast_medium_growth_FFL.csv", sep = ",")
plast_ph <- read.csv("scripts/data/ph_plast_FFL.csv", sep = ",")
plast_mg_c <- read.csv("scripts/data/mg_c_growth_FFL.csv", sep = ",")
plast_aero <- read.csv("scripts/data/plast_aero_FFL.csv", sep = ",")
plast_stringent <- read.csv("scripts/data/plast_stringent_FFL.csv", sep = ",")
plast_tempr <- read.csv("scripts/data/plast_temptr_FFL.csv", sep = ",")
plast_juice <- read.csv("scripts/data/plast_juice_FFL.csv", sep = ",")
plast_stress <- read.csv("scripts/data/plast_stress_FFL.csv", sep = ",")
plast_ox <- read.csv("scripts/data/plast_ox_FFL.csv", sep = ",")

df <- as.data.frame(rbind(colSums(rbind(non_plast,all_plast)[,3:12])*100/(nrow(non_plast)+nrow(all_plast)),
            colSums(non_plast[,3:12])*100/nrow(non_plast),
            colSums(all_plast[,3:12])*100/nrow(all_plast),
            colSums(np_medgrowthloops[,3:12])*100/nrow(np_medgrowthloops),
            colSums(plast_medgrowthloops[,3:12])*100/nrow(plast_medgrowthloops),
            colSums(plast_ph[,3:12])*100/nrow(plast_ph),
            colSums(plast_aero[,3:12])*100/nrow(plast_aero),
            colSums(plast_mg_c[,3:12])*100/nrow(plast_mg_c),
            colSums(plast_stringent[,3:12])*100/nrow(plast_stringent),
            colSums(plast_tempr[,3:12])*100/nrow(plast_tempr),
            colSums(plast_juice[,3:12])*100/nrow(plast_juice),
            colSums(plast_stress[,3:12])*100/nrow(plast_stress),
            colSums(plast_ox[,3:12])*100/nrow(plast_ox)
            ))
rownames(df) <- c("All","NP", "All_plast","NP_Medgth", "Medgth", "pH","Aero","Mg_C","Stringent","T°","Juice","Stress","Ox")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:10)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,10)]))


pdfname <- "figures/E_coli"
pdf(paste0(pdfname,"_type_FFL",".pdf"), width=14, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,2:10]), col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
       legend=c( "No_FFL","Input_Dep_Amplifying_neg","Input_Dep_Amplifying_pos","Input_Dep_Disruptive_neg","Input_Dep_Disruptive_pos","Input_Ind_Amplifying_neg","Input_Ind_Amplifying_pos","Input_Ind_Disruptive_neg","Input_Ind_Disruptive_pos") )
#Coherence
barplot(t(df[,c(2,11,12)]), col=c("grey","indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(2,13,14)]), col=c("grey","orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","orange","yellowgreen"),
       legend=c( "No_FFL","Homogenous","Heterogenous") )
dev.off()



#########################
##Diamond
non_plast <- read.csv("scripts/data/nonplast_diamond.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_diamond.csv", sep = ",")
np_medgrowthloops <- read.csv("scripts/data/np_medium_growth_diamond.csv", sep = ",")
plast_medgrowthloops <- read.csv("scripts/data/plast_medium_growth_diamond.csv", sep = ",")
plast_ph <- read.csv("scripts/data/ph_plast_diamond.csv", sep = ",")
plast_mg_c <- read.csv("scripts/data/mg_c_growth_diamond.csv", sep = ",")
plast_aero <- read.csv("scripts/data/plast_aero_diamond.csv", sep = ",")
plast_stringent <- read.csv("scripts/data/plast_stringent_diamond.csv", sep = ",")
plast_tempr <- read.csv("scripts/data/plast_temptr_diamond.csv", sep = ",")
plast_juice <- read.csv("scripts/data/plast_juice_diamond.csv", sep = ",")
plast_stress <- read.csv("scripts/data/plast_stress_diamond.csv", sep = ",")
plast_ox <- read.csv("scripts/data/plast_ox_diamond.csv", sep = ",")

df <- as.data.frame(rbind(colSums(rbind(non_plast,all_plast)[,3:12])*100/(nrow(non_plast)+nrow(all_plast)),
                          colSums(non_plast[,3:14])*100/nrow(non_plast),
                          colSums(all_plast[,3:14])*100/nrow(all_plast),
                          colSums(np_medgrowthloops[,3:14])*100/nrow(np_medgrowthloops),
                          colSums(plast_medgrowthloops[,3:14])*100/nrow(plast_medgrowthloops),
                          colSums(plast_ph[,3:14])*100/nrow(plast_ph),
                          colSums(plast_aero[,3:14])*100/nrow(plast_aero),
                          colSums(plast_mg_c[,3:14])*100/nrow(plast_mg_c),
                          colSums(plast_stringent[,3:14])*100/nrow(plast_stringent),
                          colSums(plast_tempr[,3:14])*100/nrow(plast_tempr),
                          colSums(plast_juice[,3:14])*100/nrow(plast_juice),
                          colSums(plast_stress[,3:14])*100/nrow(plast_stress),
                          colSums(plast_ox[,3:14])*100/nrow(plast_ox)
))
rownames(df) <- c("All","NP", "All_plast","NP_Medgth", "Medgth", "pH","Aero","Mg_C","Stringent","T°","Juice","Stress","Ox")
df$coherent <- rowSums2(as.matrix(df[,c(3,5,6,8,10,11)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(4,7,9,12)]))
df$homogenous_pos <- rowSums2(as.matrix(df[,c(3,6,9)]))
df$homogenous_neg <- rowSums2(as.matrix(df[,c(5,8,12)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(4,7,10,11)]))


pdfname <- "figures/E_coli"
pdf(paste0(pdfname,"_type_diamond",".pdf"), width=14, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,2:12])*100/300, col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","darkorchid1","plum1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","darkorchid1","plum1"),
       legend=c( "No_FFL","Pos_pos","Pos_mixt","Pos_neg",
                 "Neg_pos","Neg_mixt","Neg_neg",
                 "Mixt_pos","Mixt_mixt_hom","Mixt_mixt_het","Mixt_neg") )
#Coherence
barplot(t(df[,c(2,13,14)])*100/300, col=c("grey","indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(2,15:17)])*100/300, col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3"),
       legend=c( "No_FFL","Homogenous_pos","Homogenous_neg","Heterogenous") )
dev.off()

