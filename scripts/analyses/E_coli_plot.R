#Figures of E_coli Loop analyses
source("scripts/functions/functions.R")

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

df <- as.data.frame(rbind(colSums(non_plast[,3:12])*100/nrow(non_plast),
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
rownames(df) <- c("NP", "All_plast","NP_Medgth", "Medgth", "pH","Aero","Mg_C","Stringent","T°","Juice","Stress","Ox")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:10)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,10)]))


pdfname <- "figures/E_coli"
pdf(paste0(pdfname,"_type_FFL2_cutoff_3",".pdf"), width=12, height=6)
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







pdf("figures/FFL_distrib_e_coli.pdf", width=10, height=4)
layout(matrix(c(1), 1, 1, byrow = TRUE))
par(mar=c(2, 2, 2, 2), mgp = c(1.75, 0.75, 0), las=0)
barplot(t(df[,1:3]), col=c(7,3,"grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c(7,3,"grey"),
       legend=c("Coherent FFL", "Incoherent FFL", "No FFL"))
dev.off()
