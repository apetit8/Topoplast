#Figures of E_coli Loop analyses
source("scripts/functions/functions.R")

all_genes <- read.csv("scripts/data/all_genes_FFL.csv", sep = ",")[,2:5]
non_plast <- read.csv("scripts/data/nonplast_ffloops.csv", sep = ",")[,2:5]
all_plast <- read.csv("scripts/data/plast_genes_ffloops.csv", sep = ",")[,2:5]
np_medgrowthloops <- read.csv("scripts/data/np_medium_growth_ffloops.csv", sep = ",")[,2:5]
plast_medgrowthloops <- read.csv("scripts/data/plast_medium_growth_ffloops.csv", sep = ",")[,2:5]
plast_ph <- read.csv("scripts/data/ph_plast_FFL.csv", sep = ",")[,2:5]
plast_mg_c <- read.csv("scripts/data/plast_mg_c_ffloops.csv", sep = ",")[,2:5]
plast_aero <- read.csv("scripts/data/plast_aero_ffloops.csv", sep = ",")[,2:5]
plast_stringent <- read.csv("scripts/data/plast_stringent_ffloops.csv", sep = ",")[,2:5]
plast_tempr <- read.csv("scripts/data/plast_temptr_ffloops.csv", sep = ",")[,2:5]
plast_juice <- read.csv("scripts/data/plast_juice_ffloops.csv", sep = ",")[,2:5]
plast_stress <- read.csv("scripts/data/plast_stress_ffloops.csv", sep = ",")[,2:5]
plast_ox <- read.csv("scripts/data/plast_ox_ffloops.csv", sep = ",")[,2:5]

df <- rbind(colSums(all_genes)/nrow(all_genes),
            colSums(non_plast)/nrow(non_plast),
            colSums(all_plast)/nrow(all_plast),
            colSums(np_medgrowthloops)/nrow(np_medgrowthloops),
            colSums(plast_medgrowthloops)/nrow(plast_medgrowthloops),
            colSums(plast_ph)/nrow(plast_ph),
            colSums(plast_aero)/nrow(plast_aero),
            colSums(plast_mg_c)/nrow(plast_mg_c),
            colSums(plast_stringent)/nrow(plast_stringent),
            colSums(plast_tempr)/nrow(plast_tempr),
            colSums(plast_juice)/nrow(plast_juice),
            colSums(plast_stress)/nrow(plast_stress),
            colSums(plast_ox)/nrow(plast_ox)
            )
rownames(df) <- c("All","NP", "All_plast","NP_Medgth", "Medgth", "pH","Aero","Mg_C","Stringent","T°","Juice","Stress","Ox")

layout(matrix(c(1:1), 1, 1, byrow = TRUE))

barplot(t(df[,1:3]), col=c(7,3,"grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c(7,3,"grey"),
       legend=c("Coherent FFl", "Incoherent FFL","No FFL"))


pdf("figures/FFL_distrib_e_coli.pdf", width=10, height=4)
layout(matrix(c(1), 1, 1, byrow = TRUE))
par(mar=c(2, 2, 2, 2), mgp = c(1.75, 0.75, 0), las=0)
barplot(t(df[,1:3]), col=c(7,3,"grey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c(7,3,"grey"),
       legend=c("Coherent FFL", "Incoherent FFL", "No FFL"))
dev.off()
