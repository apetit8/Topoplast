#Figures of E_coli Loop analyses
source("scripts/functions/functions.R")

##FFL
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")
plast_medgrowthloops <- read.csv("scripts/data/plast_medium_E_coli_FFL.csv", sep = ",")
plast_ph <- read.csv("scripts/data/ph_plast_E_coli_FFL.csv", sep = ",")
plast_mg_c <- read.csv("scripts/data/mg_c_growth_E_coli_FFL.csv", sep = ",")
plast_aero <- read.csv("scripts/data/plast_aero_E_coli_FFL.csv", sep = ",")
plast_stringent <- read.csv("scripts/data/plast_stringent_E_coli_FFL.csv", sep = ",")
plast_tempr <- read.csv("scripts/data/plast_temptr_E_coli_FFL.csv", sep = ",")
plast_juice <- read.csv("scripts/data/plast_juice_E_coli_FFL.csv", sep = ",")
plast_stress <- read.csv("scripts/data/plast_stress_E_coli_FFL.csv", sep = ",")
plast_ox <- read.csv("scripts/data/plast_ox_E_coli_FFL.csv", sep = ",")

df1 <- as.data.frame(rbind(colSums(rbind(non_plast,all_plast)[,3:12])*100/(nrow(rbind(non_plast,all_plast))),
            colSums(non_plast[,3:12])*100/nrow(non_plast),
            colSums(all_plast[,3:12])*100/nrow(all_plast),
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
rownames(df1) <- c("All","Non-plastic\ngenes", "Plastic\ngenes","Medgth", "pH","Aero","Mg_C","Stringent","T°","Juice","Stress","Ox")




##DIAMONDs
non_plast <- read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")
plast_medgrowthloops <- read.csv("scripts/data/plast_medium_E_coli_diamond.csv", sep = ",")
plast_ph <- read.csv("scripts/data/ph_plast_E_coli_diamond.csv", sep = ",")
plast_mg_c <- read.csv("scripts/data/mg_c_growth_E_coli_diamond.csv", sep = ",")
plast_aero <- read.csv("scripts/data/plast_aero_E_coli_diamond.csv", sep = ",")
plast_stringent <- read.csv("scripts/data/plast_stringent_E_coli_diamond.csv", sep = ",")
plast_tempr <- read.csv("scripts/data/plast_temptr_E_coli_diamond.csv", sep = ",")
plast_juice <- read.csv("scripts/data/plast_juice_E_coli_diamond.csv", sep = ",")
plast_stress <- read.csv("scripts/data/plast_stress_E_coli_diamond.csv", sep = ",")
plast_ox <- read.csv("scripts/data/plast_ox_E_coli_diamond.csv", sep = ",")

df2 <- as.data.frame(rbind(colSums(rbind(non_plast,all_plast)[,3:14])*100/(nrow(rbind(non_plast,all_plast))),
                           colSums(non_plast[,3:14])*100/nrow(non_plast),
                           colSums(all_plast[,3:14])*100/nrow(all_plast),
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
rownames(df2) <- c("All","Non-plastic\ngenes", "Plastic\ngenes","Medgth", "pH","Aero","Mg_C","Stringent","T°","Juice","Stress","Ox")


barplot(t(df2[,2:12]), col=c("grey","olivedrab1","palegreen","mediumseagreen","orchid1","darkorchid1","plum1","lightsalmon1","indianred1","darkgoldenrod1","peachpuff"))





pdf("figures/Suppfig_3_type_FFL.pdf", width=14, height=10)
par(mar = c(8,2, 2,1), mfrow=c(2,1))
barplot(t(df1[,2:10]), col=c("grey","forestgreen","yellowgreen","dodgerblue3","lightskyblue","hotpink2","lightpink","orange","lightgoldenrod1"),
        legend.text = c("No FFL",colnames(df1[,3:10])), args.legend = list(ncol=4, x = "topright", inset = c(0.3, 1.22)),
        main = "Feedforward loop")

barplot(t(df2[,2:12]), col=c("grey","olivedrab1","palegreen","mediumseagreen","orchid1","darkorchid1","plum1","lightsalmon1","indianred1","darkgoldenrod1","peachpuff"),
        legend.text = c("No Diamond",colnames(df2[,3:12])), args.legend = list(ncol=4, x = "topright", inset = c(0.27, 1.22)),
        main = "Diamond loop")
dev.off()

