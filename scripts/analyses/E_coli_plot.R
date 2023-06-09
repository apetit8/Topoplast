#Figures of E_coli Loop analyses
source("scripts/functions/functions.R")

##FFL
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")
np_medgrowthloops <- read.csv("scripts/data/np_medium_growth_E_coli_FFL.csv", sep = ",")
plast_medgrowthloops <- read.csv("scripts/data/plast_medium_growth_E_coli_FFL.csv", sep = ",")
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
rownames(df1) <- c("All","Non-plastic\ngenes", "Plastic\ngenes","NP_Medgth", "Medgth", "pH","Aero","Mg_C","Stringent","T°","Juice","Stress","Ox")
df1$coherent <- rowSums2(as.matrix(df1[,c(3:6)])) #Count of 
df1$incoherent <- rowSums2(as.matrix(df1[,c(7:10)]))
df1$homogenous <- rowSums2(as.matrix(df1[,c(3,4,7,8)]))
df1$heterogenous <- rowSums2(as.matrix(df1[,c(5,6,9,10)]))

pdfname <- "figures/E_coli"
pdf(paste0(pdfname,"_type_FFL",".pdf"), width=14, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df1[,2:10]), col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
       legend=c( "No_FFL","C_Ho_neg","C_Ho_pos","C_He_neg","C_He_pos","I_Ho_neg","I_Ho_pos","I_He_neg","I_He_pos") )
#Coherence
barplot(t(df1[,c(2,11,12)]), col=c("grey","indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df1[,c(2,13,14)]), col=c("grey","orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","orange","yellowgreen"),
       legend=c( "No_FFL","Homogenous","Heterogenous") )
dev.off()

#####Motif prop
df <- as.data.frame(rbind(colSums(rbind(non_plast,all_plast)[,3:12])*100/(sum(rbind(non_plast,all_plast)[,3])),
                          colSums(non_plast[,3:12])*100/sum(non_plast[,3]),
                          colSums(all_plast[,3:12])*100/sum(all_plast[,3]),
                          colSums(np_medgrowthloops[,3:12])*100/sum(np_medgrowthloops[,3]),
                          colSums(plast_medgrowthloops[,3:12])*100/sum(plast_medgrowthloops[,3]),
                          colSums(plast_ph[,3:12])*100/sum(plast_ph[,3]),
                          colSums(plast_aero[,3:12])*100/sum(plast_aero[,3]),
                          colSums(plast_mg_c[,3:12])*100/sum(plast_mg_c[,3]),
                          colSums(plast_stringent[,3:12])*100/sum(plast_stringent[,3]),
                          colSums(plast_tempr[,3:12])*100/sum(plast_tempr[,3]),
                          colSums(plast_juice[,3:12])*100/sum(plast_juice[,3]),
                          colSums(plast_stress[,3:12])*100/sum(plast_stress[,3]),
                          colSums(plast_ox[,3:12])*100/sum(plast_ox[,3])
))
rownames(df) <- c("All","NP", "All_plast","NP_Medgth", "Medgth", "pH","Aero","Mg_C","Stringent","T°","Juice","Stress","Ox")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:10)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,10)]))


pdfname <- "figures/E_coli"
pdf(paste0(pdfname,"_percent_FFL",".pdf"), width=14, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,3:10]), col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
       legend=c( "No_FFL","C_Ho_neg","C_Ho_pos","C_He_neg","C_He_pos","I_Ho_neg","I_Ho_pos","I_He_neg","I_He_pos") )
#Coherence
barplot(t(df[,c(11,12)]), col=c("indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(13,14)]), col=c("orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("orange","yellowgreen"),
       legend=c( "No_FFL","Homogenous","Heterogenous") )
dev.off()

################################################################################
##FFL with the X gene being plastic as well
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL_from.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL_from.csv", sep = ",")
np_medgrowthloops <- read.csv("scripts/data/np_medium_growth_E_coli_FFL_from.csv", sep = ",")
plast_medgrowthloops <- read.csv("scripts/data/plast_medium_E_coli_FFL_from.csv", sep = ",")
plast_ph <- read.csv("scripts/data/ph_plast_E_coli_FFL_from.csv", sep = ",")
plast_mg_c <- read.csv("scripts/data/mg_c_growth_E_coli_FFL_from.csv", sep = ",")
plast_aero <- read.csv("scripts/data/plast_aero_E_coli_FFL_from.csv", sep = ",")
plast_stringent <- read.csv("scripts/data/plast_stringent_E_coli_FFL_from.csv", sep = ",")
plast_tempr <- read.csv("scripts/data/plast_temptr_E_coli_FFL_from.csv", sep = ",")
plast_juice <- read.csv("scripts/data/plast_juice_E_coli_FFL_from.csv", sep = ",")
plast_stress <- read.csv("scripts/data/plast_stress_E_coli_FFL_from.csv", sep = ",")
plast_ox <- read.csv("scripts/data/plast_ox_E_coli_FFL_from.csv", sep = ",")

df2 <- as.data.frame(rbind(colSums(rbind(non_plast,all_plast)[,3:13])*100/(nrow(rbind(non_plast,all_plast))),
                          colSums(non_plast[,3:13])*100/nrow(non_plast),
                          colSums(all_plast[,3:13])*100/nrow(all_plast),
                          colSums(np_medgrowthloops[,3:13])*100/nrow(np_medgrowthloops),
                          colSums(plast_medgrowthloops[,3:13])*100/nrow(plast_medgrowthloops),
                          colSums(plast_ph[,3:13])*100/nrow(plast_ph),
                          colSums(plast_aero[,3:13])*100/nrow(plast_aero),
                          colSums(plast_mg_c[,3:13])*100/nrow(plast_mg_c),
                          colSums(plast_stringent[,3:13])*100/nrow(plast_stringent),
                          colSums(plast_tempr[,3:13])*100/nrow(plast_tempr),
                          colSums(plast_juice[,3:13])*100/nrow(plast_juice),
                          colSums(plast_stress[,3:13])*100/nrow(plast_stress),
                          colSums(plast_ox[,3:13])*100/nrow(plast_ox)
))
rownames(df2) <- c("All","NP", "All_plast","NP_Medgth", "Medgth", "pH","Aero","Mg_C","Stringent","T°","Juice","Stress","Ox")
df2$coherent <- rowSums2(as.matrix(df2[,c(3:6)])) #Count of 
df2$incoherent <- rowSums2(as.matrix(df2[,c(7:10)]))
df2$homogenous <- rowSums2(as.matrix(df2[,c(3,4,7,8)]))
df2$heterogenous <- rowSums2(as.matrix(df2[,c(5,6,9,10)]))


pdfname <- "figures/E_coli"
pdf(paste0(pdfname,"_FFL_from",".pdf"), width=14, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df2[,2:11]), col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"),
       legend=c( "No_FFL","C_Ho_neg","C_Ho_pos","C_He_neg","C_He_pos","I_Ho_neg","I_Ho_pos","I_He_neg","I_He_pos","FFL_from_NP") )
#Coherence
barplot(t(df2[,c(2,12,13)]), col=c("grey","indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df2[,c(2,14,15)]), col=c("grey","orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","orange","yellowgreen"),
       legend=c( "No_FFL","Homogenous","Heterogenous") )
dev.off()

##
#FFL FROM percent
df <- as.data.frame(rbind(colSums(rbind(non_plast,all_plast)[,3:13])*100/(sum(rbind(non_plast,all_plast)[,3])-sum(rbind(non_plast,all_plast)[,13])),
                          colSums(non_plast[,3:13])*100/(sum(non_plast[,3])-sum(non_plast[,13])),
                          colSums(all_plast[,3:13])*100/(sum(all_plast[,3])-sum(all_plast[,13])),
                          colSums(np_medgrowthloops[,3:13])*100/(sum(np_medgrowthloops[,3])-sum(np_medgrowthloops[,13])),
                          colSums(plast_medgrowthloops[,3:13])*100/(sum(plast_medgrowthloops[,3])-sum(plast_medgrowthloops[,13])),
                          colSums(plast_ph[,3:13])*100/(sum(plast_ph[,3])-sum(plast_ph[,13])),
                          colSums(plast_aero[,3:13])*100/(sum(plast_aero[,3])-sum(plast_aero[,13])),
                          colSums(plast_mg_c[,3:13])*100/(sum(plast_mg_c[,3])-sum(plast_mg_c[,13])),
                          colSums(plast_stringent[,3:13])*100/(sum(plast_stringent[,3])-sum(plast_stringent[,13])),
                          colSums(plast_tempr[,3:13])*100/(sum(plast_tempr[,3])-sum(plast_tempr[,13])),
                          colSums(plast_juice[,3:13])*100/(sum(plast_juice[,3])-sum(plast_juice[,13])),
                          colSums(plast_stress[,3:13])*100/(sum(plast_stress[,3])-sum(plast_stress[,13])),
                          colSums(plast_ox[,3:13])*100/(sum(plast_ox[,3])-sum(plast_ox[,13]))
))
rownames(df) <- c("All","NP", "All_plast","NP_Medgth", "Medgth", "pH","Aero","Mg_C","Stringent","T°","Juice","Stress","Ox")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:10)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,10)]))


pdfname <- "figures/E_coli"
pdf(paste0(pdfname,"_percent_FFL_from",".pdf"), width=14, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,3:10]), col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"),
       legend=c( "C_Ho_neg","C_Ho_pos","C_He_neg","C_He_pos","I_Ho_neg","I_Ho_pos","I_He_neg","I_He_pos","NP_FFL") )
#Coherence
barplot(t(df[,c(12,13)]), col=c("indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("indianred1","dodgerblue"),
       legend=c("Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(14,15)]), col=c("orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("orange","yellowgreen"),
       legend=c("Homogenous","Heterogenous") )
dev.off()
##
#FFL FROM percent
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL_from_allplast.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL_from_allplast.csv", sep = ",")
np_medgrowthloops <- read.csv("scripts/data/np_medium_growth_E_coli_FFL_from_allplast.csv", sep = ",")
plast_medgrowthloops <- read.csv("scripts/data/plast_medium_E_coli_FFL_from_allplast.csv", sep = ",")
plast_ph <- read.csv("scripts/data/ph_plast_E_coli_FFL_from_allplast.csv", sep = ",")
plast_mg_c <- read.csv("scripts/data/mg_c_growth_E_coli_FFL_from_allplast.csv", sep = ",")
plast_aero <- read.csv("scripts/data/plast_aero_E_coli_FFL_from_allplast.csv", sep = ",")
plast_stringent <- read.csv("scripts/data/plast_stringent_E_coli_FFL_from_allplast.csv", sep = ",")
plast_tempr <- read.csv("scripts/data/plast_temptr_E_coli_FFL_from_allplast.csv", sep = ",")
plast_juice <- read.csv("scripts/data/plast_juice_E_coli_FFL_from_allplast.csv", sep = ",")
plast_stress <- read.csv("scripts/data/plast_stress_E_coli_FFL_from_allplast.csv", sep = ",")
plast_ox <- read.csv("scripts/data/plast_ox_E_coli_FFL_from_allplast.csv", sep = ",")

df <- as.data.frame(rbind(colSums(rbind(non_plast,all_plast)[,3:13]),
                          colSums(non_plast[,3:13]),
                          colSums(all_plast[,3:13]),
                          colSums(np_medgrowthloops[,3:13]),
                          colSums(plast_medgrowthloops[,3:13]),
                          colSums(plast_ph[,3:13]),
                          colSums(plast_aero[,3:13]),
                          colSums(plast_mg_c[,3:13]),
                          colSums(plast_stringent[,3:13]),
                          colSums(plast_tempr[,3:13]),
                          colSums(plast_juice[,3:13]),
                          colSums(plast_stress[,3:13]),
                          colSums(plast_ox[,3:13])
))
rownames(df) <- c("All","NP", "All_plast","NP_Medgth", "Medgth", "pH","Aero","Mg_C","Stringent","T°","Juice","Stress","Ox")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:10)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,10)]))


pdfname <- "figures/E_coli"
pdf(paste0(pdfname,"_percent_FFL_from_allplast",".pdf"), width=14, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,3:10]), col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"),
       legend=c( "C_Ho_neg","C_Ho_pos","C_He_neg","C_He_pos","I_Ho_neg","I_Ho_pos","I_He_neg","I_He_pos","NP_FFL") )
#Coherence
barplot(t(df[,c(12,13)]), col=c("indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("indianred1","dodgerblue"),
       legend=c("Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(14,15)]), col=c("orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("orange","yellowgreen"),
       legend=c("Homogenous","Heterogenous") )
dev.off()
##

#FFL FROM NP percent
################################################################################
##FFL with the X gene being plastic as well
non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL_from_NP.csv", sep = ",")
all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL_from_NP.csv", sep = ",")
np_medgrowthloops <- read.csv("scripts/data/np_medium_growth_E_coli_FFL_from_NP.csv", sep = ",")
plast_medgrowthloops <- read.csv("scripts/data/plast_medium_E_coli_FFL_from_NP.csv", sep = ",")
plast_ph <- read.csv("scripts/data/ph_plast_E_coli_FFL_from_NP.csv", sep = ",")
plast_mg_c <- read.csv("scripts/data/mg_c_growth_E_coli_FFL_from_NP.csv", sep = ",")
plast_aero <- read.csv("scripts/data/plast_aero_E_coli_FFL_from_NP.csv", sep = ",")
plast_stringent <- read.csv("scripts/data/plast_stringent_E_coli_FFL_from_NP.csv", sep = ",")
plast_tempr <- read.csv("scripts/data/plast_temptr_E_coli_FFL_from_NP.csv", sep = ",")
plast_juice <- read.csv("scripts/data/plast_juice_E_coli_FFL_from_NP.csv", sep = ",")
plast_stress <- read.csv("scripts/data/plast_stress_E_coli_FFL_from_NP.csv", sep = ",")
plast_ox <- read.csv("scripts/data/plast_ox_E_coli_FFL_from_NP.csv", sep = ",")

df <- as.data.frame(rbind(colSums(rbind(non_plast,all_plast)[,3:13])*100/(sum(rbind(non_plast,all_plast)[,3])-sum(rbind(non_plast,all_plast)[,13])),
                          colSums(non_plast[,3:13])*100/(sum(non_plast[,3])-sum(non_plast[,13])),
                          colSums(all_plast[,3:13])*100/(sum(all_plast[,3])-sum(all_plast[,13])),
                          colSums(np_medgrowthloops[,3:13])*100/(sum(np_medgrowthloops[,3])-sum(np_medgrowthloops[,13])),
                          colSums(plast_medgrowthloops[,3:13])*100/(sum(plast_medgrowthloops[,3])-sum(plast_medgrowthloops[,13])),
                          colSums(plast_ph[,3:13])*100/(sum(plast_ph[,3])-sum(plast_ph[,13])),
                          colSums(plast_aero[,3:13])*100/(sum(plast_aero[,3])-sum(plast_aero[,13])),
                          colSums(plast_mg_c[,3:13])*100/(sum(plast_mg_c[,3])-sum(plast_mg_c[,13])),
                          colSums(plast_stringent[,3:13])*100/(sum(plast_stringent[,3])-sum(plast_stringent[,13])),
                          colSums(plast_tempr[,3:13])*100/(sum(plast_tempr[,3])-sum(plast_tempr[,13])),
                          colSums(plast_juice[,3:13])*100/(sum(plast_juice[,3])-sum(plast_juice[,13])),
                          colSums(plast_stress[,3:13])*100/(sum(plast_stress[,3])-sum(plast_stress[,13])),
                          colSums(plast_ox[,3:13])*100/(sum(plast_ox[,3])-sum(plast_ox[,13]))
))
rownames(df) <- c("All","NP", "All_plast","NP_Medgth", "Medgth", "pH","Aero","Mg_C","Stringent","T°","Juice","Stress","Ox")
df$coherent <- rowSums2(as.matrix(df[,c(3:6)])) #Count of 
df$incoherent <- rowSums2(as.matrix(df[,c(7:10)]))
df$homogenous <- rowSums2(as.matrix(df[,c(3,4,7,8)]))
df$heterogenous <- rowSums2(as.matrix(df[,c(5,6,9,10)]))


pdfname <- "figures/E_coli"
pdf(paste0(pdfname,"_percent_FFL_from_NP",".pdf"), width=14, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,3:10]), col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","lightslategrey"),
       legend=c("C_Ho_neg","C_Ho_pos","C_He_neg","C_He_pos","I_Ho_neg","I_Ho_pos","I_He_neg","I_He_pos","NP_FFL") )
#Coherence
barplot(t(df[,c(12,13)]), col=c("indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("indianred1","dodgerblue"),
       legend=c("Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(14,15)]), col=c("orange","yellowgreen"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("orange","yellowgreen"),
       legend=c("Homogenous","Heterogenous") )
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

df <- as.data.frame(rbind(colSums(rbind(non_plast,all_plast)[,3:14])*100/(nrow(rbind(non_plast,all_plast))),
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
pdf(paste0(pdfname,"_type_diamond_from",".pdf"), width=14, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
#Each motif topology
barplot(t(df[,2:12]), col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","darkorchid1","plum1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1","darkorchid1","plum1"),
       legend=c( "No_FFL","Pos_pos","Pos_mixt","Pos_neg",
                 "Neg_pos","Neg_mixt","Neg_neg",
                 "Mixt_pos","Mixt_mixt_hom","Mixt_mixt_het","Mixt_neg") )
#Coherence
barplot(t(df[,c(2,13,14)]), col=c("grey","indianred1","dodgerblue"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","indianred1","dodgerblue"),
       legend=c( "No_FFL","Coherent","Incoherent") )
#Homogeneity
barplot(t(df[,c(2,15:17)]), col=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("grey","darkseagreen","yellowgreen","dodgerblue","deepskyblue3"),
       legend=c( "No_FFL","Homogenous_pos","Homogenous_neg","Heterogenous") )
dev.off()

