source("scripts/functions/functions.R")
#Plot combining both empirical data and theoretical data
ifelse(!dir.exists("figures/short"), dir.create("figures/short"), FALSE)
################################################################################
pval <- 0.05
################################################################################
for (gen in c(1000)){
  pdfname <- paste0("figures/short/fig_",gen)
    ##FFL motifs####################################################################
    non_plast <- read.csv("scripts/data/nonplast_E_coli_FFL.csv", sep = ",")
    all_plast <- read.csv("scripts/data/plast_genes_E_coli_FFL.csv", sep = ",")
    theory <- read.csv(paste0("scripts/data/short/full_1k_",gen,"_FFL.csv"), sep = ",")[,1:11]
    
    df_ffl <- as.data.frame(rbind(colSums(non_plast[,5:12])*100/sum(non_plast[,3]), colSums(all_plast[,5:12])*100/sum(all_plast[,3]),
                                  theory[2,4:11]*100/theory[2,2],theory[1,4:11]*100/theory[1,2]))
    rownames(df_ffl) <- c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory")
    
    #Stat test to know if plastic and non plastic are different from each other
    #empiric
    df_anova <- rbind(subset(non_plast, FFL==1)[,5:12], subset(all_plast, FFL==1)[,5:12])
    df_anova$Plasticity <- c(rep("No", nrow(subset(non_plast, FFL==1))),rep("Yes", nrow(subset(all_plast, FFL==1))))
    
    #colnames(df_ffl)
    # t.test(subset(non_plast, FFL==1)$I2,  subset(all_plast, FFL==1)$I2)$p.value
    legendtext <- c(paste0("C3 ", ifelse(anova(lm(C3 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("C1 ", ifelse(anova(lm(C1 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("C2 ", ifelse(anova(lm(C2 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("C4 ", ifelse(anova(lm(C4 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("I2 ", ifelse(anova(lm(I2 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("I4 ", ifelse(anova(lm(I4 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("I3 ", ifelse(anova(lm(I3 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("I1 ", ifelse(anova(lm(I1 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")))
    
    pdf(paste0(pdfname,"_FFL_motifs",".pdf"), width=5, height=6)
    par(mar = c(8,2, 2,1))
    #Each motif topology
    barplot(t(df_ffl), main="FFL motif distribution", col=c("forestgreen","yellowgreen","dodgerblue3","lightskyblue","hotpink2","lightpink","orange","lightgoldenrod1"),
            space=c(0.3,0.1,0.4,0.1), legend.text = legendtext, args.legend = list(ncol=4, x = "topright", inset = c(0.03, 1.15)))
    dev.off()
    
    ####################################
    ###Same plot but for random empirical groups
    # df <- (rbind(non_plast, all_plast))
    # 
    # df1 <- df[sample(nrow(df), nrow(all_plast)), ]
    # df2 <- subset(df, !(V1 %in% df1$V1))
    # 
    # dfrandom <- as.data.frame(rbind(colSums(df1[,5:12])*100/nrow(df1), colSums(df2[,5:12])*100/nrow(df2) ))
    # rownames(dfrandom) <- c("Random 1\nEmpiric","Random2\nEmpiric")
    # 
    # pdf(paste0(pdfname,"_random_FFL_motifs",".pdf"), width=5, height=6)
    # par(mar = c(8,2, 2,1))
    # #Each motif topology
    # barplot(t(dfrandom), main="FFL motif distribution", col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
    #         space=c(0.3,0.1), legend.text = colnames(dfrandom), args.legend = list(ncol=4, x = "topright", inset = c(0.2, 1.2)))
    # dev.off()
    
    
    ################################################################################
    ##DMD motifs####################################################################
    non_plast <- read.csv("scripts/data/nonplast_E_coli_diamond.csv", sep = ",")
    all_plast <- read.csv("scripts/data/plast_genes_E_coli_diamond.csv", sep = ",")
    theory <- read.csv(paste0("scripts/data/short/full_1k_",gen,"_DMD.csv"), sep = ",")[,1:13]
    
    df1 <- as.data.frame(rbind(colSums(non_plast[,5:14])*100/sum(non_plast[,3]), colSums(all_plast[,5:14])*100/sum(all_plast[,3]),
                               theory[2,4:13]*100/theory[2,2], theory[1,4:13]*100/theory[1,2] ))
    rownames(df1) <- c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory")
    
    #Stat test to know if plastic and non plastic are different from each other
    df_anova <- rbind(subset(non_plast, FFL==1)[,5:14], subset(all_plast, FFL==1)[,5:14])
    df_anova$Plasticity <- c(rep("No", nrow(subset(non_plast, FFL==1))),rep("Yes", nrow(subset(all_plast, FFL==1))))
    
    # colnames(df1)
    legendtext <- c(paste0("PP", ifelse(anova(lm(PP ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("PM", ifelse(anova(lm(PM ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("PN", ifelse(anova(lm(PN ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("NP", ifelse(anova(lm(NP ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("NM", ifelse(anova(lm(NM ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("NN", ifelse(anova(lm(NN ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("MP", ifelse(anova(lm(MP ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("MM2", ifelse(anova(lm(MM2 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("MM1", ifelse(anova(lm(MM1 ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")),
                    paste0("MN", ifelse(anova(lm(MN ~ Plasticity, df_anova))[1,5] <= pval, "(*)", "")))
    
    pdf(paste0(pdfname,"_DMD_motifs",".pdf"), width=5, height=6)
    par(mar = c(8,2, 2,1))
    #Each motif topology
    barplot(t(df1[1:10]), main="Diamond motif distribution", col=c("olivedrab1","palegreen","mediumseagreen","orchid1","darkorchid1","plum1","lightsalmon1","indianred1","darkgoldenrod1","peachpuff"),
            space=c(0.3,0.1,0.4,0.1), legend.text = legendtext, args.legend = list(ncol=4, x = "topright", inset = c(0.01, 1.13)))
    dev.off()
    
    ################################################################################
    ##FBL motifs####################################################################
    non_plast <- read.csv("scripts/data/nonplast_E_coli_FBL.csv", sep = ",")
    all_plast <- read.csv("scripts/data/plast_genes_E_coli_FBL.csv", sep = ",")
    theory <- read.csv(paste0("scripts/data/short/full_1k_",gen,"_FBL.csv"), sep = ",")[,]
    
    df <- as.data.frame(rbind(colSums(non_plast[,7:11])*100/sum(non_plast[,3]), colSums(all_plast[,7:11])*100/sum(all_plast[,3]), 
                              theory[2,6:10]*100/theory[2,2], theory[1,6:10]*100/theory[1,2]))
    rownames(df) <- c("Non-Plastic\nE. coli", "Plastic\nE. coli", "Non-plastic\nTheory","Plastic\nTheory")
    
    
    pdf(paste0(pdfname,"_FBL",".pdf"), width=5, height=5)
    par(mar = c(7,2, 2,1))
    #Each motif topology
    barplot(t(df), main="FBL motifs", col=c("olivedrab1","palegreen","mediumseagreen","plum","darkorchid1"),
            space=c(0.3,0.1,0.4,0.1), legend.text = colnames(df), 
            args.legend = list(ncol=2, x = "topright", inset = c(0.3, 1.2)))
    dev.off()
    
    ################################################################################
    ##Regulations###################################################################
    #Simulations
    topo_plastic <- readRDS(paste0("scripts/data/short/list_plastic_topo_full_1k_",gen,".Rds"))
    topo_nnplast <- readRDS(paste0("scripts/data/short/list_nnplast_topo_full_1k_",gen,".Rds"))
    
    THreg_sum_plast <- unlist(lapply(topo_plastic, function(i){rowSums2(abs(i))[2]} )) # [2] : the target of the network is in position #2
    THreg_sum_np <- unlist(lapply(topo_nnplast, function(i){rowSums2(abs(i))[2]} ))
    #
    length(THreg_sum_plast) <- max(length(THreg_sum_np), length(THreg_sum_plast))
    length(THreg_sum_np) <- max(length(THreg_sum_np), length(THreg_sum_plast))
    #
    pdf(paste0(pdfname,"_regnumber",".pdf"), width=5, height=4)
    par(mgp=c(2.5, 1.2, 0), mar = c(2.9,3.5, 0.1,0.1))
    boxplot(cbind(as.array(THreg_sum_np), as.array(THreg_sum_plast)), ylab = "Number of regulation per gene",
            names = c("Non-plastic\nTheory","Plastic\nTheory"))
    dev.off()
    ######################################
    #Simulations
    THreg_sum_pos_plast <- unlist(lapply(topo_plastic, function(i){ i[(i == -1)] <- 0
    return(rowSums2(abs(i))[2])} )) # [2] : the target of the network is in position #2
    THreg_sum_neg_plast <- unlist(lapply(topo_plastic, function(i){ i[(i == 1)] <- 0
    return(rowSums2(abs(i))[2])} ))
    #
    THreg_sum_pos_np <- unlist(lapply(topo_nnplast, function(i){ i[(i == -1)] <- 0
    return(rowSums2(abs(i))[2])} ))
    THreg_sum_neg_np <- unlist(lapply(topo_nnplast, function(i){ i[(i == 1)] <- 0
    return(rowSums2(abs(i))[2])} ))
    
    df <- data.frame(Pos=c(sum(THreg_sum_pos_np, na.rm=TRUE)*100/sum(THreg_sum_np, na.rm=TRUE), sum(THreg_sum_pos_plast, na.rm=TRUE)*100/sum(THreg_sum_plast, na.rm=TRUE) ),
                     Neg=c(sum(THreg_sum_neg_np, na.rm=TRUE)*100/sum(THreg_sum_np, na.rm=TRUE), sum(THreg_sum_neg_plast, na.rm=TRUE)*100/sum(THreg_sum_plast, na.rm=TRUE) ))
    row.names(df) <- c("Non-plastic\nTheory","Plastic\nTheory")
    
    pdf(paste0(pdfname,"_sign",".pdf"), width=5, height=4)
    par(mar = c(5,2, 2,1))
    barplot(t(df), main = "Regulation type", 
            legend.text = c( "Activating","Inhibiting"), args.legend = list(ncol=2, x = "topright", inset = c(0.2, 1.2)), ylim = c(0,105))
    dev.off()
}

