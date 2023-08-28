library(ggvenn)

x <- list(
  FFL = c(subset(read.csv("scripts/data/plast_genes_E_coli_nffl.csv", sep = ","), Loop_number!=0)[,2], subset(read.csv("scripts/data/nonplast_E_coli_nffl.csv", sep = ","), Loop_number!=0)[,2]),
  DMD = c(subset(read.csv("scripts/data/plast_genes_E_coli_nDMD.csv", sep = ","), Loop_number!=0)[,2], subset(read.csv("scripts/data/nonplast_E_coli_nDMD.csv", sep = ","), Loop_number!=0)[,2]), 
  FBL = c(subset(read.csv("scripts/data/plast_genes_E_coli_nFBL.csv", sep = ","), FBL_number!=0)[,2], subset(read.csv("scripts/data/nonplast_E_coli_nFBL.csv", sep = ","), FBL_number!=0)[,2])
)

pdf(paste0("figures/Loop_Venn",".pdf"), width=5, height=4)
par(mgp=c(2.5, 1.2, 0), mar = c(2.9,3.5, 0.1,0.1))
ggvenn(
  x, fill_alpha = 0.7, stroke_color = "white",
  fill_color = c("khaki1", "darkseagreen2", "darkolivegreen2"),
  stroke_size = 0.5, set_name_size = 6, text_size = 4.5
)
dev.off()