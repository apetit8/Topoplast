#Robustness analyses
source("functions/functions.R")
#####################
sims.dirs <- list.dirs("../simul/Poplastlab", recursive = TRUE)
sims.dirs <- grep(sims.dirs, pattern="sampling",  invert=TRUE, value=TRUE)
df.plastlab <- df.simul(sims.dirs, M=TRUE, all.gen = TRUE, size=600000)
df.plastlab$treatment <- str_split((str_split(df.plastlab$data.dir, "../simul/Poplastlab/", n=2, simplify = TRUE)[,2]), "/", n=2, simplify = TRUE)[,1]
df.plastlab$anc <- str_split(df.plastlab$data.dir, "/", n=6, simplify = TRUE)[,5]
#
df.slope1 <- df.plasticity(df.plastlab, treshold = 0.1, anccol=2609, envircol=2608, bygenes = TRUE)
df.slope1$envir <- df.plastlab$treatment
#####################
df.poplast <- df.simul(list.dirs("../simul/Poplast", recursive = TRUE), M=FALSE, all.gen = TRUE)
df.poplast$treatment <- str_split((str_split(df.poplast$data.dir, "../simul/Popoplast/", n=2, simplify = TRUE)[,2]), "/", n=2, simplify = TRUE)[,1]
df.poplast$anc <- str_split(df.poplast$data.dir, "/", n=6, simplify = TRUE)[,5]
#
df.slope2 <- df.plasticity(df.poplast, treshold = 0.1, anccol=2609, envircol=2608, bygenes = TRUE)
df.slope2$envir <- df.poplast$treatment
#####################


gpp1 <- ggplot(df.slope2, aes(y = Slope_ss, x = Gen, col=envir))+
  geom_point(alpha=0.2)+ geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpp1.1 <- ggplot(df.slope2, aes(y = Slope_ss_new, x = Gen, col=envir))+
  geom_point(alpha=0.2)+ geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpp2 <- ggplot(df.slope2, aes(y = Slope_tf, x = Gen, col=envir))+
  geom_point(alpha=0.2)+ geom_smooth(span=0.1, alpha=0)+  ylim(c(0,1)) +theme_bw()

gpp3 <- ggplot(df.slope2, aes(y = Slope_pp, x = Gen, col=envir))+
  geom_point(alpha=0.2)+ geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpp4 <- ggplot(df.slope2, aes(y = Slope_pn, x = Gen, col=envir))+
  geom_point(alpha=0.2)+ geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()


gpl1 <- ggplot(df.slope1, aes(y = Slope_ss, x = Gen, col=envir))+
  geom_point(alpha=0.2)+ geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpl1.1 <- ggplot(df.slope1, aes(y = Slope_ss_new, x = Gen, col=envir))+
  geom_point(alpha=0.2)+ geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpl2 <- ggplot(df.slope1, aes(y = Slope_tf, x = Gen, col=envir))+
  geom_point(alpha=0.2)+ geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpl3 <- ggplot(df.slope1, aes(y = Slope_pp, x = Gen, col=envir))+
  geom_point(alpha=0.2)+ geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpl4 <- ggplot(df.slope1, aes(y = Slope_pn, x = Gen, col=envir))+
  geom_point(alpha=0.2)+ geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()


cairo_pdf("../figures/fig_slopes.pdf", width=15, height=20)
grid.arrange(
  gpp1,gpl1, gpp1.1,gpl1.1, gpp3,gpl3,gpp4,gpl4,gpp2,gpl2,
  ncol = 2,
  nrow = 5,
  widths = c(1,1),
  clip = FALSE
)
dev.off()


#All genes######################################################################

gpp1 <- ggplot(df.slope2, aes(x = Gen, y= Slope_ss, col=envir))+
  geom_point(aes(y = abs(g32)),alpha=0.1)+ geom_point(aes(y = abs(g33)),alpha=0.1)+geom_point(aes(y = abs(g34)),alpha=0.1)+geom_point(aes(y = abs(g35)),alpha=0.1)+geom_point(aes(y = abs(g36)),alpha=0.1)+
  geom_point(aes(y = abs(g37)),alpha=0.1)+ geom_point(aes(y = abs(g38)),alpha=0.1)+geom_point(aes(y = abs(g39)),alpha=0.1)+geom_point(aes(y = abs(g40)),alpha=0.1)+geom_point(aes(y = abs(g41)),alpha=0.1)+  
  geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpp2 <- ggplot(df.slope2, aes(x = Gen, y= Slope_ss_new, col=envir))+
  geom_point(aes(y = abs(g42)),alpha=0.1)+ geom_point(aes(y = abs(g43)),alpha=0.1)+geom_point(aes(y = abs(g44)),alpha=0.1)+geom_point(aes(y = abs(g45)),alpha=0.1)+geom_point(aes(y = abs(g46)),alpha=0.1)+
  geom_point(aes(y = abs(g47)),alpha=0.1)+ geom_point(aes(y = abs(g48)),alpha=0.1)+geom_point(aes(y = abs(g49)),alpha=0.1)+geom_point(aes(y = abs(g50)),alpha=0.1)+geom_point(aes(y = abs(g51)),alpha=0.1)+  
  geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpp3 <- ggplot(df.slope2, aes(x = Gen, y= Slope_pp, col=envir))+
  geom_point(aes(y = g22),alpha=0.1)+ geom_point(aes(y = g23),alpha=0.1)+geom_point(aes(y = g24),alpha=0.1)+geom_point(aes(y = g25),alpha=0.1)+geom_point(aes(y = g26),alpha=0.1)+
  geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpp4 <- ggplot(df.slope2, aes(x = Gen, y= Slope_pn, col=envir))+
  geom_point(aes(y = abs(g27)),alpha=0.1)+ geom_point(aes(y = abs(g28)),alpha=0.1)+geom_point(aes(y = abs(g29)),alpha=0.1)+geom_point(aes(y = abs(g30)),alpha=0.1)+geom_point(aes(y = abs(g31)),alpha=0.1)+
  geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpp5 <- ggplot(df.slope2, aes(x = Gen, y= Slope_tf, col=envir))+
  geom_point(aes(y = abs(g2)),alpha=0.1)+ geom_point(aes(y = abs(g3)),alpha=0.1)+geom_point(aes(y = abs(g4)),alpha=0.1)+geom_point(aes(y = abs(g5)),alpha=0.1)+geom_point(aes(y = abs(g6)),alpha=0.1)+
  geom_point(aes(y = abs(g7)),alpha=0.1)+ geom_point(aes(y = abs(g8)),alpha=0.1)+geom_point(aes(y = abs(g9)),alpha=0.1)+geom_point(aes(y = abs(g10)),alpha=0.1)+geom_point(aes(y = abs(g11)),alpha=0.1)+ 
  geom_point(aes(y = abs(g12)),alpha=0.1)+ geom_point(aes(y = abs(g13)),alpha=0.1)+geom_point(aes(y = abs(g14)),alpha=0.1)+geom_point(aes(y = abs(g15)),alpha=0.1)+geom_point(aes(y = abs(g16)),alpha=0.1)+
  geom_point(aes(y = abs(g17)),alpha=0.1)+ geom_point(aes(y = abs(g18)),alpha=0.1)+geom_point(aes(y = abs(g19)),alpha=0.1)+geom_point(aes(y = abs(g20)),alpha=0.1)+geom_point(aes(y = abs(g21)),alpha=0.1)+  
  geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()



gpl1 <- ggplot(df.slope1, aes(x = Gen, y= Slope_ss, col=envir))+
  geom_point(aes(y = abs(g32)),alpha=0.1)+ geom_point(aes(y = abs(g33)),alpha=0.1)+geom_point(aes(y = abs(g34)),alpha=0.1)+geom_point(aes(y = abs(g35)),alpha=0.1)+geom_point(aes(y = abs(g36)),alpha=0.1)+
  geom_point(aes(y = abs(g37)),alpha=0.1)+ geom_point(aes(y = abs(g38)),alpha=0.1)+geom_point(aes(y = abs(g39)),alpha=0.1)+geom_point(aes(y = abs(g40)),alpha=0.1)+geom_point(aes(y = abs(g41)),alpha=0.1)+  
  geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpl2 <- ggplot(df.slope1, aes(x = Gen, y= Slope_ss_new, col=envir))+
  geom_point(aes(y = abs(g42)),alpha=0.1)+ geom_point(aes(y = abs(g43)),alpha=0.1)+geom_point(aes(y = abs(g44)),alpha=0.1)+geom_point(aes(y = abs(g45)),alpha=0.1)+geom_point(aes(y = abs(g46)),alpha=0.1)+
  geom_point(aes(y = abs(g47)),alpha=0.1)+ geom_point(aes(y = abs(g48)),alpha=0.1)+geom_point(aes(y = abs(g49)),alpha=0.1)+geom_point(aes(y = abs(g50)),alpha=0.1)+geom_point(aes(y = abs(g51)),alpha=0.1)+  
  geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpl3 <- ggplot(df.slope1, aes(x = Gen, y= Slope_pp, col=envir))+
  geom_point(aes(y = g22),alpha=0.1)+ geom_point(aes(y = g23),alpha=0.1)+geom_point(aes(y = g24),alpha=0.1)+geom_point(aes(y = g25),alpha=0.1)+geom_point(aes(y = g26),alpha=0.1)+
  geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpl4 <- ggplot(df.slope1, aes(x = Gen, y= Slope_pn, col=envir))+
  geom_point(aes(y = abs(g27)),alpha=0.1)+ geom_point(aes(y = abs(g28)),alpha=0.1)+geom_point(aes(y = abs(g29)),alpha=0.1)+geom_point(aes(y = abs(g30)),alpha=0.1)+geom_point(aes(y = abs(g31)),alpha=0.1)+
  geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

gpl5 <- ggplot(df.slope1, aes(x = Gen, y= Slope_tf, col=envir))+
  geom_point(aes(y = abs(g2)),alpha=0.1)+ geom_point(aes(y = abs(g3)),alpha=0.1)+geom_point(aes(y = abs(g4)),alpha=0.1)+geom_point(aes(y = abs(g5)),alpha=0.1)+geom_point(aes(y = abs(g6)),alpha=0.1)+
  geom_point(aes(y = abs(g7)),alpha=0.1)+ geom_point(aes(y = abs(g8)),alpha=0.1)+geom_point(aes(y = abs(g9)),alpha=0.1)+geom_point(aes(y = abs(g10)),alpha=0.1)+geom_point(aes(y = abs(g11)),alpha=0.1)+ 
  geom_point(aes(y = abs(g12)),alpha=0.1)+ geom_point(aes(y = abs(g13)),alpha=0.1)+geom_point(aes(y = abs(g14)),alpha=0.1)+geom_point(aes(y = abs(g15)),alpha=0.1)+geom_point(aes(y = abs(g16)),alpha=0.1)+
  geom_point(aes(y = abs(g17)),alpha=0.1)+ geom_point(aes(y = abs(g18)),alpha=0.1)+geom_point(aes(y = abs(g19)),alpha=0.1)+geom_point(aes(y = abs(g20)),alpha=0.1)+geom_point(aes(y = abs(g21)),alpha=0.1)+  
  geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()



cairo_pdf("../figures/fig_slopes_allgenes.pdf", width=15, height=20)
grid.arrange(
  gpp1,gpl1, gpp2,gpl2, gpp3,gpl3,gpp4,gpl4,gpp5,gpl5,
  ncol = 2,
  nrow = 5,
  widths = c(1,1),
  clip = FALSE
)
dev.off()


#M size######################################################################

df.slope1$M_size_ss <- df.plastlab$M_size_ss
df.slope1$M_size_ssn <- df.plastlab$M_size_ssn
df.slope1$M_size_tf <- df.plastlab$M_size_tf
df.slope1$M_size_pn <- df.plastlab$M_size_pn
df.slope1$M_size_pp <- df.plastlab$M_size_pp

ggplot(df.slope1, aes(x = Msize_ss, y= Slope_ss, col=envir))+
  geom_point(aes(y = abs(g32)))+ geom_point(aes(y = abs(g33)))+geom_point(aes(y = abs(g34)))+geom_point(aes(y = abs(g35)))+geom_point(aes(y = abs(g36)))+
  geom_point(aes(y = abs(g37)))+ geom_point(aes(y = abs(g38)))+geom_point(aes(y = abs(g39)))+geom_point(aes(y = abs(g40)))+geom_point(aes(y = abs(g41)))+  
  geom_smooth(span=0.1, alpha=0)+ ylim(c(0,1)) +theme_bw()

#####
sims.dirs <- list.dirs("../simul/Poplastlab", recursive = TRUE)
sims.dirs <- grep(sims.dirs, pattern="sampling",  invert=TRUE, value=TRUE)
df.plastlab_genes <- df.genes(sims.dirs, all.gen = TRUE, size=600000)

g1 <- ggplot(df.plastlab_genes, aes(x = M_var, y= abs(Slope), col=Treatment))+
  geom_point()+
  geom_smooth(span=0.1, alpha=0)+ 
  theme_bw()+facet_wrap(Gene_cat ~.,  ncol=3)

ggplot(df.plastlab_genes, aes(x = M_var, y= abs(Slope), col=Gen))+
  geom_point()+
  geom_smooth(span=0.1, alpha=0)+ 
  theme_bw()+facet_wrap(Gene_cat ~.,  ncol=3)

g2 <- ggplot(df.plastlab_genes, aes(x = M_var, y= abs(Slope), col=Gene_cat))+
  geom_point()+
  geom_smooth(span=0.1, alpha=0)+ 
  theme_bw()+facet_wrap(Treatment ~.,  ncol=3)


cor.test(df.plastlab_genes$M_var, abs(df.plastlab_genes$Slope))$estimate

#Do subset to see if it differss between the environments ?


cairo_pdf("../figures/fig_robustess.pdf", width=15, height=20)
grid.arrange(
  g1,g2,
  ncol = 1,
  nrow = 2,
  widths = c(1),
  clip = FALSE
)
dev.off()



sims.dirs <- list.dirs("../simul/Poplastlab", recursive = TRUE)
sims.dirs <- grep(sims.dirs, pattern="sampling",  invert=TRUE, value=TRUE)
df.plastlab_genes2 <- df.genes(sims.dirs, all.gen = FALSE, size=600000)

gg <- ggplot(df.plastlab_genes, aes(y = M_var, x= as.numeric(Gen), col=Gene_cat))+
  geom_point()+
  geom_smooth(span=0.1, alpha=0.1)+ 
  theme_bw()+facet_wrap(Treatment  ~.,  ncol=3)

cairo_pdf("../figures/fig_mvar_by_gen.pdf", width=12, height=10)
grid.arrange(
  gg,
  ncol = 1,
  nrow = 1,
  widths = c(1),
  clip = FALSE
)
dev.off()

g1 <- ggplot(df.plastlab_genes2, aes(x = M_var, y= abs(Slope), col=Treatment))+
  geom_point()+
  geom_smooth(span=0.1, alpha=0)+ 
  theme_bw()+facet_wrap(Gene_cat ~.,  ncol=3)


g2 <- ggplot(df.plastlab_genes2, aes(x = M_var, y= abs(Slope), col=Gene_cat))+
  geom_point()+
  geom_smooth(span=0.1, alpha=0)+ 
  theme_bw()+facet_wrap(Treatment ~.,  ncol=3)


subset(df.plastlab_genes, Gene_cat=="tf")

cor.test(subset(df.plastlab_genes, Gene_cat=="tf")$M_var, abs(subset(df.plastlab_genes, Gene_cat=="tf")$Slope))$estimate
cor.test(subset(df.plastlab_genes, Gene_cat=="ssn")$M_var, abs(subset(df.plastlab_genes, Gene_cat=="ssn")$Slope))$estimate
cor.test(subset(df.plastlab_genes, Gene_cat=="ss")$M_var, abs(subset(df.plastlab_genes, Gene_cat=="ss")$Slope))$estimate
cor.test(subset(df.plastlab_genes, Gene_cat=="pp")$M_var, abs(subset(df.plastlab_genes, Gene_cat=="pp")$Slope))$estimate
cor.test(subset(df.plastlab_genes, Gene_cat=="pn")$M_var, abs(subset(df.plastlab_genes, Gene_cat=="pn")$Slope))$estimate



cairo_pdf("../figures/fig_robustess2.pdf", width=15, height=20)
grid.arrange(
  g1,g2,
  ncol = 1,
  nrow = 2,
  widths = c(1),
  clip = FALSE
)
dev.off()
