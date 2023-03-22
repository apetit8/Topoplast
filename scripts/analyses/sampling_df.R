#Sampling Analyses
source("functions/functions.R")
source("functions/Analysis_Functions.R")
#####################
sims.dirs <- list.dirs("../simul/GenFluct", full.names = TRUE, recursive = TRUE)
sims.dirs <- c(grep("sampling", sims.dirs, value = TRUE))

#####################

df.smplg <- df.sampling(sims.dirs, size = 2000)

ggplot(df.smplg, aes(x = Before_Gen, y = Slope))+
  geom_point()+
  stat_smooth()+ theme_bw()

#
df <- data.frame()
for( pop in unique(df.smplg$Before_Pop)){
  var <- var(subset(df.smplg, Before_Pop == pop)$Slope)
  df <- rbind( df, c(unique(subset(df.smplg, Before_Pop == pop)$Before_Gen), var) )
}
setnames(df, 1:2, c("Before_Gen","Var_Slope"))

g1 <- ggplot(df, aes(x=Before_Gen, y=sqrt(Var_Slope)))+
  geom_point()+theme_bw()+
  geom_smooth(span=0, alpha=0)



################################################################################





cairo_pdf("../figures/fig6_repro.pdf", width=14, height=14)
grid.arrange(
  g2, g1,
  ncol = 1,
  nrow = 2,
  clip = FALSE
)
dev.off()
















