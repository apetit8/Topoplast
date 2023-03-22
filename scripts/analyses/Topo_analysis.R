#Analyses ?
source("functions/functions.R")
# source("functions/plotnet.R")
library(igraph)
library(ggstatsplot)
#####################
sims.dirs1 <- list.dirs("../simul/Evolpast/burning/Correlated", recursive = TRUE)
sims.dirs2 <- list.dirs("../simul/Evolpast/new_envir/B_Correlated", recursive = TRUE)

#####################
#DATA
##################
#Burning df
df.b <- df.simul(sims.dirs1, all.gen = TRUE)
df.b$n_envir <- str_split(df.b$data.dir, "/", n=8, simplify = TRUE)[,5]
df.b$b_envir <- str_split(df.b$data.dir, "/", n=8, simplify = TRUE)[,5]
df.b$anc_id <- str_split(str_split(df.b$data.dir, "/", n=8, simplify = TRUE)[,6], "-", simplify = TRUE)[,2]

#New envir df
df.n <- df.simul(sims.dirs2, all.gen = TRUE)
df.n$n_envir <- str_split(df.n$data.dir, "/", n=8, simplify = TRUE)[,7]
df.n$b_envir <- str_split(df.n$data.dir, "/", n=8, simplify = TRUE)[,5]
df.n$anc_id <- str_split(df.n$data.dir, "/", n=8, simplify = TRUE)[,6]
df.n$Gen <- df.n$Gen + max(df.b$Gen)

df.connect.b <- centrality.df(df.b)
df.connect.n <- centrality.df(df.n)


pc1 <- ggbetweenstats(data = df.connect.b, x = Gene_cat,  y = Sum_in,
                      centrality.plotting=FALSE,  plot.type = "box",
                      ggtheme = ggplot2::theme_bw()+
                        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)),
                      pairwise.comparisons=FALSE,  bf.message=FALSE,  results.subtitle=FALSE,
                      xlab="Gene Category", ylab="Connection Degree : in", title = NULL
)
pc2 <- ggbetweenstats(data = df.connect.b, x = Gene_cat,  y = Sum_out,
                      centrality.plotting=FALSE,  plot.type = "box",
                      ggtheme = ggplot2::theme_bw()+
                        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)),
                      pairwise.comparisons=FALSE,  bf.message=FALSE,  results.subtitle=FALSE,
                      xlab="Gene Category", ylab="Connection Degree : out", title = NULL
)

pc3 <- ggbetweenstats(data = df.connect.n, x = Cat_envir,  y = Sum_in,
                      centrality.plotting=FALSE,  plot.type = "box",
                      ggtheme = ggplot2::theme_bw()+
                        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)),
                      pairwise.comparisons=FALSE,  bf.message=FALSE,  results.subtitle=FALSE,
                      xlab="Gene Category", ylab="Connection Degree : in", title = NULL
)+
  ggplot2::scale_color_manual(values = c("aquamarine3", "deeppink", "yellow2","aquamarine3", "deeppink", "yellow2","aquamarine3", "deeppink", "yellow2","aquamarine3", "deeppink", "yellow2"))

pc4 <- ggbetweenstats(data = df.connect.n, x = Cat_envir,  y = Sum_out,
                      centrality.plotting=FALSE,  plot.type = "box",
                      ggtheme = ggplot2::theme_bw()+
                        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)),
                      pairwise.comparisons=FALSE,  bf.message=FALSE,  results.subtitle=FALSE,
                      xlab="Gene Category", ylab="Connection Degree : out", title = NULL
)+
  ggplot2::scale_color_manual(values = c("aquamarine3", "deeppink", "yellow2","aquamarine3", "deeppink", "yellow2","aquamarine3", "deeppink", "yellow2","aquamarine3", "deeppink", "yellow2"))

cairo_pdf("../figures/fig_connect_corr.pdf", width=14, height=15)
grid.arrange(
  pc1,pc3,pc2,pc4,
  ncol = 2,
  nrow = 2,
  clip = FALSE,
  widths = c(0.5,1)
)
dev.off()



#####################
sims.dirs1 <- list.dirs("../simul/burning/Noncorrelated", recursive = TRUE)
sims.dirs2 <- list.dirs("../simul/new_envir/B_Noncorrelated", recursive = TRUE)

#####################
#DATA
##################
#Burning df
df.b <- df.simul(sims.dirs1, all.gen = FALSE)
df.b$n_envir <- str_split(df.b$data.dir, "/", n=8, simplify = TRUE)[,4]
df.b$b_envir <- str_split(df.b$data.dir, "/", n=8, simplify = TRUE)[,4]
df.b$anc_id <- str_split(str_split(df.b$data.dir, "/", n=8, simplify = TRUE)[,5], "-", simplify = TRUE)[,2]

#New envir df
df.n <- df.simul(sims.dirs2, all.gen = FALSE)
df.n$n_envir <- str_split(df.n$data.dir, "/", n=8, simplify = TRUE)[,6]
df.n$b_envir <- str_split(df.n$data.dir, "/", n=8, simplify = TRUE)[,4]
df.n$anc_id <- str_split(df.n$data.dir, "/", n=8, simplify = TRUE)[,5]
df.n$Gen <- df.n$Gen + max(df.b$Gen)

df.connect.b <- centrality.df(df.b)
df.connect.n <- centrality.df(df.n)


pc1 <- ggbetweenstats(data = df.connect.b, x = Gene_cat,  y = Sum_in,
                      centrality.plotting=FALSE,  plot.type = "box",
                      ggtheme = ggplot2::theme_bw()+
                        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)),
                      pairwise.comparisons=FALSE,  bf.message=FALSE,  results.subtitle=FALSE,
                      xlab="Gene Category", ylab="Connection Degree : in", title = NULL
)
pc2 <- ggbetweenstats(data = df.connect.b, x = Gene_cat,  y = Sum_out,
                      centrality.plotting=FALSE,  plot.type = "box",
                      ggtheme = ggplot2::theme_bw()+
                        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)),
                      pairwise.comparisons=FALSE,  bf.message=FALSE,  results.subtitle=FALSE,
                      xlab="Gene Category", ylab="Connection Degree : out", title = NULL
)

pc3 <- ggbetweenstats(data = df.connect.n, x = Cat_envir,  y = Sum_in,
                      centrality.plotting=FALSE,  plot.type = "box",
                      ggtheme = ggplot2::theme_bw()+
                        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)),
                      pairwise.comparisons=FALSE,  bf.message=FALSE,  results.subtitle=FALSE,
                      xlab="Gene Category", ylab="Connection Degree : in", title = NULL
)+
  ggplot2::scale_color_manual(values = c("aquamarine3", "deeppink", "yellow2","aquamarine3", "deeppink", "yellow2","aquamarine3", "deeppink", "yellow2","aquamarine3", "deeppink", "yellow2"))

pc4 <- ggbetweenstats(data = df.connect.n, x = Cat_envir,  y = Sum_out,
                      centrality.plotting=FALSE,  plot.type = "box",
                      ggtheme = ggplot2::theme_bw()+
                        theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)),
                      pairwise.comparisons=FALSE,  bf.message=FALSE,  results.subtitle=FALSE,
                      xlab="Gene Category", ylab="Connection Degree : out", title = NULL
)+
  ggplot2::scale_color_manual(values = c("aquamarine3", "deeppink", "yellow2","aquamarine3", "deeppink", "yellow2","aquamarine3", "deeppink", "yellow2","aquamarine3", "deeppink", "yellow2"))


cairo_pdf("../figures/fig_connect_noncorr.pdf", width=14, height=15)
grid.arrange(
  pc1,pc3,pc2,pc4,
  ncol = 2,
  nrow = 2,
  clip = FALSE,
  widths = c(0.5,1)
)
dev.off()