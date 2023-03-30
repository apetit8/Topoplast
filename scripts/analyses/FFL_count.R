source("scripts/functions/functions.R")
source("scripts/functions/detectloops.R")
library(signnet)


n <- 9
w <- matrix(sample(c(-1,0,1), n*n, replace=TRUE, prob=rep(1, 3)), ncol=n, nrow=n)
diag(w) <- 0

g <- graph.adjacency(abs(t(w)), mode="directed")

fb <- feedback.from(g, from=2, cutoff=4)
# Il y a 9 feedbacks pour le gène 2, avec max 4 étapes dans le chemin

ff <- feedforward.to(g, to=2, cutoff.max=3, cutoff.min=1)

#Signs
E(g, path=fb[[1]])
E(g)$sign <- t(w)[t(w) != 0] # attention à ne pas se gourer avec la transposition de w


sapply(fb, function(x) prod(E(g, path=x)$sign) == 1)

#feedback.from(g, from=2, cutoff=4)

ff <- feedforward.to(g, to=2, cutoff.max=3, cutoff.min=1)

table(sapply(ff, function(n){
  reg1 <- E(g, path=n[[1]])$sign[length(E(g, path=n[[1]])$sign)] #Take the last sign since it's the x -> Target
  reg2 <- E(g, path=n[[2]])$sign[length(E(g, path=n[[2]])$sign)] #Take the last sign since it's the x -> Target
  coherence <- ifelse(reg1==reg2, c("Coherent"), c("Incoherent"))
  return(coherence)
}))

################################################################################

#get essential topo from simulations

anticor_UD3 <- essential.topo(df=subset(df.3, Gen==max(df.3$Gen) & envir=="Anticorrelated_UD"),
                              treshold_coeff=treshold_coeff, treshold_og=treshold_og, genes=genes, groups = list(1,2,3))


essential.graphs <- lapply(anticor_UD3, function(W){
  graph.adjacency(abs(t(W)), mode="directed")
})




eedforward.to(g, to=2, cutoff.max=3, cutoff.min=1)



coherence_count <- function(list.w, cutoff.max=3, cutoff.min=1, target=2){
  results <- sapply(list.w, function(w){
    g <- graph.adjacency(abs(t(w)), mode="directed")
    E(g)$sign <- t(w)[t(w) != 0] #Signs
    # ff <- feedforward.to(g, to=2, cutoff.max=3, cutoff.min=1)
    ff <- ifelse(!is.null(feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min)), feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min), "no_ff")
    #Coherence or incoherence of FF loops
    coherence <- sapply(ff, function(n){
      if(n[1]=="no_ff") return("No_FF") else{
        reg1 <- E(g, path=n[[1]])$sign[length(E(g, path=n[[1]])$sign)] #Take the last sign since it's the x -> Target
        reg2 <- E(g, path=n[[2]])$sign[length(E(g, path=n[[2]])$sign)] #Take the last sign since it's the x -> Target
        return(ifelse(reg1==reg2, c("FF_Coherent"), c("FF_Incoherent")))
      }})
    return(coherence)
  })
  return(table(results))
}

loop_count(anticor_UD4, cutoff.max = 3, cutoff.min = 2)
loop_count(anticor_UD3, cutoff.max = 3, cutoff.min = 1)
loop_count(corr_UD4, cutoff.max = 3, cutoff.min = 2)
loop_count(corr_UD3, cutoff.max = 3, cutoff.min = 1)



