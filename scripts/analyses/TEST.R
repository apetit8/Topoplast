#TEST DEBUG OSKOUR

igraph_options(return.vs.es=F)
#Loop and loop coherence count
c.count <- function(list.w, cutoff.max=3, cutoff.min=1, target=2, randomFF=FALSE){
  df <- data.frame(Coherent=c(rep(0, length(list.w))), Incoherent=c(rep(0, length(list.w))), No_loop=c(rep(0, length(list.w))), Loop=c(rep(0, length(list.w))))
  i <- 1
  for(i in 1:length(list.w)){
    g <- graph.adjacency(abs(t(list.w[[i]])), mode="directed") # abs(list.w[[i]]) mode="directed"
    E(g)$sign <- (list.w[[i]])[list.w[[i]] != 0] #signs ; does not work with signed=TRUE because of the negative values. 
    # OK: with directed, when reg=5, count it 5x ... Change how matrix is encoded !!
    if(is.null(feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min))){
      df[i,3] <- 1  } else{
        df[i,4] <- 1
        ff <- feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min)
        if(randomFF==TRUE){
          n <- ff[[sample(length(ff),1)]]
          reg1 <- E(g, path=n[[1]])$sign[length(E(g, path=n[[1]])$sign)] #Take the last sign since it's the x -> Target
          reg2 <- E(g, path=n[[2]])$sign[length(E(g, path=n[[2]])$sign)] #Take the last sign since it's the x -> Target
          df[i,1] <- df[i,1] + ifelse(reg1==reg2, 1, 0)
          df[i,2] <- df[i,2] + ifelse(reg1!=reg2, 1, 0)
        }else{
          for(n in ff) {
            reg1 <- E(g, path=n[[1]])$sign[length(E(g, path=n[[1]])$sign)] #Take the last sign since it's the x -> Target
            reg2 <- E(g, path=n[[2]])$sign[length(E(g, path=n[[2]])$sign)] #Take the last sign since it's the x -> Target
            df[i,1] <- df[i,1] + ifelse(reg1==reg2, 1/length(ff), 0)
            df[i,2] <- df[i,2] + ifelse(reg1!=reg2, 1/length(ff), 0)
          }}
      }}
  return(df)
}


Anticor3 <- c.count(topo.anticor3, cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE)




phloops <- lapply(phgenes[1], function(gene) {
  cc <- c.count(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min,
                      randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(cc)
})


phloops <- lapply(c("aaeR"), function(gene) {
  cc <- c.count(list(E_coli_mat[1:238,1:238]), cutoff.max = cutoff.max, cutoff.min = cutoff.min, target = which(colnames(E_coli_mat[1:100,1:100])==gene))
  # cc <- FFL.coherence(list(E_coli_mat), cutoff.max = cutoff.max, cutoff.min = cutoff.min, randomFF=FALSE, target = which(colnames(E_coli_mat)==gene))
  return(c(gene,cc))
})

#240



mat <- matrix(sample(c(-1, 1, 0), 239*239, replace=TRUE), ncol=239)
g <- graph.adjacency(abs(t(mat)), mode="directed")
E(g)$sign <- (mat)[mat != 0] #This is necessary to not attributes edge weights, and only look for weither or not nodes are connected
cc <- feedforward.to(g, to=10, cutoff.max=5, cutoff.min=1)






