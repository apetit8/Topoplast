library(igraph)
library(signnet)

my_all_simple_paths <- function(graph, from, to=igraph::V(graph), mode=c("out", "in", "all", "total"),  cutoff = -1) {
	# Technical function to deal with old versions of igraph::all_simple_paths that does not have the cutoff argument
	if (cutoff %in% names(formals((igraph::all_simple_paths)))) {
		igraph::all_simple_paths(graph=graph, from=from, to=tp, mode=mode, cutoff=cutoff)
	} else {
		a <- igraph::all_simple_paths(graph=graph, from=from, to=to, mode=mode)
		if (cutoff > 0) 
			a <- a[sapply(a, length) <= cutoff]
		a
}}

feedback.from.to <- function(graph, from, to, cutoff = -1) {
	path.from <- my_all_simple_paths(graph=graph, from=from[1], to=to[1], mode="out", cutoff=cutoff-1)
	path.to   <- my_all_simple_paths(graph=graph, from=to[1], to=from[1], mode="out", cutoff=cutoff-1)
	
	if (length(path.from) == 0 || length(path.to) == 0)
		return(NULL)
	ee <- expand.grid(seq_along(path.from), seq_along(path.to))
	ans <- lapply(1:nrow(ee), function(i) c(path.from[[ee[i,1]]], path.to[[ee[i,2]]][-1]))
	if (cutoff > 0) 
		ans <- ans[sapply(ans, length) <= cutoff]
	if (length(ans) > 0) 	# Remove answers with duplcated nodes (beyond "from")
		ans <- ans[sapply(ans, function(x) sum(duplicated(x)) == 1)]
	ans
}

feedback.from <- function(graph, from, cutoff = -1) {
	ans <- do.call(c, lapply(igraph::V(g), feedback.from.to, graph=graph, from=from, cutoff=cutoff))
		# Remove duplicated path -- this is a trick, because unique() does not work for lists
	ans <- ans[!duplicated(sapply(ans, paste0, collapse="-"))]
	ans
}

feedforward.from.to <- function(graph, from, to, cutoff.max = -1, cutoff.min = 1) {
	path.from.to <- my_all_simple_paths(graph=graph, from=from[1], to=to[1], mode="out", cutoff=cutoff.max)
	
	which.small <- which(sapply(path.from.to, length) <= cutoff.min + 1)
	ans <- unlist(lapply(which.small, function(i) {
						relevant.path <- NULL
						for (j in seq_along(path.from.to)) {
							if (j != i && (j > i || !j%in%which.small) && sum(duplicated(c(path.from.to[[i]], path.from.to[[j]]))) == 2)
								relevant.path <- append(relevant.path, j)
						}
						if (length(relevant.path) > 0) {
							return(lapply(relevant.path, function(p) list(path.from.to[[i]], path.from.to[[p]])))
						} else {
							return(NULL)
						}
					}
				 ), recursive=FALSE)
}

feedforward.to <- function(graph, to, cutoff.max = -1, cutoff.min = 1) {
	ans <- unlist(lapply(V(graph), feedforward.from.to, graph=graph, to=to, cutoff.max=cutoff.max, cutoff.min=cutoff.min), recursive=FALSE)
}



#Loop and loop coherence count
c.count <- function(list.w, cutoff.max=3, cutoff.min=1, target=2, randomFF=FALSE){
  df <- data.frame(Coherent=c(rep(0, length(list.w))), Incoherent=c(rep(0, length(list.w))), No_loop=c(rep(0, length(list.w))), Loop=c(rep(0, length(list.w))))
  for(i in 1:length(list.w)){
    g <- graph.adjacency(t(list.w[[i]]), weighted = TRUE)
    E(g)$sign <- (list.w[[i]])[list.w[[i]] != 0] #signs ; does not work with signed=TRUE because of the negative values. 
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


#Loop and loop homogeneity count
homog.count <- function(list.w, cutoff.max=3, cutoff.min=1, target=2, randomFF=FALSE){
  df <- data.frame(Heterogenous=c(rep(0, length(list.w))), Homogenous=c(rep(0, length(list.w))), No_loop=c(rep(0, length(list.w))), Loop=c(rep(0, length(list.w))))
  i <- 1
  for(i in 1:length(list.w)){
    g <- graph.adjacency(t(list.w[[i]]), weighted = TRUE)
    E(g)$sign <- (list.w[[i]])[list.w[[i]] != 0] #signs ; does not work with signed=TRUE because of the negative values. 
    if(is.null(feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min))){
      df[i,3] <- 1  } else{
        df[i,4] <- 1
        ff <- feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min)
        if(randomFF==TRUE){
          n <- ff[[sample(length(ff),1)]]
          reg1 <- prod(E(g, path=n[[1]])$sign) #When regulations are heterogeneous, the product is -1
          reg2 <- prod(E(g, path=n[[2]])$sign) 
          df[i,1] <- df[i,1] + ifelse(reg1*reg2==-1, 1, 0)
          df[i,2] <- df[i,2] + ifelse(reg1*reg2==1, 1, 0)
        }else{
          for(n in ff) {
            reg1 <- prod(E(g, path=n[[1]])$sign) #When regulations are heterogeneous, the product is -1
            reg2 <- prod(E(g, path=n[[2]])$sign) 
            df[i,1] <- df[i,1] + ifelse(reg1==-1 || reg2==-1, 1/length(ff), 0)
            df[i,2] <- df[i,2] + ifelse(reg1==1 && reg2==1, 1/length(ff), 0)
          }}
      }}
  return(df)
}



#Count number of loop for each gene
loops_n.count <- function(list.w, cutoff.max=3, cutoff.min=1, target=2){
  df <- data.frame(Loop_number=c(rep(0, length(list.w))), No_loop=c(rep(0, length(list.w))), Loop=c(rep(0, length(list.w))))
  i <- 1
  for(i in 1:length(list.w)){
    g <- graph.adjacency(t(list.w[[i]]), weighted = TRUE)
    E(g)$sign <- (list.w[[i]])[list.w[[i]] != 0] #signs ; does not work with signed=TRUE because of the negative values. 
    if(is.null(feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min))){
      df[i,2] <- 1  } else{
        df[i,3] <- 1
        ff <- feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min)
        df[i,1] <- length(ff)
      }}
  return(df)
}

