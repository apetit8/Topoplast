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


#gives a table with count of coherence and incoherence loop
coherence_count <- function(list.w, cutoff.max=3, cutoff.min=1, target=2){
  results <- sapply(list.w, function(w){
    g <- graph.adjacency(abs(t(w)), mode="directed")
    E(g)$sign <- (w)[(w) != 0] #Signs
    # ff <- feedforward.to(g, to=2, cutoff.max=3, cutoff.min=1)
    ff <- ifelse(!is.null(feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min)), feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min), "no_ff")
    #Coherence or incoherence of FF loops
    coherence <- sapply(ff, function(n){
      if(n[1]=="no_ff") return("No_FF") else{
        reg1 <- E(g, path=n[[1]])$sign[length(E(g, path=n[[1]])$sign)] #Take the last sign since it's the x -> Target
        reg2 <- E(g, path=n[[2]])$sign[length(E(g, path=n[[2]])$sign)] #Take the last sign since it's the x -> Target
        browser()
        return(ifelse(reg1==reg2, c("FF_Coherent"), c("FF_Incoherent")))
      }})
    return(coherence)
  })
  return(table(results))
}
