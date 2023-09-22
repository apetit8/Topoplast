library(igraph)
library(signnet)

.my_all_simple_paths <- function(graph, from, to=igraph::V(graph), mode=c("out", "in", "all", "total"),  cutoff = -1) {
	# Technical function to deal with old versions of igraph::all_simple_paths that do not have the cutoff argument
	if ("cutoff" %in% names(formals(igraph::all_simple_paths))) {
		igraph::all_simple_paths(graph=graph, from=from, to=to, mode=mode, cutoff=cutoff)
	} else {
		a <- igraph::all_simple_paths(graph=graph, from=from, to=to, mode=mode)
		if (cutoff > 0) 
			a <- a[sapply(a, length) <= cutoff + 1]
		a
}}

.loops.from.paths <- function(paths1, paths2 = NULL) {
	# Computes a list of loops from two lists of paths
	# Loops with duplicated nodes are removed
	# The logic of the for loops is different:
	# * both paths1 and paths2 provided: all combinations of paths1 and paths2 elements are returned
	# * only paths1 is provided: all pairs of paths from paths1 are returned
	
	if (
		( is.null(paths2) && length(paths1) < 2) || 
		(!is.null(paths2) && (length(paths1) < 1 || length(paths2) < 1)))
		return(list())
	
	ans <- list()
	
	loop1 <- if(is.null(paths2))
		     seq(1, length(paths1)-1)
		else seq_along(paths1)
		
	for (i in loop1) {
		
		loop2 <- if(is.null(paths2))
		     seq(i+1, length(paths1))
		else seq_along(paths2)
		
		for (j in loop2) {
			p1 <- paths1[[i]]
			p2 <- if(is.null(paths2))
			     paths1[[j]]
			else paths2[[j]]
			if (length(setdiff(p1, p2)) == length(p1) - 2) # no duplicated nodes
				ans[[length(ans)+1]] <- list(p1, p2)
		}
	}
	ans
}

feedback.from.to <- function(
	graph, 
	from, 
	to, 
	edges1 = c(1, igraph::gorder(graph)-1), 
	edges2 = c(1, igraph::gorder(graph)-1))
{
	# path1 is from "from" to "to", path2 from "to" to "from"
	path.from <- .my_all_simple_paths(graph=graph, from=from[1], to=to[1], mode="out", cutoff=max(edges1))
	path.to   <- .my_all_simple_paths(graph=graph, from=to[1], to=from[1], mode="out", cutoff=max(edges2))
	
	path.from.list <- lapply(seq_len(max(edges1)), function(i) path.from[sapply(path.from, length) == i + 1])
	path.to.list   <- lapply(seq_len(max(edges2)), function(i) path.to  [sapply(path.to  , length) == i + 1])
	
	seq.path1     <- seq(min(edges1), max(edges1))
	seq.path2     <- seq(min(edges2), max(edges2))
	
	ans <- list()
	
	for (i in seq.path1) {
		for (j in seq.path2) {
			ans <- c(ans, .loops.from.paths(path.from.list[[i]], path.to.list[[j]]))
		}
	}
	return(ans)
}

feedback.from <- function(
	graph, 
	from, 
	edges    = c(2, igraph::gorder(graph)),
	collapse = TRUE)
{
	ans <- do.call(c, lapply(
			               igraph::V(graph), 
			FUN          = feedback.from.to, 
			graph        = graph, 
			from         = from, 
			edges1 = c(1,1), 
			edges2 = c(min(edges)-1, max(edges)-1)
		))
			
		# Remove duplicated path -- this is a trick, because unique() does not work for lists
		# This step can hardly be avoided, as several pairs (from, to) might give the same feedback loop

	ans <- ans[!duplicated(sapply(ans, paste0, collapse="-"))]
	
	if (collapse)
		ans <- lapply(ans, function(x) c(x[[1]][1:2], x[[2]][2:length(x[[2]])]))
	
	return(ans)
}

feedforward.from.to <- function(
	graph, 
	from, 
	to, 
	edges1 = c(1, gorder(graph)-1), 
	edges2 = c(1, gorder(graph)-1)) 
{
	ans <- list()
	
	seq.path1     <- seq(min(edges1), max(edges1))
	seq.path2     <- seq(min(edges2), max(edges2))
	shortest.path <- min(seq.path1, seq.path2)
	longest.path  <- max(seq.path1, seq.path2)
	
	paths <- .my_all_simple_paths(
		graph = graph, 
		from  = from[1], 
		to    = to[1], 
		mode  = "out", 
		cutoff= longest.path)
	
	if (length(paths) < 2)
		return(list())
		
	path.list <- lapply(seq_len(longest.path), function(i) paths[sapply(paths, length) == i + 1])

	for (p1 in seq.path1) {
		for (p2 in seq.path2) {
			if (!p2 %in% seq.path1 || p2 >= p1) # Avoids mirror loops
				ans <- c(ans, 
					if (p1 == p2) .loops.from.paths(path.list[[p1]])
					else          .loops.from.paths(path.list[[p1]], path.list[[p2]])
				)
		}
	}
	ans
}

feedforward.to <- function(
	graph, 
	to, 
	edges1 = c(1, igraph::gorder(graph)-1), 
	edges2 = c(1, igraph::gorder(graph)-1)) 
{
	unlist(lapply(
			         igraph::V(graph), 
			FUN    = feedforward.from.to, 
			graph  = graph, 
			to     = to, 
			edges1 = edges1, 
			edges2 = edges2), 
		recursive=FALSE)
}

FFL.type2 <- function(list.w, edges1=2, edges2=1, target, from=(1:ncol(list.w[[1]]))){
  
  if(edges1==2 && edges2==1){
    df <- data.frame(FFL=c(rep(0, length(list.w))), No_FFL=c(rep(0, length(list.w))),
    C3=c(rep(0, length(list.w))), C1=c(rep(0, length(list.w))), C2=c(rep(0, length(list.w))), C4=c(rep(0, length(list.w))),
    I2=c(rep(0, length(list.w))), I4=c(rep(0, length(list.w))), I3=c(rep(0, length(list.w))), I1=c(rep(0, length(list.w))), NonPl_X=c(rep(0, length(list.w))))   }
 
   else if(edges1==2 && edges2==2){
    df <- data.frame(FFL=c(rep(0, length(list.w))), No_FFL=c(rep(0, length(list.w))),
    PP=c(rep(0, length(list.w))), PM=c(rep(0, length(list.w))), PN=c(rep(0, length(list.w))), NP=c(rep(0, length(list.w))),
    NM=c(rep(0, length(list.w))), NN=c(rep(0, length(list.w))), MP=c(rep(0, length(list.w))), MM2=c(rep(0, length(list.w))),
    MM1=c(rep(0, length(list.w))), MN=c(rep(0, length(list.w))), NonPl_X=c(rep(0, length(list.w))))   }
   else {
    df <- data.frame(FFL=c(rep(0, length(list.w))), No_FFL=c(rep(0, length(list.w))), NonPl_X=c(rep(0, length(list.w))))  }
  for(i in 1:length(list.w)){
    g <- graph.adjacency(abs(t(list.w[[i]])), mode="directed") 
    E(g)$sign <- (list.w[[i]])[list.w[[i]] != 0] 
    if(length(feedforward.to(g, to=target, edges1=edges1, edges2=edges2))==0){
      df[i,2] <- 1 } else{
        df[i,1] <- 1
        
        ff <- feedforward.to(g, to=target, edges1=edges1, edges2=edges2)
        if(edges1==2 && edges2==1){
        for(m in ff) {
          if(m[[1]][1] %in% from){ #If the first node is in FROM
          if(length(E(g, path=m[[1]])$sign)>1 && length(E(g, path=m[[2]])$sign)>2) stop("Wrong cut.offs")
            if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]==E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==-1) df[i,3] <- df[i,3] + 1/length(ff)
            if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]==E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==1) df[i,4] <- df[i,4] + 1/length(ff)
            if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]!=E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==-1) df[i,5] <- df[i,5] + 1/length(ff)
            if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]!=E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==1) df[i,6] <- df[i,6] + 1/length(ff)
            if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==-1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]==E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==-1) df[i,7] <- df[i,7] + 1/length(ff)
            if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==-1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]==E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==1) df[i,8] <- df[i,8] + 1/length(ff)
            if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==-1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]!=E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==-1) df[i,9] <- df[i,9] + 1/length(ff)
            if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==-1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]!=E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==1) df[i,10] <- df[i,10] + 1/length(ff)
          }else{df[i,11] <- df[i,11] + 1/length(ff)}
          }}
        else if(edges1==2 && edges2==2){
          for(m in ff) {
            if(m[[1]][1] %in% from){
              if(length(E(g, path=m[[1]])$sign)!=2 || length(E(g, path=m[[2]])$sign)!=2) stop("Wrong cut.offs")
              reg_j1 <- E(g, path=m[[1]])$sign[1]
              reg_j2 <- E(g, path=m[[1]])$sign[2]
              reg_i1 <- E(g, path=m[[2]])$sign[1]
              reg_i2 <- E(g, path=m[[2]])$sign[2]
              # browser()
              #BUG : Cutoffs are not working !!
              if(reg_j1==1 && reg_i1==1 && reg_j2==1 && reg_i2==1) df[i,3] <- df[i,3] + 1/length(ff)
              if(reg_j1==1 && reg_i1==1 && reg_j2!=reg_i2) df[i,4] <- df[i,4] + 1/length(ff)
              if(reg_j1==1 && reg_i1==1 && reg_j2==-1 && reg_i2==-1) df[i,5] <- df[i,5] + 1/length(ff)
              #
              if(reg_j1==-1 && reg_i1==-1 && reg_j2==1 && reg_i2==1) df[i,6] <- df[i,6] + 1/length(ff)
              if(reg_j1==-1 && reg_i1==-1 && reg_j2!=reg_i2) df[i,7] <- df[i,7] + 1/length(ff)
              if(reg_j1==-1 && reg_i1==-1 && reg_j2==-1 && reg_i2==-1) df[i,8] <- df[i,8] + 1/length(ff)
              #
              if(reg_j1!=reg_i1 && reg_j2==1 && reg_i2==1) df[i,9] <- df[i,9] + 1/length(ff)
              if(reg_j1!=reg_i1 && reg_j2!=reg_i2 && reg_j1==reg_j2) df[i,10] <- df[i,10] + 1/length(ff)
              if(reg_j1!=reg_i1 && reg_j2!=reg_i2 && reg_j1!=reg_j2) df[i,11] <- df[i,11] + 1/length(ff)
              if(reg_j1!=reg_i1 && reg_j2==-1 && reg_i2==-1) df[i,12] <- df[i,12] + 1/length(ff)
            }else{df[i,11] <- df[i,11] + 1/length(ff)}
          }} else         
            for(m in ff) {
            if(!(m[[1]][1] %in% from)){df[i,3] <- df[i,3] + 1/length(ff)} }
        }}
  return(df)
}


#Count number of loop for each gene
loops_n.count <- function(list.w, edges1=edges1, edges2=edges2, target=2){
  df <- data.frame(Loop_number=c(rep(0, length(list.w))), No_loop=c(rep(0, length(list.w))), Loop=c(rep(0, length(list.w))))
  for(i in 1:length(list.w)){
    g <- graph.adjacency(t(list.w[[i]]), weighted = TRUE)
    E(g)$sign <- (list.w[[i]])[list.w[[i]] != 0] #signs ; does not work with signed=TRUE because of the negative values. 
    if(is.null(feedforward.to(g, to=target, edges1=edges1, edges2=edges2))){
      df[i,2] <- 1  } else{
        df[i,3] <- 1
        ff <- feedforward.to(g, to=target, edges1=edges1, edges2=edges2)
        df[i,1] <- length(ff)
      }}
  return(df)
}


#Count number of feedbackloop for each gene
FBL_n.count <- function(list.w,  edges=2, target=2){
  df <- data.frame(FBL_number=c(rep(0, length(list.w))), No_FBL=c(rep(0, length(list.w))), FBL=c(rep(0, length(list.w))))
  for(i in 1:length(list.w)){
    g <- graph.adjacency(t(list.w[[i]]), weighted = TRUE)
    E(g)$sign <- (list.w[[i]])[list.w[[i]] != 0] #signs ; does not work with signed=TRUE because of the negative values. 
    if(is.null(feedback.from(g, from=target, edges=edges))){
      df[i,2] <- 1  } else{
        df[i,3] <- 1
        ff <- feedback.from(g, from=target, edges=edges)
        df[i,1] <- length(ff)
      }}
  return(df)
}



#FBL and FBL homogeneity count
FBL.type <- function(list.w, edges=c(2:5), target=2, randomFF=FALSE){
  df <- data.frame(FBL=c(rep(0, length(list.w))), No_FBL=c(rep(0, length(list.w))),
                   Inhibiting=c(rep(0, length(list.w))), Activating=c(rep(0, length(list.w))),
                   Size2=c(rep(0, length(list.w))),  Size3=c(rep(0, length(list.w))),  Size4=c(rep(0, length(list.w))),  Size5=c(rep(0, length(list.w))),  Size6=c(rep(0, length(list.w))) )
  for(i in 1:length(list.w)){
    g <- graph.adjacency(t(list.w[[i]]), weighted = TRUE)
    E(g)$sign <- (list.w[[i]])[list.w[[i]] != 0] #signs ; does not work with signed=TRUE because of the negative values. 
    ff <- feedback.from(g, from=target, edges=edges)
    if(length(ff)==0){
      df[i,2] <- 1  } else{
        df[i,1] <- 1
        if(randomFF==TRUE){
          nn <- ff[[sample(length(ff),1)]]
          reg1 <- E(g, path=nn)$sign #When regulations are heterogeneous, the product is -1
          df[i,3] <- df[i,3] + ifelse(reg1[1]==-1, 1, 0)
          df[i,4] <- df[i,4] + ifelse(reg1[1]==1, 1, 0)
        }else{
          for(nn in 1:length(ff)) {
            size <- length(E(g, path=ff[[nn]])$sign)
            ifelse(size==2, df[i,5] <- df[i,5] + 1/length(ff), ifelse(size==3, df[i,6] <- df[i,6] + 1/length(ff), ifelse(size==4, df[i,7] <- df[i,7] + 1/length(ff), ifelse(size==5, df[i,8] <- df[i,8] + 1/length(ff), df[i,9] <- df[i,9] + 1/length(ff))) ))
            reg1 <- E(g, path=ff[[nn]])$sign #When regulations are heterogeneous, the product is -1
            df[i,3] <- df[i,3] + ifelse(reg1[1]==-1, 1/length(ff), 0)
            df[i,4] <- df[i,4] + ifelse(reg1[1]==1, 1/length(ff), 0)
          }}
      }}
  return(df)
}


#g_mat is now useless
e_coli_prep_analyses <- function(genes_list, g, g_mat, fun="FFL", edges1=2, edges2=1, cores=2, from=FALSE){
  stopifnot(fun=="FFL" || fun=="FBL" || fun=="FFLcount"|| fun=="FBLcount" )
  if(length(edges1)==1){stopifnot(edges1 >= edges2)}
  loops <- mclapply(genes_list, function(gene) {
    #To avoid igraph::all_simple_paths to take days and weeks, we subset the regulatory networks. Only target genes and connected TFs are kept.
    gg <- induced.subgraph(g, vids = c(which(colnames(E_coli_mat)==gene), which(colnames(E_coli_mat)%in%TF_genes)) )
    E_coli_mat2 <- t((get.adjacency(gg, sparse=FALSE, attr='V3')))
    ggg <- induced.subgraph(gg, vids = as.vector(unlist(neighborhood(gg, max(edges1), nodes = which(colnames(E_coli_mat2)==gene), mode = 'all'))))
    E_coli_mat3 <- t((get.adjacency(ggg, sparse=FALSE, attr='V3')))
    E_coli_mat3 <- matrix(as.numeric(E_coli_mat3), ncol = ncol(E_coli_mat3), dimnames = dimnames(E_coli_mat3)) #convert to numeric matrix
    E_coli_mat3[is.na(E_coli_mat3)] <- 0
    if(is.character(from)) envirgenes <- which(colnames(E_coli_mat3)%in%from)
    else if(from==TRUE) envirgenes <- which(colnames(E_coli_mat3)%in%genes_list) else envirgenes <- 1:ncol(E_coli_mat3)
    print(paste0(gene," ; network size ", ncol(E_coli_mat3)))
    if(fun=="FFL") cc <- FFL.type2(list(E_coli_mat3), edges1 = edges1, edges2 = edges2, target = which(colnames(E_coli_mat3)==gene), from=envirgenes)
    if(fun=="FFLcount") cc <- loops_n.count(list(E_coli_mat3), edges1 = edges1, edges2 = edges2, target = which(colnames(E_coli_mat3)==gene))
    if(fun=="FBL") cc <- FBL.type(list(E_coli_mat3), edges=edges1, target = which(colnames(E_coli_mat3)==gene))
    if(fun=="FBLcount") cc <- FBL_n.count(list(E_coli_mat3), edges = edges1, target = which(colnames(E_coli_mat3)==gene))
    #Here, insert writing c(gene,cc) in a txt file
    return(c(gene,cc))
  }, mc.cores = cores)
  return(rbindlist(loops))
}

