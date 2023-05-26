library(igraph)
library(signnet)

my_all_simple_paths <- function(graph, from, to=igraph::V(graph), mode=c("out", "in", "all", "total"),  cutoff = -1) {
	# Technical function to deal with old versions of igraph::all_simple_paths that does not have the cutoff argument
	if ("cutoff" %in% names(formals(igraph::all_simple_paths))) {
		igraph::all_simple_paths(graph=graph, from=from, to=to, mode=mode, cutoff=cutoff)
	} else {
		a <- igraph::all_simple_paths(graph=graph, from=from, to=to, mode=mode)
		if (cutoff > 0) 
			a <- a[sapply(a, length) <= cutoff]
		a
}}

feedback.from.to <- function(graph, from, to, cutoff = -1) {
  path.from <- my_all_simple_paths(graph=graph, from=from[1], to=to[1], mode="out", cutoff=cutoff)
	path.to   <- my_all_simple_paths(graph=graph, from=to[1], to=from[1], mode="out", cutoff=cutoff)
	if (length(path.from) == 0 || length(path.to) == 0)
		return(NULL)
	ee <- expand.grid(seq_along(path.from), seq_along(path.to))
	ans <- lapply(1:nrow(ee), function(i) c(path.from[[ee[i,1]]], path.to[[ee[i,2]]][-1]))
	if (cutoff > 0) 
		ans <- ans[sapply(ans, length) <= cutoff]
	if (length(ans) > 0) 	# Remove answers with duplcated nodes (beyond "from")
		ans <- ans[sapply(ans, function(x) sum(duplicated(x)) == 1)]
	return(ans)
}

feedback.from <- function(graph, from, cutoff = -1) {
	ans <- do.call(c, lapply(igraph::V(graph), feedback.from.to, graph=graph, from=from, cutoff=cutoff))
		# Remove duplicated path -- this is a trick, because unique() does not work for lists
	ans <- ans[!duplicated(sapply(ans, paste0, collapse="-"))]
	return(ans)
}

#BUG: get 
feedforward.from.to <- function(graph, from, to, cutoff.max = -1, cutoff.min = 1) {
	path.from.to <- my_all_simple_paths(graph=graph, from=from[1], to=to[1], mode="out", cutoff=cutoff.max) 
	
	# which.small <- which(sapply(path.from.to, length) <= cutoff.min + 1 ) #Plus petit que le minimum ? Cutoff.min = longueur max de la plus petite branche ?
	which.small <- which(sapply(path.from.to, length) == cutoff.min )
	ans <- unlist(lapply(which.small, function(i) {
						relevant.path <- NULL
						for (j in seq_along(path.from.to)) {
						  #Problem: what if both branches are of the same length ? It does not count it ?
						  if (j != i && j >= i && length(path.from.to[[j]]) <= cutoff.max && sum(duplicated(c(path.from.to[[i]], path.from.to[[j]]))) == 2)
							#if (j != i && (j > i || !j%in%which.small) && sum(duplicated(c(path.from.to[[i]], path.from.to[[j]]))) == 2)
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
  #test <- feedforward.from.to(graph=graph, V(graph), to=to, cutoff.max=cutoff.max, cutoff.min=cutoff.min)
	ans <- unlist(lapply(V(graph), feedforward.from.to, graph=graph, to=to, cutoff.max=cutoff.max, cutoff.min=cutoff.min), recursive=FALSE)
}


# 
# #Loop and loop coherence count
# FFL.coherence <- function(list.w, cutoff.max=3, cutoff.min=1, target=2, randomFF=FALSE){
#   df <- data.frame(Coherent=c(rep(0, length(list.w))), Incoherent=c(rep(0, length(list.w))), No_loop=c(rep(0, length(list.w))), Loop=c(rep(0, length(list.w))))
#   for(i in 1:length(list.w)){
#     g <- graph.adjacency(t(list.w[[i]]), weighted = TRUE)
#     #E(g)$sign <- (list.w[[i]])[list.w[[i]] != 0] #signs ; does not work with signed=TRUE because of the negative values. 
#     if(is.null(feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min))){
#       df[i,3] <- 1  } else{
#         df[i,4] <- 1
#         ff <- feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min)
#         if(randomFF==TRUE){
#           n <- ff[[sample(length(ff),1)]]
#           reg1 <- E(g, path=n[[1]])$weight[length(E(g, path=n[[1]])$weight)] #Take the last weight since it's the x -> Target
#           reg2 <- E(g, path=n[[2]])$weight[length(E(g, path=n[[2]])$weight)] #Take the last weight since it's the x -> Target
#           df[i,1] <- df[i,1] + ifelse(reg1==reg2, 1, 0)
#           df[i,2] <- df[i,2] + ifelse(reg1!=reg2, 1, 0)
#         }else{
#           for(n in ff) {
#             reg1 <- E(g, path=n[[1]])$weight[length(E(g, path=n[[1]])$weight)] #Take the last weight since it's the x -> Target
#             reg2 <- E(g, path=n[[2]])$weight[length(E(g, path=n[[2]])$weight)] #Take the last weight since it's the x -> Target
#             df[i,1] <- df[i,1] + ifelse(reg1==reg2, 1/length(ff), 0)
#             df[i,2] <- df[i,2] + ifelse(reg1!=reg2, 1/length(ff), 0)
#           }}
#       }}
#   return(df)
# }

FFL.coherence <- function(list.w, cutoff.max=3, cutoff.min=1, target=2, randomFF=FALSE){
  df <- data.frame(Coherent=c(rep(0, length(list.w))), Incoherent=c(rep(0, length(list.w))), No_loop=c(rep(0, length(list.w))), Loop=c(rep(0, length(list.w))))
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

#Frequency = FALSE : the most represented motif "win"
FFL.type <- function(list.w, cutoff.max=3, cutoff.min=1, target=2, frequencies=TRUE){
  df <- data.frame(Activating=c(rep(0, length(list.w))), Inhibiting=c(rep(0, length(list.w))), Z_act=c(rep(0, length(list.w))), Z_inh=c(rep(0, length(list.w))), No_loop=c(rep(0, length(list.w))), Loop=c(rep(0, length(list.w))))
  for(i in 1:length(list.w)){
    g <- graph.adjacency(abs(t(list.w[[i]])), mode="directed") # abs(list.w[[i]]) mode="directed"
    E(g)$sign <- (list.w[[i]])[list.w[[i]] != 0] #signs ; does not work with signed=TRUE because of the negative values. 
    if(is.null(feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min))){
      df[i,5] <- 1  } else{
        df[i,6] <- 1
        ff <- feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min)
        if(frequencies==FALSE){
          act   <- 0
          inh   <- 0
          z_act <- 0
          z_inh <- 0
          #for n in ff
            for (n in ff) {
              reg1 <- E(g, path=n[[1]])$sign[length(E(g, path=n[[1]])$sign)] #Take the last sign since it's the x -> Target
              reg2 <- E(g, path=n[[2]])$sign[length(E(g, path=n[[2]])$sign)]
              act   <- act + ifelse(reg1==reg2 && reg2==1, 1, 0)
              inh   <- inh + ifelse(reg1==reg2 && reg2==-1, 1, 0)
              z_act <- z_act + ifelse(reg1!=reg2 && reg2==-1, 1, 0)
              z_inh <-  z_inh + ifelse(reg1!=reg2 && reg2==1, 1, 0)
            }
          fract <- length(which(c(act, inh, z_act, z_inh)==max(c(act, inh, z_act, z_inh))))
          df[i,1] <- df[i,1] + ifelse(act >= inh && act >= z_act && act >= z_inh, 1/fract, 0)
          df[i,2] <- df[i,2] + ifelse(inh >= act && inh >= z_act && inh >= z_inh, 1/fract, 0)
          df[i,3] <- df[i,3] + ifelse(z_act >= act && z_act >= z_inh && z_act >= inh, 1/fract, 0)
          df[i,4] <- df[i,4] + ifelse(z_inh >= act && z_inh >= z_act && z_inh >= inh, 1/fract, 0)
          
        }else{
          for(n in ff) {
            reg1 <- E(g, path=n[[1]])$sign[length(E(g, path=n[[1]])$sign)] #Take the last sign since it's the x -> Target
            reg2 <- E(g, path=n[[2]])$sign[length(E(g, path=n[[2]])$sign)] #Take the last sign since it's the x -> Target
            df[i,1] <- df[i,1] + ifelse(reg1==reg2 && reg2==1, 1/length(ff), 0)
            df[i,2] <- df[i,2] + ifelse(reg1==reg2 && reg2==-1, 1/length(ff), 0)
            df[i,3] <- df[i,3] + ifelse(reg1!=reg2 && reg2==-1, 1/length(ff), 0)
            df[i,4] <- df[i,4] + ifelse(reg1!=reg2 && reg2==1, 1/length(ff), 0)
          }}
      }}
  return(df)
}


#Custom motif categories : adapted to our simulation model with constitutive expression different from 0
#change everything to get a tab ? 5 column, "Loop","No loop", "Z_sign", "In/Coherence" and "True/False" and column of concatenate last 3 columns
FFL.type2 <- function(list.w, cutoff.max=3, cutoff.min=1, target=2){
  
  if(cutoff.max==3 && cutoff.min==1){
    df <- data.frame(FFL=c(rep(0, length(list.w))), No_FFL=c(rep(0, length(list.w))),
    ID_A_N=c(rep(0, length(list.w))), ID_A_P=c(rep(0, length(list.w))), ID_D_N=c(rep(0, length(list.w))), ID_D_N=c(rep(0, length(list.w))),
    II_A_N=c(rep(0, length(list.w))), II_A_P=c(rep(0, length(list.w))), II_D_N=c(rep(0, length(list.w))), II_D_N=c(rep(0, length(list.w))) )   }
 
   if(cutoff.max==3 && cutoff.min==3){
    df <- data.frame(FFL=c(rep(0, length(list.w))), No_FFL=c(rep(0, length(list.w))),
    Pos_pos=c(rep(0, length(list.w))), Pos_mixt=c(rep(0, length(list.w))), Pos_neg=c(rep(0, length(list.w))), Neg_pos=c(rep(0, length(list.w))),
    Neg_mixt=c(rep(0, length(list.w))), Neg_neg=c(rep(0, length(list.w))), Mixt_pos=c(rep(0, length(list.w))), Mixt_mixt_hom=c(rep(0, length(list.w))),
    Mixt_mixt_het=c(rep(0, length(list.w))), Mixt_neg=c(rep(0, length(list.w))))   }
  for(i in 1:length(list.w)){
    g <- graph.adjacency(abs(t(list.w[[i]])), mode="directed") 
    E(g)$sign <- (list.w[[i]])[list.w[[i]] != 0] 
    if(is.null(feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min))){
      df[i,2] <- 1 } else{
        df[i,1] <- 1
        
        ff <- feedforward.to(g, to=target, cutoff.max=cutoff.max, cutoff.min=cutoff.min)
        if(cutoff.max==3 && cutoff.min==1){
        for(m in ff) {
          if(length(E(g, path=m[[1]])$sign)>1 && length(E(g, path=m[[2]])$sign)>2) stop("Wrong cut.offs")
          B1 <- FFL_Z_reg(g, m, Input=1)
          B0 <- FFL_Z_reg(g, m, Input=0)
          if(B1!=B0 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]==E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==-1) df[i,3] <- df[i,3] + 1/length(ff)
          if(B1!=B0 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]==E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==1) df[i,4] <- df[i,4] + 1/length(ff)
          if(B1!=B0 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]!=E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==-1) df[i,5] <- df[i,5] + 1/length(ff)
          if(B1!=B0 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]!=E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==1) df[i,6] <- df[i,6] + 1/length(ff)
          if(B1==B0 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]==E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==-1) df[i,7] <- df[i,7] + 1/length(ff)
          if(B1==B0 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]==E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==1) df[i,8] <- df[i,8] + 1/length(ff)
          if(B1==B0 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]!=E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==-1) df[i,9] <- df[i,9] + 1/length(ff)
          if(B1==B0 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]!=E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==1) df[i,10] <- df[i,10] + 1/length(ff)
        }}
        if(cutoff.max==3 && cutoff.min==3){
          for(m in ff) {
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
          }}
      }}
  return(df)
}


#Count number of loop for each gene
loops_n.count <- function(list.w, cutoff.max=3, cutoff.min=1, target=2){
  df <- data.frame(Loop_number=c(rep(0, length(list.w))), No_loop=c(rep(0, length(list.w))), Loop=c(rep(0, length(list.w))))
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


#Count number of feedbackloop for each gene
FBL_n.count <- function(list.w, cutoff=2, target=2){
  df <- data.frame(FBL_number=c(rep(0, length(list.w))), No_FBL=c(rep(0, length(list.w))), FBL=c(rep(0, length(list.w))))
  for(i in 1:length(list.w)){
    g <- graph.adjacency(t(list.w[[i]]), weighted = TRUE)
    E(g)$sign <- (list.w[[i]])[list.w[[i]] != 0] #signs ; does not work with signed=TRUE because of the negative values. 
    if(is.null(feedback.from(g, from=target, cutoff=cutoff))){
      df[i,2] <- 1  } else{
        df[i,3] <- 1
        ff <- feedback.from(g, from=target, cutoff=cutoff)
        df[i,1] <- length(ff)
      }}
  return(df)
}



#FBL and FBL homogeneity count
FBL.type <- function(list.w, cutoff=3, target=2, randomFF=FALSE){
  df <- data.frame(Inhibiting=c(rep(0, length(list.w))), Activating=c(rep(0, length(list.w))), No_FBL=c(rep(0, length(list.w))), FBL=c(rep(0, length(list.w))))
  for(i in 1:length(list.w)){
    g <- graph.adjacency(t(list.w[[i]]), weighted = TRUE)
    E(g)$sign <- (list.w[[i]])[list.w[[i]] != 0] #signs ; does not work with signed=TRUE because of the negative values. 
    ff <- feedback.from(g, from=target, cutoff=cutoff)
    if(length(ff)==0){
      df[i,3] <- 1  } else{
        df[i,4] <- 1
        if(randomFF==TRUE){
          nn <- ff[[sample(length(ff),1)]]
          reg1 <- E(g, path=nn)$sign #When regulations are heterogeneous, the product is -1
          df[i,1] <- df[i,1] + ifelse(reg1[1]==-1, 1, 0)
          df[i,2] <- df[i,2] + ifelse(reg1[1]==1, 1, 0)
        }else{
          for(nn in 1:length(ff)) {
            reg1 <- E(g, path=ff[[nn]])$sign #When regulations are heterogeneous, the product is -1
            df[i,1] <- df[i,1] + ifelse(reg1[1]==-1, 1/length(ff), 0)
            df[i,2] <- df[i,2] + ifelse(reg1[1]==1, 1/length(ff), 0)
          }}
      }}
  return(df)
}



e_coli_prep_analyses <- function(genes_list, g, g_mat, fun="FFL", cores=50){
  stopifnot(fun=="FFL" || fun=="FBL")
  #here, create txt file, name=paste0(filename,"_",fun,".csv")
  loops <- mclapply(genes_list, function(gene) {
    #To avoid igraph::all_simple_paths to take days and weeks, we subset the regulatory networks. Only target genes and connected TFs are kept.
    gg <- induced.subgraph(g, vids = c(which(colnames(E_coli_mat)==gene), which(colnames(E_coli_mat)%in%TF_genes)) )
    E_coli_mat2 <- t((get.adjacency(gg, sparse=FALSE, attr='V3')))
    ggg <- induced.subgraph(gg, vids = as.vector(unlist(neighborhood(gg, cutoff.max-1, nodes = which(colnames(E_coli_mat2)==gene), mode = 'all'))))
    E_coli_mat3 <- t((get.adjacency(ggg, sparse=FALSE, attr='V3')))
    E_coli_mat3 <- matrix(as.numeric(E_coli_mat3), ncol = ncol(E_coli_mat3), dimnames = dimnames(E_coli_mat3)) #convert to numeric matrix
    E_coli_mat3[is.na(E_coli_mat3)] <- 0
    if(fun=="FFL") cc <- FFL.coherence(list(E_coli_mat3), cutoff.max = cutoff.max, cutoff.min = cutoff.min, target = which(colnames(E_coli_mat3)==gene))
    if(fun=="FBL") cc <- FBL.type(list(E_coli_mat3), cutoff = cutoff, target = which(colnames(E_coli_mat3)==gene))
    #Here, insert writing c(gene,cc) in a txt file
    return(c(gene,cc))
  }, mc.cores = cores)
  return(rbindlist(loops))
}



FFL_Z_reg <- function(graph, motif, a=0.5, Input=1){
  #motif must be a subgraph
  regX_Z <- half_loop_reg(E(graph, path=motif[[1]])$sign, Input=Input)
  regX_Y_Z <- half_loop_reg(E(graph, path=motif[[2]])$sign, Input=Input)
  return( a*(regX_Z+regX_Y_Z)  )
}


half_loop_reg <- function(sub_loop_sign, a=0.5, Input=1){
    if(length(sub_loop_sign)==1){
      regX_Z <- Input*sub_loop_sign[1]    }
  else if(length(sub_loop_sign)==2){
    regX_Z <- sub_loop_sign[2]*(a+(Input*sub_loop_sign[1]))  }
  else if(length(sub_loop_sign)==3){
    regX_Z <- (a+(sub_loop_sign[2]*(a+(Input*sub_loop_sign[1]))))*sub_loop_sign[3]  }
  else if(length(sub_loop_sign)==4){
    regX_Z <- (a+((a+(sub_loop_sign[2]*(a+(Input*sub_loop_sign[1]))))*sub_loop_sign[3]))*sub_loop_sign[4]  }
}
