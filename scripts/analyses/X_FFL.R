#Custom motif categories : adapted to our simulation model with constitutive expression different from 0
#change everything to get a tab ? 5 column, "Loop","No loop", "Z_sign", "In/Coherence" and "True/False" and column of concatenate last 3 columns
FFL.from <- function(list.w, edges1=2, edges2=1, target=2, from=(1:ncol(list.w[[1]]))){
  if(edges1==2 && edges2==1){
    df <- data.frame()}
    for(i in 1:length(list.w)){
    g <- graph.adjacency(abs(t(list.w[[i]])), mode="directed") 
    E(g)$sign <- (list.w[[i]])[list.w[[i]] != 0] 
    if(length(feedforward.to(g, to=target, edges1=edges1, edges2=edges2))==0){ } else{
        
        ff <- feedforward.to(g, to=target, edges1=edges1, edges2=edges2)
        if(edges1==2 && edges2==1){
          for(m in ff) {
              if(length(E(g, path=m[[1]])$sign)>1 && length(E(g, path=m[[2]])$sign)>2) stop("Wrong cut.offs")
              if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]==E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==-1) motif <- "C_Ho_neg"
              if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]==E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==1) motif <- "C_Ho_pos"
              if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]!=E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==-1) motif <- "C_He_neg"
              if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]!=E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==1) motif <- "C_He_pos"
              if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==-1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]==E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==-1) motif <- "I_Ho_neg"
              if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==-1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]==E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==1) motif <- "I_Ho_pos"
              if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==-1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]!=E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==-1) motif <- "I_He_neg"
              if(prod(c(E(g, path=m[[1]])$sign,E(g, path=m[[2]])$sign))==-1 && E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)]!=E(g, path=m[[2]])$sign[length(E(g, path=m[[2]])$sign)] && (E(g, path=m[[1]])$sign[length(E(g, path=m[[1]])$sign)])==1) motif <- "I_He_pos"
              df <- rbind(df, data.frame(V1=colnames(list.w[[i]])[m[[1]][1]], V2=motif))
              }
          }}
    }
  return(df)
}

e_coli_prep_analyses <- function(genes_list, g, g_mat, fun="FFL", edges1=2, edges2=1, cores=2, from=FALSE){
  # edges1 is supposed to be the longer edge
  stopifnot(fun=="FFL" || fun=="FFL.from" || fun=="FBL")
  #here, create txt file, name=paste0(filename,"_",fun,".csv")
  loops <- mclapply(genes_list, function(gene) {
    #To avoid igraph::all_simple_paths to take days and weeks, we subset the regulatory networks. Only target genes and connected TFs are kept.
    gg <- induced.subgraph(g, vids = c(which(colnames(E_coli_mat)==gene), which(colnames(E_coli_mat)%in%TF_genes)) )
    E_coli_mat2 <- t((get.adjacency(gg, sparse=FALSE, attr='V3')))
    ggg <- induced.subgraph(gg, vids = as.vector(unlist(neighborhood(gg, edges1, nodes = which(colnames(E_coli_mat2)==gene), mode = 'all'))))
    E_coli_mat3 <- t((get.adjacency(ggg, sparse=FALSE, attr='V3')))
    E_coli_mat3 <- matrix(as.numeric(E_coli_mat3), ncol = ncol(E_coli_mat3), dimnames = dimnames(E_coli_mat3)) #convert to numeric matrix
    E_coli_mat3[is.na(E_coli_mat3)] <- 0
    if(is.character(from)) envirgenes <- which(colnames(E_coli_mat3)%in%from)
    else if(from==TRUE) envirgenes <- which(colnames(E_coli_mat3)%in%genes_list) else envirgenes <- 1:ncol(E_coli_mat3)
    print(paste0(gene," ; network size ", ncol(E_coli_mat3)))
    if(fun=="FFL") cc <- FFL.type2(list(E_coli_mat3), edges1 = edges1, edges2 = edges2, target = which(colnames(E_coli_mat3)==gene), from=envirgenes)
    if(fun=="FFL.from") cc <- FFL.from(list(E_coli_mat3), edges1 = edges1, edges2 = edges2, target = which(colnames(E_coli_mat3)==gene), from=envirgenes)
    if(fun=="FBL") cc <- FBL.type(list(E_coli_mat3), cutoff = cutoff, target = which(colnames(E_coli_mat3)==gene))
    #Here, insert writing c(gene,cc) in a txt file
    return(cc)
  }, mc.cores = cores)
  return(rbindlist(loops))
}



X_from <- e_coli_prep_analyses(c(all_plast_genes$V1, nonplast_genes), g, E_coli_mat, fun="FFL.from", edges1=2, edges2=1, cores=1, from = from)


#Doesn't work ; problem = loose the gene number (dumb)
df <- data.frame()
for (gene in unique(X_from$V1)) {
 df <- rbind.fill(df, as.data.frame.matrix(table(subset(X_from, V1==gene)))  )
}
df$gene <- unique(X_from$V1)
rownames(df) <- unique(X_from$V1)

barplot(t(df[,1:7]), col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
       legend=c( "C_Ho_neg","C_Ho_pos","C_He_neg","C_He_pos","I_Ho_neg","I_Ho_pos","I_He_neg","I_He_pos") )


pdfname <- "figures/E_coli"
pdf(paste0(pdfname,"_X_FFL",".pdf"), width=100, height=6)
layout(matrix(c(1:1), 1, 1, byrow = TRUE))
barplot(t(df[,1:7]), col=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"))
legend("bottomleft", box.lty=0,  bg="transparent", fill=c("darkseagreen","yellowgreen","dodgerblue","deepskyblue3","indianred1","lightpink","orange","lightgoldenrod1"),
       legend=colnames(df) )
dev.off()


df[is.na(df)] <- 0
rowSums(df[,1:7])*100/sum(df[,1:7])
#Works BUT: should keep gene name instead of number.

rowSums(df[,1:7])*100/colSums(df[,1:7])


df[,1]/sum(df[,1])
df[,2]/sum(df[,2])
df[,3]/sum(df[,3])
df[,4]/sum(df[,4])
df[,5]/sum(df[,5])








