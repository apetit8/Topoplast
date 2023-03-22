source("functions/functions.R")
library(igraph)

test <- netw.group(df.b[1:5,], target=2)

plot(test[[1]], edge.width=abs(E(test[[1]])$weight)*2, edge.color=ifelse(E(test[[1]])$weight > 0, "blue","red" ))

topo_sort(test[[1]], mode = c("in"))


as_adj(test[[1]], attr="weight")

################################################################################

test <- netw.group(subset(df.b, Gen==max(df.b$Gen)), target=2)


simul.topo <- data.frame()

for(i in 1:length(test)){
untopo <-unique.topo(test[[i]], groups=list(1,2,3:ncol(test[[i]]))) #group the topology who are exancheable, start at gene 3 since i'm only looking at submodules

simul.topo[i,1] =  paste(as.character(untopo), collapse ="")
}

length(unique(simul.topo[,1])) 
#Doesn't work. Gives as many topology as I have W ...




#Keep the same gene number every time ? Does not work either ...
netw.group <- function(df, Ss=1, Pg=1, Sg=1, Tf=7, start=7, target=2){
  # Ss = number of sensor genes ; Pg = plastic genes ; Sg = stables genes ; Tf = transcription factors
  # start = column in df from which the network begun
  #Target : 2="Changing gene" ; 33="Constant" gene ; 1="Signal" gene
  file.results <- lapply(1:nrow(df), function(i) { #mclapply(files, function(ff) {
    W <- t(matrix(as.numeric(df[i,start:((Ss+Tf+Pg+Sg)*(Ss+Tf+Pg+Sg)+start-1)]), ncol = (Ss+Tf+Pg+Sg)))
    return(W)
  })
  return(file.results)
}

netw.group <- function(df, Ss=1, Pg=1, Sg=1, Tf=7, start=7, target=2){
  # Ss = number of sensor genes ; Pg = plastic genes ; Sg = stables genes ; Tf = transcription factors
  # start = column in df from which the network begun
  #Target : 2="Changing gene" ; 33="Constant" gene ; 1="Signal" gene
  file.results <- lapply(1:nrow(df), function(i) { #mclapply(files, function(ff) {
    W <- t(matrix(as.numeric(df[i,start:((Ss+Tf+Pg+Sg)*(Ss+Tf+Pg+Sg)+start-1)]), ncol = (Ss+Tf+Pg+Sg)))
    #W matrix as a graph : 
    G <- as.directed(graph.adjacency(t(W), weighted = T))
    G <- delete.edges(G, E(G)[ abs(weight) < median(abs(weight)) ]) #delete.edges(G, E(G)[ abs(weight) < 0.1 ]) ; mean(abs(weight))
    GG <- induced_subgraph(G, c(1,target, neighbors(G,target, mode="all")))
    GG <- t(as.matrix(as_adj(GG, attr="weight")))
    
    return(GG)
  })
  return(file.results)
}

test <- netw.group(subset(df.b, Gen==max(df.b$Gen)), target=2)

simul.topo <- data.frame()
for(i in 1:length(test)){
  untopo <-unique.topo(test[[i]], groups=list(1,2,3:ncol(test[[i]]))) #group the topology who are exancheable, start at gene 3 since i'm only looking at submodules
  
  simul.topo[i,1] =  paste(as.character(untopo), collapse ="")
}
length(unique(simul.topo[,1])) 
#Doesn't work either

################################################################################
df <- subset(df.b, Gen==max(df.b$Gen) & b_envir=="Correlated")

#Keep "essential" connections ? Test every connection to see impact on target gene RN, and draw from that.
#Inspired by Burda et al., 2011
min <- 0.15
max <- 0.85
target <- 3
treshold_coeff <- 0.05
treshold_og <- 0.01
#How ? 
essential_Ws <- lapply(1:nrow(df), function(i){ #df : 
  W <-  t(matrix(as.numeric(df[i,7:106]), ncol = 10))
  W2 <- W
  RN_W_coeff <- getSlope.ALR(W=W, n.env=30, target.gene=target, min=min, max=max)
  RN_W_og <- getSlope.ALR(W=W, n.env=30, target.gene=target, min=min, max=max, giveback=1)
  for(Wij in 1:length(W)){  #(but diagonal)
    W_test <- W
    W_test[Wij] <- 0
    RN_Wij_coeff <- getSlope.ALR(W=W_test, n.env=30, target.gene=target, min=min, max=max)
    RN_Wij_og <- getSlope.ALR(W=W_test, n.env=30, target.gene=target, min=min, max=max, giveback=1)
    W2[Wij] <- ifelse(RN_Wij_coeff >= (RN_W_coeff-treshold_coeff) &
                        RN_Wij_coeff <= (RN_W_coeff+treshold_coeff) &
                        RN_Wij_og >= (RN_W_og-treshold_og)&
                        RN_Wij_og <= (RN_W_og+treshold_og), 0, 1) * sign(W2[Wij])
  }
return(W2)})

simul.topo <- data.frame()
for(i in 1:length(essential_Ws)){
  untopo <-unique.topo(essential_Ws[[i]], groups=list(1,2,3,4:ncol(essential_Ws[[i]]))) #group the topology who are exancheable, start at gene 3 since i'm only looking at submodules
  
  simul.topo[i,1] =  paste(as.character(untopo), collapse =",")
}
length(unique(simul.topo[,1])) 

#Now draw the unique combinations
untopo <- unique(simul.topo[,1])
tt <- unlist(str_split(string=untopo[12], pattern= ","))
W_topo <- matrix(as.numeric(tt), ncol = 10, byrow = FALSE)
#W matrix as a graph : 
G <- as.directed(graph.adjacency(t(W_topo), weighted = T))
V(G)$color <- c("darkred","orange", "green","yellow",  "yellow", "yellow", "yellow", "yellow", "yellow", "yellow")

plot(G, edge.color=ifelse(E(G)$weight > 0, "black","red" ))


#WORKS !! Youhouuuuu
#Results as predicted :')

#Would it be possible to extract the common regulations between all networks ?




