# Script to "count" the known motifs in essentiel.topo output.
#Loop for each case of the tab ? 
#OR for each topo, go through different if to determine in which case it goes ? 


motif_tab1 <- function(list.topo){
  total <- length(list.topo)
  if(ncol(list.topo[[1]])==3){
  df.motif <- data.frame(row.names = c("FF_pos", "FF_neg", "FF_inc", "No_FF"), FB_neg=c(0,0,0,0),
                         FB_nul=c(0,0,0,0), FB_pos=c(0,0,0,0))
  for (W in list.topo){
    #look for the FB sign (W2_3 = W[3,2])
    FB <- ifelse(W[3,2]==0, 2, ifelse(W[3,2]==1, 3, 1 )) #give a number of column (only works for current df.motif !!)
    # browser()
    if(W[2,1]==0 || W[3,1]==0) df.motif[4, FB] <- as.numeric(df.motif[4, FB])+ 1/total*100 #no FF
      else if((W[2,1]==1 & W[3,1]==-1 & W[2,3]==1) || (W[2,1]==-1 & W[3,1]==-1 & W[2,3]==-1) ||
              (W[2,1]==-1 & W[3,1]==1 & W[2,3]==1) || (W[2,1]==1 & W[3,1]==1 & W[2,3]==-1)) df.motif[3, FB] <- df.motif[3, FB]+ 1/total*100 #Incoherent FF
      else if(W[2,1]==1) df.motif[1, FB] <- df.motif[1, FB]+ 1/total*100 #Positive FF
      else if(W[2,1]==-1) df.motif[2, FB] <- df.motif[2, FB]+ 1/total*100 #Negative FF
  }}
  #
  else if(ncol(list.topo[[1]])==4){
    df.motif <- data.frame(row.names = c("C_Diam_neg", "C_Diam_pos", "I_Diam_side", "I_Diam_alt", "No_diam"), FB_neg=c(0,0,0,0,0),
                           FB_nul=c(0,0,0,0,0), FB_pos=c(0,0,0,0,0))
    for (W in list.topo){
      #look for the FB sign (W2_3 = W[3,2])
      FB <- 2 #ifelse((W[3,2]==1 & W[2,3]==-1) || (W[4,2]==-1 & W[2,4]==1), 1, ifelse(W[3,2]==1 & W[3,1], 3, 1 )) #give a number of column (only works for current df.motif and topo sorting !!)
      # browser()
      if(W[3,1]==0 || W[4,1]==0 || W[2,3]==0 || W[2,4]==0) df.motif[5, FB] <- as.numeric(df.motif[5, FB])+ 1/total*100 #No Diamond
      else if(W[4,1]==-1 & W[3,1]==-1 & W[2,3]==-1 & W[2,4]==-1) df.motif[1, FB] <- df.motif[1, FB]+ 1/total*100 #Neg Diam
      else if(W[4,1]==1 & W[3,1]==1 & W[2,3]==1 & W[2,4]==1) df.motif[2, FB] <- df.motif[2, FB]+ 1/total*100 #Pos Diam
      else if(W[4,1]==1 & W[3,1]==-1 & W[2,3]==-1 & W[2,4]==1) df.motif[3, FB] <- df.motif[3, FB]+ 1/total*100 #Incoherent Diamond by side
      else if(W[4,1]==-1 & W[3,1]==1 & W[2,3]==-1 & W[2,4]==1) df.motif[4, FB] <- df.motif[4, FB]+ 1/total*100 #Incoherent Diamond, alternate sign
    }}
  return(df.motif)
}

motif_tab2 <- function(list.topo){
  #Function to count the frequency of topologies. No other choice than to do it case by case, i.e. 81 times.
  total <- length(list.topo)
  if(ncol(list.topo[[1]])!=3) stop("Not a 3 topology")
  df.motif <- data.frame(row.names = c("S_T_-1_TF_T_-1", "S_T_-1_TF_T_0","S_T_-1_TF_T_1",
                                       "S_T_0_TF_T_-1",  "S_T_0_TF_T_0","S_T_0_TF_T_1",
                                       "S_T_1_TF_T_-1","S_T_1_TF_T_0","S_T_1_TF_T_1"),
                         "S_TF_-1_T_TF_-1"=c(0,0,0,0,0,0,0,0,0), "S_TF_-1_T_TF_0"=c(0,0,0,0,0,0,0,0,0), "S_TF_-1_T_TF_1"=c(0,0,0,0,0,0,0,0,0),
                         "S_TF_0_T_TF_-1"=c(0,0,0,0,0,0,0,0,0),"S_TF_0_T_TF_0"=c(0,0,0,0,0,0,0,0,0),"S_TF_0_T_TF_1"=c(0,0,0,0,0,0,0,0,0),
                         "S_TF_1_T_TF_-1"=c(0,0,0,0,0,0,0,0,0),"S_TF_1_T_TF_0"=c(0,0,0,0,0,0,0,0,0),"S_TF_1_T_TF_1"=c(0,0,0,0,0,0,0,0,0))
  for (W in list.topo){

    if(W[2,1]==-1 & W[3,1]==-1 & W[2,3]==-1 & W[3,2]==-1) df.motif[1, 1] <- as.numeric(df.motif[1, 1])+ 1/total*100
    else if(W[2,1]==-1 & W[3,1]==-1 & W[2,3]==0 & W[3,2]==-1) df.motif[2, 1] <- as.numeric(df.motif[2, 1])+ 1/total*100 #no FF
    else if(W[2,1]==-1 & W[3,1]==-1 & W[2,3]==1 & W[3,2]==-1) df.motif[3, 1] <- as.numeric(df.motif[3, 1])+ 1/total*100 #no FF
    #
    else if(W[2,1]==0 & W[3,1]==-1 & W[2,3]==-1 & W[3,2]==-1) df.motif[4, 1] <- as.numeric(df.motif[4, 1])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==-1 & W[2,3]==0 & W[3,2]==-1) df.motif[5, 1] <- as.numeric(df.motif[5, 1])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==-1 & W[2,3]==1 & W[3,2]==-1) df.motif[6, 1] <- as.numeric(df.motif[6, 1])+ 1/total*100 #no FF
    #
    else if(W[2,1]==1 & W[3,1]==-1 & W[2,3]==-1 & W[3,2]==-1) df.motif[7, 1] <- as.numeric(df.motif[7, 1])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==-1 & W[2,3]==0 & W[3,2]==-1) df.motif[8, 1] <- as.numeric(df.motif[8, 1])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==-1 & W[2,3]==1 & W[3,2]==-1) df.motif[9, 1] <- as.numeric(df.motif[9, 1])+ 1/total*100 #no FF
    #Column 2
    else if(W[2,1]==-1 & W[3,1]==-1 & W[2,3]==-1 & W[3,2]==0) df.motif[1, 2] <- as.numeric(df.motif[1, 2])+ 1/total*100
    else if(W[2,1]==-1 & W[3,1]==-1 & W[2,3]==0 & W[3,2]==0) df.motif[2, 2] <- as.numeric(df.motif[2, 2])+ 1/total*100 #no FF
    else if(W[2,1]==-1 & W[3,1]==-1 & W[2,3]==1 & W[3,2]==0) df.motif[3, 2] <- as.numeric(df.motif[3, 2])+ 1/total*100 #no FF
    #
    else if(W[2,1]==0 & W[3,1]==-1 & W[2,3]==-1 & W[3,2]==0) df.motif[4, 2] <- as.numeric(df.motif[4, 2])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==-1 & W[2,3]==0 & W[3,2]==0) df.motif[5, 2] <- as.numeric(df.motif[5, 2])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==-1 & W[2,3]==1 & W[3,2]==0) df.motif[6, 2] <- as.numeric(df.motif[6, 2])+ 1/total*100 #no FF
    #
    else if(W[2,1]==1 & W[3,1]==-1 & W[2,3]==-1 & W[3,2]==0) df.motif[7, 2] <- as.numeric(df.motif[7, 2])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==-1 & W[2,3]==0 & W[3,2]==0) df.motif[8, 2] <- as.numeric(df.motif[8, 2])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==-1 & W[2,3]==1 & W[3,2]==0) df.motif[9, 2] <- as.numeric(df.motif[9, 2])+ 1/total*100 #no FF
    #Column 3
    else if(W[2,1]==-1 & W[3,1]==-1 & W[2,3]==-1 & W[3,2]==1) df.motif[1, 3] <- as.numeric(df.motif[1, 3])+ 1/total*100
    else if(W[2,1]==-1 & W[3,1]==-1 & W[2,3]==0 & W[3,2]==1) df.motif[2, 3] <- as.numeric(df.motif[2, 3])+ 1/total*100 #no FF
    else if(W[2,1]==-1 & W[3,1]==-1 & W[2,3]==1 & W[3,2]==1) df.motif[3, 3] <- as.numeric(df.motif[3, 3])+ 1/total*100 #no FF
    #
    else if(W[2,1]==0 & W[3,1]==-1 & W[2,3]==-1 & W[3,2]==1) df.motif[4, 3] <- as.numeric(df.motif[4, 3])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==-1 & W[2,3]==0 & W[3,2]==1) df.motif[5, 3] <- as.numeric(df.motif[5, 3])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==-1 & W[2,3]==1 & W[3,2]==1) df.motif[6, 3] <- as.numeric(df.motif[6, 3])+ 1/total*100 #no FF
    #
    else if(W[2,1]==1 & W[3,1]==-1 & W[2,3]==-1 & W[3,2]==1) df.motif[7, 3] <- as.numeric(df.motif[7, 3])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==-1 & W[2,3]==0 & W[3,2]==1) df.motif[8, 3] <- as.numeric(df.motif[8, 3])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==-1 & W[2,3]==1 & W[3,2]==1) df.motif[9, 3] <- as.numeric(df.motif[9, 3])+ 1/total*100 #no FF
    #Column 4
    if(W[2,1]==-1 & W[3,1]==0 & W[2,3]==-1 & W[3,2]==-1) df.motif[1, 4] <- as.numeric(df.motif[1, 4])+ 1/total*100
    else if(W[2,1]==-1 & W[3,1]==0 & W[2,3]==0 & W[3,2]==-1) df.motif[2, 4] <- as.numeric(df.motif[2, 4])+ 1/total*100 #no FF
    else if(W[2,1]==-1 & W[3,1]==0 & W[2,3]==1 & W[3,2]==-1) df.motif[3, 4] <- as.numeric(df.motif[3, 4])+ 1/total*100 #no FF
    #
    else if(W[2,1]==0 & W[3,1]==0 & W[2,3]==-1 & W[3,2]==-1) df.motif[4, 4] <- as.numeric(df.motif[4, 4])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==0 & W[2,3]==0 & W[3,2]==-1) df.motif[5, 4] <- as.numeric(df.motif[5, 4])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==0 & W[2,3]==1 & W[3,2]==-1) df.motif[6, 4] <- as.numeric(df.motif[6, 4])+ 1/total*100 #no FF
    #
    else if(W[2,1]==1 & W[3,1]==0 & W[2,3]==-1 & W[3,2]==-1) df.motif[7, 4] <- as.numeric(df.motif[7, 4])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==0 & W[2,3]==0 & W[3,2]==-1) df.motif[8, 4] <- as.numeric(df.motif[8, 4])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==0 & W[2,3]==1 & W[3,2]==-1) df.motif[9, 4] <- as.numeric(df.motif[9, 4])+ 1/total*100 #no FF
    #Column 5
    else if(W[2,1]==-1 & W[3,1]==0 & W[2,3]==-1 & W[3,2]==0) df.motif[1, 5] <- as.numeric(df.motif[1, 5])+ 1/total*100
    else if(W[2,1]==-1 & W[3,1]==0 & W[2,3]==0 & W[3,2]==0) df.motif[2, 5] <- as.numeric(df.motif[2, 5])+ 1/total*100 #no FF
    else if(W[2,1]==-1 & W[3,1]==0 & W[2,3]==1 & W[3,2]==0) df.motif[3, 5] <- as.numeric(df.motif[3, 5])+ 1/total*100 #no FF
    #
    else if(W[2,1]==0 & W[3,1]==0 & W[2,3]==-1 & W[3,2]==0) df.motif[4, 5] <- as.numeric(df.motif[4, 5])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==0 & W[2,3]==0 & W[3,2]==0) df.motif[5, 5] <- as.numeric(df.motif[5, 5])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==0 & W[2,3]==1 & W[3,2]==0) df.motif[6, 5] <- as.numeric(df.motif[6, 5])+ 1/total*100 #no FF
    #
    else if(W[2,1]==1 & W[3,1]==0 & W[2,3]==-1 & W[3,2]==0) df.motif[7, 5] <- as.numeric(df.motif[7, 5])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==0 & W[2,3]==0 & W[3,2]==0) df.motif[8, 5] <- as.numeric(df.motif[8, 5])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==0 & W[2,3]==1 & W[3,2]==0) df.motif[9, 5] <- as.numeric(df.motif[9, 5])+ 1/total*100 #no FF
    #Column 6
    else if(W[2,1]==-1 & W[3,1]==0 & W[2,3]==-1 & W[3,2]==1) df.motif[1, 6] <- as.numeric(df.motif[1, 6])+ 1/total*100
    else if(W[2,1]==-1 & W[3,1]==0 & W[2,3]==0 & W[3,2]==1) df.motif[2, 6] <- as.numeric(df.motif[2, 6])+ 1/total*100 #no FF
    else if(W[2,1]==-1 & W[3,1]==0 & W[2,3]==1 & W[3,2]==1) df.motif[3, 6] <- as.numeric(df.motif[3, 6])+ 1/total*100 #no FF
    #
    else if(W[2,1]==0 & W[3,1]==0 & W[2,3]==-1 & W[3,2]==1) df.motif[4, 6] <- as.numeric(df.motif[4, 6])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==0 & W[2,3]==0 & W[3,2]==1) df.motif[5, 6] <- as.numeric(df.motif[5, 6])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==0 & W[2,3]==1 & W[3,2]==1) df.motif[6, 6] <- as.numeric(df.motif[6, 6])+ 1/total*100 #no FF
    #
    else if(W[2,1]==1 & W[3,1]==0 & W[2,3]==-1 & W[3,2]==1) df.motif[7, 6] <- as.numeric(df.motif[7, 6])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==0 & W[2,3]==0 & W[3,2]==1) df.motif[8, 6] <- as.numeric(df.motif[8, 6])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==0 & W[2,3]==1 & W[3,2]==1) df.motif[9, 6] <- as.numeric(df.motif[9, 6])+ 1/total*100 #no FF
    #Column 7
    if(W[2,1]==-1 & W[3,1]==1 & W[2,3]==-1 & W[3,2]==-1) df.motif[1, 7] <- as.numeric(df.motif[1, 7])+ 1/total*100
    else if(W[2,1]==-1 & W[3,1]==1 & W[2,3]==0 & W[3,2]==-1) df.motif[2, 7] <- as.numeric(df.motif[2, 7])+ 1/total*100 #no FF
    else if(W[2,1]==-1 & W[3,1]==1 & W[2,3]==1 & W[3,2]==-1) df.motif[3, 7] <- as.numeric(df.motif[3, 7])+ 1/total*100 #no FF
    #
    else if(W[2,1]==0 & W[3,1]==1 & W[2,3]==-1 & W[3,2]==-1) df.motif[4, 7] <- as.numeric(df.motif[4, 7])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==1 & W[2,3]==0 & W[3,2]==-1) df.motif[5, 7] <- as.numeric(df.motif[5, 7])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==1 & W[2,3]==1 & W[3,2]==-1) df.motif[6, 7] <- as.numeric(df.motif[6, 7])+ 1/total*100 #no FF
    #
    else if(W[2,1]==1 & W[3,1]==1 & W[2,3]==-1 & W[3,2]==-1) df.motif[7, 7] <- as.numeric(df.motif[7, 7])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==1 & W[2,3]==0 & W[3,2]==-1) df.motif[8, 7] <- as.numeric(df.motif[8, 7])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==1 & W[2,3]==1 & W[3,2]==-1) df.motif[9, 7] <- as.numeric(df.motif[9, 7])+ 1/total*100 #no FF
    #Column 8
    else if(W[2,1]==-1 & W[3,1]==1 & W[2,3]==-1 & W[3,2]==0) df.motif[1, 8] <- as.numeric(df.motif[1, 8])+ 1/total*100
    else if(W[2,1]==-1 & W[3,1]==1 & W[2,3]==0 & W[3,2]==0) df.motif[2, 8] <- as.numeric(df.motif[2, 8])+ 1/total*100 #no FF
    else if(W[2,1]==-1 & W[3,1]==1 & W[2,3]==1 & W[3,2]==0) df.motif[3, 8] <- as.numeric(df.motif[3, 8])+ 1/total*100 #no FF
    #
    else if(W[2,1]==0 & W[3,1]==1 & W[2,3]==-1 & W[3,2]==0) df.motif[4, 8] <- as.numeric(df.motif[4, 8])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==1 & W[2,3]==0 & W[3,2]==0) df.motif[5, 8] <- as.numeric(df.motif[5, 8])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==1 & W[2,3]==1 & W[3,2]==0) df.motif[6, 8] <- as.numeric(df.motif[6, 8])+ 1/total*100 #no FF
    #
    else if(W[2,1]==1 & W[3,1]==1 & W[2,3]==-1 & W[3,2]==0) df.motif[7, 8] <- as.numeric(df.motif[7, 8])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==1 & W[2,3]==0 & W[3,2]==0) df.motif[8, 8] <- as.numeric(df.motif[8, 8])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==1 & W[2,3]==1 & W[3,2]==0) df.motif[9, 8] <- as.numeric(df.motif[9, 8])+ 1/total*100 #no FF
    #Column 9
    else if(W[2,1]==-1 & W[3,1]==1 & W[2,3]==-1 & W[3,2]==1) df.motif[1, 9] <- as.numeric(df.motif[1, 9])+ 1/total*100
    else if(W[2,1]==-1 & W[3,1]==1 & W[2,3]==0 & W[3,2]==1) df.motif[2, 9] <- as.numeric(df.motif[2, 9])+ 1/total*100 #no FF
    else if(W[2,1]==-1 & W[3,1]==1 & W[2,3]==1 & W[3,2]==1) df.motif[3, 9] <- as.numeric(df.motif[3, 9])+ 1/total*100 #no FF
    #
    else if(W[2,1]==0 & W[3,1]==1 & W[2,3]==-1 & W[3,2]==1) df.motif[4, 9] <- as.numeric(df.motif[4, 9])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==1 & W[2,3]==0 & W[3,2]==1) df.motif[5, 9] <- as.numeric(df.motif[5, 9])+ 1/total*100 #no FF
    else if(W[2,1]==0 & W[3,1]==1 & W[2,3]==1 & W[3,2]==1) df.motif[6, 9] <- as.numeric(df.motif[6, 9])+ 1/total*100 #no FF
    #
    else if(W[2,1]==1 & W[3,1]==1 & W[2,3]==-1 & W[3,2]==1) df.motif[7, 9] <- as.numeric(df.motif[7, 9])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==1 & W[2,3]==0 & W[3,2]==1) df.motif[8, 9] <- as.numeric(df.motif[8, 9])+ 1/total*100 #no FF
    else if(W[2,1]==1 & W[3,1]==1 & W[2,3]==1 & W[3,2]==1) df.motif[9, 9] <- as.numeric(df.motif[9, 9])+ 1/total*100 #no FF

  }
  return(df.motif)
}

motif_tab2(anticor_UD3)
motif_tab2(corr_UD3)

dd <- read.table("nettable.txt", header=TRUE, sep="\t")
rownames(dd) <- paste(dd[,1], dd[,2], dd[,3], dd[,4], sep="/")

mat2name <- function(mat) {
  paste(mat[2,1], mat[3,1], mat[2,3], mat[3,2], sep="/")
}


net.property <- function(corlist, property="Coherent.Net") {
  nn <- sapply(corlist, mat2name)
  dd[nn, property]
}

net.table <- function(corlist) {
  table(sapply(corlist, mat2name))
}

sapply(corr_UD3, mat2name)
corr_UD3.names <- sapply(corr_UD3, mat2name)
table(sapply(corr_UD3, mat2name))
summary(net.property(corr_UD3, property = "Coherent.Net"))
summary(net.property(corr_UD3, property = "Plastic"))
summary(net.property(corr_UD3, property = "Coherent.FB"))
summary(net.property(corr_UD3, property = "Coherent.FF"))

table(sapply(anticor_UD3, mat2name))
summary(net.property(anticor_UD3, property = "Coherent.Net"))
summary(net.property(anticor_UD3, property = "Plastic"))
summary(net.property(anticor_UD3, property = "Coherent.FB"))
summary(net.property(anticor_UD3, property = "Coherent.FF"))




pdf("figures/3g_AUD_random.pdf", width=6, height=6)
layout(matrix(c(1:9), 3, 3, byrow = TRUE))
par(mar=c(0.5, 0.5, 1, 0), mgp = c(1.75, 0.75, 0), las=0)
j<-1
for (i in unique(anticor_UD3)) {
  mainT <- length(which(sapply(1:length(anticor_UD3),function(x) length(which(paste0(anticor_UD3[[x]],
             collapse = "") == paste0(unique(anticor_UD3)[[j]], collapse = "")))) == 1))/length(anticor_UD3)*100
  #W matrix as a graph :
  G <- as.directed(graph.adjacency(t(i), weighted = T))
  V(G)$color <- c("green","orange", "yellow", "yellow")
  plot(G, layout=layout_in_circle, edge.color=ifelse(E(G)$weight > 0, "black","red" ), vertex.size=30, main=round(mainT, 3)) #, layout=layout_in_circle
  j <- j+1
  }
dev.off()



pdf("figures/3g_CUD_random.pdf", width=6, height=6)
layout(matrix(c(1:9), 3, 3, byrow = TRUE))
par(mar=c(0.5, 0.5, 1, 0), mgp = c(1.75, 0.75, 0), las=0)
j<-1
for (i in unique(corr_UD3)) {
  mainT <- length(which(sapply(1:length(corr_UD3),function(x) length(which(paste0(corr_UD3[[x]],
                                                                                     collapse = "") == paste0(unique(corr_UD3)[[j]], collapse = "")))) == 1))/length(corr_UD3)*100
  #W matrix as a graph :
  G <- as.directed(graph.adjacency(t(i), weighted = T))
  V(G)$color <- c("green","orange", "yellow", "yellow")
  plot(G, layout=layout_in_circle, edge.color=ifelse(E(G)$weight > 0, "black","red" ), vertex.size=30, main=round(mainT, 3)) #, layout=layout_in_circle
  j <- j+1
}
dev.off()


anticor_UD3[[1]]/length(anticor_UD3)

which(anticor_UD3[[1]] %in% anticor_UD3)




length(which(paste0(anticor_UD3, collapse = "") %in% paste0(anticor_UD3[[1]], collapse = "")))



