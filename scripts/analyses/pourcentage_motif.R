# Script to "count" the known motifs in essentiel.topo output.
#Loop for each case of the tab ? 
#OR for each topo, go through different if to determine in which case it goes ? 


motif_tab <- function(list.topo){
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


motif_tab(anticor_UD3)


motif_tab(corr_UD3)


motif_tab(anticor_UD4)


motif_tab(corr_UD4)








