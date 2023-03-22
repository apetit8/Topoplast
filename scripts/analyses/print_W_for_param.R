#W matrix construction site

#Network proprieties
Sensor <- 1 #Can't be regulated but can regulate all other genes
TF <- 8 #Regulated and regulate
Plastic <- 4 #Regulated
TG <- 8 #Regulated
Total <- Sensor+TF+Plastic+TG

W <- matrix(data=0, nrow=Total, ncol = Total)
cat(as.character(t(W)), file="../templates/Matrix_W.txt", sep=" ", append=FALSE)
#Columns
W[,c((Sensor+1):(Sensor+TF))] <- 1
#Rows
W[c((Sensor+1):(Sensor+TF)),c(1:Sensor)] <- 1
W[c((Sensor+TF+1):Total),c(1:(Sensor+TF))] <- 1
W[c(1:Sensor),] <- 0
#Diagonal
diag(W) <- 0

#Print
# cat(as.character(t(W)), file="../templates/Matrix_W.txt", sep=" ", append=FALSE)
W <- ifelse(W==0,"immut","normal")
cat(as.array(t(W)), file="../templates/Matrix_W.txt", sep=" ", append=TRUE)


################################################################################
#PHASE 2

#Network proprieties
Sensor <- 1 #Can't be regulated but can regulate TF
TF <- 7 #Regulated and regulate
Plastic <- 1 #Regulated and regulate TF
TG <- 1 #Regulated and regulate TF
Total <- Sensor+TF+Plastic+TG

W <- matrix(data=0, nrow=Total, ncol = Total)
cat(as.character(t(W)), file="../templates/Matrix_W.txt", sep=" ", append=FALSE)
#Columns
W[,c((Sensor+1):(Sensor+TF))] <- 1
#Rows
W[c((Sensor+1):(Sensor+TF)),c(1:Sensor)] <- 1
W[c((Sensor+TF+1):Total),c(1:(Sensor+TF))] <- 1
W[c(1:(Sensor+TF)),c((Sensor+TF+1):Total)] <- 1
W[c(1:Sensor),] <- 0
W[c((Sensor+TF):Total),1:Sensor] <- 0
#Diagonal
diag(W) <- 0

W

#Print
# cat(as.character(t(W)), file="../templates/Matrix_W.txt", sep=" ", append=FALSE)
W <- ifelse(W==0,"immut","normal")
cat(as.array(t(W)), file="../templates/Matrix_W.txt", sep=" ", append=TRUE)




