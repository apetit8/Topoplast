################################################################################
#Network proprieties
Sensor <- 1 #Can't be regulated but can regulate TF
TF <- 10 #Regulated and regulate
Plastic <- 10 #Regulated and regulate TF
TG <- 10 #Regulated and regulate TF
Control <- 5 #Regulated
Total <- Sensor+TF+Plastic+TG+Control


#A ACTUALISER
W <- matrix(data=0, nrow=Total, ncol = Total)
cat(as.character(t(W)), file="templates/Matrix_W.txt", sep=" ", append=FALSE)
#Columns
W[,c((Sensor):(Sensor+TF))] <- 1
#Rows
W[c((Sensor+1):(Sensor+TF)),c(1:Sensor)] <- 1
W[c((Sensor+TF+1):Total),c(1:(Sensor+TF))] <- 1
W[c(1:(Sensor+TF)),c((Sensor+TF+1):Total)] <- 1
W[,c((Total-Control):Total)] <- 0
W[c(1:Sensor),] <- 0
#Diagonal
diag(W) <- 0

W

#Print
# cat(as.character(t(W)), file="../templates/Matrix_W.txt", sep=" ", append=FALSE)
W <- ifelse(W==0,"immut","normal")
cat(as.array(t(W)), file="templates/Matrix_W.txt", sep=" ", append=TRUE)
