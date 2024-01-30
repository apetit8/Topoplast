ec_cyc <- read.csv("e_coli/ECOLI-regulatory-network_2024_01_29.csv", sep = ",", header = FALSE)
ec_cyc <- subset(ec_cyc, V2 != "")

for (i in 1:nrow(ec_cyc)) {
  ec_cyc[i,3] <- ifelse(grepl("+/-", ec_cyc[i,2], fixed = TRUE), 0, ifelse(grepl( "+", ec_cyc[i,2], fixed = TRUE), 1, -1))
}
#Remove the reg sign from second column
ec_cyc$V2 <- gsub("+/-","",ec_cyc$V2)
ec_cyc$V2 <- gsub("\\+","",ec_cyc$V2)
ec_cyc$V2 <- gsub("-","",ec_cyc$V2)


write.csv(ec_cyc, "e_coli/ECOLI-regulatory-network_cyc_editd_2024_01_29.csv", row.names = FALSE)



