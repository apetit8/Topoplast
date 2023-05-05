ec_cyc <- read.csv("e_coli/ECOLI-regulatory-network_cyc.csv", sep = " ", header = FALSE)
ec_cyc <- subset(ec_cyc, V2 != "")

for (i in 1:nrow(ec_cyc)) {
  ec_cyc[i,3] <- ifelse(grepl("+/-", ec_cyc[i,2], fixed = TRUE), 5, ifelse(grepl( "+", ec_cyc[i,2], fixed = TRUE), 1, -1))
}

write.csv(ec_cyc, "e_coli/ECOLI-regulatory-network_cyc_editd.csv", row.names = FALSE)
#In .csv, removed all "+" and "-" and "/" from the second column of the dataset

ec_cyc <- read.csv("e_coli/ECOLI-regulatory-network_cyc_editd.csv")





