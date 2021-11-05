library(data.table)

sample_data <- fread("data/sample_key.csv")

chisq_table <- sample_data[,c(4,5,7)]
chisq_table <- chisq_table[chisq_table$Behaviour=="E",]
chisq_table <- chisq_table[,-c(2)]

##If I filter for only Evasive samples and
##test for an interaction between location and viral expression - does sig result mean
##interaction between location and parasitism attack success??
##Has lincoln ever looked at what % of attacked ASW are parasitised?
chi_sq_res <- list(chisq.test(table(chisq_table)))
