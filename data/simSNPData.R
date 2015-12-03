nSnp <- 2000
p <- cbind(rbeta(nSnp , 50, 50), rbeta(nSnp, 5, 5))
pops <- rep(c("Control", "Treat"), times = c(53, 51))
nSamp <- data.frame(Control = table(pops)[["Control"]],
                    Treat = table(pops)[["Treat"]])


mySnpSampler <- function(p, na.prob=0.01, .nSamp = nSamp){
  alleles <- rbind(
    matrix(sample(c("A","B"), 
                  2*nSamp$Control,
                  replace = TRUE,
                  c(p[1], 1-p[1])),
           ncol=2),
    matrix(sample(c("A","B"), 
                  2*nSamp$Treat,
                  replace = TRUE,
                  c(p[2], 1-p[2])),
           ncol=2)
  )
  alleles <- apply(alleles, MARGIN=1, paste, collapse="") 
  naRep <- sample(c(TRUE, FALSE), length(alleles), replace=TRUE, prob = c(na.prob, 1-na.prob))
  alleles[naRep] <- NA
  alleles[alleles == "BA"] <- "AB"
  alleles
}

genotypes <- apply(p, MARGIN=1, mySnpSampler)
colnames(genotypes) <- paste0("SNP", 1:nSnp)
genotypes <- data.frame(Population = pops, genotypes)
write.csv(genotypes, 
          file.path("Day 2", "session 4 Digging Deeper in R", "data", "snps.csv"),
          row.names = FALSE)
