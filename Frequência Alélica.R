library(adegenet)
library(vcfR)

#importar vcf pro R
vcf <- read.vcfR(file.choose())

#meta data, fixed data e genotype data
#genotype é o que tá importando aqui, fixed é informação do cromossomo e meta é meta

#extrair genotype data em dataframe
gt <- extract.gt(vcf, element = "GT")


#transpor para depois adicionar populações
transposto <- as.data.frame(t(gt))
transposto <- tibble::rownames_to_column(transposto,"ID")

View(GENEINPUT)

  # Subset the data for the current popCode value
GENEINPUT[GENEINPUT == 'IBS,IBS'] <- "IBS"
GENEINPUT[GENEINPUT == 'IBS,MSL'] <- "IBS"
GENEINPUT[GENEINPUT == '1|0'] <- "0 1"
GENEINPUT[GENEINPUT == '1|1'] <- "1 1"
GENEINPUT[GENEINPUT == '0|1'] <- "0 1"
GENEINPUT[GENEINPUT == '0|0'] <- "0 0"

  popData <- subset(GENEINPUT, popCode == popCode)
  View(popData)
  df <- arrange(popData, AFR, EUR, AMR, EAS, SAS)
  
  # Convert the data to genind format
  popGenind <- df2genind(X = popData[, -c(1, 2, 3)], pop= GENEINPUT$popCode, sep = " ", ploidy = 2, ncode = 2)
  
  # Convert the genind object to a genpop object
  popGenpop <- genind2genpop(popGenind)
  
  # Compute the frequency table
  popFreq <- makefreq(popGenpop, quiet = FALSE, missing = NA, truenames = TRUE)
  View(popFreq)
  # Print the frequency table to the console
  cat("\n\n")
  print(paste("Frequency table for population", popCode))
  print(popFreq)
  
  
  #tentar separar com linha nula as superpopulações
  popFreq2 <- popFreq
 popFreq2 <- rbind(popFreq2[1:26, ], X, popFreq2[- (1:), ])
 View(popFreq2)

 data_new <- rbind(data[1:3, ],            # new_row,
                   data[- (1:3), ])
 data_new 
 View(popFreq2)

X <- c()
populations_order <- c("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI", X, "CLM", "MXL", "PEL",
                       "PUR", X, "CDX", "CHB", "CHS", "JPT", "KHV", X, "CEU", "FIN", "GBR",
                       "IBS", "TSI", X, "BEB", "GIH", "ITU", "PJL", "STU")
  popFreq2 <- popFreq2[populations_order, ]
View(popFreq2)
    fileName <- paste0(popCode, "_freq.csv", sep = c())
  write.table(popFreq, file = fileName, sep = ",", col.names = TRUE, quote = FALSE)
  
  # Write the frequency table to an Excel file
  library(openxlsx)
  excelFileName <- paste0(popCode, "_freq.xlsx", sep = c())
  write.xlsx(popFreq, file = excelFileName, col.names = TRUE)
  View(popData)