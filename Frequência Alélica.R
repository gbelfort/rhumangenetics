library('vcfR')
library('adegenet')
library('tibble')
library('dplyr')

#importar vcf pro R
vcf <- read.vcfR(file.choose())

#meta data, fixed data e genotype data
#genotype é o que tá importando aqui, fixed é informação do cromossomo e meta é meta

#extrair genotype data em dataframe
gt <- extract.gt(vcf, element = "GT")


#transpor para depois adicionar populações
transposto <- as.data.frame(t(gt))
transposto <- tibble::rownames_to_column(transposto,"ID")

#importar correspondência indivíduo-população
dados_pop_1kgp <- read.csv(file.choose(), sep = "\t")
dados_pop_1kgp <- dados_pop_1kgp %>%
  rename(ID = Sample.name)
#organizar por ID
Pops1KG <- dados_pop_1kgp[order(dados_pop_1kgp$ID),]


merge <- left_join(transposto, Pops1KG, by= "ID")

GENEDADO <- merge%>%
  select(ID, Population.code, Population.name, everything()) %>%
  select(-c("Sex","Biosample.ID","Superpopulation.code","Superpopulation.name",
            "Population.elastic.ID", "Data.collections"))

# renomear Population.code e Population.name para popCode e popName. adegenet não funciona se tiver nome de coluna com '.'
names(GENEDADO)[names(GENEDADO) == "Population.code"] <- "popCode"
names(GENEDADO)[names(GENEDADO) == "Population.name"] <- "popName"
# Um indivíduo com pop.code IBS,IBS. Isso seria interpretado como uma população diferente já que uso popCode como base, mudei para IBS.
GENEDADO[GENEDADO == 'IBS,IBS'] <- "IBS"
GENEDADO[GENEDADO == 'IBS,MSL'] <- "IBS"
# Trocar "|" por " "
GENEDADO[] <- lapply(GENEDADO, function(x) gsub("\\|", " ", x))
#ignorar isso abaixo
GENEDADO$Phenotype2 <- GENEDADO$Phenotype
View(GENEDADO)

  popData <- subset(GENEDADO, popCode == popCode)
  
  # Converter dados para genind
  popGenind <- df2genind(X = popData[, -c(1, 2, 3)], pop= GENEDADO$popCode, sep = " ", ploidy = 2, ncode = 2)
  
  # Converter genind em genpop
  popGenpop <- genind2genpop(popGenind)
  
  # Fazer a frequência da tabela
  popFreq <- makefreq(popGenpop, quiet = FALSE, missing = NA, truenames = TRUE)
  popFreq<- as.data.frame(popFreq)
  
  # Mostrar a frequência no Console
  cat("\n\n")
  print(paste("Frequency table for population", popCode))
  print(popFreq)
  
  # Colocar a tabela de frequência num arquivo .CSV
  fileName <- paste0("_freq.csv", sep = "")
  write.table(popFreq, file = fileName, sep = ",", col.names = TRUE, quote = FALSE)

  library(writexl)
  write_xlsx(popFreq, "popFreq.xlsx")
  
  library(writexl)
  # Colocar a tabela de frequência num arquivo Excel
  library(openxlsx)
  excelFileName <- paste0(popFreq, "_freq.xlsx", sep = "")
  write.xlsx(popFreq, file = excelFileName, colNames = TRUE)
  View(popData)
  
  write.xlsx(popFreq, "C:/Users/guilh/OneDrive/Área de Trabalho/freqSLC1kg.xlsx", rowNames = TRUE)
  SLCOr<- as.data.frame(popFreq)
  
