#======================================================
# Calculo da Frequência Alélica de VCF
#=======================================================
#
# (C) Copyright 2024, by GP-PGx-UFTM and Contributors.
#
# 
#-----------------
#  
#-----------------
#
# Original Author: Guilherme Belfort Almeida, Caique Manochio
# Contributor(s):
# Updated by (and date): Guilherme Belfort Almeida 28/02/2024
#
# Dependencies: R
#
# Command line:



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
GENEDADO <- select(ID, Population.code, Population.name, everything()) %>%
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

View(GENEDADO)

  popData <- subset(GENEDADO, popCode == popCode)
  
  # Converter dados para genind
  popGenind <- df2genind(X = popData[, -c(1, 2, 3)], pop= GENEDADO$popCode, sep = " ", ploidy = 2, ncode = 2)
  
  # Converter genind em genpop
  popGenpop <- genind2genpop(popGenind)
  
  # Fazer a frequência da tabela
  popFreq <- makefreq(popGenpop, quiet = FALSE, missing = NA, truenames = TRUE)
  
  #transformar em dataframe próprio
  df <- as.data.frame(popFreq)
  

  
  # Colocar a tabela de frequência num arquivo excel
  library(writexl)
  write.xlsx(df, "GENE.xlsx", col.names = TRUE, row.names = TRUE)
  
