#======================================================
# Calculo de correlação entre frequencia alélica e media de ancestralidade 
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
# Updated by (and date): Guilherme Belfort Almeida 05/02/2024
#
# Dependencies: R
#
# Command line:

library(clipr)
library(readxl)

#Lê o Excel da Frequência Alélica e cria uma lista com os SNPs
df <- read.xlsx(file.choose())
lista_snps <- unique(grep("^.+(.1)$",names(df), value=TRUE, perl = TRUE))



calcular_correlacoes <- function(arquivo_excel, lista_snps) {
  

  # cria uma lista vazia para armazenar os resultados
  resultados <- list()
  
  # loop pelos SNPs na lista
  for (snp in lista_snps) {
    # converte a coluna do SNP para numérico
    df[[snp]] <- as.numeric(df[[snp]])
    
    # realiza as correlações com cada ancestralidade
    EAS <- cor.test(df[[snp]], df$EAS, method = "pearson")
    EUR <- cor.test(df[[snp]], df$EUR, method = "pearson")
    NAT <- cor.test(df[[snp]], df$NAT, method = "pearson")
    AFRL <- cor.test(df[[snp]], df$AFRL, method = "pearson")
    EAS2 <- cor.test(df[[snp]], df$EAS2, method = "pearson")
    SAS <- cor.test(df[[snp]], df$SAS, method = "pearson")
    AFR <- cor.test(df[[snp]], df$AFR, method = "pearson")
    
    # armazena os valores de p e R em um data frame
    valores_p_r <- data.frame(
      SNP = lista_snps <- gsub('\\.1','',snp),
      Ancestralidade = c("EAS", "EUR", "NAT", "AFRL", "EAS2", "SAS", "AFR"),
      Valor_p = c(EAS$p.value, EUR$p.value, NAT$p.value, AFRL$p.value, EAS2$p.value, SAS$p.value, AFR$p.value),
      Valor_r = c(EAS$estimate, EUR$estimate, NAT$estimate, AFRL$estimate, EAS2$estimate, SAS$estimate, AFR$estimate)
    )
    
    # adiciona o data frame à lista de resultados
    resultados[[snp]] <- valores_p_r
  }
  
  # retorna a lista de resultados
  return(resultados)
}

resultados <- calcular_correlacoes(#"Caminho+nomeoutputaqui",
  lista_snps)

CorLista <- dplyr::bind_rows(resultados)
# exporta os resultados para um arquivo Excel
write.xlsx(CorLista, #"Caminho+nomeoutputaqui2",
           rowNames = FALSE)
