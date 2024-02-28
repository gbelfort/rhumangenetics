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
# Original Author: Guilherme Belfort Almeida
# Contributor(s):  
# Updated by (and date): Guilherme Belfort Almeida 28/02/2024
#
# Dependencies: R
#
# Command line:

library(openxlsx)
library(clipr)

#Lê o arquivo da Frequência Alélica e cria uma lista com os SNPs
df <- read.xlsx(file.choose())
lista_snps.1 <- unique(grep("^.+(.1)$",names(df), value=TRUE, perl = TRUE))
lista_snps.2 <- unique(grep("^.+(.2)$",names(df), value=TRUE, perl = TRUE))
lista_snps.3 <- unique(grep("^.+(.3)$",names(df), value=TRUE, perl = TRUE))
lista_snps.4 <- unique(grep("^.+(.4)$",names(df), value=TRUE, perl = TRUE))
lista_snps <- c(lista_snps.1, lista_snps.2, lista_snps.3, lista_snps.4)

calcular_correlacoes <- function(arquivo_excel, lista_snps) {
  

  # cria uma lista vazia para armazenar os resultados
  resultados <- list()
  
  # loop pelos SNPs na lista
  for (snp in lista_snps) {
    # converte a coluna do SNP para numérico
    df[[snp]] <- as.numeric(df[[snp]])
    
    # realiza as correlações com cada ancestralidade
    N_EUR <- cor.test(df[[snp]], df$N_EUR, method = "pearson")
    W_AFR <- cor.test(df[[snp]], df$W_AFR, method = "pearson")
    NAT <- cor.test(df[[snp]], df$NAT, method = "pearson")
    EAS_W <- cor.test(df[[snp]], df$EAS_W, method = "pearson")
    S_EUR <- cor.test(df[[snp]], df$S_EUR, method = "pearson")
    E_AFR <- cor.test(df[[snp]], df$E_AFR, method = "pearson")
    EAS_E <- cor.test(df[[snp]], df$EAS_E, method = "pearson")
    SAS <- cor.test(df[[snp]], df$SAS, method = "pearson")
    
    # armazena os valores de p e R em um data frame
    valores_p_r <- data.frame(
      SNP = lista_snps <- gsub('\\.1','',snp),
      Ancestralidade = c("N_EUR", "W_AFR", "NAT", "EAS_W", "S_EUR", "E_AFR", "EAS_E","SAS"),
      Valor_p = c(N_EUR$p.value, W_AFR$p.value, NAT$p.value, EAS_W$p.value,
                  S_EUR$p.value, E_AFR$p.value, EAS_E$p.value, SAS$p.value),
      Valor_r = c(N_EUR$estimate, W_AFR$estimate, NAT$estimate, EAS_W$estimate,
                  S_EUR$estimate, E_AFR$estimate, EAS_E$estimate, SAS$estimate)
    )
    
    # adiciona o data frame à lista de resultados
    resultados[[snp]] <- valores_p_r
  }
  
  # retorna a lista de resultados
  return(resultados)
}

CorLista <- dplyr::bind_rows(resultados)
# exporta os resultados para um arquivo Excel
write.xlsx(CorLista, "caminho",
           rowNames = FALSE)
