#=========================================================================
# Calculo de regressao entre frequencia alelica e media de ancestralidade 
#=========================================================================
#
# (C) Copyright 2024, by GP-PGx-UFTM and Contributors.
#
# 
#-----------------
#  
#-----------------
#
# Original Author: Guilherme Belfort-Almeida and Fernanda Rodrigues-Soares
# Contributor(s):  
# Updated by (and D/M/Y): Guilherme Belfort Almeida 07/12/2024
#
# Dependencies: R
#
# Command line:

library(openxlsx)
library(dplyr)

# Lê o arquivo da Frequência Alélica e cria uma lista com os SNPs
df <- read.xlsx(file.choose())

# Limpeza das colunas no dataframe
nomes_colunas <- names(df)
nomes_colunas <- gsub("\\.1$", "", nomes_colunas)
colunas_para_excluir <- grep("\\.0$", nomes_colunas, value = TRUE)
df <- df[, !names(df) %in% colunas_para_excluir]
nomes_colunas <- gsub("\\.1$", "", nomes_colunas)
names(df) <- nomes_colunas

# Lista de SNPs
lista_snps <- setdiff(names(df), c("POP", "NAT2", "AFRL", "EASL", "SAS", 
                                   "NAT", "EURN", "AFRO", "EURS", "EASO"))


  
calcular_correlacoes <- function(arquivo_excel, lista_snps) {
  resultados <- list()
  
  for (snp in lista_snps) {
    # Verifica se o SNP é numérico e converte se necessário
    if (!is.numeric(df[[snp]])) {
      df[[snp]] <- as.numeric(df[[snp]])
    }
    
    # Cria data frame das ancestralidades e assegura tipos numéricos
    media_k9 <- data.frame(NAT2 = as.numeric(df$NAT2), 
                           AFRL = as.numeric(df$AFRL), 
                           EASL = as.numeric(df$EASL), 
                           SAS = as.numeric(df$SAS),
                           NAT = as.numeric(df$NAT), 
                           EURN = as.numeric(df$EURN), 
                           AFRO = as.numeric(df$AFRO), 
                           EURS = as.numeric(df$EURS), 
                           EASO = as.numeric(df$EASO))
    lista_ances <- names(media_k9)
    
    resultados_df <- data.frame(SNP = character(), 
                                Beta = numeric(), 
                                p_regr = numeric(),
                                p_adj = numeric(),
                                R_sqr_ajustado = numeric(), 
                                Ancestralidade = character(),
                                stringsAsFactors = FALSE)
    
    p_values <- numeric() # Lista para armazenar p-values antes da correção
    
    for (nome in lista_ances) {
      reg_model <- lm(df[[snp]] ~ media_k9[[nome]])
      summary_reg <- summary(reg_model)
      
      # Verifica se há coeficientes válidos
      if (nrow(summary_reg$coefficients) > 1) {
        beta <- summary_reg$coefficients[2, "Estimate"]
        p_value <- summary_reg$coefficients[2, "Pr(>|t|)"]
        adjusted_r_squared <- summary_reg$adj.r.squared
      } else {
        beta <- NA
        p_value <- NA
        adjusted_r_squared <- NA
      }
      
      p_values <- c(p_values, p_value) # Armazena p-values
      
      resultados_df <- rbind(resultados_df, 
                             data.frame(SNP = snp, 
                                        Ancestralidade = nome,
                                        Beta = beta,
                                        p_adj = NA, # Placeholder
                                        R_sqr_ajustado = adjusted_r_squared,
                                        p_regr = p_value))
    }
    
    # Bonferroni correction
    p_adj_values <- p.adjust(p_values, method = "bonferroni")
    resultados_df$p_adj <- p_adj_values
    
    resultados[[snp]] <- resultados_df
  }
  
  return(resultados)
}

# Chama a função e armazena o resultado
resultados <- calcular_correlacoes(arquivo_excel = "caminho_para_o_arquivo.xlsx", lista_snps = lista_snps)

# Combina os dataframes individuais em um único dataframe
Regress <- dplyr::bind_rows(resultados)

# Exporta os resultados para um arquivo Excel
write.xlsx(Regress, "C:/path", rowNames = FALSE)

