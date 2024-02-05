library(readxl) # Carrega o pacote readxl para ler arquivos do Excel
library(dplyr) # Carrega o pacote dplyr para manipular dados
library(tools) # Carrega o pacote tools para trabalhar com nomes de arquivos

#Input
InputGene <- "C:/Users/guilh/OneDrive/Área de Trabalho/Tuberculose/XPO1.xlsx" # Define o caminho do arquivo de entrada
#Ou escolher na pasta InputGene <- file.choose() # Ou permite ao usuário escolher o arquivo de entrada

# Diretório Output (se não existir o caminho ele é criado sozinho aqui)
directory <- "C:/Users/guilh/OneDrive/Área de Trabalho/Tuberculose/pvalor" # Define o caminho do diretório de saída

nome_gene <- file_path_sans_ext(basename(InputGene)) # Extrai o nome do gene do nome do arquivo de entrada
# Função para fazer os testes de forma automática
fun_qui_quadrado <- function(planilha, sheet_name) { # Define uma função que recebe uma planilha e um nome de planilha como argumentos
  colunas <- c(1, 4, 7, 8) # Define um vetor com os índices das colunas de interesse
  nomes_colunas <- colnames(planilha)[colunas] # Extrai os nomes das colunas de interesse
  
  planilha <- planilha[, colunas] # Filtra a planilha para manter apenas as colunas de interesse
  nome_pop <- colnames(planilha)[1] # Extrai o nome da primeira coluna
  col_indices <- c(1, 2, 3, 4) # Define um vetor com os índices das colunas a serem renomeadas
  new_names <- c("...1", "...4", "...7", "...8") # Define um vetor com os novos nomes das colunas
  
  for (j in seq_along(col_indices)) { # Inicia um loop para renomear as colunas
    col_index <- col_indices[j] # Obtém o índice da coluna atual
    new_name <- new_names[j] # Obtém o novo nome da coluna atual
    colnames(planilha)[col_index] <- new_name # Atribui o novo nome à coluna atual
  }
  
  planilha[, 2:4] <- lapply(planilha[, 2:4], gsub, pattern = ",", replacement = ".") # Substitui as vírgulas por pontos nas colunas 2 a 4
  planilha$...4 <- as.numeric(planilha$...4) # Converte a coluna ...4 em numérica
  planilha$...8 <- as.numeric(planilha$...8) # Converte a coluna ...8 em numérica
  planilha$ind <- planilha$...4 - planilha$...8 # Cria uma nova coluna com a diferença entre as colunas ...4 e ...8
  
  # Extrai as SNP names da coluna 1 da planilha
  snps <- unique(grep("^rs(?!.*_(CONTROLE|CASO)$).*", planilha$...1, value = TRUE, perl = TRUE)) # Usa uma expressão regular para encontrar os nomes dos SNPs na coluna ...1 e remove os duplicados
  
  # Data frame só com caso/controle
  rstotais <- planilha[grepl("^rs.*_(CASO|CONTROLE)$", planilha[[1]]), ] # Filtra a planilha para manter apenas as linhas que contêm CASO ou CONTROLE na coluna ...1
  rstotais$...7 <- as.numeric(rstotais$...7) # Converte a coluna ...7 em numérica
  
  
  # Cria uma lista para armazenar os resultados
  resultados <- list() # Cria uma lista vazia para armazenar os resultados dos testes qui-quadrado
  
  # Corre um loop para realizar o teste qui-quadrado para cada SNP
  for (snp in snps) { # Inicia um loop para cada SNP
    filtered_df <- rstotais %>%
      filter(grepl(paste0("^", snp, "_(CASO|CONTROLE)$"), ...1)) # Filtra o data frame rstotais para manter apenas as linhas que contêm o SNP atual e CASO ou CONTROLE na coluna ...1
    
    if (nrow(filtered_df) > 0) { # Verifica se o data frame filtrado tem alguma linha
      # tryCatch prevenir erro = ou < que 0
      tryCatch(
        {
          resultados[[snp]] <- chisq.test(filtered_df[, c('ind', '...8')])  # Realiza o teste qui-quadrado entre as colunas ind e ...8 do data frame filtrado e armazena o resultado na lista resultados
        },
        error = function(e) {
          cat("Erro ao realizar o teste qui-quadrado para o SNP", snp, "\n") # Em caso de erro, imprime uma mensagem com o nome do SNP
        }
      )
      
      # Cria um data frame vazio
      df <- data.frame(SNP = character(), Freq_Caso = numeric(), N_Caso = numeric(),
                       Freq_Controle = numeric(), N_Controle = numeric(), p_valor = numeric(), stringsAsFactors = FALSE) # Cria um data frame vazio com as colunas SNP, Freq_Caso, N_Caso, Freq_Controle, N_Controle e p_valor
      
      
      # Preenche o data frame com os valores do resultado do teste qui-quadrado e os valores de N Caso e N Controle
      for (snp in snps) { # Inicia um loop para cada SNP
        if (snp %in% names(resultados)) { # Verifica se o SNP atual está na lista de resultados
          filtro_caso <- grepl(paste0("^", snp, "_CASO$"), rstotais$...1) # Cria um vetor lógico para filtrar as linhas que contêm o SNP atual e CASO na coluna ...1
          filtro_controle <- grepl(paste0("^", snp, "_CONTROLE$"), rstotais$...1) # Cria um vetor lógico para filtrar as linhas que contêm o SNP atual e CONTROLE na coluna ...1
          Freq_Caso <- rstotais[filtro_caso, "...7"] # Extrai o valor da frequência do caso para o SNP atual
          N_Caso <- rstotais[filtro_caso, "...4"] # Extrai o valor de N Caso para o SNP atual
          Freq_Controle<- rstotais[filtro_controle, "...7"] # Extrai o valor da frequência do controle para o SNP atual
          N_Controle <- rstotais[filtro_controle, "...4"] # Extrai o valor de N Controle para o SNP atual
          p_valor <- resultados[[snp]]$p.value # Extrai o valor de p do resultado do teste qui-quadrado para o SNP atual
          df <- rbind(df, data.frame(SNP = snp, Freq_Caso = Freq_Caso,
                                     N_Caso = N_Caso, Freq_Controle = Freq_Controle, N_Controle = N_Controle, p_valor = p_valor)) # Adiciona uma nova linha ao data frame com os valores extraídos
        } else {
          cat("Não foram encontrados resultados para o SNP", snp, "\n") # Em caso de ausência de resultados, imprime uma mensagem com o nome do SNP
        }
      }
      
      colnames(df)[2:5] <- c("Freq_Caso", "N_Caso", "Freq_Controle", "N_Controle") # Renomeia as colunas 2 a 5 do data frame
      
      
      if (!dir.exists(directory)) { # Verifica se o diretório de saída existe
        dir.create(directory, recursive = TRUE) # Se não existir, cria o diretório de forma recursiva
      }
      
      file_path <- file.path(directory, paste0(nome_gene, "_", sheet_name, "_pvalor", ".csv")) # Define o caminho completo do arquivo CSV de saída
      
      write.csv2(df, file = file_path, row.names = FALSE) # Salva o data frame em um arquivo CSV com ponto e vírgula como separador e sem os nomes das linhas
    } else {
      cat("Não foram encontrados resultados para o SNP", snp, "\n") # Em caso de ausência de resultados, imprime uma mensagem com o nome do SNP
    }
  }
  
  return(resultados) # Retorna a lista de resultados dos testes qui-quadrado
}

# Usa a função lapply para ler todas as planilhas do arquivo Excel de entrada e armazená-las em uma lista
planilhas <- lapply(excel_sheets(InputGene), read_excel, path = InputGene, col_names = TRUE)

# Itera em loop cada planilha para fazer o teste qui-quadrado
for (i in seq_along(planilhas)) {
  planilha <- planilhas[[i]]  # Obtém o data frame da planilha atual
  sheet_name <- excel_sheets(InputGene)[i]  # Obtém o nome da planilha atual
  resultado <- fun_qui_quadrado(planilha, sheet_name) # Chama a função fun_qui_quadrado com a planilha e o nome da planilha como argumentos e armazena o resultado em uma variável
  # Process or save the results as needed
}