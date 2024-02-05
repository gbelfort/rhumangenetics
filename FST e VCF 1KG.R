#======================================================
# Calculo do Fst
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
# Contributor(s):  Caique Manochio
# Updated by (and date): Guilherme Belfort Almeida 05/02/2024
#
# Dependencies: R
#
# Command line:

# RStudio -> Clicar na lupa acima e substituir GENEDADO pelo nome do gene a ser usado

library('vcfR') 
library('adegenet')
library('tibble')
library('dplyr')

# importar vcf pro R

vcfGENEDADO <- read.vcfR(file.choose())

# meta data, fixed data e genotype data

# montar genotype data em dataframe
gt <- extract.gt(vcfGENEDADO, element = "GT")


# transpor para depois adicionar populações na vertical
transposto <- as.data.frame(t(gt))
transposto <- tibble::rownames_to_column(transposto,"ID")


# importar correspondência indivíduo-população direto da fase correspondente do 1KG
dados_pop_1kgp <- read.csv(file.choose(), sep = "\t")
# renomear para igualar ao nome do gt
dados_pop_1kgp <- dados_pop_1kgp %>%
  rename(ID = Sample.name)
# organizar por ID
Pops1KG <- dados_pop_1kgp[order(dados_pop_1kgp$ID),]

# left_join atribuindo indivíduo filtrado (2501) a população, tirando indivíduos não estudados (uns 600)
merge <- left_join(transposto, Pops1KG, by= "ID")

# Ao juntar, sobram algumas colunas desnecessárias como Sex, Biosample.ID, Superpopulation.code, Superpopulation.name, Population.elastic.ID e Data.collections
# Elas ficam no final do data frame. Tirei manualmente, prestando atenção para não tirar o Population.code e o Population.name
# Sobram ID, Population.code, Population.name e os 'rs'.
View(merge)
# Encontrar o número das colunas a retirar, criar novo dataframe sem elas e colocar a ordem (ID, code, name e o resto everything)
# Exemplo abaixo

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
GENEDADO[GENEDADO == '1|0'] <- "1 0"
GENEDADO[GENEDADO == '1|1'] <- "1 1"
GENEDADO[GENEDADO == '0|1'] <- "0 1"
GENEDADO[GENEDADO == '0|0'] <- "0 0"
GENEDADO[GENEDADO == '2|0'] <- "2 0"
GENEDADO[GENEDADO == '0|2'] <- "0 2"
GENEDADO[GENEDADO == '1|2'] <- "1 2"
GENEDADO[GENEDADO == '2|1'] <- "2 1"
GENEDADO[GENEDADO == '2|2'] <- "2 2"
#ad infinitum "X|X" <- "X X"

# Fst:

library(poppr)
library(hierfstat)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(scales)

# Criar genind a partir de data frame. Separador padrão '|', coluna indivíduos GENEDADO$ID, populações GENEDADO$popCode
# Ploidia 2|2

GENEDADO_gen <- df2genind(
  GENEDADOpolimorf,
  sep = '|',
  ind.names = GENEDADO$ID,
  pop = GENEDADO$popCode,
  ploidy = 2|2,
  type = c("codom", "PA"),
  check.ploidy = getOption("adegenet.check.ploidy"))

# Opcional: Selecionar apenas alguns "rs" dentro do VCF, por exemplo para usar apenas os polimorfismos. 
# Manual, fazer uma lista com todos os polimorfismos antes e colar abaixo após 'popName', seguindo o formato
GENEDADOpolimorf <- select(GENEDADO, 'ID', 'popCode','popName', "rs74837985", "rs140584594",
                           "rs143315534","rs147668562","rs449856")
)
# Esse GENEDADOpolimorf substitui o banco GENEDADO no comando df2genind() acima
                           



# Criar subconjunto dos dados para diminuir tempo de computação
GENEDADO_gen_sub = popsub(GENEDADO_gen, sublist = c("ACB","ASW","ESN","GWD","LWK","MSL","YRI","CLM","MXL","PEL","PUR","CDX","CHB","CHB","CHS","JPT","KHV","CEU","FIN","GBR","IBS","TSI","BEB","GIH","ITU","PJL","STU"))

# FST teste de todos os pares. Método WC84 (estimativa FSTs seguindo Nei (1987) e Weir & Cockerham (1984))
# 3 casas decimais
GENEDADO_fst = genet.dist(GENEDADO_gen_sub, method = "WC84") %>% round(digits = 3)




# Escolher ordem das populações
lab_order = c("ACB","ASW","ESN","GWD","LWK","MSL","YRI","CLM","MXL","PEL","PUR","CDX","CHB","CHB","CHS","JPT","KHV","CEU","FIN","GBR","IBS","TSI","BEB","GIH","ITU","PJL","STU")

# Arrumar ordem de colunas e linhas para a acima
fst.mat = as.matrix(GENEDADO_fst)
fst.mat1 = fst.mat[lab_order, ]
fst.mat2 = fst.mat1[, lab_order]

# Criar data frame
ind = which(upper.tri(fst.mat2), arr.ind = TRUE)
fst.df = data.frame(Site1 = dimnames(fst.mat2)[[2]][ind[,2]],
                    Site2 = dimnames(fst.mat2)[[1]][ind[,1]],
                    Fst = fst.mat2[ ind ])
# No fim, isso faz uma tabela em que a coluna 91 compara CHB com CHB, achando Fst 0. Então retirar a coluna 91 para não incluir isso
fst.df <- fst.df[-c(91),]

# Manter a ordem para o gráfico
fst.df$Site1 = factor(fst.df$Site1, levels = unique(fst.df$Site1))
fst.df$Site2 = factor(fst.df$Site2, levels = unique(fst.df$Site2))

# ?? Convert minus values to zero
fst.df$Fst[fst.df$Fst < 0] = 0

# Sumário do data frame
fst.df %>% str

# Legenda com F em itálico e ST diminuído
fst.label = expression(italic("F")[ST])

# Extract middle Fst value for gradient argument
mid = max(fst.df$Fst) / 2

# Gráfico heatmap
#3 cores
ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 3)+
  scale_fill_gradient2(low = "blue", mid = "pink", high = "red", midpoint = mid, name = fst.label, limits = c(0, max(fst.df$Fst)), breaks = c(0.05, 0.15))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10)
  )


#2 cores
ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 3)+
  scale_fill_gradient(low = "white", high = "red", name = fst.label, limits = c(min(fst.df$Fst), max(fst.df$Fst)))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=10)
  )



#alternativa

library(ggplot2)
library(reshape2)

# Transformar dados no formato longo
dados_long <- melt(fst.mat2)

# Criar gráfico heatmap com geom_tile
ggplot(data = dados_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(colour = "black") +
  geom_text(aes(label = value), color="black", size = 3) +
  scale_fill_gradient(low = "white", high = "red", name = fst.label, limits = c(min(fst.df$Fst), max(fst.df$Fst)))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=10)
  )
