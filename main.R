#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Master's Thesis - Remote Sensing                      %
# Environmental Engineering - ISA/UL - Lisbon, Portugal %
# (c) 2014 by Jonas Schmedtmann                         %
#                                                       %
# MAIN SCRIPT                                           %
#                                                       %
# Implements the program logic using functions from     %
# functions.R. All steps typically found in remote      %
# sensing can be found here. Documented in portuguese.  %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Indice de passos seguidos:
#   01. AQUISICAO DE DADOS E PROCESSAMENTO
#   02. SELECCAO DE DADOS
#   03. ANALISE EXPLORATORIA DOS DADOS
#   04. SELECCAO DE VARIAVEIS
#   05. 06. 07. ESCOLHA/TREINO/VALIDACAO DOS CLASSIFICADORES


#CMD+CTRL+S: save all
#CTRL+ALT+->: go to other file
#CMD+SHIFT+P: run previous selection


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. AQUISICAO DE DADOS E PROCESSAMENTO                                               ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Iniciar
source('functions.R')
init()

#Definir area de estudo, carregar dados, recortar e projectar
coordsArea <- cbind(AREA.X,AREA.Y)
cdg <- carregaRecortaDadosEspaciais(coordsArea, PROJ4.UTM)
plot(cdg$area);plot(cdg$parc2005, add=T)

#Obter e corrigir todas as imagens
rm(todasImagens)
todasImagens <- constroiListaImagensLandsat(landsatPath=CAMINHO.LANDSAT,
                                            areaEstudo=cdg$area,
                                            prefixo="CORR_14.08",
                                            ano=2005,
                                            corrige=TRUE)

#Obter lista de todos os dados
rm(listaDados)
listaDados <- constroiListaDados(ano=2005)

#Obter conteudo relevante para TODAS as parcelas, em forma de data.frame
listaTodasParcelasIniciais <- constroiTodasParcelas()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. SELECCAO DE DADOS                                                                ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Excluir a partida codigos que fazem sentido considerar
codExclusao  <-  c(87,88,666)
listaTodasParcelas <- listaTodasParcelasIniciais[!(listaTodasParcelasIniciais$cultura %in% codExclusao),]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Determinar as culturas que ocupam maior parte da area (90% ou 95%)

#Todas as culturas com respectiva area total
dados <- data.table(listaTodasParcelas)
areasCulturas <- dados[,list(area=sum(area),numParc=length(area)),by=cultura]
areasCulturas <- areasCulturas[order(areasCulturas$area, decreasing = TRUE),]
areasCulturas <- cbind(areasCulturas,cumsum(areasCulturas$area))
areaTotal <- sum(areasCulturas$area)

#Culturas acima e abaixo dos limiares
plot(1:nrow(areasCulturas),areasCulturas$V2)
abline(h=areaTotal*0.9,col='blue')
abline(h=areaTotal*0.95,col='orange')
abline(h=areaTotal*0.98,col='red')

limite<-0.95
cultInfl <- areasCulturas[areasCulturas$V2 < areaTotal*limite,]
cultInfl <- cultInfl[!cultInfl$cultura == 27]
fraccaoInfl <- sum(cultInfl$area)/areaTotal
cultInfluentes <- cultInfl$cultura
length(cultInfluentes)

#Numero de culturas restantes
nrow(areasCulturas)-13

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Seleccionar apenas aquelas que representam 95% da area e reclassificar
listaTodasParcelas <- listaTodasParcelas[listaTodasParcelas$cultura %in% cultInfluentes,]

novasClasses <- as.data.frame(cbind(cultInfluentes, c(1,2,3,4,5,6,7,8,5,9,10,11,12)))
colnames(novasClasses) <- c('cultura','novaClasse')
nClasses <- length(table(novasClasses$novaClasse))

for(i in 1:length(listaTodasParcelas$cultura))
  listaTodasParcelas$cultura[i] <- novasClasses$novaClasse[which(novasClasses$cultura==listaTodasParcelas$cultura[i])]
#Este objecto esta guardado como ficheiro .Robj


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. ANALISE EXPLORATORIA DOS DADOS COMPLETOS                                         ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3.1 ANOVA para averiguar se cultura e parcela afectam reflectancias ----


#%%%%%%%%%%%%%%%%%%%%%%%%%%
#Primeira abordagem
dadosANOVA <- constroiDadosANOVA(data=5, banda=4, dimAmostra=200)
table(dadosANOVA$cultura)
plot(dadosANOVA$reflectancias, dadosANOVA$cultura)

#Modelo da ANOVA
dados.aov <- aov(reflectancias ~ cultura/parcela, data=dadosANOVA)
summary(dados.aov)



#fazer analise de residuos(normalidade etc.)
par(mfrow=c(2,2))
plot(dados.aov,which=c(1,2,4,5))



#%%%%%%%%%%%%%%%%%%%%%%%%%%
#Nova abordagem - averigua so se a cultura tem efeitos sobre cada uma das combinacoes banda/data e NDVI's
FValues <- c()
pValues <- c()
efeito <- c()
for(i in 7:length(listaTodasParcelas))
{
  parcelas.aov <- aov(listaTodasParcelas[,i] ~ as.factor(listaTodasParcelas$cultura))
  FV <- summary(parcelas.aov)[[1]][["F value"]][[1]]
  pV <- summary(parcelas.aov)[[1]][["Pr(>F)"]][[1]]
  
  FValues <- c(FValues,FV)
  pValues <- c(pValues,pV)
  if(pV <= 0.05) ef <- 1 else ef <- 0
  efeito <- c(efeito,ef)
}
nn <- colnames(listaTodasParcelas[,7:length(listaTodasParcelas)])
resultadoANOVA <- data.frame(cbind(nn,FValues,pValues,efeito))
resultadoANOVA$FValues <- as.numeric(as.character(resultadoANOVA$FValues))
resultadoANOVA <- resultadoANOVA[order(-resultadoANOVA$FValues),]; resultadoANOVA




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. SELECCAO DE VARIAVEIS                                                            ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4.1 Conjunto completo ----

#Preparacao dos dados para os classificadores - Apenas assinaturas espectrais das primeiras 6 datas
dadosClassificadores <- listaTodasParcelas[,c(5,7:42)]
dadosClassificadores$cultura <- as.factor(dadosClassificadores$cultura)
dadosClassificadores[,-1][dadosClassificadores[,-1] > 1] <- 1


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4.2 Sub-conjunto do conjunto completo ----

#Usando a funcao trim.matrix do package subselect
dadosTrim <- dadosClassificadores
criterioExclusaoVars <- 0.02
classTrim <- trim.matrix(cor(dadosTrim[,-1]), criterioExclusaoVars);classTrim

classTrim$names.discarded == classTrimOLD$names.discarded

varsRetirar <- classTrim$numbers.discarded+1  #+1 porque na data.frame original, a coluna 1 ?? a cultura
varsSobram <- 1:length(dadosTrim)
varsSobram <- varsSobram[! varsSobram %in% varsRetirar]

#Preparacao dos NOVOS DADOS APOS SELECCAO DE VARIAVEIS
dadosClassificadoresSub <- dadosTrim[,varsSobram]
dadosClassificadoresSub$cultura <- as.factor(dadosClassificadoresSub$cultura)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. 6. 7. ESCOLHA/TREINO/VALIDACAO DOS CLASSIFICADORES                               ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5.1 KNN - K vizinhos mais proximos ----

#%%%%%%%%%%%%%%%%%%%%%%%%
#Dataset 1: Todos os dados ---- PRONTO!!!

#Tuning
resTune.comp.knn <- matrix(nrow=21, ncol=0)
for(i in 1:10)
{
  print(i)
  KNN.tune.comp <- tune.knn(dadosClassificadores[,-1], dadosClassificadores[,1], k=1:20)
  k <- KNN.tune.comp[1][[1]][1,1]
  erros <- KNN.tune.comp$performances$error
  resTune.comp.knn <- cbind(resTune.comp.knn, c(k, erros))
}
#REsultado: usar k=7, porque ?? o mais frequente e tem o erro mais baixo


#Validacao cruzada
KNN.comp.cruz <- validacaoCruzada(n = 10, tipo = 'KNN', dados = dadosClassificadores, lambda = 0.8, k = 7)
KNN.comp.cruz$result$correcTot



#%%%%%%%%%%%%%%%%%%%%%%%%
#Dataset 2: Apenas sub-conjunto de dados ---- PRONTO!!!

#Tuning
resTune.sub.knn <- matrix(nrow=21, ncol=0)
for(i in 1:10)
{
  print(i)
  KNN.tune.sub <- tune.knn(dadosClassificadoresSub[,-1], dadosClassificadoresSub[,1], k=1:20)
  k <- KNN.tune.sub[1][[1]][1,1]
  erros <- KNN.tune.sub$performances$error
  resTune.sub.knn <- cbind(resTune.sub.knn, c(k, erros))
}
#REsultado: usar k=7, porque ?? o mais frequente e tem o erro mais baixo

#Validacao cruzada
KNN.sub.cruz <- validacaoCruzada(n = 10, tipo = 'KNN', dados = dadosClassificadoresSub, lambda = 0.8, k = 7)

KNN.sub.cruz$result$correcTot

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5.2 SVM - Maquinas de vectores de suporte ----

#%%%%%%%%%%%%%%%%%%%%%%%%
#Dataset 1: Todos os dados ---- PRONTO!!!

#Tuning
rGamma <- c(0.4, 0.6, 0.8, 1, 1.2, 1.4)
rCost <- c(0.5, 1, 1.5, 2, 2.5, 3)
resTune.comp.v1 <- matrix(ncol = length(rCost), nrow = length(rGamma))

for(i in 1:length(rGamma))
  for(j in 1:length(rCost))
  {
    print(paste0("i = ", i, " and j = ", j))
    res <- SVM.tune.comp <- tune.svm(dadosClassificadores[,-1], dadosClassificadores[,1], gamma = rGamma[i], cost = rCost[j])
    
    resTune.comp.v1[i,j] <- res$best.performance
  }
colnames(resTune.comp.v1) <- rCost
rownames(resTune.comp.v1) <- rGamma
min(resTune.comp.v1)


#Validacao cruzada
SVM.comp.cruz <- validacaoCruzada(n = 10, tipo = 'SVM', dados = dadosClassificadores, lambda = 0.8, gamma = 0.4, cost = 1.5)

SVM.comp.cruz$result$correcTot


#%%%%%%%%%%%%%%%%%%%%%%%%
#Dataset 2: Apenas sub-conjunto de dados ---- PRONTO!!!

#Tuning
rGamma <- c(0.4, 0.6, 0.8, 1, 1.2, 1.4)
rCost <- c(0.5, 1, 1.5, 2, 2.5, 3)
resTune.sub.v1 <- matrix(ncol = length(rCost), nrow = length(rGamma))

for(i in 1:length(rGamma))
  for(j in 1:length(rCost))
  {
    print(paste0("i = ", i, " and j = ", j))
    res <- SVM.tune.sub <- tune.svm(dadosClassificadoresSub[,-1], dadosClassificadoresSub[,1], gamma = rGamma[i], cost = rCost[j])

    resTune.sub.v1[i,j] <- res$best.performance
  }
colnames(resTune.sub.v1) <- rCost
rownames(resTune.sub.v1) <- rGamma
min(resTune.sub.v1)


#Validacao cruzada
SVM.sub.cruz <- validacaoCruzada(n = 10, tipo = 'SVM', dados = dadosClassificadoresSub, lambda = 0.8, gamma = 0.4, cost = 2.5)


#Estimacao dos q_j
SVM.sub.qi <- estimaQ(classificador = SVM.sub.cruz$classificador, lambdas = c(0.6, 0.7, 0.8, 0.9, 0.95))
SVM.sub.qi.NOVO <- NOVOestimaQ(classificador = SVM.sub.cruz$classificador, lambdas = c(0.8))



#Assessment do metodo
SVM.sub.assess <- assessmentMetodo(classificador = SVM.sub.cruz$classificador,  qis = SVM.sub.qi)
SVM.sub.assess.NOVO <- assessmentMetodo(classificador = SVM.sub.cruz$classificador,  qis = SVM.sub.qi.NOVO)







#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Geração de estruturas de dados finais ----

#%%%%%%%%%%%%%%%%%%%%%%%%
# Data.frame com todas as parcelas e informacoes provenientes do metodo
comp1 <- SVM.sub.analise$resultClass
comp1$parcela <- as.character(comp1$parcela)
comp1 <- comp1[order(comp1$parcela, decreasing = TRUE),]

comp2 <- listaTodasParcelas
comp2$id <- as.character(comp2$id)
comp2 <- comp2[order(comp2$id, decreasing = TRUE),]

igual <- comp1$verdade == comp1$class1
igual[igual == TRUE]  <- 'Igual'
igual[igual == FALSE]  <- 'Diferente'

qj <- c()
for(i in 1:nrow(comp1)) qj[i] <- SVM.sub.analise$l0.8$qi[comp1$class1[i]]

resultado <- comp1$prob1 >= qj
resultado[resultado == TRUE]  <- 'Aceite'
resultado[resultado == FALSE]  <- 'Rejeitado'

dadosFinais <- cbind(comp1, igual, qj, resultado, area=comp2$area)


#%%%%%%%%%%%%%%%%%%%%%%%%
#Lista de dados como a original, mas agora so com as parcelas que foram efectivamente usadas e com estrutura simplificada

finalListaDados <- list()
for(i in 1:length(listaDados))
  for(j in 1:length(listaDados[[i]][[1]]))
  {
    parc <- names(listaDados[[i]][[1]])[j]
    todasParc <- as.vector(listaTodasParcelas$id)
    if(parc %in% todasParc)
    {
      finalListaDados[[parc]] <- listaDados[[i]][[1]][[j]]
      cultura <- finalListaDados[[parc]][['atributos']][['CULTURA']]
      finalListaDados[[parc]][['atributos']][['CULTURA']] <- novasClasses$novaClasse[which(novasClasses$cultura == cultura)]
    }
  }







