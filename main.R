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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. SELECCAO DE DADOS                                                                ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Excluir a partida codigos que fazem sentido considerar
codExclusao  <-  c(87,88,666)
listaTodasParcelas <- listaTodasParcelasIniciais[!(listaTodasParcelasIniciais$cultura %in% codExclusao),]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Determinar as culturas que ocupam maior parte da ?rea (90% ou 95%)

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
dadosANOVA <- constroiDadosANOVA(ano=2005,data=5,banda=1,dimAmostra=15)
table(dadosANOVA$cultura)
plot(dadosANOVA)

#Modelo da ANOVA
dados.aov <- aov(reflectancias ~ cultura/parcela, data=dadosANOVA)
summary(dados.aov)

coef(dados.aov)
model.tables(dados.aov, type="means")

#fazer analise de residuos(normalidade etc.)
plot(dados.aov,which=1)

#Teste para avergiuar se tem distrubuicao normal
shapiro.test(dadosANOVA$reflectancias)
#CONCLUSAO: e tudo menos normal...

plot(dadosANOVA$reflectancias)


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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3.2 Assinaturas espectrais POR CULTURA ----

#Assinaturas todas
assinaturas <- matrix(nrow=0,ncol=length(BANDAS.LANDSAT))
for(i in 1:nrow(listaTodasParcelas))
{
  assinatura <- c()
  for(j in seq(0,6*7,8)+7)
    assinatura <- c(assinatura,listaTodasParcelas[i,j])
  
  assinaturas <- rbind(assinaturas,assinatura)
}

plot(BANDAS.LANDSAT,assinaturas[1,],type='l',ylim=c(0,0.35))
for(i in 1:length(assinaturas))
  lines(BANDAS.LANDSAT,assinaturas[i,],col=i)

#Assinaturas medias por ocupacao de solo
dados <- data.table(listaTodasParcelas)
assinaturasPorCulturaD3 <- dados[,list(area=sum(area),
                                     B1_d3=mean(B1_d3),
                                     B2_d3=mean(B2_d3),
                                     B3_d3=mean(B3_d3),
                                     B4_d3=mean(B4_d3),
                                     B5_d3=mean(B5_d3),
                                     B6_d3=mean(B6_d3)),by=cultura]
assinaturasPorCulturaD3 <- as.data.frame(assinaturasPorCulturaD3)
assinaturasPorCulturaD3 <- assinaturasPorCulturaD3[order(assinaturasPorCulturaD3$cultura, decreasing = F),]


plot(BANDAS.LANDSAT,assinaturasPorCulturaD3[1,3:8],type='l',ylim=c(min(assinaturasPorCulturaD3$B6_d3),max(assinaturasPorCulturaD3$B4_d3)))
text(5,assinaturasPorCulturaD3[1,7],assinaturasPorCulturaD3[1,1])
for(i in 1:nrow(assinaturasPorCulturaD3))
{
  lines(BANDAS.LANDSAT,assinaturasPorCulturaD3[i,3:8],col=i)
  text(5,assinaturasPorCulturaD3[i,7],assinaturasPorCulturaD3[i,1],col=i)
}

#Destacar pastagem permanente
lines(BANDAS.LANDSAT,assinaturasPorCulturaD3[1,3:8],col='red',lwd=5)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3.3 NDVI's por cultura ----
NDVIPorCulturaD4 <- dados[,list(area=sum(area),
                              NDVI_d1=mean(NDVI_d1),
                              NDVI_d2=mean(NDVI_d2),
                              NDVI_d3=mean(NDVI_d3),
                              NDVI_d4=mean(NDVI_d4),
                              NDVI_d5=mean(NDVI_d5),
                              NDVI_d6=mean(NDVI_d6),
                              NDVI_d7=mean(NDVI_d7),
                              NDVI_d8=mean(NDVI_d8)),by=cultura]
NDVIPorCulturaD4 <- as.data.frame(NDVIPorCulturaD4)

datasGrafico <- datas <- names(listaDados[[1]][[1]][[1]]$DR)
substr(datasGrafico,0,1) <- "0"
datasGrafico <- as.numeric(datasGrafico)

plot(datasGrafico,NDVIPorCulturaD4[1,3:10],type='l',ylim=c(0,1))
text(5,NDVIPorCulturaD4[1,7],NDVIPorCulturaD4[1,1])
for(i in 1:nrow(NDVIPorCulturaD4))
{
  lines(datasGrafico,NDVIPorCulturaD4[i,3:10],col=i)
  text(5,NDVIPorCulturaD4[i,7],NDVIPorCulturaD4[i,1],col=i)
}

lines(datasGrafico,NDVIPorCulturaD4[1,3:10],col='red',lwd=5)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3.4 Historgramas/densidades ----

#Histogramas
h <- hist(listaTodasParcelas$B4_d1,breaks=50,prob=T)
lines(density(listaTodasParcelas$B4_d1),col='red')

#Historgramas/densidades ao longo do TEMPO com BANDA fixa
plot(density(listaTodasParcelas$B4_d1),col=1,xlim=c(0,0.5))
lines(density(listaTodasParcelas$B4_d2),col=2)
lines(density(listaTodasParcelas$B4_d3),col=3)
lines(density(listaTodasParcelas$B4_d4),col=4)
lines(density(listaTodasParcelas$B4_d5),col=5)
lines(density(listaTodasParcelas$B4_d6),col=6)

#Historgramas/densidades ao longo da BANDA com TEMPO fixo
plot(density(listaTodasParcelas$B1_d4),col=1,xlim=c(0,0.5))
lines(density(listaTodasParcelas$B2_d4),col=2)
lines(density(listaTodasParcelas$B3_d4),col=3)
lines(density(listaTodasParcelas$B4_d4),col=4)
lines(density(listaTodasParcelas$B5_d4),col=5)
lines(density(listaTodasParcelas$B6_d4),col=6)

#Normalidade das reflectancias - Rejeita se a hipotese de virem duma populacao com distribuicao normal
hist(listaTodasParcelas$B4_d1,breaks=20)
qqnorm(listaTodasParcelas$B4_d1)
shapiro.test(listaTodasParcelas$B4_d1[1:5000])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3.5 Vizualizacao 3D dos dados nas 3 variaveis mais importantes ----

#Seleccao das 3 var's mais importantes
trim.matrix(cor(dadosTreino[,-1]),0.5)

#Graficos 3D
library(rgl)
cultGrafico <- c(1,3,9,13)
plot3d(x=dadosTreino$B3_d3[!dadosTreino$cultura %in% cultGrafico],
       y=dadosTreino$B4_d3[!dadosTreino$cultura %in% cultGrafico],
       z=dadosTreino$B4_d6[!dadosTreino$cultura %in% cultGrafico],
       xlab="B3_d3",
       ylab="B4_d3",
       zlab="B4_d6",
       size=3)
points3d(x=dadosTreino$B3_d3[dadosTreino$cultura == 1],
         y=dadosTreino$B4_d3[dadosTreino$cultura == 1],
         z=dadosTreino$B4_d6[dadosTreino$cultura == 1],
         col="red") #Pastagem
points3d(x=dadosTreino$B3_d3[dadosTreino$cultura == 3],
         y=dadosTreino$B4_d3[dadosTreino$cultura == 3],
         z=dadosTreino$B4_d6[dadosTreino$cultura == 3],
         col="blue")  #Arroz
points3d(x=dadosTreino$B3_d3[dadosTreino$cultura == 9],
         y=dadosTreino$B4_d3[dadosTreino$cultura == 9],
         z=dadosTreino$B4_d6[dadosTreino$cultura == 9],
         col="green")   #Aveia
points3d(x=dadosTreino$B3_d3[dadosTreino$cultura == 13],
         y=dadosTreino$B4_d3[dadosTreino$cultura == 13],
         z=dadosTreino$B4_d6[dadosTreino$cultura == 13],
         col="orange")  #Beterraba



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. SELECCAO DE VARIAVEIS                                                            ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4.1 Conjunto completo ----

#Preparacao dos dados para os classificadores - Apenas assinaturas espectrais das primeiras 6 datas
#dadosClassificadores <- listaTodasParcelas[,c(5,7:12,15:20,23:28,31:36,39:44,47:52)]
dadosClassificadores <- listaTodasParcelas[,c(5,7:42)]
dadosClassificadores$cultura <- as.factor(dadosClassificadores$cultura)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4.2 Sub-conjunto do conjunto completo ----

#Usando a funcao trim.matrix do package subselect
dadosTrim <- dadosClassificadores
criterioExclusaoVars <- 0.02
classTrim <- trim.matrix(cor(dadosTrim[,-1]), criterioExclusaoVars);classTrim

classTrim$names.discarded == classTrimOLD$names.discarded

varsRetirar <- classTrim$numbers.discarded+1  #+1 porque na data.frame original, a coluna 1 é a cultura
varsSobram <- 1:length(dadosTrim)
varsSobram <- varsSobram[! varsSobram %in% varsRetirar]

#Preparacao dos NOVOS DADOS APOS SELECCAO DE VARIAVEIS
dadosClassificadoresSub <- dadosTrim[,varsSobram]
dadosClassificadoresSub$cultura <- as.factor(dadosClassificadoresSub$cultura)

#Seleccao dos dados
dadosTreinoSub <- dadosClassificadoresSub[amostra,]
dadosValidacaoSub <- dadosClassificadoresSub[-amostra,]

#So usar SD's da resultados muito maus, e mesmo misturando SD's com medias da resultados maus
#NDVI's sao melhores mas mesmo assim piores que reflectancias
#dadosValidacaoSubAreas <- listaTodasParcelas$area[-amostra]

x <- cor(dadosTrim[,-1])
anneal(cor(dadosTrim[,-1]), kmin=6, kmax=12, crit="rm")
rm.coef(cor(dadosTrim[,-1]), ind=varsSobram)

pca <- prcomp(dadosTrim[,-1], cor = T)
summary(pca)
head(dadosTrim[,-1])

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
#REsultado: usar k=7, porque é o mais frequente e tem o erro mais baixo


#Validacao cruzada
KNN.comp.cruz <- validacaoCruzada(n = 10, tipo = 'KNN', dados = dadosClassificadores, lambda = 0.8, k = 7)




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
#REsultado: usar k=7, porque é o mais frequente e tem o erro mais baixo

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

SVM.sub.cruz$result$correcTot

#Old tuning version
# SVM.tune.sub <- tune.svm(dadosClassificadoresSub[,-1], dadosClassificadoresSub[,1], gamma=c(0.5,0.6), cost=c(3.2, 3.5))
# gamma.sub <- SVM.tune.sub[1][[1]][1,1];gamma.sub
# cost.sub <- SVM.tune.sub[1][[1]][1,2];cost.sub



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Fazer a validacao ao contrario para ter a certeza que esta tudo bem

#resultCV e o resultado da validacao cruzada efectuada
obtemPercentagemAceites <- function(resultCV, lambda)
{
  #Juncao das classificacoes em cada fold
  class <- NA
  for(i in 1:10)
    class <- rbind(class, resultCV$resultTodos[[i]]$result)
  class <- class[-1,]
  
  #Seleccao apenas das classificacoes aceites
  qi <- resultCV$result$qi
  classAceites <- class[class$prob1 >= qi[class$class1],]
  
  #Matriz de erro e accuracy ststistics
  matErro <- as.data.frame.matrix(table(class = factor(classAceites$class1, levels=levels(classAceites$verdade)), verdade = classAceites$verdade))
  correcTot <- cc(matErro);correcTot
  userA <- diag(as.matrix(matErro))/rowSums(matErro);userA
  prodA <- diag(as.matrix(matErro))/colSums(matErro);prodA
  
  #Percentagem de parcelas com classificacao aceites
  percClass <- nrow(classAceites)/nrow(dadosClassificadores);percClass
  
  #QUESTAO: como medir a percentagem de parcelas com classificacao aceite?
  #Isto e a versao com a verdade, que faz mais sentido
  percClassCult <- table(factor(classAceites$class1, levels=levels(classAceites$verdade)))/table(dadosClassificadores$cultura);percClassCult
  
  #Isto e a versao com a classificacao
  percClassCult <- table(factor(classAceites$class1, levels=levels(classAceites$verdade)))/table(class$class1);percClassCult
}

#TODO:
#Fazer funcao com isto
#Resultado da validacao cruzada tem que ser os qi para todos os lambdas e nao a percentagem de parcelas aceites, isso faz a nova funçao.





#Explicacao para resultados acima nao fazerem sentido nenhum -> o numero de parcelas que podem ser classificadas devia ser calculada DEPOIS da validaçao cruzada!
SVM.sub.cruz$resultTodos$valid3$result[SVM.sub.cruz$resultTodos$valid3$result$class1 == 7,]
SVM.sub.cruz$resultTodos$valid3$qi
SVM.sub.cruz$resultTodos$valid3$percsDecs[7+1,]














class <- NA
for(i in 1:10)
  class <- rbind(class, SVM.sub.cruz$resultTodos[[i]]$result)
class <- class[-1,]
qi <- SVM.sub.cruz$result$qi
classAceites <- class[class$prob1 >= qi[class$class1],]


matErro <- as.data.frame.matrix(table(class = factor(classAceites$class1, levels=levels(classAceites$verdade)), verdade = classAceites$verdade))
correcTot <- cc(matErro);correcTot
userA <- diag(as.matrix(matErro))/rowSums(matErro);userA
prodA <- diag(as.matrix(matErro))/colSums(matErro);prodA


percClass <- nrow(classAceites)/nrow(dadosClassificadores);percClass

#QUESTAO: como medir a percentagem de parcelas com classificacao aceite?
#Isto e a versao com a verdade, que faz mais sentido
percClassCult <- table(factor(classAceites$class1, levels=levels(classAceites$verdade)))/table(dadosClassificadores$cultura);percClassCult

#Isto e a versao com a classificacao
percClassCult <- table(factor(classAceites$class1, levels=levels(classAceites$verdade)))/table(class$class1);percClassCult

