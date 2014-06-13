#########################################################
# Master's Thesis - Remote Sensing                      #
# Environmental Engineering - ISA/UL - Lisbon, Portugal #
# (c) 2014 by Jonas Schmedtmann                         #
#                                                       #
# MAIN SCRIPT                                           #
#                                                       #
# Implements the program logic using functions from     #
# functions.R. All steps typically found in remote      #
# sensing can be found here. Documented in portuguese.  #
#########################################################


#Indice de passos seguidos:
#   01. AQUISICAO DE DADOS E PROCESSAMENTO
#   02. SELECCAO DE DADOS
#   03. ANALISE EXPLORATORIA DOS DADOS
#   04. SELECCAO DE VARIAVEIS
#   05. 06. 07. ESCOLHA/TREINO/VALIDACAO DOS CLASSIFICADORES


##########################################################################################
# 1. AQUISICAO DE DADOS E PROCESSAMENTO                                                  #
##########################################################################################

#Clear console
cat("\014")

#Iniciar
source('functions.R')
init()

#Definir area de estudo, carregar dados, recortar e projectar
areaX <- c(508262,528870,511269,490710,508262)
areaY <- c(4365569,4362462,4288360,4291713,4365569)
coordsArea <- cbind(areaX,areaY)
cdg <- carregaRecortaDadosEspaciais(coordsArea,PROJ4.UTM)

#Obter todas as imagens
rm(todasImagens)
todasImagens<- constroiListaImagensLandsat(CAMINHO.LANDSAT, cdg$area, "CORR_12.05", anos=2005, corrige=TRUE)

#Obter lista de todos os dados
rm(listaDados)
listaDados <- constroiListaDados(anos=2005)

#Obter conteudo relevante para TODAS as parcelas, em forma de data.frame
listaTodasParcelasIniciais <- constroiTodasParcelas()



##########################################################################################
# 2. SELECCAO DE DADOS                                                                   #
##########################################################################################

##########################################################################
#Excluir a partida codigos que fazem sentido considerar
codExclusao  <-  c(87,88,666)
listaTodasParcelas <- listaTodasParcelasIniciais[!(listaTodasParcelasIniciais$cultura %in% codExclusao),]


##########################################################################
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
fraccaoInfl <- sum(cultInfl$area)/areaTotal
cultInfluentes <- cultInfl$cultura


##########################################################################
#Seleccionar apenas aquelas que representam 95% da area e reclassificar
listaTodasParcelas <- listaTodasParcelas[listaTodasParcelas$cultura %in% cultInfluentes,]

novasClasses98 <- as.data.frame(cbind(cultInfluentes,c(1,2,3,4,5,6,7,6,8,9,10,11,12,13,14,5,15,4,16,17,8)))
novasClasses <- as.data.frame(cbind(cultInfluentes,c(1,2,3,4,5,6,7,6,8,9,10,11,12,13,14,5)))
colnames(novasClasses) <- c('cultura','novaClasse')
nClasses <- length(table(novasClasses$novaClasse))

for(i in 1:length(listaTodasParcelas$cultura))
  listaTodasParcelas$cultura[i] <- novasClasses$novaClasse[which(novasClasses$cultura==listaTodasParcelas$cultura[i])]




##########################################################################################
# 3. ANALISE EXPLORATORIA DOS DADOS COMPLETOS                                            #
##########################################################################################

##########################################################################
#ANOVA para averiguar se cultura e parcela afectam reflectancias

###################
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


###############
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
resultadoANOVA <- resultadoANOVA[order(-resultadoANOVA$FValues),]


##########################################################################
#Assinaturas espectrais POR CULTURA

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

plot(BANDAS.LANDSAT,assinaturasPorCulturaD3[1,3:8],type='l',ylim=c(min(assinaturasPorCulturaD3$B6_d3),max(assinaturasPorCulturaD3$B4_d3)))
text(5,assinaturasPorCulturaD3[1,7],assinaturasPorCulturaD3[1,1])
for(i in 1:nrow(assinaturasPorCulturaD3))
{
  lines(BANDAS.LANDSAT,assinaturasPorCulturaD3[i,3:8],col=i)
  text(5,assinaturasPorCulturaD3[i,7],assinaturasPorCulturaD3[i,1],col=i)
}

#Destacar pastagem permanente
lines(BANDAS.LANDSAT,assinaturasPorCulturaD3[1,3:8],col='red',lwd=5)


##########################################################################
#NDVI's por cultura
NDVIPorCulturaD4 <- dados[,list(area=sum(area),
                              NDVI_1=mean(NDVI_1),
                              NDVI_2=mean(NDVI_2),
                              NDVI_3=mean(NDVI_3),
                              NDVI_4=mean(NDVI_4),
                              NDVI_5=mean(NDVI_5),
                              NDVI_6=mean(NDVI_6),
                              NDVI_7=mean(NDVI_7),
                              NDVI_8=mean(NDVI_8)),by=cultura]
NDVIPorCulturaD4 <- as.data.frame(NDVIPorCulturaD4)

datasGrafico <- datas
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



##########################################################################
#Historgramas/densidades

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

##########################################################################
#Vizualizacao 3D dos dados nas 3 variaveis mais importantes

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



##########################################################################################
# 4. SELECCAO DE VARIAVEIS                                                               #
##########################################################################################


##########################################################################
#CONJUNTO COMPLETO

#Preparacao dos dados para os classificadores
dadosClassificadores <- listaTodasParcelas[,c(5,7:12,15:20,23:28,31:36,39:44,47:52)]
dadosClassificadores$cultura <- as.factor(dadosClassificadores$cultura)

#Amostra a usar para o modelo - metade para treino e metade para validacao
set.seed(23)
amostra <- sample(1:nrow(dadosClassificadores), ceiling(nrow(dadosClassificadores)/2))

#Seleccao dos dados
dadosTreino <- dadosClassificadores[amostra,]
dadosValidacao <- dadosClassificadores[-amostra,]


##########################################################################
#SUB-CONJUNTO

#Usando a funcao trim.matrix do package subselect
dadosTrim <- listaTodasParcelas[,c(5,(7:54))]
dadosTrim <- dadosClassificadores
criterioExclusaoVars <- 0.02
classTrim <- trim.matrix(cor(dadosTrim[,-1]),criterioExclusaoVars);classTrim

varsRetirar <- classTrim$numbers.discarded+1  #+1 porque na data.frame original, a coluna 1 é a cultura
varsSobram <- 1:length(dadosTrim)
varsSobram <- varsSobram[! varsSobram %in% varsRetirar]

#Preparacao dos NOVOS DADOS APOS SELECCAO DE VARIAVEIS
dadosClassificadoresSub <- dadosTrim[,varsSobram]
dadosClassificadoresSub$cultura <- as.factor(dadosClassificadoresSub$cultura)

#Amostra a usar para o modelo - metade para treino e metade para validacao
set.seed(23)
amostra <- sample(1:nrow(dadosClassificadoresSub), ceiling(nrow(dadosClassificadoresSub)/2))

#Seleccao dos dados
dadosTreinoSub <- dadosClassificadoresSub[amostra,]
dadosValidacaoSub <- dadosClassificadoresSub[-amostra,]

#So usar SD's da resultados muito maus, e mesmo misturando SD's com medias da resultados maus
#NDVI's sao melhores mas mesmo assim piores que reflectancias


##########################################################################################
# 5. 6. 7. ESCOLHA/TREINO/VALIDACAO DOS CLASSIFICADORES                                  #
##########################################################################################

#QUESTOES:
# - Como vamos poder validar isto sem as declaracoes? 


##########################################################################
#1. KNN - K VIZINHOS MAIS PROXIMOS

###################
#TODOS OS DADOS
classCult.knn <- knn(dadosTreino[,-1], dadosValidacao[,-1], dadosTreino[,1], prob=TRUE, k=10)
classCult.knn.tune <- tune.knn(dadosTreino[,-1], dadosTreino[,1], k=1:20);classCult.knn.tune
melhorK <- classCult.knn.tune[1][[1]][1,1]

#Validacao
resultadoCult.knn <- table(pred=classCult.knn, true=dadosValidacao[,1]);resultadoCult.knn
resultadoCult.knn <- as.data.frame.matrix(resultadoCult.knn)
classDiag.knn <- diag(as.matrix(resultadoCult.knn));classDiag.knn
certos.knn <- sum(classDiag.knn)/sum(resultadoCult.knn);certos.knn

#Probabilidades
probs.knn <- unlist(attributes(classCult.knn),use.names=FALSE)
probs.knn <- as.numeric(probs.knn[-(1:(nClasses+1))])
probs.round.knn <- round(probs.knn,1)
hist(probs.knn,col='blue')
table(pred=classCult.knn, true=dadosValidacao[,1], probs=probs.round.knn)


###################
#APENAS SUB-CONJUNTO DE DADOS
classCultSub.knn <- knn(dadosTreinoSub[,-1], dadosValidacaoSub[,-1], dadosTreinoSub[,1], prob=TRUE, k=10)
classCultSub.tune.knn <- tune.knn(dadosTreinoSub[,-1], dadosTreinoSub[,1], k=1:20);classCult.knn.sub.tune

#Validacao
resultadoCultSub.knn <- table(pred=classCultSub.knn, true=dadosValidacaoSub[,1]);resultadoCultSub.knn
resultadoCultSub.knn <- as.data.frame.matrix(resultadoCultSub.knn)
classDiagSub.knn <- diag(as.matrix(resultadoCultSub.knn));classDiagSub.knn
certosSub.knn <- sum(classDiagSub.knn)/sum(resultadoCultSub.knn);certosSub.knn

#Probabilidades
probsSub.knn <- unlist(attributes(classCultSub.knn),use.names=FALSE)
probsSub.knn <- as.numeric(probsSub.knn[-(1:(nClasses+1))])
probsSub.round.knn <- round(probsSub.knn,1)
hist(probsSub.knn,col='green')
table(pred=classCultSub.knn, true=dadosValidacaoSub[,1], probs=probsSub.round.knn)

#So com P >= 60%
jjj <- table(pred=classCultSub.knn, true=dadosValidacaoSub[,1], probs=probsSub.round.knn)
jjj <- as.matrix(jjj[,,5]+jjj[,,6]+jjj[,,7]+jjj[,,8]+jjj[,,9])
sum(diag(jjj))/sum(jjj)


###################
#ANALISE GRAFICA
#Percentagem de classificacoes correctas em cada cultura
percCertos <- c()
for(i in 1:ncol(resultadoCult.knn))
  percCertos[i] <- resultadoCult.knn[i,i]/sum(resultadoCult.knn[,i])
percCertos <- 100*percCertos

percCertosSub <- c()
for(i in 1:ncol(resultadoCult.knn))
  percCertosSub[i] <- resultadoCultSub.knn[i,i]/sum(resultadoCultSub.knn[,i])
percCertosSub <- 100*percCertosSub

par(las=2, mar=c(8,4.5,1,1))
barplot(rbind(percCertos,percCertosSub), beside=T, space=c(0,2.5), border=NA, names.arg=codigosNovos2005$codNomesAbrev[1:14], ylab='Percentagem de classificacoes correctas (%)', ylim=c(0,100), col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()
legend(x=30, y=98, legend=c("Completo","Sub-conjunto"), fill=c(rgb(0,0,0),rgb(0.7,0.7,0.7)))


#Percentagem de classificacoes correctas em cada valor de probabilidade
classProb <- table(pred=classCult.knn, true=dadosValidacao[,1], probs=probs.round.knn)
percCertos <- c()
for (i in 1:length(classProb[1,1,]))
{
  diag <- diag(as.matrix(classProb[,,i]))
  percCertos[i] <- sum(diag)/sum(classProb[,,i])
}
percCertos <- 100*percCertos

classProbSub <- table(pred=classCultSub.knn, true=dadosValidacaoSub[,1], probs=probsSub.round.knn)
percCertoSubs <- c()
for (i in 1:length(classProbSub[1,1,]))
{
  diag <- diag(as.matrix(classProbSub[,,i]))
  percCertosSub[i] <- sum(diag)/sum(classProbSub[,,i])
}
percCertosSub <- 100*percCertosSub

par(las=1, mar=c(4.5,4.5,1,1))
barplot(rbind(percCertos,percCertosSub), beside=T, space=c(0,2.5), border=NA, ylab='Percentagem de classificacoes correctas (%)', xlab='Probabilidade associada a classificao (%)', ylim=c(0,100), names.arg=seq(20,100,10),  col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()
legend(x=2, y=98, legend=c("Completo","Sub-conjunto"), fill=c(rgb(0,0,0),rgb(0.7,0.7,0.7)))

plot(seq(0.2,1,0.1),percCertos/100)

#Frequ?ncia da ocorrencia de cada valor de probabilidade
compProbs.knn <- rbind(table(probs.round.knn),table(probsSub.round.knn))
barplot(compProbs.knn, beside=T, space=c(0,4), legend=c("Completo","Sub-conjunto"), border=NA, col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()





##########################################################################
#2. SVM - MAQUINAS DE VECTORES DE SUPORTE

###################
#TODOS OS DADOS
classCult.svm <- svm(cultura ~ ., data=dadosTreino, probability=TRUE, kernel='radial', gamma=0.1, cost=3)
plot(classCult.svm, dadosTreino, B3_d3 ~ B4_d3)

#Tuning
classCult.tune.svm <- tune.svm(dadosTreino[,-1], dadosTreino[,1], gamma=seq(0, 1, by = .1));classCult.tune.svm
classCult.tune.svm <- tune.svm(dadosTreino[,-1], dadosTreino[,1], gamma=0.1, cost=seq(1, 10, by = 1));classCult.tune.svm
classCult.tune.svm <- tune.svm(dadosTreino[,-1], dadosTreino[,1], gamma=seq(0, 0.2, by=0.1), cost=seq(2, 4, by=1));classCult.tune.svm


#Validacao
validacao.svm <- predict(classCult.svm, dadosValidacao[,-1])
resultadoCult.svm <- table(pred = validacao.svm, true = dadosValidacao[,1]);resultadoCult.svm
resultadoCult.svm <- as.data.frame.matrix(resultadoCult.svm)
classDiag.svm <- diag(as.matrix(resultadoCult.svm));classDiag.svm
certos.svm <- sum(classDiag.svm)/sum(resultadoCult.svm);certos.svm

#Probabilidades
probs.svm <- predict(classCult.svm, dadosValidacao[,-1], probability=T)
probs.svm <- attr(probs.svm, "prob")
probs.svm <- as.vector(apply(probs.svm,1,max))
probs.round.svm <- round(probs.svm,1)
hist(probs.round.svm, breaks=14, col='blue')
table(pred=validacao.svm, true=dadosValidacao[,1], probs=probs.round.svm)


###################
#APENAS SUB-CONJUNTO DE DADOS
classCultSub.svm <- svm(cultura ~ ., data=dadosTreinoSub, probability=TRUE, kernel='radial', gamma=0.2, cost=2)
plot(classCultSub.svm, dadosTreinoSub, B3_d3 ~ B4_d3)

#Tuning
classCultSub.tune.svm <- tune.svm(dadosTreinoSub[,-1], dadosTreinoSub[,1], gamma=seq(0, 0.5, by = .1));classCultSub.tune.svm
classCultSub.tune.svm <- tune.svm(dadosTreinoSub[,-1], dadosTreinoSub[,1], gamma=0.2, cost=seq(1, 8, by = 1));classCultSub.tune.svm

#Validacao
validacaoSub.svm <- predict(classCultSub.svm, dadosValidacaoSub[,-1])
resultadoCultSub.svm <- table(pred = validacaoSub.svm, true = dadosValidacaoSub[,1]);resultadoCultSub.svm
resultadoCultSub.svm <- as.data.frame.matrix(resultadoCultSub.svm)
classDiagSub.svm <- diag(as.matrix(resultadoCultSub.svm));classDiagSub.svm
certosSub.svm <- sum(classDiagSub.svm)/sum(resultadoCultSub.svm);certosSub.svm

#Probabilidades
probsSub.svm <- predict(classCultSub.svm, dadosValidacaoSub[,-1], probability=T)
probsSub.svm <- attr(probsSub.svm, "prob")
probsSub.svm <- as.vector(apply(probsSub.svm,1,max))
probsSub.round.svm <- round(probsSub.svm,1)
hist(probsSub.round.svm, breaks=14, col='green')
table(pred=validacaoSub.svm, true=dadosValidacaoSub[,1], probs=probsSub.round.svm)

#So com P >= 60%
jjj <- table(pred=validacaoSub.svm, true=dadosValidacaoSub[,1], probs=probsSub.round.svm)
jjj <- as.matrix(jjj[,,5]+jjj[,,6]+jjj[,,7]+jjj[,,8]+jjj[,,9])
sum(diag(jjj))/sum(jjj)


###################
#Comparacao das probabilidades
compProbs.svm <- rbind(table(probs.round.svm),table(probsSub.round.svm))
barplot(compProbs.svm, beside=T, space=c(0,5), legend=c("Completo","Sub-conjunto"), xlab="Probabilidade", ylab="Numero de parcelas", border=NA, col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()
#Conclusao: com dados nao-completos, existem mais classificacaes com 100% de certeza



