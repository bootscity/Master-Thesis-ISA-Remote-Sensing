#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Master's Thesis - Remote Sensing                      %
# Environmental Engineering - ISA/UL - Lisbon, Portugal %
# (c) 2014 by Jonas Schmedtmann                         %
#                                                       %
# "STAFF" SCRIPT                                        %
#                                                       %
# This script contains random chunks of code like       %
# experiments, no longer used code or just informations %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   ABORDAGEM INICIAL:    CRITERIO 1 --- UTILIZACAOO DE DIFERENTES CLASSIFICADORES    ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Classificadores utilizados:
# 1. LDA - An??lise discriminante linear
# 2. KNN - k vizinhos mais pr??ximos
# 3. SVM - M??quinas de vectores de suporte
# 4. ??rvores de decis??o


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#0. PREPARA????O E VISUALIZA????O DOS DADOS

#Cria????o dos dados com as seguintes caracteristicas:
#  - pastagem ?? factor que ?? 0 se nao for pastagem, 1 se for pastagem permanente e 2 se for pastagem pobre
#  - usam-se dados de todas as bandas nas primeiras 6 datas do ano
#  - metade dos dados vai ser para treino dos classificadores e metade para valida????o

pastagem <- c()
seq <- c()
dadosClassificadores <- listaTodasParcelas[,c(5,7:12,15:20,23:28,31:36,39:44,47:52)]

for(i in 1:length(dadosClassificadores$cultura))
{
  if(dadosClassificadores$cultura[i] == 1)
    pastagem[i] <- 1
  else
    pastagem[i] <- 0
}
dadosClassificadoresPast <- NA
dadosClassificadoresPast <- cbind(pastagem,dadosClassificadores)
dadosClassificadoresPast <- dadosClassificadoresPast[,-2]
dadosClassificadoresPast$pastagem <- as.factor(as.numeric(as.vector(dadosClassificadoresPast$pastagem)))

#Amostra a usar para o modelo - metade para treino e metade para valida????o
set.seed(23)
amostra <- sample(1:nrow(dadosClassificadoresPast), ceiling(nrow(dadosClassificadoresPast)/2))

#Seleccao dos dados
dadosTreino <- dadosClassificadoresPast[amostra,]
dadosValidacao <- dadosClassificadoresPast[-amostra,]

#Visualiza????o dos dados 2D
plot(c(),xlim=range(dadosTreino$B4_d3),ylim=range(dadosTreino$B4_d5),xlab="NDVI DOY 125",ylab="NDVI DOY 189")
text(dadosTreino$B4_d3[dadosTreino$pastagem == 0],dadosTreino$B4_d5[dadosTreino$pastagem == 0],dadosTreino$pastagem[dadosTreino$pastagem == 0])
text(dadosTreino$B4_d3[dadosTreino$pastagem == 1],dadosTreino$B4_d5[dadosTreino$pastagem == 1],dadosTreino$pastagem[dadosTreino$pastagem == 1], col="red")

#Visualiza????o dos dados 3D
library(rgl)
plot3d(x=dadosTreino$B4_d2[dadosTreino$pastagem == 0],
       y=dadosTreino$B4_d3[dadosTreino$pastagem == 0],
       z=dadosTreino$B4_d5[dadosTreino$pastagem == 0],
       xlab="B4 DOY 61",
       ylab="B4 DOY 125",
       zlab="B4 DOY 189",
       size=3)
points3d(x=dadosTreino$B4_d2[dadosTreino$pastagem == 1],
         y=dadosTreino$B4_d3[dadosTreino$pastagem == 1],
         z=dadosTreino$B4_d5[dadosTreino$pastagem == 1],
         col="red")

#Resultados da optimiza????o da LDA
plot3d(x=dadosTreino$B3_d5[dadosTreino$pastagem == 0],
       y=dadosTreino$B5_d5[dadosTreino$pastagem == 0],
       z=dadosTreino$B6_d5[dadosTreino$pastagem == 0],
       xlab="B3_d5",
       ylab="B5_d5",
       zlab="B6_d5",
       size=3)
points3d(x=dadosTreino$B3_d5[dadosTreino$pastagem == 1],
         y=dadosTreino$B5_d5[dadosTreino$pastagem == 1],
         z=dadosTreino$B6_d5[dadosTreino$pastagem == 1],
         col="red")

#Matriz de gr??ficos
pairs(~B4_d2+B4_d3+B4_d5, data=dadosTreino)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#1. LDA - ANALISE DISCRIMINANTE LINEAR 

#Treino do classificador
class.lda <- lda(dadosTreino[,-1],grouping=dadosTreino$pastagem);class.lda
plot(class.lda,col=(as.numeric(dadosTreino[,1])))

#Valida????o
validacao.lda <- predict(class.lda, dadosValidacao[,-1])
resultado.lda <- table(pred = validacao.lda$class, true = dadosValidacao[,1]);resultado.lda
classAgreement(resultado.lda)



#Optimiza????o com package subselect
class.lda.Hmat <- ldaHmat(dadosTreino[,-1], dadosTreino$pastagem)
class.lda.opt <- anneal(class.lda.Hmat$mat,kmin=3,kmax=5,H=class.lda.Hmat$H,crit="RM")
class.lda.opt$subsets

#nomes das 3 variaveis mais importantes
colnames(dadosTreino[,class.lda.opt$subsets])


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#2. SVM - M??QUINAS DE VECTORES DE SUPORTE

#Treino do classificador
class.svm <- svm(pastagem ~ ., data=dadosTreino, cost=100, gamma=0.9)
plot(class.svm, dadosTreino, B1_d3 ~ B5_d6)

#Optimiza????o - tentar com MAIS hipoteses
class.svm.tune <- tune.svm(pastagem ~ ., data=dadosTreino, gamma=seq(.1, 1, by = .1), cost = seq(100,1000, by = 100));class.svm.tune

#resultados: gamma=0.9 e cost=100

#Valida????o
validacao.svm <- predict(class.svm, dadosValidacao[,-1])
resultado.svm <- table(pred = validacao.svm, true = dadosValidacao[,1]);resultado.svm
classAgreement(resultado.svm)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#3. KNN - K VIZINHOS MAIS PR??XIMOS

#Treino e valida????o simultanea
class.knn <- knn(dadosTreino[,-1], dadosValidacao[,-1], dadosTreino[,1], k=16)
resultado.knn <- table(pred=class.knn, true=dadosValidacao[,1]);resultado.knn
classAgreement(resultado.knn)

#Optimiza????o
class.knn.tune <- tune.knn(dadosTreino[,-1], dadosTreino[,1], k = 1:20);class.knn.tune


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#4. ??RVORES DE DECIS??O

#Ajuste inicial da ??rvore
class.tree <- tree(pastagem ~ ., data=dadosTreino, control= tree.control(nobs=nrow(dadosTreino), minsize=1))
plot(class.tree);text(class.tree,cex=.8)

#Valida????o
validacao.tree <- predict(class.tree,type="class")
resultado.tree <- table(pred = validacao.tree, true = dadosTreino[,1]);resultado.svm
classAgreement(resultado.tree)

#Determina????o da melhor ??rvore (com menor estimativa de erro de classifica????o)
cvDev <- cv.tree(class.tree)$dev
opt <- which(cvDev == min(cvDev))
class.tree.prune <- prune.tree(class.tree,best=opt)
plot(class.tree.prune);text(class.tree.prune,cex=.8)



#OUTRA ABORDAGEM COM OUTRO PACKAGE QUE PARECE SER MELHOR
#Treino do classificador
class.tree2 <- ctree(pastagem ~ ., data=dadosTreino)
plot(class.tree2)

#Valida????o
validacao.tree2 <- predict(class.tree2, dadosValidacao[,-1])
resultado.tree2 <- table(pred = validacao.tree2, true = dadosValidacao[,1]);resultado.tree2
classAgreement(resultado.tree2)


#mais uma abordagem
class.tree3 <- rpart(pastagem ~ ., data=dadosTreino)
plot(class.tree3)
text(class.tree3,cex=.8)

#Valida????o
validacao.tree3 <- predict(class.tree3, dadosValidacao[,-1],type="class")
resultado.tree3 <- table(pred = validacao.tree3, true = dadosValidacao[,1]);resultado.tree3
classAgreement(resultado.tree3)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#5. REGRESS??O LOG??STICA MULTINOMIAL

#Com package NNET
library(nnet)
class.multi <- multinom(pastagem ~ ., data=dadosTreino)
summary(class.multi)

validacao.multi <- predict(class.multi, dadosValidacao[,-1])
resultado.multi <- table(pred = validacao.multi, true = dadosValidacao[,1]);resultado.multi





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   CENAS DE "ANALISE EXPLORATORIA"                                                   ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
#   EXPERIENCIAS                                                                      ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Experiencia de estimacao das densidades de probabilidade de cada classe de
#culturas com metodo KNN, do tipo P(x|C_i), ao contrario do tradicional P(C_i|x)

#QUESTOES:
# - Como vamos poder validar isto?
# - Como fazer uma verdadeira 'funcao' com isto para podermos fazer testes de hipoteses?

classe <- 1
k <- 10      #Numero de vizinhos a considerar - resultado do KNN.TUNE
d <- ncol(dadosTreinoSub[,-1])      #Numero de dimensoes
dadosTreinoClasse <- dadosTreinoSub[dadosTreinoSub$cultura == classe,]
dadosValidacaoClasse <- dadosValidacaoSub[dadosValidacaoSub$cultura == classe,]
n <- nrow(dadosTreinoClasse)     #Numero de observacoes na amostra de treino???

probs <- c()
for (i in 1:nrow(dadosValidacaoClasse))
{
  observacao <- dadosValidacaoClasse[i,-1]
  distancias <- sort(as.matrix(dist(rbind(observacao,dadosTreinoClasse[,-1])))[-1,1])
  distKVizinho <- distancias[k]
  
  #Calculo da probabilidade
  if(d%%2 == 0)
    cD <- pi^(d/2)/factorial(d/2)
  else
    cD <- pi^(d/2)/gamma(d/2+1)
  
  p <- k/(n*distKVizinho*cD)
  probs[i] <- p
}

plot(sort(probs))
hist(probs,breaks=20)



#Experiencia com 2 dimensoes e pontos aleatorios
x <- sample(seq(1,20,0.1),10,replace=T)
y <- sample(seq(1,20,0.1),10,replace=T)
novo.x <- 9
novo.y <- 7
plot(x,y)
text(x,y+0.4,2:11)
points(novo.x,novo.y,col='red')
text(novo.x,novo.y+0.4,1,col='red')

distancias <- sort(as.matrix(dist(rbind(cbind(novo.x,novo.y),cbind(x,y))))[-1,1])
distKVizinho <- distancias[5]

draw.circle(novo.x,novo.y,distKVizinho,border='red')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Experiencias com NDVI's (apenas para ano 2005)

#Fazer plot
datasGrafico <- datas
substr(datasGrafico,0,1) <- "0"
datasGrafico <- as.numeric(datasGrafico)

#plot(datasGrafico,NDVIs[1,],ylim=c(-1,1),col=100,type='l',lty=1,pch=16,ylab="NDVI", xlab="DOY 2005")
plot(datasGrafico,listaTodasParcelas[1,55:62],ylim=c(-1,1),col=1,type='l',lty=1,ylab="NDVI", xlab="DOY 2005")
for(i in 1:nrow(listaTodasParcelas))
  lines(datasGrafico,listaTodasParcelas[i,55:62],ylim=c(-1,1),col=i,lty=1,pch=16)


#distinguir entre sequeiro e regadio
todasParcelasSEQ <- listaTodasParcelas[listaTodasParcelas$seq == TRUE,]
todasParcelasREG <- listaTodasParcelas[listaTodasParcelas$seq == FALSE,]

plot(datasGrafico,todasParcelasSEQ[1,55:62],ylim=c(-1,1),col='burlywood2',type='l',lty=1,ylab="NDVI", xlab="DOY 2005")
for(i in 1:100)
  lines(datasGrafico,todasParcelasSEQ[i,55:62],ylim=c(-1,1),col='burlywood2',lty=1,pch=16)
for(i in 1:100)
  lines(datasGrafico,todasParcelasREG[i,55:62],ylim=c(-1,1),col='darkolivegreen3',lty=1,pch=16)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Experi??ncias com outros indices de vegeta??ao

listaTodasParcelas$NDVI_1[50]
(listaTodasParcelas$B4_d1[50]-listaTodasParcelas$B3_d1[50])/(listaTodasParcelas$B4_d1[50]+listaTodasParcelas$B3_d1[50])

nir <- listaTodasParcelas$B4_d4
red <- listaTodasParcelas$B3_d4
NDVId4 <- (nir-red)/(nir+red)
NDVId4.orig <- listaTodasParcelas$NDVI_4

plot(NDVId4,NDVId4.orig)
abline(lm(NDVId4.orig~NDVId4),col='red')


MSAVId4 <- (2*nir+1-sqrt((2*nir+1)^2-8*(nir-red)))/2
plot(NDVId4,MSAVId4)
abline(lm(MSAVId4~NDVId4),col='red')
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Experimentar diferen??a entre usar apenas os interiores das parcelas e usar parcelas inteiras
testeTotInt <- matrix(nrow=0,ncol=5)
for(i in 1:length(listaDados))
{
  for(j in 1:length(listaDados[[i]][[1]]))
  {
    parc <- listaDados[[i]][[1]][[j]]
    conteudoTotal <- obtemConteudoParcela(parc,apenasInterior=FALSE)
    conteudoInterior <- obtemConteudoParcela(parc,apenasInterior=TRUE)
    nTotal <- length(conteudoTotal[[1]][,1])
    nInterior <- length(conteudoInterior[[1]][,1])
    sdTotal <- sd(conteudoTotal[[1]][,1])
    sdInterior <- sd(conteudoInterior[[1]][,1])
    area <- parc$atributos$AREA
    linha <- c(nTotal,sdTotal,nInterior,sdInterior,area)
    testeTotInt <- rbind(testeTotInt,linha)
  } 
  print(i)
}

mean(testeTotInt[,1])   #M??dia do numero de pixeis quando se considera tambem a fronteira
mean(testeTotInt[,3])   #M??dia do numero de pixeis quando se considera apenas o interior

mean(testeTotInt[,1])/mean(testeTotInt[,5]/10000)   #M??dia de pixeis/ha -- fronteira
mean(testeTotInt[,3])/mean(testeTotInt[,5]/10000)   #M??dia de pixeis/ha -- interior

mean(testeTotInt[,2],na.rm=T)   #Media dos desvios padrao -- fronteira
mean(testeTotInt[,4],na.rm=T)   #Media dos desvios padrao -- interior

length(testeTotInt[testeTotInt[,1] == 0,1])/length(testeTotInt[,1])  #% de parcelas sem pixeis -- fronteira
length(testeTotInt[testeTotInt[,3] == 0,3])/length(testeTotInt[,3])  #% de parcelas sem pixeis -- interior

table(as.vector(testeTotInt[,2] > testeTotInt[,4]))   #Numero de vezes que sd ?? maior utilizando a fronteira
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Experiencias para encontrar bom algoritmo para toma decisão final
xx <- SVM.sub$tabClassProb[,,9]+SVM.sub$tabClassProb[,,8]+SVM.sub$tabClassProb[,,7]
diag(as.matrix(xx))/rowSums(xx)   #Esta é a user's accuracy (linhas)

#grafico de confiança vs. num ou % de parcelas em que se toma uma decisao 

#Juntar codigos 1,2,6 na tabela de erro
SVM.comp$tabClass
condensaMatriz(SVM.comp$tabClass, c(2,5,8,9))

SVM.comp$result[1:10,]


#Tentativa de escolha da P ideal: classificacoes correctas se escolhermos apenas classificaoes com P >= x
correcCum <- c()
correcCumCult <- matrix(data=0, nrow=length(SVM.sub$tabClassProb[1,1,]),ncol=dim(SVM.sub$tabClassProb[,,1])[2])
for (i in 1:length(SVM.sub$tabClassProb[1,1,]))
{
  m <- matrix(data=0, nrow=dim(SVM.sub$tabClassProb[,,1])[1], ncol=dim(SVM.sub$tabClassProb[,,1])[2])
  for (j in i:length(SVM.sub$tabClassProb[1,1,]))
    m <- m + SVM.sub$tabClassProb[,,j]
  
  #Correccao total
  correcCum[i] <- cc(m)
  
  #Correccao para cada cultura
  for (j in 1:dim(SVM.sub$tabClassProb[,,1])[2])
    correcCumCult[i,j] <- m[j,j]/sum(m[,j])
}

#Correccao acumulada total
plot(seq(0.2,1,0.1),correcCum,xlim=c(1,0.2),type='l',ylim=c(0.2,1)xlab='Probabilidade maxima',ylab='Classificacoes correctas')
abline(h=0.8, col='red')

#Correccao acumulada para cada cultura
plot(seq(0.2,1,0.1),correcCumCult[,1],xlim=c(1,0.2),type='l',ylim=c(0,1),xlab='Probabilidade maxima',ylab='Classificacoes correctas')
abline(h=0.8, col='red')
for(i in 2:ncol(correcCumCult))
  lines(seq(0.2,1,0.1),correcCumCult[,i],xlim=c(1,0.2),col=i)


#SUGESTAO:
#Podemos fixar um nivel de qualidade desejado (eixo do y, exemplo de 80% de classificacoes correctas) e fazer uma curva deste genero para cada classe e ver que valor de probabilidade é o minimo aceitavel para assegurar o nivel desejado.
#Desta forma fica fixo um valor de P=x
#Talvez perguntar ao IFAP qual o nivel que eles pretendem?
#Depois dá jeito tambem ver a quantidade de parcelas que podemos de facto classsifcar fixando um valor de probabilidade x.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Dados SHP para professor
x <- c(494261, 524261, 524261, 494261, 494261)
y <- c(4301509, 4301509, 4331509, 4331509, 4301509)
cords <- cbind(x,y)
cdgProf <- carregaRecortaDadosEspaciais(cords, PROJ4.UTM)

plot(cdgProf$area, add=T)
plot(cdgProf$parc2005, add=T)
plot(maio$b3)

#Dados RASTER para professor
maio <- corrigeLeConjuntoImagensLandsat(conjuntoPath = paste0(CAMINHO.LANDSAT,"/LE72040332005125ASN00"), areaEstudo = cdgProf$area, prefixo = 'CORR_15.06', corrige = FALSE)

for(i in 1:length(maio))
{
  maio[[i]] <- crop(maio[[i]], cords)
  maio[[i]][maio[[i]] == 0] <- NA
}

junho <- corrigeLeConjuntoImagensLandsat(conjuntoPath = paste0(CAMINHO.LANDSAT,"/LE72040332005157EDC00"), areaEstudo = cdgProf$area, prefixo = 'CORR_15.06', corrige = FALSE)

for(i in 1:length(maio))
{
  junho[[i]] <- crop(junho[[i]], cords)
  junho[[i]][junho[[i]] == 0] <- NA
}

save(cdgProf, file='cdgProf')
save(maio, file='maio')
save(junho, file='junho')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   CAIXOTE DO LIXO                                                                   ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Bom livro em PT
#http://cran.r-project.org/doc/contrib/Torgo-ProgrammingIntro.pdf
#http://www.statmethods.net/advgraphs/parameters.html


#Obter conteudo de uma parcela para diferentes datas e bandas
parc  <-  listaDados$n1541$a2005$p1432302680003.301.00
conteudoParcIn <- obtemConteudoParcela(parc)
conteudoParcInSemViz <- obtemConteudoParcela(parc,excluiVizinhos=TRUE)
conteudoParc <- obtemConteudoParcela(parc,apenasInterior=FALSE)

plot(parc$DR$d45[,,7], parc$DR$d45[,,8])
lines(parc$coordsPoly, col="red")


#Package kknn
library(kknn)
classTodCult.kknn <- kknn(dadosTreino2[,-1], dadosValidacao2[,-1], dadosTreino2[,1], k=10)


#Areas das parcelas
hist(listaTodasParcelas$area,breaks=20)
boxplot(listaTodasParcelas$area)
plot(sort(listaTodasParcelas$area))
summary(listaTodasParcelas$area)

which(listaTodasParcelas$area==max(listaTodasParcelas$area)) #Maior parcela de todas
listaTodasParcelas$cultura[which(listaTodasParcelas$area==max(listaTodasParcelas$area))]  #Ocupa????o da maior parcela de todas

par(col="black")
par(opar)

plot(listaDados$n197$a2003$p1162347095001.301.00$DR$d29$X, listaDados$n197$a2003$p1162347095001.301.00$DR$d29$Y)
lines(listaDados$n197$a2003$p1162347095001.301.00$coordsPoly, col="red")


#Frequencias para assinaturas espectrais
frequencias <- c(0.485, 0.56, 0.66, 0.835, 1.65, 2.22)


#Aceder a dados geograficos importados
cdg$parc2003@polygons[[1]]@Polygons[[1]]@area

cdg$area@polygons[[1]]@Polygons[[1]]@coords

plot(cdg$parc2003, col="red")
plot(cdg$conc, col=rgb(1,0,0,alpha=0), add=TRUE)


#Informaçoes de desempenho de tarefas
proc.time()
system.time()


#Usar isto para percorrer todas os NIFAPs
for(i in 1:length(listaDados))
  for(j in 1:length(listaDados[[i]][[1]]))
  {
    parc <- listaDados[[i]][[1]][[j]]
    .
    .
    .
  } 



#Subselect mais complicado para excluir variaveis
classHmat <- ldaHmat(dadosTreino[,-1], dadosTreino$cultura)
classSubselect <- eleaps(classHmat$mat,kmin=10,kmax=10,H=classHmat$H,r=classHmat$r,crit="Tau2",timelimit=600)
classSubselect



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   TESTES  DE APRENDIZAGEM DE NOVAS FUNÇÕES DO R                                     ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%
#USING APPLY
j <- c(3,4,5,-2,4,-3)
j[j<0]

m <- matrix(data=cbind(rnorm(20, 0), rnorm(20, 2), rnorm(20, 5)), nrow=20, ncol=3)
apply(m, 1, mean) #applying across rows of matrix
apply(m, 2, mean) #applying across cols of matrix
apply(m, 2, function(x) length(x[x<0]))
apply(m, 2, is.vector)
apply(m, 2, function(j) mean(j[j>0]))

#USING LAPPLY
sapply(1:3, function(x) x^2)  #Returns vector
lapply(1:3, function(x) x^2)  #Returns lsit

#%%%%%%%%%%%%%%%%%%%%%%%%%%
#more stuff
colnames(m) <- c('method1', 'method2', 'method3')
head(m)
m[,'method1']
df <- as.data.frame(m);df
df.stack <- stack(df);df.stack
unstack(df.stack)


x <- 1:10
names(x) <- letters[1:10]
x['f']
x[x > 5] <- 0;x

y <- seq(1, 30, by=3);y
which(y %% 2 == 0)




sales <- expand.grid(country = c('USA', 'UK', 'FR'), product = c(1, 2, 3))
sales$revenue <- rnorm(dim(sales)[1], mean=100, sd=10)
usd2eur <- 1.434
transform(sales, euro = revenue * usd2eur)
subset(sales, 
       country == 'USA' & product %in% c(1, 2), 
       select = c('product', 'revenue'))

#DO.CALL
foo <- list()
foo[[1]] <- data.frame(a=1:5, b=11:15)
foo[[2]] <- data.frame(a=101:105, b=111:115)
foo[[3]] <- data.frame(a=200:210, b=300:310)
do.call(rbind, foo)



#%%%%%%%%%%%%%%%%%%%%%%%%%%
#using apply for lists
mylist <- list()
for(i in 1:5) {
  mylist[[i]] <- list()
  for(j in 1:5) {
    mylist[[i]][[j]] <- i*j
  }
}

sapply(mylist[[1]], identity)
sapply(mylist, function(x) lapply(x, identity))
str(mylist)

mylist<-list(1:3, "dog", TRUE)
sapply(mylist, identity)
unlist(sapply(mylist, identity))












#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#DataFrame com os valores de V, P0, P1, P2
solidos <- read.xlsx2("s.xlsx", sheetIndex=2, colClasses=rep('numeric',4))


SST <- 10^6*(solidos$P1 - solidos$P0)/solidos$V
SSF <- 10^6*(solidos$P2 - solidos$P0)/solidos$V;SSF
SSV <- SST - SSF
solidos <- cbind(solidos, SST, SSF, SSV)

medSST <- c()
for(i in seq(1, nrow(solidos), 2))
  medSST <- cbind(medSST, mean(c(SST[i], SST[i+1]), na.rm = T))

date <- c('A', 'A', 'B', 'B', 'B', 'C', 'C', 'D', 'A', 'A', 'A')
SST <- round(runif(11, 0, 3), 1)
j <- data.frame(date,SST);j
id <- rle(j$date)
medias <- tapply(j$SST, rep(1:length(id$lengths), id$lengths), mean) 
data.frame(date = id$values, medias)

