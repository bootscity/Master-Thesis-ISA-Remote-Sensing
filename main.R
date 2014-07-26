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
                                            prefixo="CORR_15.06",
                                            anos=2005,
                                            corrige=FALSE)

#Obter lista de todos os dados
rm(listaDados)
listaDados <- constroiListaDados(anos=2005)

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Seleccionar apenas aquelas que representam 95% da area e reclassificar
listaTodasParcelas <- listaTodasParcelas[listaTodasParcelas$cultura %in% cultInfluentes,]

novasClasses <- as.data.frame(cbind(cultInfluentes,c(1,2,3,4,5,6,7,8,5,9,10,11,12)))
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
dadosClassificadores <- listaTodasParcelas[,c(5,7:12,15:20,23:28,31:36,39:44,47:52)]
dadosClassificadores$cultura <- as.factor(dadosClassificadores$cultura)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4.2 Sub-conjunto do conjunto completo ----

#Usando a funcao trim.matrix do package subselect
dadosTrim <- dadosClassificadores
criterioExclusaoVars <- 0.02
classTrim <- trim.matrix(cor(dadosTrim[,-1]), criterioExclusaoVars);classTrim

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. 6. 7. ESCOLHA/TREINO/VALIDACAO DOS CLASSIFICADORES                               ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5.1 KNN - K vizinhos mais proximos ----

#%%%%%%%%%%%%%%%%%%%%%%%%
#Todos os dados ----

#Tuning
KNN.tune.comp <- tune.knn(dadosClassificadores[,-1], dadosClassificadores[,1], k=c(1,20))
k.comp <- KNN.tune[1][[1]][1,1];k.comp

#Validacao cruzada
KNN.comp.cruz <- validacaoCruzada(n = 10, tipo = 'KNN', dados = dadosClassificadores, lambda = 0.8, k = k.comp)



#%%%%%%%%%%%%%%%%%%%%%%%%
#Apenas sub-conjunto de dados ----

#Tuning
KNN.tune.sub <- tune.knn(dadosClassificadoresSub[,-1], dadosClassificadoresSub[,1], k=c(1,20))
k.sub <- KNN.tune[1][[1]][1,1];k.sub

#Validacao cruzada
KNN.sub.cruz <- validacaoCruzada(n = 10, tipo = 'KNN', dados = dadosClassificadoresSub, lambda = 0.8, k = k.sub)




#KNN.comp <- treinaValida(tipo = 'KNN', treino = dadosTreino, valid = dadosValidacao, k = 10)
#qi.KNN.comp <- fixaLimiteQ(classificador = KNN.comp, lambda = 0.8)
#KNN.sub  <- treinaValida(tipo = 'KNN', treino = dadosTreinoSub, valid = dadosValidacaoSub, k = 10)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5.2 SVM - Maquinas de vectores de suporte ----

#%%%%%%%%%%%%%%%%%%%%%%%%
#Todos os dados ----

#Tuning: usa de o do sub-conjunto porque este nao faz...

#Validacao cruzada
SVM.comp.cruz <- validacaoCruzada(n = 10, tipo = 'SVM', dados = dadosClassificadores, lambda = 0.8, gamma = gamma.sub, cost = cost.sub)



#%%%%%%%%%%%%%%%%%%%%%%%%
#Apenas sub-conjunto de dados ----

#Tuning
SVM.tune.sub <- tune.svm(dadosClassificadoresSub[,-1], dadosClassificadoresSub[,1], gamma=1:3)  #1
SVM.tune.sub <- tune.svm(dadosClassificadoresSub[,-1], dadosClassificadoresSub[,1], gamma=c(0.2, 0.6, 1))  #0.6
SVM.tune.sub <- tune.svm(dadosClassificadoresSub[,-1], dadosClassificadoresSub[,1], gamma=c(0.4, 0.6, 0.8)) #0.4
SVM.tune.sub <- tune.svm(dadosClassificadoresSub[,-1], dadosClassificadoresSub[,1], gamma=c(0.3, 0.4, 0.5)) #0.5
SVM.tune.sub <- tune.svm(dadosClassificadoresSub[,-1], dadosClassificadoresSub[,1], gamma=0.5, cost=1:3) #3
SVM.tune.sub <- tune.svm(dadosClassificadoresSub[,-1], dadosClassificadoresSub[,1], gamma=0.5, cost=c(2.5, 3, 3.5)) #3
SVM.tune.sub <- tune.svm(dadosClassificadoresSub[,-1], dadosClassificadoresSub[,1], gamma=0.5, cost=c(2.8, 3, 3.2)) #3.2
SVM.tune.sub <- tune.svm(dadosClassificadoresSub[,-1], dadosClassificadoresSub[,1], gamma=c(0.5,0.6), cost=c(3.2, 3.5))

gamma.sub <- SVM.tune.sub[1][[1]][1,1];gamma.sub
cost.sub <- SVM.tune.sub[1][[1]][1,2];cost.sub

#Validacao cruzada
SVM.sub.cruz <- validacaoCruzada(n = 10, tipo = 'SVM', dados = dadosClassificadoresSub, lambda = 0.8, gamma = gamma.sub, cost = cost.sub)



#VERIFICACOES
# as.vector(SVM.sub.cruz$result$percsDecs[,'0.8'])
# as.vector(c(SVM.sub.cruz$result$percTotDec, SVM.sub.cruz$result$percParcDec))
# 
# SVM.sub.cruz$resultTodos$valid1$nParcDec[4]/SVM.sub.cruz$resultTodos$valid1$nParcTot[4]
# SVM.sub.cruz$resultTodos$valid2$nParcDec[4]/SVM.sub.cruz$resultTodos$valid2$nParcTot[4]
# SVM.sub.cruz$resultTodos$valid3$nParcDec[4]/SVM.sub.cruz$resultTodos$valid3$nParcTot[4]
# SVM.sub.cruz$resultTodos$valid4$nParcDec[4]/SVM.sub.cruz$resultTodos$valid4$nParcTot[4]
# SVM.sub.cruz$resultTodos$valid5$nParcDec[4]/SVM.sub.cruz$resultTodos$valid5$nParcTot[4]
# SVM.sub.cruz$resultTodos$valid6$nParcDec[4]/SVM.sub.cruz$resultTodos$valid6$nParcTot[4]
# SVM.sub.cruz$resultTodos$valid7$nParcDec[4]/SVM.sub.cruz$resultTodos$valid7$nParcTot[4]
# SVM.sub.cruz$resultTodos$valid8$nParcDec[4]/SVM.sub.cruz$resultTodos$valid8$nParcTot[4]
# SVM.sub.cruz$resultTodos$valid9$nParcDec[4]/SVM.sub.cruz$resultTodos$valid9$nParcTot[4]
# SVM.sub.cruz$resultTodos$valid10$nParcDec[4]/SVM.sub.cruz$resultTodos$valid10$nParcTot[4]
# 
# SVM.sub.cruz$resultTodos$valid7$qi[4]
# SVM.sub.cruz$resultTodos$valid7$todasConfiancas
# 
# (SVM.sub.cruz$resultTodos$valid1$nParcDec/SVM.sub.cruz$resultTodos$valid1$nParcTot)==as.vector(SVM.sub.cruz$resultTodos$valid1$percsDecs[,'0.8'])[-1]



#resultValidCruz <- list(KNN.comp.cruz, KNN.sub.cruz, SVM.comp.cruz, SVM.sub.cruz)
#SVM.comp <- treinaValida(tipo = 'SVM', treino = dadosTreino, valid = dadosValidacao, gamma = 0.1, cost = 3)
#qi.SVM.comp <- fixaLimiteQ(classificador = SVM.comp, lambda = 0.8)
#SVM.sub  <- treinaValida(tipo = 'SVM', treino = dadosTreinoSub, valid = dadosValidacaoSub, gamma = 0.2, cost = 2)
#qi.SVM.sub <- fixaLimiteQ(classificador = SVM.sub, lambda = 0.8)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Fazer a validacao ao contrario para ter a certeza que esta tudo bem
set.seed(1)
amostra <- sample(1:nrow(dadosClassificadores), ceiling(nrow(dadosClassificadores)/2))
dadosTreinoSub <- dadosClassificadoresSub[amostra,]
dadosValidacaoSub <- dadosClassificadoresSub[-amostra,]

validValid <- treinaValida(tipo = 'SVM', treino = dadosValidacaoSub, valid = dadosTreinoSub, gamma = gamma.sub, cost = cost.sub)
vv <- c()
for (i in 1:12)
{
  Rq <- validValid$result[validValid$result$prob1 >= SVM.sub.cruz$result$qi[i],]
  matriz <- as.data.frame.matrix(squareTable(x = Rq$class1, y = Rq$verdade))  #para garantir que e quadrada
  confiancas <- diag(as.matrix(matriz))/rowSums(matriz);
  vv[i] <- confiancas[i]
}
round(vv,3)