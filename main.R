#########################################################
# Master Thesis - Remote Sensing                        #
# Environmental Engineering - ISA/UL - Lisbon, Portugal #
# (c) 2014 by Jonas Schmedtmann                         #
#########################################################



##########################################################################################
# 1. AQUISICAO DE DADOS E PROCESSAMENTO                                                  #
##########################################################################################

#Clear console
cat("\014")

#Iniciar
source('functions.R')
init()

#Definir area de estudo, carregar dados, recortar e projectar
areaX<-c(508262,528870,511269,490710,508262)
areaY<-c(4365569,4362462,4288360,4291713,4365569)
coordsArea<-cbind(areaX,areaY)
cdg<-carregaRecortaDadosEspaciais(coordsArea,PROJ4.UTM)

#Obter todas as imagens
rm(todasImagens)
todasImagens<<-constroiListaImagensLandsat(CAMINHO.LANDSAT, cdg$area, "CORR_12.05", anos=2005, corrige=TRUE)
#todasImagens<<-constroiListaImagensLandsat(CAMINHO.LANDSAT, cdg$area, "CORR_12.05", anos=2005, corrige=FALSE)

#Obter lista de todos os dados
rm(listaDados)
listaDados<-constroiListaDados(anos=2005)

#Obter conteudo relevante para TODAS as parcelas, em forma de data.frame
listaTodasParcelasIniciais<-constroiTodasParcelas()
listaTodasParcelasIniciaisSD<-constroiTodasParcelas(funcao=sd)



##########################################################################################
# 2. SELECCAO DE DADOS                                                                   #
##########################################################################################

##########################################################################
#Excluir a partida codigos que fazem sentido considerar
codExclusao<-c(87,88,666)
listaTodasParcelas<-listaTodasParcelasIniciais[!(listaTodasParcelasIniciais$cultura %in% codExclusao),]


##########################################################################
#Determinar as culturas que ocupam maior parte da ?rea (90% ou 95%)

#Todas as culturas com respectiva area total
dados<-data.table(listaTodasParcelas)
areasCulturas<-dados[,list(area=sum(area),numParc=length(area)),by=cultura]
areasCulturas<-areasCulturas[order(areasCulturas$area, decreasing = TRUE),]
areasCulturas<-cbind(areasCulturas,cumsum(areasCulturas$area))
areaTotal<-sum(areasCulturas$area)

#Culturas acima e abaixo dos limiares
plot(1:nrow(areasCulturas),areasCulturas$V2)
abline(h=areaTotal*0.9,col='blue')
abline(h=areaTotal*0.95,col='orange')
abline(h=areaTotal*0.98,col='red')

cultInfl<-areasCulturas[areasCulturas$V2 < areaTotal*0.95,]
fraccaoInfl<-sum(cultInfl$area)/areaTotal
cultInfluentes<-cultInfl$cultura


##########################################################################
#Seleccionar apenas aquelas que representam 95% da ?rea e reclassificar
listaTodasParcelas<-listaTodasParcelas[listaTodasParcelas$cultura %in% cultInfluentes,]

novasClasses98<-as.data.frame(cbind(cultInfluentes,c(1,2,3,4,5,6,7,6,8,9,10,11,12,13,14,5,15,4,16,17,8)))
novasClasses<-as.data.frame(cbind(cultInfluentes,c(1,2,3,4,5,6,7,6,8,9,10,11,12,13,14,5)))
colnames(novasClasses)<-c('cultura','novaClasse')
nClasses<-length(table(novasClasses$novaClasse))

for(i in 1:length(listaTodasParcelas$cultura))
  listaTodasParcelas$cultura[i]<-novasClasses$novaClasse[which(novasClasses$cultura==listaTodasParcelas$cultura[i])]


##########################################################################
#Preparacao dos dados para os classificadores

dadosClassificadores<-listaTodasParcelas[,c(5,7:12,15:20,23:28,31:36,39:44,47:52)]
dadosClassificadores$cultura<-as.factor(dadosClassificadores$cultura)

#Amostra a usar para o modelo - metade para treino e metade para validacao
set.seed(23)
amostra<-sample(1:nrow(dadosClassificadores), ceiling(nrow(dadosClassificadores)/2))

#Seleccao dos dados
dadosTreino<-dadosClassificadores[amostra,]
dadosValidacao<-dadosClassificadores[-amostra,]



##########################################################################################
# 3. ANALISE EXPLORATORIA DOS DADOS                                                      #
##########################################################################################

##########################################################################
#ANOVA para averiguar se cultura e parcela afectam reflectancias

###################
#Primeira abordagem
dadosANOVA<-constroiDadosANOVA(ano=2005,data=5,banda=1,dimAmostra=15)
table(dadosANOVA$cultura)
plot(dadosANOVA)

#Modelo da ANOVA
dados.aov<-aov(reflectancias ~ cultura/parcela, data=dadosANOVA)
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
FValues<-c()
pValues<-c()
efeito<-c()
for(i in 7:length(listaTodasParcelas))
{
  parcelas.aov<-aov(listaTodasParcelas[,i] ~ as.factor(listaTodasParcelas$cultura))
  FV<-summary(parcelas.aov)[[1]][["F value"]][[1]]
  pV<-summary(parcelas.aov)[[1]][["Pr(>F)"]][[1]]
  
  FValues<-c(FValues,FV)
  pValues<-c(pValues,pV)
  if(pV <= 0.05) ef<-1 else ef<-0
  efeito<-c(efeito,ef)
}
nn<-colnames(listaTodasParcelas[,7:length(listaTodasParcelas)])
resultadoANOVA<-data.frame(cbind(nn,FValues,pValues,efeito))
resultadoANOVA$FValues<-as.numeric(as.character(resultadoANOVA$FValues))
resultadoANOVA<-resultadoANOVA[order(-resultadoANOVA$FValues),]


##########################################################################
#Assinaturas espectrais POR CULTURA

#Assinaturas todas
assinaturas<-matrix(nrow=0,ncol=length(BANDAS.LANDSAT))
for(i in 1:nrow(listaTodasParcelas))
{
  assinatura<-c()
  for(j in seq(0,6*7,8)+7)
    assinatura<-c(assinatura,listaTodasParcelas[i,j])
  
  assinaturas<-rbind(assinaturas,assinatura)
}

plot(BANDAS.LANDSAT,assinaturas[1,],type='l',ylim=c(0,0.35))
for(i in 1:length(assinaturas))
  lines(BANDAS.LANDSAT,assinaturas[i,],col=i)

#Assinaturas medias por ocupacao de solo
dados<-data.table(listaTodasParcelas)
assinaturasPorCulturaD3<-dados[,list(area=sum(area),
                                     B1_d3=mean(B1_d3),
                                     B2_d3=mean(B2_d3),
                                     B3_d3=mean(B3_d3),
                                     B4_d3=mean(B4_d3),
                                     B5_d3=mean(B5_d3),
                                     B6_d3=mean(B6_d3)),by=cultura]
assinaturasPorCulturaD3<-as.data.frame(assinaturasPorCulturaD3)

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
NDVIPorCulturaD4<-dados[,list(area=sum(area),
                              NDVI_1=mean(NDVI_1),
                              NDVI_2=mean(NDVI_2),
                              NDVI_3=mean(NDVI_3),
                              NDVI_4=mean(NDVI_4),
                              NDVI_5=mean(NDVI_5),
                              NDVI_6=mean(NDVI_6),
                              NDVI_7=mean(NDVI_7),
                              NDVI_8=mean(NDVI_8)),by=cultura]
NDVIPorCulturaD4<-as.data.frame(NDVIPorCulturaD4)

datasGrafico<-datas
substr(datasGrafico,0,1)<-"0"
datasGrafico<-as.numeric(datasGrafico)

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
h<-hist(listaTodasParcelas$B4_d1,breaks=50,prob=T)
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
cultGrafico<-c(1,3,9,13)
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
#Seleccao das variaveis
#Matriz de variancias-covariancias das assinaturas espectrais
corAssinaturas<-cor(listaTodasParcelas[,c(7:12,15:20,23:28,31:36,39:44,47:52)])

#Usando apenas o trim.matrix do package subselect
#Questao: Como escolher valor de tolval? Experimentar classificacoes com varios?
criterioExclusaoVars<-0.02
classTrim<-trim.matrix(cor(dadosTreino[,-1]),criterioExclusaoVars);classTrim

varsRetirar<-classTrim$numbers.discarded+1  #+1 porque na data.frame original, a coluna 1 ? a cultura
varsSobram<-1:length(dadosTreino)
varsSobram<-varsSobram[! varsSobram %in% varsRetirar]



#Nao vale a pena usar isto...
classHmat<-ldaHmat(dadosTreino[,-1], dadosTreino$cultura)
classSubselect<-eleaps(classHmat$mat,kmin=10,kmax=10,H=classHmat$H,r=classHmat$r,crit="Tau2",timelimit=600)
classSubselect


##########################################################################
#Preparacao dos NOVOS DADOS APOS SELECCAO DE VARIAVEIS
dadosClassificadoresSub<-dadosClassificadores[,varsSobram]

#Amostra a usar para o modelo - metade para treino e metade para validacao
set.seed(23)
amostra<-sample(1:nrow(dadosClassificadoresSub), ceiling(nrow(dadosClassificadoresSub)/2))

#Seleccao dos dados
dadosTreinoSub<-dadosClassificadoresSub[amostra,]
dadosValidacaoSub<-dadosClassificadoresSub[-amostra,]



##########################################################################################
# 5. 6. 7. ESCOLHA/TREINO/VALIDACAO DOS CLASSIFICADORES                                  #
##########################################################################################

#QUESTOES:
# - Como vamos poder validar isto sem as declaracoes? 


##########################################################################
#1. KNN - K VIZINHOS MAIS PROXIMOS

###################
#TODOS OS DADOS
classCult.knn<-knn(dadosTreino[,-1], dadosValidacao[,-1], dadosTreino[,1], prob=TRUE, k=10)
classCult.knn.tune<-tune.knn(dadosTreino[,-1], dadosTreino[,1], k=1:20);classCult.knn.tune

#Validacao
resultadoCult.knn<-table(pred=classCult.knn, true=dadosValidacao[,1]);resultadoCult.knn
resultadoCult.knn<-as.data.frame.matrix(resultadoCult.knn)
classDiag.knn<-diag(as.matrix(resultadoCult.knn));classDiag.knn
certos.knn<-sum(classDiag.knn)/sum(resultadoCult.knn);certos.knn

#Probabilidades
probs.knn<-unlist(attributes(classCult.knn),use.names=FALSE)
probs.knn<-as.numeric(probs.knn[-(1:(nClasses+1))])
probs.round.knn<-round(probs.knn,1)
hist(probs.knn,col='blue')
table(pred=classCult.knn, true=dadosValidacao[,1], probs=probs.round.knn)


###################
#APENAS SUB-CONJUNTO DE DADOS
classCultSub.knn<-knn(dadosTreinoSub[,-1], dadosValidacaoSub[,-1], dadosTreinoSub[,1], prob=TRUE, k=10)
classCultSub.tune.knn<-tune.knn(dadosTreinoSub[,-1], dadosTreinoSub[,1], k=1:20);classCult.knn.sub.tune

#Validacao
resultadoCultSub.knn<-table(pred=classCultSub.knn, true=dadosValidacaoSub[,1]);resultadoCultSub.knn
resultadoCultSub.knn<-as.data.frame.matrix(resultadoCultSub.knn)
classDiagSub.knn<-diag(as.matrix(resultadoCult.knn.sub));classDiagSub.knn
certosSub.knn<-sum(classDiagSub.knn)/sum(resultadoCultSub.knn);certosSub.knn

#Probabilidades
probsSub.knn<-unlist(attributes(classCultSub.knn),use.names=FALSE)
probsSub.knn<-as.numeric(probsSub.knn[-(1:(nClasses+1))])
probsSub.round.knn<-round(probsSub.knn,1)
hist(probsSub.knn,col='green')
table(pred=classCultSub.knn, true=dadosValidacaoSub[,1], probs=probsSub.round.knn)

#So com P >= 60%
jjj<-table(pred=classCultSub.knn, true=dadosValidacaoSub[,1], probs=probsSub.round.knn)
jjj<-as.matrix(jjj[,,5]+jjj[,,6]+jjj[,,7]+jjj[,,8]+jjj[,,9])
sum(diag(jjj))/sum(jjj)


###################
#ANALISE GRAFICA
#Percentagem de classificacoes correctas em cada cultura
percCertos<-c()
for(i in 1:ncol(resultadoCult.knn))
  percCertos[i]<-resultadoCult.knn[i,i]/sum(resultadoCult.knn[,i])
percCertos<-100*percCertos

percCertosSub<-c()
for(i in 1:ncol(resultadoCult.knn))
  percCertosSub[i]<-resultadoCultSub.knn[i,i]/sum(resultadoCultSub.knn[,i])
percCertosSub<-100*percCertosSub

par(las=2, mar=c(8,4.5,1,1))
barplot(rbind(percCertos,percCertosSub), beside=T, space=c(0,2.5), border=NA, names.arg=codigosNovos2005$codNomesAbrev[1:14], ylab='Percentagem de classificacoes correctas (%)', ylim=c(0,100), col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()
legend(x=30, y=98, legend=c("Completo","Sub-conjunto"), fill=c(rgb(0,0,0),rgb(0.7,0.7,0.7)))


#Percentagem de classificacoes correctas em cada valor de probabilidade
classProb<-table(pred=classCult.knn, true=dadosValidacao[,1], probs=probs.round.knn)
percCertos<-c()
for (i in 1:length(classProb[1,1,]))
{
  diag<-diag(as.matrix(classProb[,,i]))
  percCertos[i]<-sum(diag)/sum(classProb[,,i])
}
percCertos<-100*percCertos

classProbSub<-table(pred=classCultSub.knn, true=dadosValidacaoSub[,1], probs=probsSub.round.knn)
percCertoSubs<-c()
for (i in 1:length(classProbSub[1,1,]))
{
  diag<-diag(as.matrix(classProbSub[,,i]))
  percCertosSub[i]<-sum(diag)/sum(classProbSub[,,i])
}
percCertosSub<-100*percCertosSub

par(las=1, mar=c(4.5,4.5,1,1))
barplot(rbind(percCertos,percCertosSub), beside=T, space=c(0,2.5), border=NA, ylab='Percentagem de classificacoes correctas (%)', xlab='Probabilidade associada a classificao (%)', ylim=c(0,100), names.arg=seq(20,100,10),  col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()
legend(x=2, y=98, legend=c("Completo","Sub-conjunto"), fill=c(rgb(0,0,0),rgb(0.7,0.7,0.7)))

plot(seq(0.2,1,0.1),percCertos/100)

#Frequ?ncia da ocorrencia de cada valor de probabilidade
compProbs.knn<-rbind(table(probs.round.knn),table(probsSub.round.knn))
barplot(compProbs.knn, beside=T, space=c(0,4), legend=c("Completo","Sub-conjunto"), border=NA, col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()





##########################################################################
#2. SVM - MAQUINAS DE VECTORES DE SUPORTE

###################
#TODOS OS DADOS
classCult.svm<-svm(cultura ~ ., data=dadosTreino, probability=TRUE, kernel='radial', gamma=0.1, cost=3)
plot(classCult.svm, dadosTreino, B3_d3 ~ B4_d3)

#Tuning
classCult.tune.svm<-tune.svm(dadosTreino[,-1], dadosTreino[,1], gamma=seq(0, 1, by = .1));classCult.tune.svm
classCult.tune.svm<-tune.svm(dadosTreino[,-1], dadosTreino[,1], gamma=0.1, cost=seq(1, 10, by = 1));classCult.tune.svm
classCult.tune.svm<-tune.svm(dadosTreino[,-1], dadosTreino[,1], gamma=seq(0, 0.2, by=0.1), cost=seq(2, 4, by=1));classCult.tune.svm


#Validacao
validacao.svm<-predict(classCult.svm, dadosValidacao[,-1])
resultadoCult.svm<-table(pred = validacao.svm, true = dadosValidacao[,1]);resultadoCult.svm
resultadoCult.svm<-as.data.frame.matrix(resultadoCult.svm)
classDiag.svm<-diag(as.matrix(resultadoCult.svm));classDiag.svm
certos.svm<-sum(classDiag.svm)/sum(resultadoCult.svm);certos.svm

#Probabilidades
probs.svm<-predict(classCult.svm, dadosValidacao[,-1], probability=T)
probs.svm<-attr(probs.svm, "prob")
probs.svm<-as.vector(apply(probs.svm,1,max))
probs.round.svm<-round(probs.svm,1)
hist(probs.round.svm, breaks=14, col='blue')
table(pred=validacao.svm, true=dadosValidacao[,1], probs=probs.round.svm)


###################
#APENAS SUB-CONJUNTO DE DADOS
classCultSub.svm<-svm(cultura ~ ., data=dadosTreinoSub, probability=TRUE, kernel='radial', gamma=0.2, cost=2)
plot(classCultSub.svm, dadosTreinoSub, B3_d3 ~ B4_d3)

#Tuning
classCultSub.tune.svm<-tune.svm(dadosTreinoSub[,-1], dadosTreinoSub[,1], gamma=seq(0, 0.5, by = .1));classCultSub.tune.svm
classCultSub.tune.svm<-tune.svm(dadosTreinoSub[,-1], dadosTreinoSub[,1], gamma=0.2, cost=seq(1, 8, by = 1));classCultSub.tune.svm

#Validacao
validacaoSub.svm<-predict(classCultSub.svm, dadosValidacaoSub[,-1])
resultadoCultSub.svm<-table(pred = validacaoSub.svm, true = dadosValidacaoSub[,1]);resultadoCultSub.svm
resultadoCultSub.svm<-as.data.frame.matrix(resultadoCultSub.svm)
classDiagSub.svm<-diag(as.matrix(resultadoCultSub.svm));classDiagSub.svm
certosSub.svm<-sum(classDiagSub.svm)/sum(resultadoCultSub.svm);certosSub.svm

#Probabilidades
probsSub.svm<-predict(classCultSub.svm, dadosValidacaoSub[,-1], probability=T)
probsSub.svm<-attr(probsSub.svm, "prob")
probsSub.svm<-as.vector(apply(probsSub.svm,1,max))
probsSub.round.svm<-round(probsSub.svm,1)
hist(probsSub.round.svm, breaks=14, col='green')
table(pred=validacaoSub.svm, true=dadosValidacaoSub[,1], probs=probsSub.round.svm)

#So com P >= 60%
jjj<-table(pred=validacaoSub.svm, true=dadosValidacaoSub[,1], probs=probsSub.round.svm)
jjj<-as.matrix(jjj[,,5]+jjj[,,6]+jjj[,,7]+jjj[,,8]+jjj[,,9])
sum(diag(jjj))/sum(jjj)


###################
#Comparacao das probabilidades
compProbs.svm<-rbind(table(probs.round.svm),table(probsSub.round.svm))
barplot(compProbs.svm, beside=T, space=c(0,5), legend=c("Completo","Sub-conjunto"), xlab="Probabilidade", ylab="Numero de parcelas", border=NA, col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()
#Conclusao: com dados nao-completos, existem mais classificacaes com 100% de certeza







##########################################################################
#Abordagem diferente - Estimacao das densidades de probabilidade de cada classe de culturas

#QUESTOES:
# - Como vamos poder validar isto?
# - Como fazer uma verdadeira 'funcao' com isto para podermos fazer testes de hipoteses?

classe<-1
k<-10      #Numero de vizinhos a considerar - resultado do KNN.TUNE
d<-ncol(dadosTreinoSub[,-1])      #Numero de dimensoes
dadosTreinoClasse<-dadosTreinoSub[dadosTreinoSub$cultura == classe,]
dadosValidacaoClasse<-dadosValidacaoSub[dadosValidacaoSub$cultura == classe,]
n<-nrow(dadosTreinoClasse)     #Numero de observacoes na amostra de treino???

probs<-c()
for (i in 1:nrow(dadosValidacaoClasse))
{
  observacao<-dadosValidacaoClasse[i,-1]
  distancias<-sort(as.matrix(dist(rbind(observacao,dadosTreinoClasse[,-1])))[-1,1])
  distKVizinho<-distancias[k]
  
  #Calculo da probabilidade
  if(d%%2 == 0)
    cD<-pi^(d/2)/factorial(d/2)
  else
    cD<-pi^(d/2)/gamma(d/2+1)
  
  p<-k/(n*distKVizinho*cD)
  probs[i]<-p
}

plot(sort(probs))
hist(probs,breaks=20)


#########
#Experiencia com 2 dimensoes e pontos aleatorios
x<-sample(seq(1,20,0.1),10,replace=T)
y<-sample(seq(1,20,0.1),10,replace=T)
novo.x<-9
novo.y<-7
plot(x,y)
text(x,y+0.4,2:11)
points(novo.x,novo.y,col='red')
text(novo.x,novo.y+0.4,1,col='red')

distancias<-sort(as.matrix(dist(rbind(cbind(novo.x,novo.y),cbind(x,y))))[-1,1])
distKVizinho<-distancias[5]

draw.circle(novo.x,novo.y,distKVizinho,border='red')


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################







##########################################################################################
#   CRITERIO 1 --- UTILIZA????O DE DIFERENTES CLASSIFICADORES - TREINO E VALIDA????O         #
##########################################################################################

#Classificadores utilizados:
# 1. LDA - An??lise discriminante linear
# 2. KNN - k vizinhos mais pr??ximos
# 3. SVM - M??quinas de vectores de suporte
# 4. ??rvores de decis??o


##########################################################################
#0. PREPARA????O E VISUALIZA????O DOS DADOS

#Cria????o dos dados com as seguintes caracteristicas:
#  - pastagem ?? factor que ?? 0 se nao for pastagem, 1 se for pastagem permanente e 2 se for pastagem pobre
#  - usam-se dados de todas as bandas nas primeiras 6 datas do ano
#  - metade dos dados vai ser para treino dos classificadores e metade para valida????o

pastagem<-c()
seq<-c()
dadosClassificadores<-listaTodasParcelas[,c(5,7:12,15:20,23:28,31:36,39:44,47:52)]

for(i in 1:length(dadosClassificadores$cultura))
{
  if(dadosClassificadores$cultura[i] == 1)
    pastagem[i]<-1
  else
    pastagem[i]<-0
}
dadosClassificadoresPast<-NA
dadosClassificadoresPast<-cbind(pastagem,dadosClassificadores)
dadosClassificadoresPast<-dadosClassificadoresPast[,-2]
dadosClassificadoresPast$pastagem<-as.factor(as.numeric(as.vector(dadosClassificadoresPast$pastagem)))

#Amostra a usar para o modelo - metade para treino e metade para valida????o
set.seed(23)
amostra<-sample(1:nrow(dadosClassificadoresPast), ceiling(nrow(dadosClassificadoresPast)/2))

#Seleccao dos dados
dadosTreino<-dadosClassificadoresPast[amostra,]
dadosValidacao<-dadosClassificadoresPast[-amostra,]

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


##########################################################################
#1. LDA - ANALISE DISCRIMINANTE LINEAR 

#Treino do classificador
class.lda<-lda(dadosTreino[,-1],grouping=dadosTreino$pastagem);class.lda
plot(class.lda,col=(as.numeric(dadosTreino[,1])))

#Valida????o
validacao.lda<-predict(class.lda, dadosValidacao[,-1])
resultado.lda<-table(pred = validacao.lda$class, true = dadosValidacao[,1]);resultado.lda
classAgreement(resultado.lda)



#Optimiza????o com package subselect
class.lda.Hmat<-ldaHmat(dadosTreino[,-1], dadosTreino$pastagem)
class.lda.opt<-anneal(class.lda.Hmat$mat,kmin=3,kmax=5,H=class.lda.Hmat$H,crit="RM")
class.lda.opt$subsets

#nomes das 3 variaveis mais importantes
colnames(dadosTreino[,class.lda.opt$subsets])


##########################################################################
#2. SVM - M??QUINAS DE VECTORES DE SUPORTE

#Treino do classificador
class.svm<-svm(pastagem ~ ., data=dadosTreino, cost=100, gamma=0.9)
plot(class.svm, dadosTreino, B1_d3 ~ B5_d6)

#Optimiza????o - tentar com MAIS hipoteses
class.svm.tune<-tune.svm(pastagem ~ ., data=dadosTreino, gamma=seq(.1, 1, by = .1), cost = seq(100,1000, by = 100));class.svm.tune

#resultados: gamma=0.9 e cost=100

#Valida????o
validacao.svm<-predict(class.svm, dadosValidacao[,-1])
resultado.svm<-table(pred = validacao.svm, true = dadosValidacao[,1]);resultado.svm
classAgreement(resultado.svm)


##########################################################################
#3. KNN - K VIZINHOS MAIS PR??XIMOS

#Treino e valida????o simultanea
class.knn<-knn(dadosTreino[,-1], dadosValidacao[,-1], dadosTreino[,1], k=16)
resultado.knn<-table(pred=class.knn, true=dadosValidacao[,1]);resultado.knn
classAgreement(resultado.knn)

#Optimiza????o
class.knn.tune<-tune.knn(dadosTreino[,-1], dadosTreino[,1], k = 1:20);class.knn.tune


##########################################################################
#4. ??RVORES DE DECIS??O

#Ajuste inicial da ??rvore
class.tree<-tree(pastagem ~ ., data=dadosTreino, control= tree.control(nobs=nrow(dadosTreino), minsize=1))
plot(class.tree);text(class.tree,cex=.8)

#Valida????o
validacao.tree<-predict(class.tree,type="class")
resultado.tree<-table(pred = validacao.tree, true = dadosTreino[,1]);resultado.svm
classAgreement(resultado.tree)

#Determina????o da melhor ??rvore (com menor estimativa de erro de classifica????o)
cvDev<-cv.tree(class.tree)$dev
opt<-which(cvDev == min(cvDev))
class.tree.prune<-prune.tree(class.tree,best=opt)
plot(class.tree.prune);text(class.tree.prune,cex=.8)


#####
#OUTRA ABORDAGEM COM OUTRO PACKAGE QUE PARECE SER MELHOR
#Treino do classificador
class.tree2<-ctree(pastagem ~ ., data=dadosTreino)
plot(class.tree2)

#Valida????o
validacao.tree2<-predict(class.tree2, dadosValidacao[,-1])
resultado.tree2<-table(pred = validacao.tree2, true = dadosValidacao[,1]);resultado.tree2
classAgreement(resultado.tree2)

#####
#mais uma abordagem
class.tree3<-rpart(pastagem ~ ., data=dadosTreino)
plot(class.tree3)
text(class.tree3,cex=.8)

#Valida????o
validacao.tree3<-predict(class.tree3, dadosValidacao[,-1],type="class")
resultado.tree3<-table(pred = validacao.tree3, true = dadosValidacao[,1]);resultado.tree3
classAgreement(resultado.tree3)

##########################################################################
#5. REGRESS??O LOG??STICA MULTINOMIAL

#Com package NNET
library(nnet)
class.multi<-multinom(pastagem ~ ., data=dadosTreino)
summary(class.multi)

validacao.multi<-predict(class.multi, dadosValidacao[,-1])
resultado.multi<-table(pred = validacao.multi, true = dadosValidacao[,1]);resultado.multi


##########################################################################################
#   CRITERIO 2                                                                           #
##########################################################################################








##########################################################################################
#   EXPERIENCIAS                                                                         #
##########################################################################################



##########################################################################
#Experiencias com NDVI's (apenas para ano 2005)

#Fazer plot
datasGrafico<-datas
substr(datasGrafico,0,1)<-"0"
datasGrafico<-as.numeric(datasGrafico)

#plot(datasGrafico,NDVIs[1,],ylim=c(-1,1),col=100,type='l',lty=1,pch=16,ylab="NDVI", xlab="DOY 2005")
plot(datasGrafico,listaTodasParcelas[1,55:62],ylim=c(-1,1),col=1,type='l',lty=1,ylab="NDVI", xlab="DOY 2005")
for(i in 1:nrow(listaTodasParcelas))
  lines(datasGrafico,listaTodasParcelas[i,55:62],ylim=c(-1,1),col=i,lty=1,pch=16)


#distinguir entre sequeiro e regadio
todasParcelasSEQ<-listaTodasParcelas[listaTodasParcelas$seq == TRUE,]
todasParcelasREG<-listaTodasParcelas[listaTodasParcelas$seq == FALSE,]

plot(datasGrafico,todasParcelasSEQ[1,55:62],ylim=c(-1,1),col='burlywood2',type='l',lty=1,ylab="NDVI", xlab="DOY 2005")
for(i in 1:100)
  lines(datasGrafico,todasParcelasSEQ[i,55:62],ylim=c(-1,1),col='burlywood2',lty=1,pch=16)
for(i in 1:100)
  lines(datasGrafico,todasParcelasREG[i,55:62],ylim=c(-1,1),col='darkolivegreen3',lty=1,pch=16)

##########################################################################






##########################################################################
#Experi??ncias com outros indices de vegeta??ao

listaTodasParcelas$NDVI_1[50]
(listaTodasParcelas$B4_d1[50]-listaTodasParcelas$B3_d1[50])/(listaTodasParcelas$B4_d1[50]+listaTodasParcelas$B3_d1[50])

nir<-listaTodasParcelas$B4_d4
red<-listaTodasParcelas$B3_d4
NDVId4<-(nir-red)/(nir+red)
NDVId4.orig<-listaTodasParcelas$NDVI_4

plot(NDVId4,NDVId4.orig)
abline(lm(NDVId4.orig~NDVId4),col='red')


MSAVId4<-(2*nir+1-sqrt((2*nir+1)^2-8*(nir-red)))/2
plot(NDVId4,MSAVId4)
abline(lm(MSAVId4~NDVId4),col='red')


##########################################################################






##########################################################################
#Experimentar diferen??a entre usar apenas os interiores das parcelas e usar parcelas inteiras
testeTotInt<-matrix(nrow=0,ncol=5)
for(i in 1:length(listaDados))
{
  for(j in 1:length(listaDados[[i]][[1]]))
  {
    parc<-listaDados[[i]][[1]][[j]]
    conteudoTotal<-obtemConteudoParcela(parc,apenasInterior=FALSE)
    conteudoInterior<-obtemConteudoParcela(parc,apenasInterior=TRUE)
    nTotal<-length(conteudoTotal[[1]][,1])
    nInterior<-length(conteudoInterior[[1]][,1])
    sdTotal<-sd(conteudoTotal[[1]][,1])
    sdInterior<-sd(conteudoInterior[[1]][,1])
    area<-parc$atributos$AREA
    linha<-c(nTotal,sdTotal,nInterior,sdInterior,area)
    testeTotInt<-rbind(testeTotInt,linha)
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
##########################################################################




##########################################################################################
#   EXPORTAR GRAFICOS E QUADROS PARA LATEX                                               #
##########################################################################################


#Graficos
opar<-par()
par(opar)

#Graficos
library(tikzDevice)
tikz('../05 Escrito/TEX/test.tex',width=5,height=4)
par(mar=c(3.5,3.5,0.5,0.5), mgp=c(2,0.6,0))   #bottom, left, top, right
barplot(compProbs.svm, beside=T, space=c(0,5), legend=c("Completo","Sub-conjunto"), xlab="Probabilidade", ylab="Numero de parcelas", border=NA, col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()
dev.off()


#Quadros
library(xtable)
dados.xtab<-xtable(dadosTreino[1:10,1:5])
align(dados.xtab)<-rep("c",6)
digits(dados.xtab)<-3
print.xtable(dados.xtab, booktabs=T, include.rownames=F)



##########################################################################################
#   EXPORTAR E IMPORTAR GRANDES OBJECTOS                                                 #
##########################################################################################

#Gravar objectos grandes como ficheiros
save(cdg, file="cdg2005.Robj")
save(listaDados, file="listaDados2005.Robj")
save(todasImagens, file="todasImagens2005.Robj")
save(listaTodasParcelasIniciais, file="listaTodasParcelasIniciais2005.Robj")

rm(listaDados)
rm(cdg)
rm(todasImagens)

#Carregar objectos grandes
load("listaDados2005.Robj")
load("cdg2005.Robj")
load("todasImagens2005.Robj")


