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


rowMeans(cbind(x, y))

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

#Amostra a usar para o modelo - metade para treino e metade para validacao
set.seed(23)
amostra <- sample(1:nrow(dadosClassificadores), ceiling(nrow(dadosClassificadores)/2))

#Seleccao dos dados
dadosTreino <- dadosClassificadores[amostra,]
dadosValidacao <- dadosClassificadores[-amostra,]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4.2 Sub-conjunto do conjunto completo ----

#Usando a funcao trim.matrix do package subselect
dadosTrim <- dadosClassificadores
criterioExclusaoVars <- 0.02
classTrim <- trim.matrix(cor(dadosTrim[,-1]),criterioExclusaoVars);classTrim

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
dadosValidacaoSubAreas <- listaTodasParcelas$area[-amostra]





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. 6. 7. ESCOLHA/TREINO/VALIDACAO DOS CLASSIFICADORES                               ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5.1 KNN - K vizinhos mais proximos ----

#Todos os dados
KNN.comp <- treinaValida(tipo = 'KNN', treino = dadosTreino, valid = dadosValidacao, k = 10)
qi.KNN.comp <- fixaLimiteQ(classificador = KNN.comp, lambda = 0.8)


#Apenas sub-conjunto de dados
KNN.sub  <- treinaValida(tipo = 'KNN', treino = dadosTreinoSub, valid = dadosValidacaoSub, k = 10)

#Tuning: kRange=c(1,10)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5.2 SVM - Maquinas de vectores de suporte ----

#Todos os dados
SVM.comp <- treinaValida(tipo = 'SVM', treino = dadosTreino, valid = dadosValidacao, gamma = 0.1, cost = 3)
qi.SVM.comp <- fixaLimiteQ(classificador = SVM.comp, lambda = 0.8)

#Tuning: gamma=seq(0, 1, by = .1)) e gamma=0.1, cost=seq(1, 10, by = 1))


#Apenas sub-conjunto de dados
SVM.sub  <- treinaValida(tipo = 'SVM', treino = dadosTreinoSub, valid = dadosValidacaoSub, gamma = 0.2, cost = 2)
qi.SVM.sub <- fixaLimiteQ(classificador = SVM.sub, lambda = 0.8)

#Tuning: gamma=seq(0, 0.5, by = .1) e gamma=0.2, cost=seq(1, 8, by = 1)





SVM.comp.cruz <- validacaoCruzada(n = 10, tipo = 'SVM', dados = dadosClassificadores, lambda = 0.8, gamma = 0.2, cost = 2)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Graficos interessantes
plot(SVM.sub$result$class1, SVM.sub$result$prob1,ylim=c(0,1))
points(qi,col='red', bg='red', pch=21)  #Limites qi
points(as.vector(percParcDec/100), col='blue', bg='blue', pch=21)   #Percentagem em que é tomada uma decisao 

hist(SVM.sub$result$prob1,prob=T)
lines(density(SVM.sub$result$prob1))

plot(SVM.sub$userA,qi.SVM.sub$qi)
plot(qi.SVM.sub$percParcDec)


#Plot das confianças em funçao de q
plot(seq(0, 1, by = 0.01),todasConfiancas[,1], type='l',ylim=c(0,1))
for(i in 2:12)
  lines(seq(0, 1, by = 0.01),todasConfiancas[,i], type='l',col=i)
abline(h=lambda,col='red')


#Percentagem de parcelas em que se toma uma decisao em funcao de lambda
lambdas <- seq(0.5,1,by=0.05)
percsDecs.KNN.comp <- percsDecs.KNN.sub <- percsDecs.SVM.comp <- percsDecs.SVM.sub <- c()
for(k in 1:length(lambdas))
{
  percsDecs.KNN.comp[k] <- fixaLimiteQ(classificador = KNN.comp, lambda = lambdas[k])$percTotDec
  percsDecs.KNN.sub[k] <- fixaLimiteQ(classificador = KNN.sub, lambda = lambdas[k])$percTotDec
  percsDecs.SVM.comp[k] <- fixaLimiteQ(classificador = SVM.comp, lambda = lambdas[k])$percTotDec
  percsDecs.SVM.sub[k] <- fixaLimiteQ(classificador = SVM.sub, lambda = lambdas[k])$percTotDec  
}

#KNN
plot(lambdas*100, percsDecs.KNN.comp*100, xaxt='n', yaxt='n', xlab='Confianca', ylab='Percentagem de parcelas com decisao', ylim=c(0,100), col='green', bg='green', pch=21)
lines(lambdas*100, percsDecs.KNN.comp*100, col='green')
#text(lambdas*100, percsDecs*100-6, round(percsDecs*100,1), col='blue', cex=0.75)
axis(1, at=lambdas*100, labels=lambdas*100)
axis(2, at=seq(0,100,10), labels=seq(0,100,10))

points(lambdas*100, percsDecs.KNN.sub*100, col='yellow', bg='yellow', pch=21)
lines(lambdas*100, percsDecs.KNN.sub*100, col='yellow')

#SVM
points(lambdas*100, percsDecs.SVM.comp*100, col='blue', bg='blue', pch=21)
lines(lambdas*100, percsDecs.SVM.comp*100, col='blue')

points(lambdas*100, percsDecs.SVM.sub*100, col='red', bg='red', pch=21)
lines(lambdas*100, percsDecs.SVM.sub*100, col='red')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5.3 Comparacao de KNN e SVM para conjunto total e sub-conjunto ----

#User's accuracies (ALTERAR NO KNN!!!!)
par(las=2, mar=c(10,4.5,1,1))
barplot(rbind(KNN.comp$userA,KNN.sub$userA,SVM.comp$userA,SVM.sub$userA), 
        beside=T, 
        space=c(0,2), 
        border=NA, 
        names.arg=codigosNovos2005$NOVO_NOME[1:12], 
        ylab='Classificacoes correctas', 
        ylim=c(0,1), 
        col=c(rgb(0,0,1,1),rgb(0,0,1,0.5),rgb(1,0,0,1),rgb(1,0,0,0.5)));box()
legend(x=42.5, y=0.98,
       legend=c("KNN Completo","KNN Sub-conjunto","SVM Completo","SVM Sub-conjunto"),
       fill=c(rgb(0,0,1,1),rgb(0,0,1,0.5),rgb(1,0,0,1),rgb(1,0,0,0.5)))


#Percentagem de classificacoes correctas em cada valor de probabilidade
par(las=1, mar=c(4.5,4.5,1,1))
barplot(rbind(KNN.comp$correcProb,KNN.sub$correcProb,SVM.comp$correcProb,SVM.sub$correcProb),
        beside=T,
        space=c(0,2.5),
        border=NA,
        ylab='Classificacoes correctas',
        xlab='Probabilidade associada a classificao',
        ylim=c(0,1),
        col=c(rgb(0,0,1,1),rgb(0,0,1,0.5),rgb(1,0,0,1),rgb(1,0,0,0.5)));box()
legend(x=2, y=0.98,
       legend=c("KNN Completo","KNN Sub-conjunto","SVM Completo","SVM Sub-conjunto"),
       fill=c(rgb(0,0,1,1),rgb(0,0,1,0.5),rgb(1,0,0,1),rgb(1,0,0,0.5)))


#Frequencia da ocorrencia de cada valor de probabilidade
barplot(rbind(table(KNN.comp$result$prob1R),table(KNN.sub$result$prob1R),table(SVM.comp$result$prob1R),table(SVM.sub$result$prob1R)), 
        beside=T, 
        ylim=c(0,1400), 
        space=c(0,2.5),
        border=NA,
        ylab='Frequencia',
        xlab='Probabilidade associada a classificao',
        col=c(rgb(0,0,1,1),rgb(0,0,1,0.5),rgb(1,0,0,1),rgb(1,0,0,0.5)));box()
legend(x=1.3, y=1380,
       legend=c("KNN Completo","KNN Sub-conjunto","SVM Completo","SVM Sub-conjunto"),
       fill=c(rgb(0,0,1,1),rgb(0,0,1,0.5),rgb(1,0,0,1),rgb(1,0,0,0.5)))




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Cenas para apresentacao IFAP ----

jjj <- cbind(SVM.sub$result,area=dadosValidacaoSubAreas,SVM.sub$probs)
jjj <- jjj[order(jjj$area, decreasing = TRUE),]

#Exemplo com ALTA probabilidade
jjj[4,]
parcNome <- as.character(jjj$parcela[4])
b2 <- listaDados[[parcNIFAP(parcNome)]][['a2005']][[parcNome]][['DR']][['d189']][,,2]
b3 <- listaDados[[parcNIFAP(parcNome)]][['a2005']][[parcNome]][['DR']][['d189']][,,3]
b4 <- listaDados[[parcNIFAP(parcNome)]][['a2005']][[parcNome]][['DR']][['d189']][,,4]
bX <- listaDados[[parcNIFAP(parcNome)]][['a2005']][[parcNome]][['DR']][['d189']][,,7]
bY <- listaDados[[parcNIFAP(parcNome)]][['a2005']][[parcNome]][['DR']][['d189']][,,8]
poly <- listaDados[[parcNIFAP(parcNome)]][['a2005']][[parcNome]][['coordsPoly']]
b2R <- raster(b2, xmn=min(bX), xmx=max(bX), ymn=min(bY), ymx=max(bY), crs=PROJ4.UTM)
b3R <- raster(b3, xmn=min(bX), xmx=max(bX), ymn=min(bY), ymx=max(bY), crs=PROJ4.UTM)
b4R <- raster(b4, xmn=min(bX), xmx=max(bX), ymn=min(bY), ymx=max(bY), crs=PROJ4.UTM)

plot(s$layer.1)
lines(poly,type='l',lwd=4,col='yellow')

s <- stack(b2R,b3R,b4R)
plotRGB(s, 3, 2, 1, stretch='lin',colNA='black')

out <- cbind(jjj[4,c(3,4)],0.6384283,jjj[4,9])
colnames(out) <- c('class', 'prob', 'probCult', 'area')


#%%%%%%%%%%%%
#Exemplo com BAIXA probabilidade
jjj[204,]   #204
parcNome <- as.character(jjj$parcela[204])
b2 <- listaDados[[parcNIFAP(parcNome)]][['a2005']][[parcNome]][['DR']][['d189']][,,2]
b3 <- listaDados[[parcNIFAP(parcNome)]][['a2005']][[parcNome]][['DR']][['d189']][,,3]
b4 <- listaDados[[parcNIFAP(parcNome)]][['a2005']][[parcNome]][['DR']][['d189']][,,4]
bX <- listaDados[[parcNIFAP(parcNome)]][['a2005']][[parcNome]][['DR']][['d189']][,,7]
bY <- listaDados[[parcNIFAP(parcNome)]][['a2005']][[parcNome]][['DR']][['d189']][,,8]
poly <- listaDados[[parcNIFAP(parcNome)]][['a2005']][[parcNome]][['coordsPoly']]
b2R <- raster(b2, xmn=min(bX), xmx=max(bX), ymn=min(bY), ymx=max(bY), crs=PROJ4.UTM)
b3R <- raster(b3, xmn=min(bX), xmx=max(bX), ymn=min(bY), ymx=max(bY), crs=PROJ4.UTM)
b4R <- raster(b4, xmn=min(bX), xmx=max(bX), ymn=min(bY), ymx=max(bY), crs=PROJ4.UTM)

plot(s$layer.1)
lines(poly,type='l',lwd=4,col='yellow')

s <- stack(b2R,b3R,b4R)
plotRGB(s, 3, 2, 1, stretch='lin',colNA='black')

