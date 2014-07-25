#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Master's Thesis - Remote Sensing                      %
# Environmental Engineering - ISA/UL - Lisbon, Portugal %
# (c) 2014 by Jonas Schmedtmann                         %
#                                                       %
# TABLES AND FIGURES OUTPUT                             %
#                                                       %
# This script contains code to produce figures and      %
# tables for the final written document                 %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# OUTPUTS DE FIGURAS                                                                  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Lista de possiveis figuras
# - Mapa com area de estudo
# - Assinaturas espectrais medias por cultura (todas ou so algumas? meter desvios padrao?)
# - User's accuracies globais (faz sentido incluir erro sob forma de desvio padrao nos graficos? e faz sentido estes graficos serem de barras?)
# - Limites qi por cultura
# - Percentagem de parcelas em que e tomada uma decisao, por parcela (fixando confianca)
# - Percentagem de parcelas com decisao em funcao do nivel de confianca pretendido


#Lista de possiveis quadros
# - Lista de imagens Landsat usadas
# - Lista de culturas e numero de parcelas e area correspondente, pre e pos seleccao de dados para o estudo (pode se tambem atribuir um codigo a cada cultura e depois meter esse codigo nos graficos)
# - Resultados da correccao de imagens
# - Exclusao de variaveis (como na apresentacao - incluir aqui as datas usadas)
# - Comparacao entre classificadores (usando que medidas? so a percentagem de boas correccoes? total ou por cultura? e o kappa?)
# - Matriz de erro global de um dos classificadores (tinha que ser a soma das matrizes de erro de cada uma das 10 validacoes na validacao cruzada)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#GRAFICOS ACTUALIZADOS COM VALIDACAO CRUZADA ----

#User's accuracies globais
par(las=2, mar=c(9.5,3,0.5,0.1), mgp=c(2,0.5,0), tck = 0.012)
barplot(rbind(KNN.comp.cruz$result$userA, KNN.sub.cruz$result$userA, SVM.comp.cruz$result$userA, SVM.sub.cruz$result$userA), 
        beside=T, 
        space=c(0,2), 
        border=NA, 
        names.arg=codigosNovos2005$NOVO_NOME[1:12], 
        ylab='User\'s accuracy', 
        ylim=c(0,1), 
        col=c(rgb(0,0,1,1),rgb(0,0,1,0.5),rgb(1,0,0,1),rgb(1,0,0,0.5)));box()
legend(x=42.5, y=0.98,
       legend=c("KNN Completo","KNN Sub-conjunto","SVM Completo","SVM Sub-conjunto"),
       fill=c(rgb(0,0,1,1),rgb(0,0,1,0.5),rgb(1,0,0,1),rgb(1,0,0,0.5)))


#Limites qi por cultura
barplot(rbind(KNN.comp.cruz$result$qi, KNN.sub.cruz$result$qi, SVM.comp.cruz$result$qi, SVM.sub.cruz$result$qi), 
        beside=T, 
        space=c(0,2), 
        border=NA, 
        names.arg=codigosNovos2005$NOVO_NOME[1:12], 
        ylab='Limite qi', 
        ylim=c(0,1), 
        col=c(rgb(0,0,1,1),rgb(0,0,1,0.5),rgb(1,0,0,1),rgb(1,0,0,0.5)));box()
legend(x=42.5, y=0.98,
       legend=c("KNN Completo","KNN Sub-conjunto","SVM Completo","SVM Sub-conjunto"),
       fill=c(rgb(0,0,1,1),rgb(0,0,1,0.5),rgb(1,0,0,1),rgb(1,0,0,0.5)))


#Percentagem de parcelas em que e tomada uma decisao, por parcela
barplot(rbind(KNN.comp.cruz$result$percParcDec*100, KNN.sub.cruz$result$percParcDec*100, SVM.comp.cruz$result$percParcDec*100, SVM.sub.cruz$result$percParcDec*100), 
        beside=T, 
        space=c(0,2), 
        border=NA, 
        names.arg=codigosNovos2005$NOVO_NOME[1:12], 
        ylab='Percentagem de parcelas com decisao', 
        ylim=c(0,100), 
        col=c(rgb(0,0,1,1),rgb(0,0,1,0.5),rgb(1,0,0,1),rgb(1,0,0,0.5)));box()
legend(x=42.5, y=98,
       legend=c("KNN Completo","KNN Sub-conjunto","SVM Completo","SVM Sub-conjunto"),
       fill=c(rgb(0,0,1,1),rgb(0,0,1,0.5),rgb(1,0,0,1),rgb(1,0,0,0.5)))

#Percentagem de parcelas com decisao em funcao do nivel de confianca pretendido
opar <- par()
par <- opar
par(las=0, mar=c(3.2,3.2,0.5,0.5), mgp=c(2,0.5,0), tck = 0.012)
lambdas <- seq(0.5,1,by=0.05)
plot(lambdas*100,
     SVM.comp.cruz$result$percsDecs*100,
     xaxt='n',
     yaxt='n',
     xlab='Confianca',
     ylab='Percentagem de parcelas com decisao',
     ylim=c(0,100),
     col=g[1],
     bg=g[1],
     pch=21,
     type='o')
axis(1, at=lambdas*100, labels=lambdas*100)
axis(2, at=seq(0,100,10), labels=seq(0,100,10))
axis(3, at=lambdas*100)
axis(4, at=seq(0,100,10))

points(lambdas*100, SVM.sub.cruz$result$percsDecs*100, col=g[2], bg=g[2], pch=22, type='o')
points(lambdas*100, KNN.comp.cruz$result$percsDecs*100, col=g[3], bg=g[3], pch=23, type='o')
points(lambdas*100, KNN.sub.cruz$result$percsDecs*100, col=g[4], bg=g[4], pch=24, type='o')


g <- gray.colors(4, start = 0.0, end = 0.6, gamma = 1.5)

barplot(KNN.comp.cruz$result$percParcDec*100, space = 0.3, density = 20, angle=45, ylim=c(0,100))
box()

#line: type = "l"
#points: type = "p"
#both: type = "o"


plot(SVM.comp.cruz$result$userA, SVM.comp.cruz$result$qi)
abline(lm(SVM.comp.cruz$result$qi ~ SVM.comp.cruz$result$userA), col='red')

plot(SVM.comp.cruz$result$userA, SVM.comp.cruz$result$percParcDec)
abline(lm(SVM.comp.cruz$result$percParcDec ~ SVM.comp.cruz$result$userA), col='red')




#Visualizacao de 'assinatura espectral' total
dados <- data.table(dadosClassificadoresSub)
assinaturasSub <- dados[,list(B1_d2=mean(B1_d2),
                              B1_d5=mean(B1_d5),
                              B1_d6=mean(B1_d6),
                              B3_d1=mean(B3_d1),
                              B4_d1=mean(B4_d1),
                              B4_d3=mean(B4_d3),
                              B4_d4=mean(B4_d4),
                              B4_d5=mean(B4_d5),
                              B4_d6=mean(B4_d6),
                              B5_d4=mean(B5_d4),
                              B5_d6=mean(B5_d6),
                              B6_d3=mean(B6_d3),
                              B1_d2sd=sd(B1_d2),
                              B1_d5sd=sd(B1_d5),
                              B1_d6sd=sd(B1_d6),
                              B3_d1sd=sd(B3_d1),
                              B4_d1sd=sd(B4_d1),
                              B4_d3sd=sd(B4_d3),
                              B4_d4sd=sd(B4_d4),
                              B4_d5sd=sd(B4_d5),
                              B4_d6sd=sd(B4_d6),
                              B5_d4sd=sd(B5_d4),
                              B5_d6sd=sd(B5_d6),
                              B6_d3sd=sd(B6_d3)),by=cultura]
assinaturasSub <- as.data.frame(assinaturasSub)
assinaturasSub <- assinaturasSub[order(assinaturasSub$cultura, decreasing = F),]

d <- c('14 Fev', '2 Mar', '5 Mai', '6 Jun', '8 Jul', '25 Ago')

par(mfrow=c(6,2), las=2, mar=c(1.2,1.2,0,0), oma=c(3,1.5,0.3,0.3), mgp=c(2,0.5,0), tck = 0.04)
#Todos
for(i in 1:nrow(assinaturasSub))
{
  plot(1:12, assinaturasSub[i,2:13], type='l', pch=20, ylim=c(0,0.4), axes=F, xlab='', ylab='')
  box()
  
  if(i == 11 | i == 12)
    axis(1, at=1:12, labels=c(d[2], d[5], d[6], d[1], d[1], d[3], d[4], d[5], d[6], d[4], d[6], d[3]))
  else
    axis(1, at=1:12, labels=rep('',12))
  
  if(i %% 2 == 1)
    axis(2, at=seq(0, 0.4, 0.1), labels=seq(0, 0.4, 0.1))
  else
    axis(2, at=seq(0, 0.4, 0.1), labels=rep('',5))
  
  axis(3, at=1:12, labels=rep('',12))
  axis(4, at=seq(0, 0.4, 0.1), labels=rep('',5))
  
  abline(v=c(3.5,4.5, 9.5, 11.5), lty=3)
  arrows(x0 = 1:12, x1 = 1:12, y0 = as.numeric(assinaturasSub[i,2:13]), y1 = as.numeric(assinaturasSub[i,2:13]+assinaturasSub[i,14:25]), length = 0.04, angle = 90)
  arrows(x0 = 1:12, x1 = 1:12, y0 = as.numeric(assinaturasSub[i,2:13]), y1 = as.numeric(assinaturasSub[i,2:13]-assinaturasSub[i,14:25]), length = 0.04, angle = 90)
  
  text(0.6, 0.34, codigosNovos2005$NOVO_NOME[i], pos = 4)
  
  if(i == 11 | i == 12)
  {
    text(1.45, 0.03, 'B1', pos = 4)
    text(3.45, 0.03, 'B3', pos = 4)
    text(6.45, 0.03, 'B4', pos = 4)
    text(9.95, 0.03, 'B5', pos = 4)
    text(11.45, 0.03, 'B7', pos = 4)
  }
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Graficos antes de validacao cruzada
plot(SVM.sub$result$class1, SVM.sub$result$prob1,ylim=c(0,1))
points(qi,col='red', bg='red', pch=21)  #Limites qi
points(as.vector(percParcDec/100), col='blue', bg='blue', pch=21)   #Percentagem em que e tomada uma decisao 

hist(SVM.sub$result$prob1,prob=T)
lines(density(SVM.sub$result$prob1))

plot(SVM.sub$userA,qi.SVM.sub$qi)
plot(qi.SVM.sub$percParcDec)


#Plot das confiancas em funcao de q
plot(seq(0, 1, by = 0.01),SVM.sub.cruz$resultTodos$valid1$todasConfiancas[,4], type='l',ylim=c(0,1))
for(i in 2:12)
  lines(seq(0, 1, by = 0.01),SVM.sub.cruz$resultTodos$valid1$todasConfiancas[,i], type='l',col=i)
abline(h=lambda,col='red')







#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Mapa com area de estudo ----

outX <- c( 490733,  573315,  573315,  490733,  490733)
outY <- c(4368546, 4368546, 4286706, 4286706, 4368546)

parcs <- cdg$parc2005
parcsTrans <- spTransform(parcs, PROJ4.UTM.LL)
pt <- readShapePoly(paste0(CAMINHO.SHP,"/Enquadramentos PT IGEO/Distritos"), proj4string=PROJ4.ETRS)
ptTrans <- spTransform(pt, PROJ4.UTM.LL)
areaEstudo <- SpatialPolygons(list(Polygons(list(Polygon(cbind(outX, outY))), "1")), proj4string=PROJ4.UTM)
areaEstudoTrans <- spTransform(areaEstudo, PROJ4.UTM.LL)
img <- corrigeLeConjuntoImagensLandsat(conjuntoPath = paste0(CAMINHO.LANDSAT,"/LE72040332005157EDC00"), prefixo = 'MAPA', corrige = F)

b2 <- projectRaster(img$b2,  crs = PROJ4.UTM.LL)
b3 <- projectRaster(img$b3,  crs = PROJ4.UTM.LL)
b4 <- projectRaster(img$b4,  crs = PROJ4.UTM.LL)
s <- stack(b2, b3, b4)

par <- opar

par(las=1, oma=c(1.7,1.7,0.5,0.5), mgp=c(2,0.5,0), tck = 0.012, las=0)
plotRGB(s, 3, 2, 1, stretch='lin', colNA = 'transparent')
axis(1, at=c(-10, -9.5, -9, -8.5, -8, -7.5), labels = c(-10, -9.5, -9, -8.5, -8, -7.5))
axis(2, at=c(38, 38.5, 39, 39.5), labels=c(38, 38.5, 39, 39.5))
axis(3)
axis(4)
box()
plot(ptTrans, add=T, border='black')
plot(parcsTrans, add=T, col='green', border='green')
plot(areaEstudoTrans, add=T, border='red', lwd=4)


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



#zoom
#drawExtent

