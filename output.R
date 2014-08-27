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
# 1) Mapa com area de estudo
# 2) Assinaturas espectrais 'modificadas' medias por cultura
# 3) User's accuracies globais
# 4) Limites qi por cultura
# 5) Percentagem de parcelas em que e tomada uma decisao, por parcela e por nivel de confianca pretendido
# 6) Percentagem de parcelas em que e tomada uma decisao, GLOBAL, por nivel de confianca pretendido


g <- gray.colors(4, start = 0.0, end = 0.6, gamma = 1.5)
opar <- par()
par <- opar

#Graficos de linhas:
#line: type = "l"
#points: type = "p"
#both: type = "o"


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1) Mapa com area de estudo ---- PRONTO!!!


#Seleccao das parcelas que foram usadas
contem <- paste0('p',cdg$parc2005$COD_CRUZA) %in% listaTodasParcelas$id
table(contem)
parcsUsadas <- cdg$parc2005[contem,]
parcsUsadasTrans <- spTransform(parcsUsadas, PROJ4.UTM.LL)
plot(parcsUsadasTrans)
#writeSpatialShape(x = parcsUsadas, fn = 'PARC_USADAS_2005')

#Area de estudo
outX <- c( 490733,  573315,  573315,  490733,  490733)
outY <- c(4368546, 4368546, 4286706, 4286706, 4368546)
areaEstudo <- SpatialPolygons(list(Polygons(list(Polygon(cbind(outX, outY))), "1")), proj4string=PROJ4.UTM)
areaEstudoTrans <- spTransform(areaEstudo, PROJ4.UTM.LL)
plot(areaEstudoTrans, add=T)
writeSpatialShape(x = areaEstudo, fn = 'AREA_GLOBAL')

areaEstudoTrans
convertCoord(8.8)

#Enquadramento
pt <- readShapePoly(paste0(CAMINHO.SHP,"/Enquadramentos PT IGEO/Distritos"), proj4string=PROJ4.ETRS)
ptTrans <- spTransform(pt, PROJ4.UTM.LL)
plot(ptTrans, add=T)

pais <- readShapePoly(paste0(CAMINHO.SHP,"/Enquadramentos PT IGEO/portugal"), proj4string=PROJ4.ETRS)
paisTrans <- spTransform(pais, PROJ4.UTM.LL)

scene <- readShapePoly(paste0(CAMINHO.SHP,"/204-33-lisboa"), proj4string=PROJ4.UTM)
sceneTrans <- spTransform(scene, PROJ4.UTM.LL)

#Imagem
img <- corrigeLeConjuntoImagensLandsat(conjuntoPath = paste0(CAMINHO.LANDSAT,"/LE72040332005189EDC00"), prefixo = 'MAPA_06.08', corrige = T, areaEstudo = areaEstudo)

b1 <- projectRaster(img$b1,  crs = PROJ4.UTM.LL) #Blue
b2 <- projectRaster(img$b2,  crs = PROJ4.UTM.LL) #Green
b3 <- projectRaster(img$b3,  crs = PROJ4.UTM.LL) #Red
b4 <- projectRaster(img$b4,  crs = PROJ4.UTM.LL) #NIR
s <- stack(b2, b3, b4)
sTrue <- stack(b1, b2, b3)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Mapa de toda a area
par(las=1, oma=c(2,2,1,1), mgp=c(5,0.75,0), tck = -0.012, las=0, cex=2)
#plotRGB(sTrue, 3, 2, 1, stretch='lin', colNA = 'black')
#extDraw <- drawExtent(show=TRUE, col="red")
jExt <- extent(-9.1, -8.16, 38.73, 39.46)
plotRGB(sTrue, 3, 2, 1, stretch='lin', colNA = 'black', ext = jExt)


#Axes
#plot(seq(-9, -8, 0.2), c(38.8, 39, 39.1, 39.2, 39.3, 39.4), xlim = range(at), axes = F);box()
axis(1, at = c(-9, -8.8, -8.6, -8.4, -8.2), labels = expression(paste(9*degree, " 0", 0*minute, " W"),
                                                                    paste(8*degree, " ", 48*minute, " W"),
                                                                    paste(8*degree, " ", 36*minute, " W"),
                                                                    paste(8*degree, " ", 24*minute, " W"),
                                                                    paste(8*degree, " ", 12*minute, " W")))
axis(2, at = c(38.8, 39, 39.2, 39.4), labels = expression(paste(38*degree, " ", 48*minute, " N"),
                                                          paste(39*degree, " 0", 0*minute, " N"),
                                                          paste(39*degree, " ", 12*minute, " N"),
                                                          paste(39*degree, " ", 24*minute, " N")))
box()
plot(parcsUsadasTrans, add=T, col='yellow', border='transparent', lwd=3)
plot(smallEstudo, border='red', add=T, lwd=5)
#Guardar com 2000x2002 px

convertCoord(38.73)
convertCoord(39.46)
convertCoord(9.1)
convertCoord(8.16)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Mapa de uma pequena zona para visuzlizacao de culturas em falsa cor
#smallerDraw <- drawExtent(show=TRUE, col="red")
smallerDraw <- smallerDraw

smallX <- c( -9.022535,  -8.854367,  -8.854367,  -9.022535,  -9.022535)
smallY <- c(38.99166, 38.99166, 38.8821, 38.8821, 38.99166)
smallEstudo <- SpatialPolygons(list(Polygons(list(Polygon(cbind(smallX, smallY))), "1")), proj4string=PROJ4.UTM.LL)
within <- gWithin(parcsUsadasTrans, smallEstudo, byid=TRUE)
parcsClip <- subset(parcsUsadasTrans,as.vector(within))

trigo <- parcsClip[parcsClip$COD_CULTUR == 1 & parcsClip$AREA > 5,]
milho <- parcsClip[parcsClip$COD_CULTUR == 6 & parcsClip$AREA > 5,]
arroz <- parcsClip[parcsClip$COD_CULTUR == 24 & parcsClip$AREA > 5,]
forrageiras <- parcsClip[parcsClip$COD_CULTUR == 142 & parcsClip$AREA > 5,]
pastagens <- parcsClip[parcsClip$COD_CULTUR == 143 & parcsClip$AREA > 5,]

#Plot it
par(las=1, oma=c(0,0,0,0), mgp=c(5,0.75,0), tck = -0.012, las=0, cex=1)
plotRGB(s, 3, 2, 1, stretch='lin', colNA = 'black', ext = smallerDraw)
plot(pastagens, add=T, border=cores2[4], lwd=5)   #Azul claro
plot(forrageiras, add=T, border=cores1[3], lwd=5)   #Cor de laranja estranho
plot(milho, add=T, border=cores1[2], lwd=5)   #Amarelo
plot(arroz, add=T, border=cores2[2], lwd=5)   #Azul
plot(trigo, add=T, border=cores1[1], lwd=5)   #Verde
legend(x=-9.02, y=38.99, border='white', box.col='white',
       legend=c("PGL","FOR","MAI","RIC","WHE"),
       fill=c(cores2[4], cores1[3], cores1[2], cores2[2], cores1[1]))


cores1 <- terrain.colors(5)
cores2 <- topo.colors(8)
pie(rep(1, 5), col = terrain.colors(5))
pie(rep(1, 8), col = topo.colors(8))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Mapa de enquadramento (exportar com 1000px)
plot(ptTrans, col="#DDDDDD", border='white')
plot(paisTrans, border="black", add=T)
plot(areaEstudoTrans, border='red', add=T, lwd=2)
plot(scene, border="black", add=T, lwd=2)

text(-8.975, 39.05, "1", cex=1.3)
text(-8.6, 39.05, "2", cex=1.3)
text(-10.275, 38.05, "Landsat 7 scene", pos=4, cex=1.3)
text(-10.275, 37.9, "Path 204, row 33", pos=4, cex=1.3)
text(-8.625, 39.35, "Study area", cex=1.3)


#Calculos de areas
newX <- c( -9.1,  -8.16,  -8.16,  -9.1,  -9.1)
newY <- c(39.46, 39.46, 38.73, 38.73, 39.46)
newEstudo <- SpatialPolygons(list(Polygons(list(Polygon(cbind(newX, newY))), "1")), proj4string=PROJ4.UTM.LL)
newEstudoTrans <- spTransform(newEstudo, PROJ4.UTM)
#Resultado: 6390km^2

studyAreaKm2 <-6390
areaParcKm2 <- areaTot/100
areaParcKm2/studyAreaKm2


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2) Assinaturas espectrais 'modificadas' medias por cultura ---- PRONTO!!!
dados <- data.table(dadosClassificadoresSub)
assinaturasSub <- dados[,list(B1_d4=mean(B1_d4),
                              B3_d1=mean(B3_d1),
                              B3_d3=mean(B3_d3),
                              B3_d6=mean(B3_d6),
                              B4_d1=mean(B4_d1),
                              B4_d2=mean(B4_d2),
                              B4_d3=mean(B4_d3),
                              B4_d4=mean(B4_d4),
                              B4_d5=mean(B4_d5),
                              B4_d6=mean(B4_d6),
                              B5_d2=mean(B5_d2),
                              B6_d5=mean(B6_d5),
                              B1_d4=sd(B1_d4),
                              B3_d1=sd(B3_d1),
                              B3_d3=sd(B3_d3),
                              B3_d6=sd(B3_d6),
                              B4_d1=sd(B4_d1),
                              B4_d2=sd(B4_d2),
                              B4_d3=sd(B4_d3),
                              B4_d4=sd(B4_d4),
                              B4_d5=sd(B4_d5),
                              B4_d6=sd(B4_d6),
                              B5_d2=sd(B5_d2),
                              B6_d5=sd(B6_d5)),by=cultura]
assinaturasSub <- as.data.frame(assinaturasSub)
assinaturasSub <- assinaturasSub[order(assinaturasSub$cultura, decreasing = F),]

d <- c('10 Nov', '14 Feb', '5 May', '6 Jun', '8 Jul', '25 Aug')

tikz('../05 Escrito/TEX/signatures.tex',width=6,height=8)
par(mfrow=c(6,2), las=2, mar=c(1.2,1.2,0,0), oma=c(3,1.5,0.3,0.3), mgp=c(3,0.8,0), tck = 0.03)
for(i in 1:nrow(assinaturasSub))
{
  #Grafico em si
  plot(1:12, assinaturasSub[i,2:13], type='l', pch=20, ylim=c(-0.04,0.44), axes=F, xlab='', ylab='')
  box()
  
  #Eixo das datas na parte de baixo
  if(i == 11 | i == 12)
    axis(1, at=1:12, labels=c(d[4], d[1], d[3], d[6], d[1], d[2], d[3], d[4], d[5], d[6], d[2], d[5]))
  else
    axis(1, at=1:12, labels=rep('',12))
  
  #Eixo das reflectancias no lado esquerdo
  if(i %% 2 == 1)
    axis(2, at=seq(0, 0.4, 0.1), labels=seq(0, 0.4, 0.1))
  else
    axis(2, at=seq(0, 0.4, 0.1), labels=rep('',5))
  
  #Restantes eixos
  axis(3, at=1:12, labels=rep('',12))
  axis(4, at=seq(0, 0.4, 0.1), labels=rep('',5))
  
  #Tracejado e barras de erros
  abline(v=c(1.5,4.5, 10.5, 11.5), lty=3)
  arrows(x0 = 1:12, x1 = 1:12, y0 = as.numeric(assinaturasSub[i,2:13]), y1 = as.numeric(assinaturasSub[i,2:13]+assinaturasSub[i,14:25]), length = 0.03, angle = 90)
  arrows(x0 = 1:12, x1 = 1:12, y0 = as.numeric(assinaturasSub[i,2:13]), y1 = as.numeric(assinaturasSub[i,2:13]-assinaturasSub[i,14:25]), length = 0.03, angle = 90)
  
  #Nomes das culturas
  text(0.6, 0.39, NOMES_CULT[i], pos = 4)
  
  #Numeros das bandas
  if(i == 11 | i == 12)
  {
    text(0.5, -0.01, 'B1', pos = 4)
    text(2.5, -0.01, 'B3', pos = 4)
    text(7, -0.01, 'B4', pos = 4)
    text(10.5, -0.01, 'B5', pos = 4)
    text(11.5, -0.01, 'B7', pos = 4)
  }
}
dev.off()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3) User's accuracies globais ----
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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4) Limites qi por cultura ----
par(las=2, mar=c(9.5,3,0.5,0.1), mgp=c(2,0.5,0), tck = 0.012)
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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5) Percentagem de parcelas em que e tomada uma decisao, por cultura e por nivel de confianca pretendido ----
g5 <- gray.colors(5, start = 0, end = 1, gamma = 1.5)

tikz('../05 Escrito/TEX/confidenceDecesion.tex', width=6, height=4, sanitize = T)
par(las=2, mar=c(2.5,3,0.5,0.1), mgp=c(2,0.5,0), tck = 0.01, cex=1)
barplot(t(SVM.sub.cruz$result$percsDecs)[c(3,5,7,9,11),2:13]*100, 
        beside=T, 
        space=c(0,2), 
        border='black', 
        names.arg=COD_CULT, 
        ylab='Percentage of parcels classified automatically', 
        ylim=c(0,100), #113
        axes=F,
        col=g5);box()
axis(2, at=seq(0, 100, 10), labels=seq(0, 100, 10))
axis(4, at=seq(0, 100, 10), labels=seq(0, 100, 10))
legend(x=75, y=101, ncol=1, bty='n', border='black',
       legend=c("60%","70%","80%","90%", "100%"),
       fill=g5)     #x=70
dev.off()


#Para abstract
tikz('../05 Escrito/TEX/confidenceDecesionABSTRACT.tex', width=6, height=4.5, sanitize = T)
par(las=2, mar=c(6.5,3,3.0,0.1), mgp=c(2,0.5,0), tck = 0.01, cex=1)
barplot(t(SVM.sub.cruz$result$percsDecs)[c(3,5,7,9,11),2:13]*100, 
        beside=T, 
        space=c(0,2), 
        border='black', 
        names.arg=ABRV_CULT, 
        ylab='Percentage of parcels classified automatically', 
        ylim=c(0,100), #113
        axes=F,
        col=g5);box()
title('Control with Remote Sensing: proportion of parcels\nautomatically classified at different confidence levels', line = 0.7)
axis(2, at=seq(0, 100, 10), labels=seq(0, 100, 10))
axis(4, at=seq(0, 100, 10), labels=seq(0, 100, 10))
legend(x=75, y=101, ncol=1, bty='n', border='black',
       legend=c("60%","70%","80%","90%", "100%"),
       fill=g5)     #x=70
dev.off()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 6) Percentagem de parcelas em que e tomada uma decisao, GLOBAL, por nivel de confianca pretendido ----
par(las=0, mar=c(3.2,3.2,0.5,0.5), mgp=c(2,0.5,0), tck = 0.012)
lambdas <- seq(0.5,1,by=0.05)
plot(lambdas*100,
     SVM.comp.cruz$result$percsDecs[1,]*100,
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

points(lambdas*100, SVM.sub.cruz$result$percsDecs[1,]*100, col=g[2], bg=g[2], pch=22, type='o')
points(lambdas*100, KNN.comp.cruz$result$percsDecs[1,]*100, col=g[3], bg=g[3], pch=23, type='o')
points(lambdas*100, KNN.sub.cruz$result$percsDecs[1,]*100, col=g[4], bg=g[4], pch=24, type='o')




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Outros (por exemplo antes de validacao cruzada)
plot(SVM.sub$result$class1, SVM.sub$result$prob1,ylim=c(0,1))
points(qi,col='red', bg='red', pch=21)  #Limites qi
points(as.vector(percParcDec/100), col='blue', bg='blue', pch=21)   #Percentagem em que e tomada uma decisao 

hist(SVM.sub$result$prob1,prob=T)
lines(density(SVM.sub$result$prob1))

plot(SVM.sub$userA,qi.SVM.sub$qi)
plot(qi.SVM.sub$percParcDec)


#Plot das confiancas em funcao de q
plot(seq(0, 1, by = 0.001),SVM.sub.cruz$resultTodos$valid7$todasConfiancas[,4], type='l',ylim=c(0,1))
for(i in 2:12)
  lines(seq(0, 1, by = 0.001),SVM.sub.cruz$resultTodos$valid7$todasConfiancas[,i], type='l',col=i)
abline(h=lambda,col='red')


plot(SVM.comp.cruz$result$userA, SVM.comp.cruz$result$qi)
abline(lm(SVM.comp.cruz$result$qi ~ SVM.comp.cruz$result$userA), col='red')

plot(SVM.comp.cruz$result$userA, SVM.comp.cruz$result$percParcDec)
abline(lm(SVM.comp.cruz$result$percParcDec ~ SVM.comp.cruz$result$userA), col='red')






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#Cenas para apresentacao IFAP

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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# OUTPUTS DE QUADROS                                                                  ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Lista de possiveis quadros
# 1) Lista de culturas e numero de parcelas e area correspondente, pre e pos seleccao de dados para o estudo (pode se tambem atribuir um codigo a cada cultura e depois meter esse codigo nos graficos)
# - Lista de imagens Landsat usadas
# 2) Resultados da correccao de imagens
# - Comparacao entre classificadores (usando que medidas? so a percentagem de boas correccoes? total ou por cultura? e o kappa?)
# - Matriz de erro global de um dos classificadores (tinha que ser a soma das matrizes de erro de cada uma das 10 validacoes na validacao cruzada)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1) Lista com culturas, parcelas e areas

#Numero de parcelas
nParc <- as.vector(table(listaTodasParcelas$cultura))

#Areas
dados <- data.table(listaTodasParcelas)
areasParc <- dados[,list(cultura=cultura,
                         area=sum(area),
                         areaMed=mean(area),
                         areaSD=sd(area)),by=cultura]
areasParc <- as.data.frame(areasParc)
areasParc <- areasParc[order(areasParc$cultura, decreasing = F),]
areaTot <- sum(areasParc$area)
relArea <- (areasParc$area/areaTot)*100

#Criacao da tabela
sumarioParc <- data.frame(NOMES_CULT, 
                          COD_CULT,
                          nParc,
                          round(areasParc$area),
                          round(relArea, 1),
                          paste(round(areasParc$areaMed, 1), 'pm', round(areasParc$areaSD, 1)))
colnames(sumarioParc) <- c('Land cover class', 'Class label', 'N', 'Total area (ha)', 'Relative area (%)', 'Mean area (ha)')

#LaTeX output
sumarioParc.xtab <- xtable(sumarioParc)
align(sumarioParc.xtab) <- rep("l", 5)
digits(sumarioParc.xtab) <- 1
print.xtable(sumarioParc.xtab, booktabs=T, include.rownames=F)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2) Resultados da correccao de imagens

resCorrec <- read.xlsx('sumario.xlsx', 1, header = T, , colClasses=c('character', 'character' ,rep('numeric',5)))
resCorrec <- resCorrec[,c(-8,-9,-10)]


resCorrec.xtab <- xtable(resCorrec)
align(resCorrec.xtab) <- rep("l", 8)
digits(resCorrec.xtab) <- 3
print.xtable(resCorrec.xtab, booktabs=T, include.rownames=F)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3) Resultados do accuracy assessment

KNN.1.user <- KNN.comp.cruz$result$userA
KNN.1.prod <- KNN.comp.cruz$result$prodA
KNN.2.user <- KNN.sub.cruz$result$userA
KNN.2.prod <- KNN.sub.cruz$result$prodA
SVM.1.user <- SVM.comp.cruz$result$userA
SVM.1.prod <- SVM.comp.cruz$result$prodA
SVM.2.user <- SVM.sub.cruz$result$userA
SVM.2.prod <- SVM.sub.cruz$result$prodA

assessment <- data.frame(KNN.1.user, KNN.1.prod, KNN.2.user, KNN.2.prod, SVM.1.user, SVM.1.prod, SVM.2.user, SVM.2.prod)
assessment <- rbind(assessment, c(KNN.comp.cruz$result$correcTot*100,
                                  NA,
                                  KNN.sub.cruz$result$correcTot*100,
                                  NA,
                                  SVM.comp.cruz$result$correcTot*100,
                                  NA,
                                  SVM.sub.cruz$result$correcTot*100,
                                  NA))

assessment <- cbind(c(COD_CULT, 'Overall accuracy'), assessment)


assessment.xtab <- xtable(assessment)
align(assessment.xtab) <- rep("r", 10)
digits(assessment.xtab) <- 2
print.xtable(assessment.xtab, booktabs=T, include.rownames=F)

