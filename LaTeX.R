#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Master's Thesis - Remote Sensing                      %
# Environmental Engineering - ISA/UL - Lisbon, Portugal %
# (c) 2014 by Jonas Schmedtmann                         %
#                                                       %
# LATEX SCRIPT                                          %
#                                                       %
# This script takes care of exporting tables and plots  %
# to LaTeX using tikzDevice and xtable                  %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#COLOR PALETES
#rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n, alpha = 1)
#heat.colors(n, alpha = 1)
#terrain.colors(n, alpha = 1)
#topo.colors(n, alpha = 1)
#cm.colors(n, alpha = 1)



#Init
opar <- par()
par(opar)


#Graficos
tikz('../05 Escrito/TEX/areasCulturas.tex',width=5,height=4)
par(mar=c(3.5,3.5,0.5,0.5), mgp=c(2,0.6,0))   #bottom, left, top, right
barplot(compProbs.svm, beside=T, space=c(0,5), legend=c("Completo","Sub-conjunto"), xlab="Probabilidade", ylab="Numero de parcelas", border=NA, col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()
dev.off()






tikz('../05 Escrito/TEX/assinaturas.tex',width=5,height=4)
par(las=1, mar=c(4.5,4.5,1,1), mgp=c(2.7,0.8,0))

#Pastagem
plot(BANDAS.LANDSAT,assinaturasPorCulturaD3[1,3:8],type='l',ylim=c(0.05,0.3),col=2,xlab='Banda Landsat',ylab='Reflectancia')
text(5,assinaturasPorCulturaD3[1,7]+0.01, as.character(codigosNovos2005$NOVO_NOME[assinaturasPorCulturaD3[1,1]]),col=2)
#Sup for
lines(BANDAS.LANDSAT,assinaturasPorCulturaD3[2,3:8],col=3)
text(5,assinaturasPorCulturaD3[2,7]-0.005, as.character(codigosNovos2005$NOVO_NOME[assinaturasPorCulturaD3[2,1]]),col=3)
#Milho
lines(BANDAS.LANDSAT,assinaturasPorCulturaD3[3,3:8],col=7)
text(5,assinaturasPorCulturaD3[3,7]+0.01, as.character(codigosNovos2005$NOVO_NOME[assinaturasPorCulturaD3[3,1]]),col=7)
#Pousio
lines(BANDAS.LANDSAT,assinaturasPorCulturaD3[5,3:8],col=4)
text(5.4,assinaturasPorCulturaD3[5,7]+0.007, as.character(codigosNovos2005$NOVO_NOME[assinaturasPorCulturaD3[5,1]]),col=4)
#Arroz
lines(BANDAS.LANDSAT,assinaturasPorCulturaD3[4,3:8],col=6)
text(5,assinaturasPorCulturaD3[4,7]+0.01, as.character(codigosNovos2005$NOVO_NOME[assinaturasPorCulturaD3[4,1]]),col=6)
#cevada
lines(BANDAS.LANDSAT,assinaturasPorCulturaD3[10,3:8],col=8)
text(5,assinaturasPorCulturaD3[10,7]+0.01, as.character(codigosNovos2005$NOVO_NOME[assinaturasPorCulturaD3[10,1]]),col=8)

dev.off()



tikz('../05 Escrito/TEX/signatures.tex',width=6,height=8)

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


dev.off()



#Quadros
dados.xtab <- xtable(dadosTreino[1:10,1:5])
align(dados.xtab) <- rep("c",6)
digits(dados.xtab) <- 3
print.xtable(dados.xtab, booktabs=T, include.rownames=F)




