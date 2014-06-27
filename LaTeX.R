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



#Quadros
dados.xtab <- xtable(dadosTreino[1:10,1:5])
align(dados.xtab) <- rep("c",6)
digits(dados.xtab) <- 3
print.xtable(dados.xtab, booktabs=T, include.rownames=F)




