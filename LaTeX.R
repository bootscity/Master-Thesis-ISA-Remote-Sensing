#########################################################
# Master's Thesis - Remote Sensing                      #
# Environmental Engineering - ISA/UL - Lisbon, Portugal #
# (c) 2014 by Jonas Schmedtmann                         #
#                                                       #
# LATEX SCRIPT                                          #
#                                                       #
# This script takes care of exporting tables and plots  #
# to LaTeX using tikzDevice and xtable                  #
#########################################################


#Init
opar <- par()
par(opar)
library(tikzDevice)
library(xtable)


#Graficos
tikz('../05 Escrito/TEX/test.tex',width=5,height=4)
par(mar=c(3.5,3.5,0.5,0.5), mgp=c(2,0.6,0))   #bottom, left, top, right
barplot(compProbs.svm, beside=T, space=c(0,5), legend=c("Completo","Sub-conjunto"), xlab="Probabilidade", ylab="Numero de parcelas", border=NA, col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()
dev.off()


#Quadros
dados.xtab <- xtable(dadosTreino[1:10,1:5])
align(dados.xtab) <- rep("c",6)
digits(dados.xtab) <- 3
print.xtable(dados.xtab, booktabs=T, include.rownames=F)


