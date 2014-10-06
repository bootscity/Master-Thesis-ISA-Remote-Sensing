#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Master's Thesis - Remote Sensing                      %
# Environmental Engineering - ISA/UL - Lisbon, Portugal %
# (c) 2014 by Jonas Schmedtmann                         %
#                                                       %
# LARGE OBJECT MANAGEMENT SCRIPT                        %
#                                                       %
# This handles the writing and reading of large objects %
# into and from .Robj files.                            %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Carregar objectos grandes

load("cdg2005_29.07.Robj")
load("todasImagens2005_15.06.Robj")
load('listaDados2005_13.09.Robj')
load("listaTodasParcelas2005_13.09.Robj")
load('listaTodasParcelasIniciais2005_29.07.Robj')

#Gravar objectos grandes como ficheiros

save(cdg, file='cdg2005_29.07.Robj')
save(todasImagens, file="todasImagens2005_13.09.Robj")
save(listaDados, file='listaDados2005_13.09.Robj')
save(finalListaDados, file='finalListaDados2005_28.09.Robj')
save(listaTodasParcelasIniciais, file="listaTodasParcelasIniciais2005_13.09.Robj")
save(listaTodasParcelas, file="listaTodasParcelas2005_13.09.Robj")
save(s, file="sRasterStackMapa.Robj")
save(sTrue, file="sTrueRasterStackMapa.Robj")


