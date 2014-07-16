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

load("todasImagens2005_15.06.Robj")
load('listaDados2005_15.06.Robj')
load("listaTodasParcelas2005_15.06.Robj")
load('listaTodasParcelasIniciais2005_15.06.Robj')

#Gravar objectos grandes como ficheiros

save(cdg, file='cdg2005_15.06.Robj')
save(todasImagens, file="todasImagens2005_15.06.Robj")
save(listaDados, file='listaDados2005_15.06.Robj')
save(listaTodasParcelasIniciais, file="listaTodasParcelasIniciais2005_15.06.Robj")
save(listaTodasParcelas, file="listaTodasParcelas2005_15.06.Robj")
save(resultValidCruz, file="resultValidCruz_01.07.Robj")
save(s, file="stackLE72040332005157EDC00longLat.Robj")

#Eliminar grandes objectos depois de escritos
rm(listaDados)
rm(cdg)
rm(todasImagens)
rm(resultValidCruz)
