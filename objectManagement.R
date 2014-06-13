#########################################################
# Master's Thesis - Remote Sensing                      #
# Environmental Engineering - ISA/UL - Lisbon, Portugal #
# (c) 2014 by Jonas Schmedtmann                         #
#                                                       #
# LARGE OBJECT MANAGEMENT SCRIPT                        #
#                                                       #
# This handles the writing and reading of large objects #
# into and from .Robj files.                            #
#########################################################


#Carregar objectos grandes
load("cdg2005.Robj")
load("todasImagens2005.Robj")
load("listaDados2005.Robj")
load("listaTodasParcelasIniciais2005.Robj")


#Gravar objectos grandes como ficheiros
save(cdg, file="cdg2005.Robj")
save(todasImagens, file="todasImagens2005.Robj")
save(listaDados, file="listaDados2005.Robj")
save(listaTodasParcelasIniciais, file="listaTodasParcelasIniciais2005.Robj")


#Eliminar grandes objectos depois de escritos
rm(listaDados)
rm(cdg)
rm(todasImagens)

