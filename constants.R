#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Master's Thesis - Remote Sensing                      %
# Environmental Engineering - ISA/UL - Lisbon, Portugal %
# (c) 2014 by Jonas Schmedtmann                         %
#                                                       %
# CONSTANTS SCRIPT                                      %
#                                                       %
# This script contains the constants used in this       %
# project and will be called by the main() function.    %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#Caminhos para shapefiles e imagens de satelite
CAMINHO.SHP <<- "../01 Dados/02 Dados espaciais"
CAMINHO.LANDSAT <<- "../01 Dados/03 Imagens satelite/LANDSAT"

#Strings de projeccao de informacao geografica
PROJ4.IGEOE <<- CRS("+proj=tmerc +lat_0=39.66666666666666 +lon_0=-8.131906111111112 +k=1 +x_0=200000 +y_0=300000 +ellps=intl +units=m +nadgrids=ptLX_e89.gsb +wktext +no_defs")
PROJ4.ETRS <<- CRS("+proj=tmerc +lat_0=39.66825833333333 +lon_0=-8.133108333333334 +k=1 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")
PROJ4.UTM <<- CRS("+proj=utm +zone=29 +ellps=WGS84 +datum=WGS84 +units=m")
ANOS <<- c(2003,2004,2005)

#Constantes relativas a correccao de imagens de Landsat
BANDAS.LANDSAT <<- c(1,2,3,4,5,7)
ESOLAR.TM5 <<- c(1983,1796,1536,1031,220,83.44)
ESOLAR.ETM7 <<- c(1997,1812,1533,1039,230.8,84.9)
TAMANHO.CELULA <<- 30
#http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20090027884.pdf



