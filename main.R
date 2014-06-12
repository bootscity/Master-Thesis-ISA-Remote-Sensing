##########################################################################################
# 0. FUNCOES                                                                             #
##########################################################################################

# Indice de funcoes:
#   01. init
#   02. carregaRecortaDadosEspaciais
#   04. corrigeLeConjuntoImagensLandsat
#   05. constroiListaImagensLandsat
#   06. constroiListaDados
#   XX. obtemConteudoParcela
#   XX. parcNDVI
#   XX. parcNIFAP
#   XX. constroiDadosANOVA
#   XX. verificaAprovacaoFinalValidacao
#   XX. constroiTodasParcelas


#1. Funcao init:
#   funcao: inicializa o programa: carrega packages, define todas as constantes superglobais
#   input:  NA
#   output: NA
init<-function()
{
  #Carregamento de packages
  library(maptools)
  library(raster)
  library(sp)
  library(rgeos)
  library(rgdal)
  library(proj4)
  library(xlsx)
  library(MASS)
  library(data.table)
  library(e1071)
  library(class)
  library(tree)
  library(party)
  library(rpart)
  library(subselect)
  library(plotrix)
  
  #Definicao de codigos
  codigos2005<-read.delim("culturas2005.txt",header=T)
  codigosNovos2005<-read.delim("culturasNovas2005.txt",header=T)
  codNomesAbrev<-c('Past. perm.', 'Sup. forrageira','Arroz','Milho','Trigo','Pousio','Cevada','Vinha','Aveia','Lolium','Sup. nao usada', 'Past. pobre','Betteraba sac.','Olival','Sorgo','Girassol','Horticolas')
  codigosNovos2005<-cbind(codigosNovos2005,codNomesAbrev)
  
  #Definicao de constantes
  CAMINHO.SHP<<-"../01 Dados/02 Dados espaciais"
  CAMINHO.LANDSAT<<-"../01 Dados/03 Imagens satelite/LANDSAT"
  
  PROJ4.IGEOE<<-CRS("+proj=tmerc +lat_0=39.66666666666666 +lon_0=-8.131906111111112 +k=1 +x_0=200000 +y_0=300000 +ellps=intl +units=m +nadgrids=ptLX_e89.gsb +wktext +no_defs")
  PROJ4.ETRS<<-CRS("+proj=tmerc +lat_0=39.66825833333333 +lon_0=-8.133108333333334 +k=1 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")
  PROJ4.UTM<<-CRS("+proj=utm +zone=29 +ellps=WGS84 +datum=WGS84 +units=m")
  ANOS<<-c(2003,2004,2005)
  
  #http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20090027884.pdf
  BANDAS.LANDSAT<<-c(1,2,3,4,5,7)
  ESOLAR.TM5<<-c(1983,1796,1536,1031,220,83.44)
  ESOLAR.ETM7<<-c(1997,1812,1533,1039,230.8,84.9)
  
  TAMANHO.CELULA<<-30
  
}


#2. Funcao carregaRecortaDadosEspaciais:
#   funcao: importar todas as shapefiles necessarias ao trabalho, transformar para sistema de coordenadas
#           indicado e recortar parcelas pela area de estudo
#   input:  coordsAreaEstudo: matriz de coordenadas X,Y que vai ser a area de estudo
#           crs: string de projeccao da classe CRS
#   output: lista com os diferentes cdg's transformados e recortados (apenas as parelas)
carregaRecortaDadosEspaciais<-function(coordsAreaEstudo,crs)
{
  #Ler os dados dos ficheiros shp
  #p2003<-readShapePoly(paste0(CAMINHO.SHP,"/2003/2003_PARC"), proj4string=PROJ4.IGEOE)
  #p2004<-readShapePoly(paste0(CAMINHO.SHP,"/2004/2004_PARC"), proj4string=PROJ4.IGEOE)
  p2005<-readShapePoly(paste0(CAMINHO.SHP,"/2005/2005_PARC"), proj4string=PROJ4.IGEOE)
  #c<-readShapePoly(paste0(CAMINHO.SHP,"/Enquadramentos PT IGEO/CONCELHOS"), proj4string=PROJ4.ETRS)
  
  #Transformar cdg's para sistema de coordenadas dado no argumento crs
  #p2003Trans<-spTransform(p2003, crs)
  #p2004Trans<-spTransform(p2004, crs)
  p2005Trans<-spTransform(p2005, crs)
  #cTrans<-spTransform(c, crs)
  
  #Construir area de estudo a partir da matriz de entrada
  areaEstudo<-SpatialPolygons(list(Polygons(list(Polygon(coordsAreaEstudo)), "1")),proj4string=crs)
  
  #Recortar as parcelas pela area de estudo
  #p2003inter<-gIntersects(p2003Trans, areaEstudo, byid=TRUE)
  #p2003clip<-subset(p2003Trans,as.vector(p2003inter))
  #p2003clip@data$ID<-1:length(p2003clip@data$ID)  
  
  p2005inter<-gIntersects(p2005Trans, areaEstudo, byid=TRUE)
  p2005clip<-subset(p2005Trans,as.vector(p2005inter))
  p2005clip@data$ID<-1:length(p2005clip@data$ID)  
  
  #out<-list(parc2003=p2003clip, parc2004=p2004clip, parc2005=p2005clip, area=areaEstudo, conc=c)
  out<-list(parc2005=p2005clip, area=areaEstudo)
  return(out)
}


#4. Funcao corrigeLeConjuntoImagensLandsat:
#   funcao: fazer correccao radiometrica de um conjunto de imagens landsat em formato raster com metodo
#           DOS4, ou apenas ler imagens ja corrigidas
#   input:  conjuntoPath: caminho para o conjunto de imagens em formato raster a corrigir
#           areaEstudo:   cdg com poligono pelo qual as imagens serao cortadas
#           prefixo:      string que e juntada ao nome original da imagem, na imagem corrigida
#           corrige:      FALSE se nao for para corrigir as imagens, mas sim ler imagens ja corrigidas
#   output: lista de imagens corrigidas em formato raster, essas imagens tambem sao escritas em ficheiros
corrigeLeConjuntoImagensLandsat<-function(conjuntoPath, areaEstudo, prefixo, corrige=TRUE)
{
  out<-list()
  splitConjuntoPath<-strsplit(conjuntoPath, "/")
  codigoImagens<-splitConjuntoPath[[1]][length(splitConjuntoPath[[1]])]
  
  #O QUE FAZER SE FOR PARA CORRIGIR OS DADOS (DEFAULT)
  if (corrige == TRUE)
  {
    #Metadados
    metadataFile<-paste0(conjuntoPath,"/",codigoImagens,"_MTL.txt")
    metadata<-read.table(metadataFile, nrows = 186)
    
    resultadosDOS4<-data.frame(Imagem=NA,k=NA,tau=NA,T_v=NA,T_z=NA,E_down=NA,L_p=NA,ro_min=NA,ro_max=NA,ro_mean=NA)
    
    #Corrige cada imagem do conjunto dado
    for(j in 1:length(BANDAS.LANDSAT))
    {
      #Abrir GeoTIFF
      nomeImagem<-paste0(codigoImagens,"_B",BANDAS.LANDSAT[j],".TIF")
      ficheiroImagem<-paste0(conjuntoPath,"/",nomeImagem)
      imagem<-raster(ficheiroImagem)
      
      #Fazer clip do GEOTIFF pela area de estudo
      imagemArea<-crop(imagem, areaEstudo)
      
      #Mascarar valores que estao fora da area de estudo como NA
      imagemArea<-mask(x=imagemArea, mask=areaEstudo, updatevalue=NA)
      
      #1) Radiancia no sensor
      #A banda 7 e diferente porque existem duas bandas 6
      if (BANDAS.LANDSAT[j] == 7)
      {
        Lmax<-as.numeric(as.character(metadata$V3[83]))
        Lmin<-as.numeric(as.character(metadata$V3[84]))
        Qmax<-as.numeric(as.character(metadata$V3[103]))
        Qmin<-as.numeric(as.character(metadata$V3[104]))
      } else {
        Lmax<-as.numeric(as.character(metadata$V3[69+2*BANDAS.LANDSAT[j]-2]))
        Lmin<-as.numeric(as.character(metadata$V3[69+2*BANDAS.LANDSAT[j]-1]))
        Qmax<-as.numeric(as.character(metadata$V3[89+2*BANDAS.LANDSAT[j]-2]))
        Qmin<-as.numeric(as.character(metadata$V3[89+2*BANDAS.LANDSAT[j]-1]))
      }
      
      LSat<-Lmin+(Lmax-Lmin)*((imagemArea-Qmin)/(Qmax-Qmin))
      
      #2) Correccao atmosferica - DOS4
      #Estimacao iterativa de T_v, T_z, E_down, L_p
      LSatMin<-cellStats(LSat, min)
      doyStr<-substr(as.character(metadata$V3[5]),14,16)
      doy<-as.numeric(doyStr)
      d<-(1-0.017^2)/(1+0.017*cos(365.25/360*doy*(pi/180)))
      E0<-ESOLAR.ETM7[j]/d^2
      tetaElev<-as.numeric(as.character(metadata$V3[62]))
      tetav<-2*pi/360*(90-tetaElev)
      tetaz<-0
      
      Tv<-1
      Tz<-1
      EDown<-0
      tau<-0
      k<-0
      repeat
      {
        Lp<-LSatMin-0.01*(E0*cos(tetaz)*Tz+EDown)*Tv/pi
        EDown<-pi*Lp
        oldTau<-tau
        tau<--cos(tetaz)*log(1-((4*pi*(LSatMin-0.01*(E0*cos(tetaz)*Tz+EDown)*Tv/pi))/(E0*cos(tetaz))))
        Tv<-exp(-tau/cos(tetav))
        Tz<-exp(-tau/cos(tetaz))
        k<-k+1
        if(abs(oldTau-tau)<0.000001){
          break
        }
      }
      if (Tv > 1) Tv<-1
      if (Tz > 1) Tz<-1
      if (EDown < 0) EDown<-0
      if (Lp < 0) Lp<-0
      
      #3) Reflectancia (objecto da classe Raster contendo imagem com reflectancia)
      reflect<-(pi*(LSat-Lp))/(Tv*(E0*cos(tetaz)*Tz+EDown))   
      reflect[reflect<0]<-0
      reflect[reflect>1]<-1
      
      #Resultados para data.frame resultadosDOS4
      roMin<-cellStats(reflect, min)
      roMax<-cellStats(reflect, max)
      roMean<-cellStats(reflect, mean)
      resultadosDOS4<-rbind(resultadosDOS4,c(nomeImagem,k,tau,Tv,Tz,EDown,Lp,roMin,roMax,roMean))
      #resultadosDOS4<-resultadosDOS4[-1,]
      
      #Fazer output do raster
      writeRaster(reflect,paste0(conjuntoPath,"/",prefixo,"_",nomeImagem), format="GTiff", overwrite=T)
      outItemName<-paste0("b",BANDAS.LANDSAT[j])
      out[[outItemName]]<-reflect
    }
    
    #Devolver resultados
    write.xlsx(resultadosDOS4, paste0(conjuntoPath,"/",prefixo,"_output.xlsx"), "Resultados DOS4", showNA=FALSE)
    
  }
  else
  {
    #O QUE FAZER SE FOR APENAS PARA LER OS DADOS *JA CORRIGIDOS* E NAO PARA CORRIGIR
    for(j in 1:length(BANDAS.LANDSAT))
    {
      ficheiroImagem<-paste0(conjuntoPath,"/",prefixo,"_",codigoImagens,"_B",BANDAS.LANDSAT[j],".TIF")
      reflect<-raster(ficheiroImagem)
      outItemName<-paste0("b",BANDAS.LANDSAT[j])
      out[[outItemName]]<-reflect
    }
  }
  
  return(out)
}


#5. Funcao constroiListaImagensLandsat:
#   funcao: construir uma lista com todas as imagens em forma de matriz de todas as datas Landsat que estao
#           numa dada directoria
#   input:  landsatPath: caminho para a directoria que contem todas as datas de Landsat
#           areaEstudo:  cdg com poligono pelo qual as imagens serao cortadas
#           prefixo:     string que e juntada ao nome original da imagem, na imagem corrigida
#           anos:        vector com os anos para os quais se pretende obter a lista de imagens
#           corrige: FALSE se nao for para corrigir as imagens, mas sim ler imagens ja corrigidas
#   output: lista com todas as datas e imagens Landsat em forma de matriz, incluindo detalhes da imagem
constroiListaImagensLandsat<-function(landsatPath, areaEstudo, prefixo, anos=ANOS, corrige=FALSE)
{
  conjuntos<-list.files(landsatPath, full.names = FALSE)
  todosAnos<-0;doys<-0
  out<-list()
  
  #Leitura dos dados a partir dos nomes das directorias
  for(i in 1:length(conjuntos))
  {
    todosAnos[i]<-as.numeric(substr(conjuntos[i], 10, 13))
    doys[i]<-as.numeric(substr(conjuntos[i], 14, 16))
  }
  
  for(i in 1:length(anos))
  {
    #ANO
    ano<-anos[i]
    strAno<-paste0("a", ano)
    out[[strAno]]<-list()
    
    #DATAS E IMAGENS
    for(j in 1:length(doys))
    {
      if(ano != todosAnos[j]) next
      doy<-doys[j]
      strDoy<-paste0("d", doy)
      out[[strAno]][[strDoy]][["DR"]]<-corrigeLeConjuntoImagensLandsat(paste0(landsatPath,"/",conjuntos[j]), areaEstudo, prefixo, corrige)
      
      #Obter dados de coordenadas min e max
      img<-out[[strAno]][[strDoy]][["DR"]][["b1"]]
      xmin<-xmin(img)
      ymax<-ymax(img)
      cellsize<-res(img)[1]
      TAMANHO.CELULA<<-cellsize
      out[[strAno]][[strDoy]][["info"]]<-list("xmin"=xmin, "ymax"=ymax,"ymax2"=ymax, "cellsize"=cellsize)
      
      #obter coordenadas de todos os pontos
      coordRas<-coordinates(img)
      lin<-nrow(img)
      col<-ncol(img)
      xCoord<-matrix(coordRas[,1],nrow=lin,ncol=col,byrow=T)
      yCoord<-matrix(coordRas[,2],nrow=lin,ncol=col,byrow=T)
      out[[strAno]][[strDoy]][["X"]]<-xCoord
      out[[strAno]][[strDoy]][["Y"]]<-yCoord
      
      #Finalmente, converter as imagens a matrizes para mais tarde facilitar o processamento
      for (k in 1:length(out[[strAno]][[strDoy]][["DR"]]))
        out[[strAno]][[strDoy]][["DR"]][[k]]<-as.matrix(out[[strAno]][[strDoy]][["DR"]][[k]])
    }
  }
  
  return(out)
}


#6. Funcao constroiListaDados:
#   funcao: construir uma lista com toda a informacao disponivel sobre parcelas e dados de reflectancias 
#   input:  
#   output: lista com toda a informacao disponivel sobre parcelas e dados de reflectancias, tendo esta forma:
#   lista$nNIFAP$aANO$pPARCELA$DR$dDATA$bBANDA
#                             $atributos
#                             $coordsPoly
constroiListaDados<-function(anos=ANOS)
{
  #Geracao dos dados de base para a criacao da lista, chamando outras funcoes
  #cdg<-leDadosEspaciais()
  #cdg<-transformaDadosEspaciais(cdg, PROJ4.UTM)
  todosNIFAP<-c(cdg$parc2005@data$NIFAP)#, cdg$parc2004@data$NIFAP, cdg$parc2005@data$NIFAP)
  unicosNIFAP<-sort(unique(todosNIFAP))
  listaGlobal<-list()
  
  #Aqui e feita a correccao de todas as imagens
  #todasImagens<<-constroiListaImagensLandsat(CAMINHO.LANDSAT, cdg$area, "CORR_12.05", anos=2005, corrige=F)
  
  #Construcao da lista
  for(i in 1:length(unicosNIFAP))
  {
    print(paste0("i = ",i))
    print(proc.time())
    #NIFAP
    NIFAP<-unicosNIFAP[i]
    strNIFAP<-paste0("n",NIFAP)
    listaGlobal[[strNIFAP]]<-list()
    
    #for(j in 1:1)
    for(j in 1:length(anos))
    {
      #ANO
      ano<-anos[j]
      strAno<-paste0("a", ano)
      listaGlobal[[strNIFAP]][[strAno]]<-list()
      
      #PARCELA
      esteCDG<-cdg[[j]]
      select<-esteCDG@data[esteCDG$NIFAP==NIFAP,]
      if (nrow(select) == 0) next
      IDs<-as.vector(as.numeric(select$ID))
      parcelas<-as.vector(select$COD_CRUZAM)
      culturas<-as.vector(select$COD_CULTUR)
      sequeiro<-as.vector(select$COD_SEQ_RE)
      
      for(k in 1:length(parcelas))
      {
        #print(paste0("k = ",k))
        #print(proc.time())
        parcela<-parcelas[k]
        strParcela<-paste0("p", parcela)
        listaGlobal[[strNIFAP]][[strAno]][[strParcela]]<-list()
        esteID<-IDs[k]
        if (sequeiro[k] == "S"){
          esteSeq<-TRUE
        } else {
          esteSeq<-FALSE}
        
        #ATRIBUTOS
        estaArea<-esteCDG@polygons[[esteID]]@Polygons[[1]]@area
        listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["atributos"]]<-list(ID=esteID,
                                                                             CULTURA=culturas[k],
                                                                             AREA=estaArea,
                                                                             SEQEIRO=esteSeq)
        
        #COORDENADAS DAS PARCELAS
        estasCoordenadas<-coordinates(esteCDG@polygons[[esteID]]@Polygons[[1]])
        listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["coordsPoly"]]<-estasCoordenadas
        
        #DR EM DIFERENTES DATAS
        listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["DR"]]<-list()
        
        #Definicao do rectangulo para extrair com base nas coordenadas da parcela
        minX<-min(estasCoordenadas[,1])
        maxX<-max(estasCoordenadas[,1])
        minY<-min(estasCoordenadas[,2])
        maxY<-max(estasCoordenadas[,2])
        
        for(l in 1:length(todasImagens[[j]]))
        {
          #DATAS
          datas<-names(todasImagens[[j]])
          estaData<-datas[[l]]
          numBandas<-length(todasImagens[[j]][[l]][["DR"]])
          #listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["DR"]][[estaData]]<-array(dim=c(1,1,numBandas))
          
          #BANDAS
          for(m in 1:numBandas)
          {
            #print(m);print(proc.time())
            estaImg<-todasImagens[[j]][[l]][["DR"]][[m]]
            
            #Extrair valores de interesse da matriz
            minXglobal<-todasImagens[[j]][[l]][["info"]][["xmin"]]
            maxYglobal<-todasImagens[[j]][[l]][["info"]][["ymax"]]
            cellsize<-todasImagens[[j]][[l]][["info"]][["cellsize"]]
            
            selectMinX<<-floor((minX-minXglobal)/cellsize)
            selectMaxX<<-ceiling((maxX-minXglobal)/cellsize)
            selectMinY<<-floor((maxYglobal-maxY)/cellsize)
            selectMaxY<<-ceiling((maxYglobal-minY)/cellsize)
            
            if (selectMaxX >= 0 && selectMaxY >= 0 && selectMinX >= 0 && selectMinY >= 0)
              imgContent<-estaImg[selectMinY:selectMaxY,selectMinX:selectMaxX]
            
            #Criar um novo array na primeira iteracao
            if (m==1)
              listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["DR"]][[estaData]]<-array(dim=c(nrow(imgContent),ncol(imgContent),numBandas+2))
            
            listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["DR"]][[estaData]][,,m]<-imgContent            
          }
          
          if (selectMaxX >= 0 && selectMaxY >= 0 && selectMinX >= 0 && selectMinY >= 0)
          {
            listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["DR"]][[estaData]][,,m+1]<-todasImagens[[j]][[l]][["X"]][selectMinY:selectMaxY,selectMinX:selectMaxX]
            listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["DR"]][[estaData]][,,m+2]<-todasImagens[[j]][[l]][["Y"]][selectMinY:selectMaxY,selectMinX:selectMaxX]
          }
        }
      }
    }
  }
  return(listaGlobal)
}


#XX. Funcao obtemConteudoParcela
#   funcao: 
#   input:  
#   output: lista que para cada data tem uma matriz que tem uma coluna para cada banda, contendo todas as
#           reflectancias
obtemConteudoParcela<-function(parcela,apenasInterior=TRUE,excluiVizinhos=FALSE)
{
  out<-list()
  poly.x<-parcela$coordsPoly[,1]
  poly.y<-parcela$coordsPoly[,2]
  
  datas<-names(parcela$DR)
  for(i in 1:length(datas))
  {
    #out[[datas[i]]]<-list()
    #bandas<-names(parcela$DR[[datas[i]]])
    numBandas<-length(parcela$DR[[datas[i]]][1,1,])
    for(j in 1:(numBandas-2))  #-2 para nao incluir as "bandas" X e Y
    {
      #O que fazer se queremos excluir os vizinhos
      if(excluiVizinhos==FALSE)
      {
        #Obter coordenadas dos CENTROS dos pixeis
        img.x<-as.vector(parcela$DR[[i]][,,numBandas-1])
        img.y<-as.vector(parcela$DR[[i]][,,numBandas])
        
        #Obter a banda em questao
        bandaVec<-as.vector(parcela$DR[[i]][,,j])
      }
      else
      {
        #Obter coordenadas dos CENTROS dos pixeis
        mat.x<-parcela$DR[[i]][,,numBandas-1]
        img.x<-as.vector(mat.x[seq(1,nrow(mat.x),2),seq(1,ncol(mat.x),2)])
        
        mat.y<-parcela$DR[[i]][,,numBandas]
        img.y<-as.vector(mat.y[seq(1,nrow(mat.y),2),seq(1,ncol(mat.y),2)])
        
        #Obter a banda em questao
        mat<-parcela$DR[[i]][,,j]
        matSemViz<-mat[seq(1,nrow(mat),2),seq(1,ncol(mat),2)]
        bandaVec<-as.vector(matSemViz)
      }
      
      #O que fazer quando nao interessa se o pixel esta todo dentro da parcela ou nao
      if(apenasInterior==FALSE)
      {
        dentro<-point.in.polygon(img.x, img.y, poly.x, poly.y)
        conteudo<-bandaVec[dentro>0]
        
        #Criar uma nova matriz na primeira iteracao
        if(j==1)
          out[[datas[i]]]<-matrix(nrow=length(conteudo),ncol=numBandas-2)
        
        out[[datas[i]]][,j]<-conteudo
      }
      #O que fazer quando queremos apenas os pixeis com 100% dentro da parcela
      else
      {
        img.x.1<-img.x-TAMANHO.CELULA/2
        img.y.1<-img.y+TAMANHO.CELULA/2
        dentro.1<-point.in.polygon(img.x.1, img.y.1, poly.x, poly.y)
        img.x.2<-img.x+TAMANHO.CELULA/2
        img.y.2<-img.y+TAMANHO.CELULA/2
        dentro.2<-point.in.polygon(img.x.2, img.y.2, poly.x, poly.y)
        img.x.3<-img.x+TAMANHO.CELULA/2
        img.y.3<-img.y-TAMANHO.CELULA/2
        dentro.3<-point.in.polygon(img.x.3, img.y.3, poly.x, poly.y)
        img.x.4<-img.x-TAMANHO.CELULA/2
        img.y.4<-img.y-TAMANHO.CELULA/2
        dentro.4<-point.in.polygon(img.x.4, img.y.4, poly.x, poly.y)
        dentroFinal<-dentro.1*dentro.2*dentro.3*dentro.4
        conteudo<-bandaVec[dentroFinal>0]
        
        #Criar uma nova matriz na primeira iteracao
        if(j==1)
          out[[datas[i]]]<-matrix(nrow=length(conteudo),ncol=numBandas-2)
        
        out[[datas[i]]][,j]<-conteudo
      }
    }
  }
  return(out)
}


#XX. Funcao parcNDVI
#   funcao: 
#   input:  
#   output: lista que para cada data tem o NDVI de cada pixel na parcela
parcNDVI<-function(parcela,apenasInterior=TRUE)
{
  out<-list()
  conteudo<-obtemConteudoParcela(parcela,apenasInterior)
  datas<-names(conteudo)
  
  for(i in 1:length(datas))
  {
    NDVI<-(conteudo[[i]][,4]-conteudo[[i]][,3])/(conteudo[[i]][,4]+conteudo[[i]][,3])
    out[[datas[i]]]<-NDVI
  }
  return(out)
}


#XX. Funcao parcNIFAP
#   funcao: 
#   input:  
#   output: NIFAP ao qual pertence a parcela introduzida
parcNIFAP<-function(nomeParc)
{
  for(i in 1:length(listaDados))
    for(j in 1:length(listaDados[[i]]))
      for(k in 1:length(listaDados[[i]][[j]]))
      {
        if (names(listaDados[[i]][[j]])[[k]] == nomeParc)
          return(names(listaDados)[[i]])
      }
}


#XX. Funcao constroiDadosANOVA
#   funcao: 
#   input:  
#   output: data.frame com dados organizados para se efectuar uma ANOVA hierarquizada
constroiDadosANOVA<-function(ano,data,banda,dimAmostra)
{
  #ano<-which(ANOS==ano)    Isto so vai funcionar quando a listaDados tiver todos os anos
  numParcANOVA<<-1
  dadosANOVA<-matrix(nrow=0,ncol=3)
  set.seed(0)    #Depois e para tirar
  amostraNIFAPs<-sort(sample(1:length(listaDados),dimAmostra))
  
  for(i in amostraNIFAPs)
    #for(k in 1:length(listaDados[[i]][[ano]]))Isto so vai funcionar quando a listadados tiver todos os anos
    for(k in 1:length(listaDados[[i]][[1]]))
    {
      parc<-listaDados[[i]][[1]][[k]]
      reflectancias<-obtemConteudoParcela(parc,apenasInterior=TRUE,excluiVizinhos=TRUE)[[data]][,banda]
      
      #Passar ao ciclo seguinte se houver um unico NA na parcela
      if(any(is.na(reflectancias))) next
      
      cultura<-rep(parc$atributos$CULTURA,length(reflectancias))
      parcela<-rep(names(listaDados[[i]][[1]])[[k]],length(reflectancias))
      novaParc<-cbind(reflectancias,cultura,parcela)
      dadosANOVA<-rbind(dadosANOVA,novaParc)
      numParcANOVA<<-numParcANOVA+1
    }
  
  #Converter resultado para data.frame com formatos correctos
  dadosANOVA<-as.data.frame(dadosANOVA)
  dadosANOVA$reflectancias<-as.numeric(as.vector(dadosANOVA$reflectancias))
  dadosANOVA$cultura<-as.factor(dadosANOVA$cultura)
  dadosANOVA$parcela<-as.factor(dadosANOVA$parcela)
  
  return(dadosANOVA)
}


#XX. Funcao verificaAprovacaoFinalValidacao
#   funcao: 
#   input:  NIFAP: string na forma 'nXXXXXX' com o numero do processo que se pretende verificar. Para
#                  efeitos de VALIDACAO, o processo de entrada e composto por parcelas cuja ocupacao
#                  do solo e CONHECIDA
#           ano:   ano para o qual foi feita a verificacao
#   output: 0/1 se o processo obtem aprovacao de financiamento - VALIDACAO
verificaAprovacaoFinalValidacao<-function(NIFAP,ano)
{
  processo<-listaDados[[NIFAP]][[paste0("a",ano)]]
  
  #Calculo da area do processo e obtencao das culturas e areas de cada parcela
  areas<-c()
  culturas<-c()
  for(i in 1:length(processo))
  {
    areas[i]<-processo[[i]]$atributos$AREA/10000
    culturas[i]<-processo[[i]]$atributos$CULTURA
  }
  somaAreas<-sum(areas)
  dados<-data.table(culturas,areas)
  dadosPorCultura<-dados[,list(areasCult=sum(areas)),by=culturas]
  dadosPorCultura<-dadosPorCultura[order(dadosPorCultura$areasCult, decreasing = TRUE),]
  
  ##########################################################################
  #Criterio 1
  dadosC1<-rbind(dadosPorCultura[dadosPorCultura$culturas == 143,], dadosPorCultura[dadosPorCultura$culturas == 144,])
  somaAreasC1<-sum(dadosC1$areasCult)
  percentagemC1<-somaAreasC1/somaAreas
  if(percentagemC1 >= 0.05) C1<-1 else C1<-0
  
  ##########################################################################
  #Criterio 2
  
  #isto ainda nao esta bem, porque considera cada COD_CULTURA como uma cultura diferente. E preciso juntar aquelas que sao culturas iguais (p.e. trigo duro e trigo mole?) e retirar coisas como inclutos ou pousios
  
  #retirar da data.frame pastagens e pousios e incultos (ocupacoes que nao sao culturas)
  codNaoCulturas<-c(21,26,31,52,87,88,89,121,125,138,666)
  codPastagem<-c(143,144)
  codExclusao<-c(codNaoCulturas,codPastagem)
  dadosPorCulturaC2<-dadosPorCultura[!(dadosPorCultura$culturas %in% codExclusao),]
  
  #Criterio 2.1
  if(somaAreas < 10)
  {
    C2<-1
    numCulturas<--1
    percentagemC2<--1
  }
  
  #Criterio 2.2
  else if(somaAreas >= 10 && somaAreas < 30)
  {
    numCulturas<-length(dadosPorCulturaC2$culturas)
    areaCultPrinc<-dadosPorCulturaC2$areasCult[1]
    percentagemC2<-areaCultPrinc/somaAreas
    
    if(numCulturas >= 2 && percentagemC2 <= 0.75) C2<-1 else C2<-0
  }
  
  #Criterio 2.3
  else if(somaAreas > 30)
  {
    numCulturas<-length(dadosPorCulturaC2$culturas)
    if(numCulturas >= 2)
      areaCultPrinc<-sum(dadosPorCulturaC2$areasCult[c(1,2)])
    else
      areaCultPrinc<-sum(dadosPorCulturaC2$areasCult[1])
    
    percentagemC2<-areaCultPrinc/somaAreas
    
    if(numCulturas >= 3 & percentagemC2 <= 0.95) C2<-1 else C2<-0
  }
  
  ##########################################################################
  #Crierio 3 nao pode ser implementado
  
  ##########################################################################
  #Decisao final
  final<-C1*C2
  return(data.frame(NIFAP=NIFAP,somaAreas=somaAreas, percentagemC1=percentagemC1, C1=C1, numCulturas=numCulturas, percentagemC2=percentagemC2, C2=C2, final=final))
}


#XX. Funcao constroiTodasParcelas
#   funcao: 
#   input:  
#   output: data.frame com informacoes relevantes sobre todas as parcelas existentes na listaDados
constroiTodasParcelas<-function(excluiVazias=TRUE, funcao=mean)
{
  count<-1
  for(i in 1:length(listaDados))
    #for(i in 1:10)
    for(j in 1:length(listaDados[[i]]))
      for(k in 1:length(listaDados[[i]][[j]]))
      {
        print(paste0("count = ",count))
        count<-count+1
        
        #Obter dados base para parcela
        parc<-listaDados[[i]][[j]][[k]]
        id<-names(listaDados[[i]][[j]])[k]
        nifap<-names(listaDados)[i]
        ano<-names(listaDados[[i]])[j]
        area<-parc$atributos$AREA/10000
        cultura<-parc$atributos$CULTURA
        seq<-parc$atributos$SEQEIRO
        
        #Obter reflectancias para diferentes datas
        reflectancias<-obtemConteudoParcela(parc,apenasInterior=TRUE)
        matRef<-matRefNomes<-matrix(nrow=length(reflectancias),ncol=ncol(reflectancias[[1]]))
        for(m in 1:length(reflectancias))
          for(n in 1:ncol(reflectancias[[1]]))
          {  
            matRef[m,n]<-funcao(reflectancias[[m]][,n])
            matRefNomes[m,n]<-paste0("B",n,"_d",m)
          }
        
        vecRef<-as.vector(matRef)
        vecRefNomes<-as.vector(matRefNomes)
        
        #Obter NDVI's para diferentes datas
        NDVIs<-parcNDVI(parc)
        vecNDVI<-c()
        for(l in 1:length(NDVIs))
          vecNDVI[l]<-mean(NDVIs[[l]])
        
        #Saltar esta parcela se houver algum NA e se isso for desejado
        if(excluiVazias == TRUE && (any(is.na(vecRef)) || any(is.na(vecNDVI)) || 0 %in% vecRef || -1 %in% vecNDVI)) next
        
        #Construir linha para data.frame
        if(i==1 && j==1 && k==1)
        {
          nDatas<-length(names(NDVIs))
          nomesCol<-c("id","NIFAP","ano","area","cultura","seq",vecRefNomes,paste0("NDVI_d",1:nDatas))
          out<-matrix(nrow=0,ncol=length(nomesCol))
          colnames(out)<-nomesCol
        }
        linha<-c(id,nifap,ano,as.numeric(area),cultura,seq,vecRef,vecNDVI)
        
        out<-rbind(out,linha)
        
      }
  
  #Construir data.frame com numeros
  out<-as.data.frame(out)
  out$area<-as.numeric(as.character(out$area))
  out$cultura<-as.numeric(as.character(out$cultura))
  
  for(i in 7:length(out))
    out[,i]<-as.numeric(as.character(out[,i]))
  
  return(out)
}



##########################################################################################
# 1. AQUISICAO DE DADOS E PROCESSAMENTO                                                  #
##########################################################################################

#Clear console
cat("\014")

#Iniciar
init()

#Definir area de estudo, carregar dados, recortar e projectar
areaX<-c(508262,528870,511269,490710,508262)
areaY<-c(4365569,4362462,4288360,4291713,4365569)
coordsArea<-cbind(areaX,areaY)
cdg<-carregaRecortaDadosEspaciais(coordsArea,PROJ4.UTM)

#Obter todas as imagens
rm(todasImagens)
todasImagens<<-constroiListaImagensLandsat(CAMINHO.LANDSAT, cdg$area, "CORR_12.05", anos=2005, corrige=TRUE)
#todasImagens<<-constroiListaImagensLandsat(CAMINHO.LANDSAT, cdg$area, "CORR_12.05", anos=2005, corrige=FALSE)

#Obter lista de todos os dados
rm(listaDados)
listaDados<-constroiListaDados(anos=2005)

#Obter conteudo relevante para TODAS as parcelas, em forma de data.frame
listaTodasParcelasIniciais<-constroiTodasParcelas()
listaTodasParcelasIniciaisSD<-constroiTodasParcelas(funcao=sd)



##########################################################################################
# 2. SELECCAO DE DADOS                                                                   #
##########################################################################################

##########################################################################
#Excluir a partida codigos que fazem sentido considerar
codExclusao<-c(87,88,666)
listaTodasParcelas<-listaTodasParcelasIniciais[!(listaTodasParcelasIniciais$cultura %in% codExclusao),]


##########################################################################
#Determinar as culturas que ocupam maior parte da ?rea (90% ou 95%)

#Todas as culturas com respectiva area total
dados<-data.table(listaTodasParcelas)
areasCulturas<-dados[,list(area=sum(area),numParc=length(area)),by=cultura]
areasCulturas<-areasCulturas[order(areasCulturas$area, decreasing = TRUE),]
areasCulturas<-cbind(areasCulturas,cumsum(areasCulturas$area))
areaTotal<-sum(areasCulturas$area)

#Culturas acima e abaixo dos limiares
plot(1:nrow(areasCulturas),areasCulturas$V2)
abline(h=areaTotal*0.9,col='blue')
abline(h=areaTotal*0.95,col='orange')
abline(h=areaTotal*0.98,col='red')

cultInfl<-areasCulturas[areasCulturas$V2 < areaTotal*0.95,]
fraccaoInfl<-sum(cultInfl$area)/areaTotal
cultInfluentes<-cultInfl$cultura


##########################################################################
#Seleccionar apenas aquelas que representam 95% da ?rea e reclassificar
listaTodasParcelas<-listaTodasParcelas[listaTodasParcelas$cultura %in% cultInfluentes,]

novasClasses98<-as.data.frame(cbind(cultInfluentes,c(1,2,3,4,5,6,7,6,8,9,10,11,12,13,14,5,15,4,16,17,8)))
novasClasses<-as.data.frame(cbind(cultInfluentes,c(1,2,3,4,5,6,7,6,8,9,10,11,12,13,14,5)))
colnames(novasClasses)<-c('cultura','novaClasse')
nClasses<-length(table(novasClasses$novaClasse))

for(i in 1:length(listaTodasParcelas$cultura))
  listaTodasParcelas$cultura[i]<-novasClasses$novaClasse[which(novasClasses$cultura==listaTodasParcelas$cultura[i])]


##########################################################################
#Preparacao dos dados para os classificadores

dadosClassificadores<-listaTodasParcelas[,c(5,7:12,15:20,23:28,31:36,39:44,47:52)]
dadosClassificadores$cultura<-as.factor(dadosClassificadores$cultura)

#Amostra a usar para o modelo - metade para treino e metade para validacao
set.seed(23)
amostra<-sample(1:nrow(dadosClassificadores), ceiling(nrow(dadosClassificadores)/2))

#Seleccao dos dados
dadosTreino<-dadosClassificadores[amostra,]
dadosValidacao<-dadosClassificadores[-amostra,]



##########################################################################################
# 3. ANALISE EXPLORATORIA DOS DADOS                                                      #
##########################################################################################

##########################################################################
#ANOVA para averiguar se cultura e parcela afectam reflectancias

###################
#Primeira abordagem
dadosANOVA<-constroiDadosANOVA(ano=2005,data=5,banda=1,dimAmostra=15)
table(dadosANOVA$cultura)
plot(dadosANOVA)

#Modelo da ANOVA
dados.aov<-aov(reflectancias ~ cultura/parcela, data=dadosANOVA)
summary(dados.aov)

coef(dados.aov)
model.tables(dados.aov, type="means")

#fazer analise de residuos(normalidade etc.)
plot(dados.aov,which=1)

#Teste para avergiuar se tem distrubuicao normal
shapiro.test(dadosANOVA$reflectancias)
#CONCLUSAO: e tudo menos normal...

plot(dadosANOVA$reflectancias)


###############
#Nova abordagem - averigua so se a cultura tem efeitos sobre cada uma das combinacoes banda/data e NDVI's
FValues<-c()
pValues<-c()
efeito<-c()
for(i in 7:length(listaTodasParcelas))
{
  parcelas.aov<-aov(listaTodasParcelas[,i] ~ as.factor(listaTodasParcelas$cultura))
  FV<-summary(parcelas.aov)[[1]][["F value"]][[1]]
  pV<-summary(parcelas.aov)[[1]][["Pr(>F)"]][[1]]
  
  FValues<-c(FValues,FV)
  pValues<-c(pValues,pV)
  if(pV <= 0.05) ef<-1 else ef<-0
  efeito<-c(efeito,ef)
}
nn<-colnames(listaTodasParcelas[,7:length(listaTodasParcelas)])
resultadoANOVA<-data.frame(cbind(nn,FValues,pValues,efeito))
resultadoANOVA$FValues<-as.numeric(as.character(resultadoANOVA$FValues))
resultadoANOVA<-resultadoANOVA[order(-resultadoANOVA$FValues),]


##########################################################################
#Assinaturas espectrais POR CULTURA

#Assinaturas todas
assinaturas<-matrix(nrow=0,ncol=length(BANDAS.LANDSAT))
for(i in 1:nrow(listaTodasParcelas))
{
  assinatura<-c()
  for(j in seq(0,6*7,8)+7)
    assinatura<-c(assinatura,listaTodasParcelas[i,j])
  
  assinaturas<-rbind(assinaturas,assinatura)
}

plot(BANDAS.LANDSAT,assinaturas[1,],type='l',ylim=c(0,0.35))
for(i in 1:length(assinaturas))
  lines(BANDAS.LANDSAT,assinaturas[i,],col=i)

#Assinaturas medias por ocupacao de solo
dados<-data.table(listaTodasParcelas)
assinaturasPorCulturaD3<-dados[,list(area=sum(area),
                                     B1_d3=mean(B1_d3),
                                     B2_d3=mean(B2_d3),
                                     B3_d3=mean(B3_d3),
                                     B4_d3=mean(B4_d3),
                                     B5_d3=mean(B5_d3),
                                     B6_d3=mean(B6_d3)),by=cultura]
assinaturasPorCulturaD3<-as.data.frame(assinaturasPorCulturaD3)

plot(BANDAS.LANDSAT,assinaturasPorCulturaD3[1,3:8],type='l',ylim=c(min(assinaturasPorCulturaD3$B6_d3),max(assinaturasPorCulturaD3$B4_d3)))
text(5,assinaturasPorCulturaD3[1,7],assinaturasPorCulturaD3[1,1])
for(i in 1:nrow(assinaturasPorCulturaD3))
{
  lines(BANDAS.LANDSAT,assinaturasPorCulturaD3[i,3:8],col=i)
  text(5,assinaturasPorCulturaD3[i,7],assinaturasPorCulturaD3[i,1],col=i)
}

#Destacar pastagem permanente
lines(BANDAS.LANDSAT,assinaturasPorCulturaD3[1,3:8],col='red',lwd=5)


##########################################################################
#NDVI's por cultura
NDVIPorCulturaD4<-dados[,list(area=sum(area),
                              NDVI_1=mean(NDVI_1),
                              NDVI_2=mean(NDVI_2),
                              NDVI_3=mean(NDVI_3),
                              NDVI_4=mean(NDVI_4),
                              NDVI_5=mean(NDVI_5),
                              NDVI_6=mean(NDVI_6),
                              NDVI_7=mean(NDVI_7),
                              NDVI_8=mean(NDVI_8)),by=cultura]
NDVIPorCulturaD4<-as.data.frame(NDVIPorCulturaD4)

datasGrafico<-datas
substr(datasGrafico,0,1)<-"0"
datasGrafico<-as.numeric(datasGrafico)

plot(datasGrafico,NDVIPorCulturaD4[1,3:10],type='l',ylim=c(0,1))
text(5,NDVIPorCulturaD4[1,7],NDVIPorCulturaD4[1,1])
for(i in 1:nrow(NDVIPorCulturaD4))
{
  lines(datasGrafico,NDVIPorCulturaD4[i,3:10],col=i)
  text(5,NDVIPorCulturaD4[i,7],NDVIPorCulturaD4[i,1],col=i)
}

lines(datasGrafico,NDVIPorCulturaD4[1,3:10],col='red',lwd=5)


##########################################################################
#Historgramas/densidades

#Histogramas
h<-hist(listaTodasParcelas$B4_d1,breaks=50,prob=T)
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

##########################################################################
#Vizualizacao 3D dos dados nas 3 variaveis mais importantes

#Seleccao das 3 var's mais importantes
trim.matrix(cor(dadosTreino[,-1]),0.5)

#Graficos 3D
library(rgl)
cultGrafico<-c(1,3,9,13)
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



##########################################################################################
# 4. SELECCAO DE VARIAVEIS                                                               #
##########################################################################################


##########################################################################
#Seleccao das variaveis
#Matriz de variancias-covariancias das assinaturas espectrais
corAssinaturas<-cor(listaTodasParcelas[,c(7:12,15:20,23:28,31:36,39:44,47:52)])

#Usando apenas o trim.matrix do package subselect
#Questao: Como escolher valor de tolval? Experimentar classificacoes com varios?
criterioExclusaoVars<-0.02
classTrim<-trim.matrix(cor(dadosTreino[,-1]),criterioExclusaoVars);classTrim

varsRetirar<-classTrim$numbers.discarded+1  #+1 porque na data.frame original, a coluna 1 ? a cultura
varsSobram<-1:length(dadosTreino)
varsSobram<-varsSobram[! varsSobram %in% varsRetirar]



#Nao vale a pena usar isto...
classHmat<-ldaHmat(dadosTreino[,-1], dadosTreino$cultura)
classSubselect<-eleaps(classHmat$mat,kmin=10,kmax=10,H=classHmat$H,r=classHmat$r,crit="Tau2",timelimit=600)
classSubselect


##########################################################################
#Preparacao dos NOVOS DADOS APOS SELECCAO DE VARIAVEIS
dadosClassificadoresSub<-dadosClassificadores[,varsSobram]

#Amostra a usar para o modelo - metade para treino e metade para validacao
set.seed(23)
amostra<-sample(1:nrow(dadosClassificadoresSub), ceiling(nrow(dadosClassificadoresSub)/2))

#Seleccao dos dados
dadosTreinoSub<-dadosClassificadoresSub[amostra,]
dadosValidacaoSub<-dadosClassificadoresSub[-amostra,]



##########################################################################################
# 5. 6. 7. ESCOLHA/TREINO/VALIDACAO DOS CLASSIFICADORES                                  #
##########################################################################################

#QUESTOES:
# - Como vamos poder validar isto sem as declaracoes? 


##########################################################################
#1. KNN - K VIZINHOS MAIS PROXIMOS

###################
#TODOS OS DADOS
classCult.knn<-knn(dadosTreino[,-1], dadosValidacao[,-1], dadosTreino[,1], prob=TRUE, k=10)
classCult.knn.tune<-tune.knn(dadosTreino[,-1], dadosTreino[,1], k=1:20);classCult.knn.tune

#Validacao
resultadoCult.knn<-table(pred=classCult.knn, true=dadosValidacao[,1]);resultadoCult.knn
resultadoCult.knn<-as.data.frame.matrix(resultadoCult.knn)
classDiag.knn<-diag(as.matrix(resultadoCult.knn));classDiag.knn
certos.knn<-sum(classDiag.knn)/sum(resultadoCult.knn);certos.knn

#Probabilidades
probs.knn<-unlist(attributes(classCult.knn),use.names=FALSE)
probs.knn<-as.numeric(probs.knn[-(1:(nClasses+1))])
probs.round.knn<-round(probs.knn,1)
hist(probs.knn,col='blue')
table(pred=classCult.knn, true=dadosValidacao[,1], probs=probs.round.knn)


###################
#APENAS SUB-CONJUNTO DE DADOS
classCultSub.knn<-knn(dadosTreinoSub[,-1], dadosValidacaoSub[,-1], dadosTreinoSub[,1], prob=TRUE, k=10)
classCultSub.tune.knn<-tune.knn(dadosTreinoSub[,-1], dadosTreinoSub[,1], k=1:20);classCult.knn.sub.tune

#Validacao
resultadoCultSub.knn<-table(pred=classCultSub.knn, true=dadosValidacaoSub[,1]);resultadoCultSub.knn
resultadoCultSub.knn<-as.data.frame.matrix(resultadoCultSub.knn)
classDiagSub.knn<-diag(as.matrix(resultadoCult.knn.sub));classDiagSub.knn
certosSub.knn<-sum(classDiagSub.knn)/sum(resultadoCultSub.knn);certosSub.knn

#Probabilidades
probsSub.knn<-unlist(attributes(classCultSub.knn),use.names=FALSE)
probsSub.knn<-as.numeric(probsSub.knn[-(1:(nClasses+1))])
probsSub.round.knn<-round(probsSub.knn,1)
hist(probsSub.knn,col='green')
table(pred=classCultSub.knn, true=dadosValidacaoSub[,1], probs=probsSub.round.knn)

#So com P >= 60%
jjj<-table(pred=classCultSub.knn, true=dadosValidacaoSub[,1], probs=probsSub.round.knn)
jjj<-as.matrix(jjj[,,5]+jjj[,,6]+jjj[,,7]+jjj[,,8]+jjj[,,9])
sum(diag(jjj))/sum(jjj)


###################
#ANALISE GRAFICA
#Percentagem de classificacoes correctas em cada cultura
percCertos<-c()
for(i in 1:ncol(resultadoCult.knn))
  percCertos[i]<-resultadoCult.knn[i,i]/sum(resultadoCult.knn[,i])
percCertos<-100*percCertos

percCertosSub<-c()
for(i in 1:ncol(resultadoCult.knn))
  percCertosSub[i]<-resultadoCultSub.knn[i,i]/sum(resultadoCultSub.knn[,i])
percCertosSub<-100*percCertosSub

par(las=2, mar=c(8,4.5,1,1))
barplot(rbind(percCertos,percCertosSub), beside=T, space=c(0,2.5), border=NA, names.arg=codigosNovos2005$codNomesAbrev[1:14], ylab='Percentagem de classificacoes correctas (%)', ylim=c(0,100), col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()
legend(x=30, y=98, legend=c("Completo","Sub-conjunto"), fill=c(rgb(0,0,0),rgb(0.7,0.7,0.7)))


#Percentagem de classificacoes correctas em cada valor de probabilidade
classProb<-table(pred=classCult.knn, true=dadosValidacao[,1], probs=probs.round.knn)
percCertos<-c()
for (i in 1:length(classProb[1,1,]))
{
  diag<-diag(as.matrix(classProb[,,i]))
  percCertos[i]<-sum(diag)/sum(classProb[,,i])
}
percCertos<-100*percCertos

classProbSub<-table(pred=classCultSub.knn, true=dadosValidacaoSub[,1], probs=probsSub.round.knn)
percCertoSubs<-c()
for (i in 1:length(classProbSub[1,1,]))
{
  diag<-diag(as.matrix(classProbSub[,,i]))
  percCertosSub[i]<-sum(diag)/sum(classProbSub[,,i])
}
percCertosSub<-100*percCertosSub

par(las=1, mar=c(4.5,4.5,1,1))
barplot(rbind(percCertos,percCertosSub), beside=T, space=c(0,2.5), border=NA, ylab='Percentagem de classificacoes correctas (%)', xlab='Probabilidade associada a classificao (%)', ylim=c(0,100), names.arg=seq(20,100,10),  col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()
legend(x=2, y=98, legend=c("Completo","Sub-conjunto"), fill=c(rgb(0,0,0),rgb(0.7,0.7,0.7)))

plot(seq(0.2,1,0.1),percCertos/100)

#Frequ?ncia da ocorrencia de cada valor de probabilidade
compProbs.knn<-rbind(table(probs.round.knn),table(probsSub.round.knn))
barplot(compProbs.knn, beside=T, space=c(0,4), legend=c("Completo","Sub-conjunto"), border=NA, col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()





##########################################################################
#2. SVM - MAQUINAS DE VECTORES DE SUPORTE

###################
#TODOS OS DADOS
classCult.svm<-svm(cultura ~ ., data=dadosTreino, probability=TRUE, kernel='radial', gamma=0.1, cost=3)
plot(classCult.svm, dadosTreino, B3_d3 ~ B4_d3)

#Tuning
classCult.tune.svm<-tune.svm(dadosTreino[,-1], dadosTreino[,1], gamma=seq(0, 1, by = .1));classCult.tune.svm
classCult.tune.svm<-tune.svm(dadosTreino[,-1], dadosTreino[,1], gamma=0.1, cost=seq(1, 10, by = 1));classCult.tune.svm
classCult.tune.svm<-tune.svm(dadosTreino[,-1], dadosTreino[,1], gamma=seq(0, 0.2, by=0.1), cost=seq(2, 4, by=1));classCult.tune.svm


#Validacao
validacao.svm<-predict(classCult.svm, dadosValidacao[,-1])
resultadoCult.svm<-table(pred = validacao.svm, true = dadosValidacao[,1]);resultadoCult.svm
resultadoCult.svm<-as.data.frame.matrix(resultadoCult.svm)
classDiag.svm<-diag(as.matrix(resultadoCult.svm));classDiag.svm
certos.svm<-sum(classDiag.svm)/sum(resultadoCult.svm);certos.svm

#Probabilidades
probs.svm<-predict(classCult.svm, dadosValidacao[,-1], probability=T)
probs.svm<-attr(probs.svm, "prob")
probs.svm<-as.vector(apply(probs.svm,1,max))
probs.round.svm<-round(probs.svm,1)
hist(probs.round.svm, breaks=14, col='blue')
table(pred=validacao.svm, true=dadosValidacao[,1], probs=probs.round.svm)


###################
#APENAS SUB-CONJUNTO DE DADOS
classCultSub.svm<-svm(cultura ~ ., data=dadosTreinoSub, probability=TRUE, kernel='radial', gamma=0.2, cost=2)
plot(classCultSub.svm, dadosTreinoSub, B3_d3 ~ B4_d3)

#Tuning
classCultSub.tune.svm<-tune.svm(dadosTreinoSub[,-1], dadosTreinoSub[,1], gamma=seq(0, 0.5, by = .1));classCultSub.tune.svm
classCultSub.tune.svm<-tune.svm(dadosTreinoSub[,-1], dadosTreinoSub[,1], gamma=0.2, cost=seq(1, 8, by = 1));classCultSub.tune.svm

#Validacao
validacaoSub.svm<-predict(classCultSub.svm, dadosValidacaoSub[,-1])
resultadoCultSub.svm<-table(pred = validacaoSub.svm, true = dadosValidacaoSub[,1]);resultadoCultSub.svm
resultadoCultSub.svm<-as.data.frame.matrix(resultadoCultSub.svm)
classDiagSub.svm<-diag(as.matrix(resultadoCultSub.svm));classDiagSub.svm
certosSub.svm<-sum(classDiagSub.svm)/sum(resultadoCultSub.svm);certosSub.svm

#Probabilidades
probsSub.svm<-predict(classCultSub.svm, dadosValidacaoSub[,-1], probability=T)
probsSub.svm<-attr(probsSub.svm, "prob")
probsSub.svm<-as.vector(apply(probsSub.svm,1,max))
probsSub.round.svm<-round(probsSub.svm,1)
hist(probsSub.round.svm, breaks=14, col='green')
table(pred=validacaoSub.svm, true=dadosValidacaoSub[,1], probs=probsSub.round.svm)

#So com P >= 60%
jjj<-table(pred=validacaoSub.svm, true=dadosValidacaoSub[,1], probs=probsSub.round.svm)
jjj<-as.matrix(jjj[,,5]+jjj[,,6]+jjj[,,7]+jjj[,,8]+jjj[,,9])
sum(diag(jjj))/sum(jjj)


###################
#Comparacao das probabilidades
compProbs.svm<-rbind(table(probs.round.svm),table(probsSub.round.svm))
barplot(compProbs.svm, beside=T, space=c(0,5), legend=c("Completo","Sub-conjunto"), xlab="Probabilidade", ylab="Numero de parcelas", border=NA, col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()
#Conclusao: com dados nao-completos, existem mais classificacaes com 100% de certeza







##########################################################################
#Abordagem diferente - Estimacao das densidades de probabilidade de cada classe de culturas

#QUESTOES:
# - Como vamos poder validar isto?
# - Como fazer uma verdadeira 'funcao' com isto para podermos fazer testes de hipoteses?

classe<-1
k<-10      #Numero de vizinhos a considerar - resultado do KNN.TUNE
d<-ncol(dadosTreinoSub[,-1])      #Numero de dimensoes
dadosTreinoClasse<-dadosTreinoSub[dadosTreinoSub$cultura == classe,]
dadosValidacaoClasse<-dadosValidacaoSub[dadosValidacaoSub$cultura == classe,]
n<-nrow(dadosTreinoClasse)     #Numero de observacoes na amostra de treino???

probs<-c()
for (i in 1:nrow(dadosValidacaoClasse))
{
  observacao<-dadosValidacaoClasse[i,-1]
  distancias<-sort(as.matrix(dist(rbind(observacao,dadosTreinoClasse[,-1])))[-1,1])
  distKVizinho<-distancias[k]
  
  #Calculo da probabilidade
  if(d%%2 == 0)
    cD<-pi^(d/2)/factorial(d/2)
  else
    cD<-pi^(d/2)/gamma(d/2+1)
  
  p<-k/(n*distKVizinho*cD)
  probs[i]<-p
}

plot(sort(probs))
hist(probs,breaks=20)


#########
#Experiencia com 2 dimensoes e pontos aleatorios
x<-sample(seq(1,20,0.1),10,replace=T)
y<-sample(seq(1,20,0.1),10,replace=T)
novo.x<-9
novo.y<-7
plot(x,y)
text(x,y+0.4,2:11)
points(novo.x,novo.y,col='red')
text(novo.x,novo.y+0.4,1,col='red')

distancias<-sort(as.matrix(dist(rbind(cbind(novo.x,novo.y),cbind(x,y))))[-1,1])
distKVizinho<-distancias[5]

draw.circle(novo.x,novo.y,distKVizinho,border='red')


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################







##########################################################################################
#   CRITERIO 1 --- UTILIZA????O DE DIFERENTES CLASSIFICADORES - TREINO E VALIDA????O         #
##########################################################################################

#Classificadores utilizados:
# 1. LDA - An??lise discriminante linear
# 2. KNN - k vizinhos mais pr??ximos
# 3. SVM - M??quinas de vectores de suporte
# 4. ??rvores de decis??o


##########################################################################
#0. PREPARA????O E VISUALIZA????O DOS DADOS

#Cria????o dos dados com as seguintes caracteristicas:
#  - pastagem ?? factor que ?? 0 se nao for pastagem, 1 se for pastagem permanente e 2 se for pastagem pobre
#  - usam-se dados de todas as bandas nas primeiras 6 datas do ano
#  - metade dos dados vai ser para treino dos classificadores e metade para valida????o

pastagem<-c()
seq<-c()
dadosClassificadores<-listaTodasParcelas[,c(5,7:12,15:20,23:28,31:36,39:44,47:52)]

for(i in 1:length(dadosClassificadores$cultura))
{
  if(dadosClassificadores$cultura[i] == 1)
    pastagem[i]<-1
  else
    pastagem[i]<-0
}
dadosClassificadoresPast<-NA
dadosClassificadoresPast<-cbind(pastagem,dadosClassificadores)
dadosClassificadoresPast<-dadosClassificadoresPast[,-2]
dadosClassificadoresPast$pastagem<-as.factor(as.numeric(as.vector(dadosClassificadoresPast$pastagem)))

#Amostra a usar para o modelo - metade para treino e metade para valida????o
set.seed(23)
amostra<-sample(1:nrow(dadosClassificadoresPast), ceiling(nrow(dadosClassificadoresPast)/2))

#Seleccao dos dados
dadosTreino<-dadosClassificadoresPast[amostra,]
dadosValidacao<-dadosClassificadoresPast[-amostra,]

#Visualiza????o dos dados 2D
plot(c(),xlim=range(dadosTreino$B4_d3),ylim=range(dadosTreino$B4_d5),xlab="NDVI DOY 125",ylab="NDVI DOY 189")
text(dadosTreino$B4_d3[dadosTreino$pastagem == 0],dadosTreino$B4_d5[dadosTreino$pastagem == 0],dadosTreino$pastagem[dadosTreino$pastagem == 0])
text(dadosTreino$B4_d3[dadosTreino$pastagem == 1],dadosTreino$B4_d5[dadosTreino$pastagem == 1],dadosTreino$pastagem[dadosTreino$pastagem == 1], col="red")

#Visualiza????o dos dados 3D
library(rgl)
plot3d(x=dadosTreino$B4_d2[dadosTreino$pastagem == 0],
       y=dadosTreino$B4_d3[dadosTreino$pastagem == 0],
       z=dadosTreino$B4_d5[dadosTreino$pastagem == 0],
       xlab="B4 DOY 61",
       ylab="B4 DOY 125",
       zlab="B4 DOY 189",
       size=3)
points3d(x=dadosTreino$B4_d2[dadosTreino$pastagem == 1],
         y=dadosTreino$B4_d3[dadosTreino$pastagem == 1],
         z=dadosTreino$B4_d5[dadosTreino$pastagem == 1],
         col="red")

#Resultados da optimiza????o da LDA
plot3d(x=dadosTreino$B3_d5[dadosTreino$pastagem == 0],
       y=dadosTreino$B5_d5[dadosTreino$pastagem == 0],
       z=dadosTreino$B6_d5[dadosTreino$pastagem == 0],
       xlab="B3_d5",
       ylab="B5_d5",
       zlab="B6_d5",
       size=3)
points3d(x=dadosTreino$B3_d5[dadosTreino$pastagem == 1],
         y=dadosTreino$B5_d5[dadosTreino$pastagem == 1],
         z=dadosTreino$B6_d5[dadosTreino$pastagem == 1],
         col="red")

#Matriz de gr??ficos
pairs(~B4_d2+B4_d3+B4_d5, data=dadosTreino)


##########################################################################
#1. LDA - ANALISE DISCRIMINANTE LINEAR 

#Treino do classificador
class.lda<-lda(dadosTreino[,-1],grouping=dadosTreino$pastagem);class.lda
plot(class.lda,col=(as.numeric(dadosTreino[,1])))

#Valida????o
validacao.lda<-predict(class.lda, dadosValidacao[,-1])
resultado.lda<-table(pred = validacao.lda$class, true = dadosValidacao[,1]);resultado.lda
classAgreement(resultado.lda)



#Optimiza????o com package subselect
class.lda.Hmat<-ldaHmat(dadosTreino[,-1], dadosTreino$pastagem)
class.lda.opt<-anneal(class.lda.Hmat$mat,kmin=3,kmax=5,H=class.lda.Hmat$H,crit="RM")
class.lda.opt$subsets

#nomes das 3 variaveis mais importantes
colnames(dadosTreino[,class.lda.opt$subsets])


##########################################################################
#2. SVM - M??QUINAS DE VECTORES DE SUPORTE

#Treino do classificador
class.svm<-svm(pastagem ~ ., data=dadosTreino, cost=100, gamma=0.9)
plot(class.svm, dadosTreino, B1_d3 ~ B5_d6)

#Optimiza????o - tentar com MAIS hipoteses
class.svm.tune<-tune.svm(pastagem ~ ., data=dadosTreino, gamma=seq(.1, 1, by = .1), cost = seq(100,1000, by = 100));class.svm.tune

#resultados: gamma=0.9 e cost=100

#Valida????o
validacao.svm<-predict(class.svm, dadosValidacao[,-1])
resultado.svm<-table(pred = validacao.svm, true = dadosValidacao[,1]);resultado.svm
classAgreement(resultado.svm)


##########################################################################
#3. KNN - K VIZINHOS MAIS PR??XIMOS

#Treino e valida????o simultanea
class.knn<-knn(dadosTreino[,-1], dadosValidacao[,-1], dadosTreino[,1], k=16)
resultado.knn<-table(pred=class.knn, true=dadosValidacao[,1]);resultado.knn
classAgreement(resultado.knn)

#Optimiza????o
class.knn.tune<-tune.knn(dadosTreino[,-1], dadosTreino[,1], k = 1:20);class.knn.tune


##########################################################################
#4. ??RVORES DE DECIS??O

#Ajuste inicial da ??rvore
class.tree<-tree(pastagem ~ ., data=dadosTreino, control= tree.control(nobs=nrow(dadosTreino), minsize=1))
plot(class.tree);text(class.tree,cex=.8)

#Valida????o
validacao.tree<-predict(class.tree,type="class")
resultado.tree<-table(pred = validacao.tree, true = dadosTreino[,1]);resultado.svm
classAgreement(resultado.tree)

#Determina????o da melhor ??rvore (com menor estimativa de erro de classifica????o)
cvDev<-cv.tree(class.tree)$dev
opt<-which(cvDev == min(cvDev))
class.tree.prune<-prune.tree(class.tree,best=opt)
plot(class.tree.prune);text(class.tree.prune,cex=.8)


#####
#OUTRA ABORDAGEM COM OUTRO PACKAGE QUE PARECE SER MELHOR
#Treino do classificador
class.tree2<-ctree(pastagem ~ ., data=dadosTreino)
plot(class.tree2)

#Valida????o
validacao.tree2<-predict(class.tree2, dadosValidacao[,-1])
resultado.tree2<-table(pred = validacao.tree2, true = dadosValidacao[,1]);resultado.tree2
classAgreement(resultado.tree2)

#####
#mais uma abordagem
class.tree3<-rpart(pastagem ~ ., data=dadosTreino)
plot(class.tree3)
text(class.tree3,cex=.8)

#Valida????o
validacao.tree3<-predict(class.tree3, dadosValidacao[,-1],type="class")
resultado.tree3<-table(pred = validacao.tree3, true = dadosValidacao[,1]);resultado.tree3
classAgreement(resultado.tree3)

##########################################################################
#5. REGRESS??O LOG??STICA MULTINOMIAL

#Com package NNET
library(nnet)
class.multi<-multinom(pastagem ~ ., data=dadosTreino)
summary(class.multi)

validacao.multi<-predict(class.multi, dadosValidacao[,-1])
resultado.multi<-table(pred = validacao.multi, true = dadosValidacao[,1]);resultado.multi


##########################################################################################
#   CRITERIO 2                                                                           #
##########################################################################################








##########################################################################################
#   EXPERIENCIAS                                                                         #
##########################################################################################



##########################################################################
#Experiencias com NDVI's (apenas para ano 2005)

#Fazer plot
datasGrafico<-datas
substr(datasGrafico,0,1)<-"0"
datasGrafico<-as.numeric(datasGrafico)

#plot(datasGrafico,NDVIs[1,],ylim=c(-1,1),col=100,type='l',lty=1,pch=16,ylab="NDVI", xlab="DOY 2005")
plot(datasGrafico,listaTodasParcelas[1,55:62],ylim=c(-1,1),col=1,type='l',lty=1,ylab="NDVI", xlab="DOY 2005")
for(i in 1:nrow(listaTodasParcelas))
  lines(datasGrafico,listaTodasParcelas[i,55:62],ylim=c(-1,1),col=i,lty=1,pch=16)


#distinguir entre sequeiro e regadio
todasParcelasSEQ<-listaTodasParcelas[listaTodasParcelas$seq == TRUE,]
todasParcelasREG<-listaTodasParcelas[listaTodasParcelas$seq == FALSE,]

plot(datasGrafico,todasParcelasSEQ[1,55:62],ylim=c(-1,1),col='burlywood2',type='l',lty=1,ylab="NDVI", xlab="DOY 2005")
for(i in 1:100)
  lines(datasGrafico,todasParcelasSEQ[i,55:62],ylim=c(-1,1),col='burlywood2',lty=1,pch=16)
for(i in 1:100)
  lines(datasGrafico,todasParcelasREG[i,55:62],ylim=c(-1,1),col='darkolivegreen3',lty=1,pch=16)

##########################################################################






##########################################################################
#Experi??ncias com outros indices de vegeta??ao

listaTodasParcelas$NDVI_1[50]
(listaTodasParcelas$B4_d1[50]-listaTodasParcelas$B3_d1[50])/(listaTodasParcelas$B4_d1[50]+listaTodasParcelas$B3_d1[50])

nir<-listaTodasParcelas$B4_d4
red<-listaTodasParcelas$B3_d4
NDVId4<-(nir-red)/(nir+red)
NDVId4.orig<-listaTodasParcelas$NDVI_4

plot(NDVId4,NDVId4.orig)
abline(lm(NDVId4.orig~NDVId4),col='red')


MSAVId4<-(2*nir+1-sqrt((2*nir+1)^2-8*(nir-red)))/2
plot(NDVId4,MSAVId4)
abline(lm(MSAVId4~NDVId4),col='red')


##########################################################################






##########################################################################
#Experimentar diferen??a entre usar apenas os interiores das parcelas e usar parcelas inteiras
testeTotInt<-matrix(nrow=0,ncol=5)
for(i in 1:length(listaDados))
{
  for(j in 1:length(listaDados[[i]][[1]]))
  {
    parc<-listaDados[[i]][[1]][[j]]
    conteudoTotal<-obtemConteudoParcela(parc,apenasInterior=FALSE)
    conteudoInterior<-obtemConteudoParcela(parc,apenasInterior=TRUE)
    nTotal<-length(conteudoTotal[[1]][,1])
    nInterior<-length(conteudoInterior[[1]][,1])
    sdTotal<-sd(conteudoTotal[[1]][,1])
    sdInterior<-sd(conteudoInterior[[1]][,1])
    area<-parc$atributos$AREA
    linha<-c(nTotal,sdTotal,nInterior,sdInterior,area)
    testeTotInt<-rbind(testeTotInt,linha)
  } 
  print(i)
}

mean(testeTotInt[,1])   #M??dia do numero de pixeis quando se considera tambem a fronteira
mean(testeTotInt[,3])   #M??dia do numero de pixeis quando se considera apenas o interior

mean(testeTotInt[,1])/mean(testeTotInt[,5]/10000)   #M??dia de pixeis/ha -- fronteira
mean(testeTotInt[,3])/mean(testeTotInt[,5]/10000)   #M??dia de pixeis/ha -- interior

mean(testeTotInt[,2],na.rm=T)   #Media dos desvios padrao -- fronteira
mean(testeTotInt[,4],na.rm=T)   #Media dos desvios padrao -- interior

length(testeTotInt[testeTotInt[,1] == 0,1])/length(testeTotInt[,1])  #% de parcelas sem pixeis -- fronteira
length(testeTotInt[testeTotInt[,3] == 0,3])/length(testeTotInt[,3])  #% de parcelas sem pixeis -- interior

table(as.vector(testeTotInt[,2] > testeTotInt[,4]))   #Numero de vezes que sd ?? maior utilizando a fronteira
##########################################################################




##########################################################################################
#   EXPORTAR GRAFICOS E QUADROS PARA LATEX                                               #
##########################################################################################


#Graficos
opar<-par()
par(opar)

#Graficos
library(tikzDevice)
tikz('../05 Escrito/TEX/test.tex',width=5,height=4)
par(mar=c(3.5,3.5,0.5,0.5), mgp=c(2,0.6,0))   #bottom, left, top, right
barplot(compProbs.svm, beside=T, space=c(0,5), legend=c("Completo","Sub-conjunto"), xlab="Probabilidade", ylab="Numero de parcelas", border=NA, col=c(rgb(0,0,0),rgb(0.7,0.7,0.7)));box()
dev.off()


#Quadros
library(xtable)
dados.xtab<-xtable(dadosTreino[1:10,1:5])
align(dados.xtab)<-rep("c",6)
digits(dados.xtab)<-3
print.xtable(dados.xtab, booktabs=T, include.rownames=F)



##########################################################################################
#   EXPORTAR E IMPORTAR GRANDES OBJECTOS                                                 #
##########################################################################################

#Gravar objectos grandes como ficheiros
save(cdg, file="cdg2005.Robj")
save(listaDados, file="listaDados2005.Robj")
save(todasImagens, file="todasImagens2005.Robj")
save(listaTodasParcelasIniciais, file="listaTodasParcelasIniciais2005.Robj")

rm(listaDados)
rm(cdg)
rm(todasImagens)

#Carregar objectos grandes
load("listaDados2005.Robj")
load("cdg2005.Robj")
load("todasImagens2005.Robj")


##########################################################################################
#   CAIXOTE DO LIXO                                                                      #
##########################################################################################


#Bom livro em PT
#http://cran.r-project.org/doc/contrib/Torgo-ProgrammingIntro.pdf
#http://www.statmethods.net/advgraphs/parameters.html


#Obter conteudo de uma parcela para diferentes datas e bandas
parc<-listaDados$n1541$a2005$p1432302680003.301.00
conteudoParcIn<-obtemConteudoParcela(parc)
conteudoParcInSemViz<-obtemConteudoParcela(parc,excluiVizinhos=TRUE)
conteudoParc<-obtemConteudoParcela(parc,apenasInterior=FALSE)

plot(parc$DR$d45[,,7], parc$DR$d45[,,8])
lines(parc$coordsPoly, col="red")


#Package kknn
library(kknn)
classTodCult.kknn<-kknn(dadosTreino2[,-1], dadosValidacao2[,-1], dadosTreino2[,1], k=10)




#Areas das parcelas
hist(listaTodasParcelas$area,breaks=20)
boxplot(listaTodasParcelas$area)
plot(sort(listaTodasParcelas$area))
summary(listaTodasParcelas$area)

which(listaTodasParcelas$area==max(listaTodasParcelas$area)) #Maior parcela de todas
listaTodasParcelas$cultura[which(listaTodasParcelas$area==max(listaTodasParcelas$area))]  #Ocupa????o da maior parcela de todas


par(col="black")
par(opar)

plot(listaDados$n197$a2003$p1162347095001.301.00$DR$d29$X, listaDados$n197$a2003$p1162347095001.301.00$DR$d29$Y)
lines(listaDados$n197$a2003$p1162347095001.301.00$coordsPoly, col="red")

#Assinaturas espectrais ao longo do ano
frequencias<-c(0.485, 0.56, 0.66, 0.835, 1.65, 2.22)

todasMedias<-matrix(nrow=length(conteudoParc), ncol=length(conteudoParc[[1]]))
for(i in 1:length(conteudoParc))
{
  assinatura<-c()
  for(j in 1:length(conteudoParc[[i]]))
    assinatura[j]<-conteudoParc[[i]][[j]][["media"]]
  todasMedias[i,]<-assinatura
}

for(i in 1:nrow(todasMedias))
{
  if(i==1)
    plot(frequencias, todasMedias[i,], type="o", lty=1, col=i, ylim=c(min(todasMedias),max(todasMedias)))
  else
    lines(frequencias, todasMedias[i,], type="o", lty=1, col=i)
}

#Aceder a dados geograficos importados
cdg$parc2003@polygons[[1]]@Polygons[[1]]@area

cdg$area@polygons[[1]]@Polygons[[1]]@coords

plot(cdg$parc2003, col="red")
plot(cdg$conc, col=rgb(1,0,0,alpha=0), add=TRUE)

proc.time()
system.time()

##########################################################################
#Usar isto para percorrer todas os NIFAPs
for(i in 1:length(listaDados))
  for(j in 1:length(listaDados[[i]][[1]]))
  {
    parc<-listaDados[[i]][[1]][[j]]
    .
    .
    .
  } 
##########################################################################
