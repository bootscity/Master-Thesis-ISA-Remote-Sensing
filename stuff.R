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


##########################################################################################
#   TESTES                                                                               #
##########################################################################################

###############################
#USING APPLY
j<-c(3,4,5,-2,4,-3)
j[j<0]

m<-matrix(data=cbind(rnorm(20, 0), rnorm(20, 2), rnorm(20, 5)), nrow=20, ncol=3)
apply(m, 1, mean) #applying across rows of matrix
apply(m, 2, mean) #applying across cols of matrix
apply(m, 2, function(x) length(x[x<0]))
apply(m, 2, is.vector)
apply(m, 2, function(j) mean(j[j>0]))

#USING LAPPLY
sapply(1:3, function(x) x^2)  #Returns vector
lapply(1:3, function(x) x^2)  #Returns lsit

########################
#more stuff
colnames(m)<-c('method1', 'method2', 'method3')
head(m)
m[,'method1']
df<-as.data.frame(m);df
df.stack<-stack(df);df.stack
unstack(df.stack)


x<-1:10
names(x)<-letters[1:10]
x['f']
x[x > 5]<-0;x

y<-seq(1, 30, by=3);y
which(y %% 2 == 0)




sales <- expand.grid(country = c('USA', 'UK', 'FR'), product = c(1, 2, 3))
sales$revenue <- rnorm(dim(sales)[1], mean=100, sd=10)
usd2eur <- 1.434
transform(sales, euro = revenue * usd2eur)
subset(sales, 
       country == 'USA' & product %in% c(1, 2), 
       select = c('product', 'revenue'))

#DO.CALL
foo <- list()
foo[[1]] <- data.frame(a=1:5, b=11:15)
foo[[2]] <- data.frame(a=101:105, b=111:115)
foo[[3]] <- data.frame(a=200:210, b=300:310)
do.call(rbind, foo)



#########################
#using apply for lists
mylist<-list()
for(i in 1:5) {
  mylist[[i]]<-list()
  for(j in 1:5) {
    mylist[[i]][[j]]<-i*j
  }
}

sapply(mylist[[1]], identity)
sapply(mylist, function(x) lapply(x, identity))
str(mylist)

mylist<-list(1:3, "dog", TRUE)
sapply(mylist, identity)
unlist(sapply(mylist, identity))