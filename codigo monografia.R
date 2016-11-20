# Codigo de la monografia
library(gamlss)
library(moments)
data=read.csv(file.choose(),header = TRUE,sep = ",")
attach(data)
str(data)

x=perdida/1000000  #la variable perdida la estoy redefiniendo y volviendo en Mlls
summary(x)    #resumen estadistico express

medidas<-function(x) {
Media<-round(mean(x),3)
Min<-round(min(x),3)
Máx<-round(max(x),3)
Desviacion<-round(sd(x),3)
P25<-round(quantile(x,0.25),3)
P50<-round(quantile(x,0.5),3)
P75<-round(quantile(x,0.75),3)
P95<-round(quantile(x,0.95),3)
P99<-round(quantile(x,0.99),3)

return(data.frame(Media,Min,Máx,Desviacion,P25,P50,P75,P95,P99))

}

medidas(x)


medidas2<-function(x) {
  
skwP95<-skewness(x[x<quantile(x,0.95)],na.rm = FALSE)
kurtP95<-kurtosis(x[x<quantile(x,0.95)],na.rm = FALSE)
skwP99<-skewness(x[x<quantile(x,0.99)],na.rm = FALSE)
kurtP99<-kurtosis(x[x<quantile(x,0.99)],na.rm = FALSE)
skwP999<-skewness(x[x<quantile(x,0.999)],na.rm = FALSE)
kurtP999<-kurtosis(x[x<quantile(x,0.999)],na.rm = FALSE)

return(data.frame(skwP95,kurtP95,skwP99,kurtP99,skwP999,kurtP999))
}

medidas2(x)


#Gráfico de perdidas totales a diferentes percentiles

win.graph()
par(mfrow = c(1, 3))
boxplot(x[x<quantile(x,0.95)],main="Boxplot Pérdidas al P95",
        ylab = "Montos Variable pérdida")
boxplot(x[x<quantile(x,0.99)],main="Boxplot Pérdidas al P99",
        ylab = "Montos Variable pérdida")
boxplot(x[x<quantile(x,0.999)],main="Boxplot Pérdidas al P99.9",
        ylab = "Montos Variable pérdida")

win.graph()
qqnorm(x,main="Prueba de normalidad variable Pérdida",
       xlab="Cuantiles teóricos",ylab="Cuantiles muestrales")
qqline(x,col="blue",lwd=2)
shapiro.test(x)$p.value
legend("topleft",y=null,title="Valor P Shapiro Test",shapiro.test(x)$p.value)


#Definición perdidas a diferentes percentiles

y=x[x<quantile(x,0.95)]
z=x[x<quantile(x,0.99)]
a=x[x<quantile(x,0.999)]


medidas(y)
medidas(z)
medidas(a)

#Histograma de pérdidas a diferentes percentiles

win.graph()
par(mfrow = c(3, 1))
hist(y,freq = FALSE,breaks = 250,main="Histograma Pérdidas al P95",
     ylab = "Frecuencias",xlab = "Variable pérdida")
lines(density(y),na.rm=TRUE,col="blue",lwd=2)
hist(z,freq = FALSE,breaks = 250,main="Histograma Pérdidas al P99",
     ylab = "Frecuencias",xlab = "Variable pérdida")
lines(density(a),na.rm=TRUE,col="blue",lwd=2)
hist(a,freq = FALSE,breaks = 250,main="Histograma Pérdidas al P99.9",
     ylab = "Frecuencias",xlab = "Variable pérdida")
lines(density(a),na.rm=TRUE,col="blue",lwd=2)

#Ajustes de la variable pérdidas a diferentes percentiles
?fitDist

ajustey=fitDist(y,type = "realplus")
ajustey$fits
ajustey

ajustez=fitDist(z,type = "realplus")
ajustez$fits
ajustez

ajustea=fitDist(a,type = "realplus")
ajustea$fits
ajustea

ajustey$fits
ajustez$fits
ajustea$fits

#Tabla de Número de eventos en cada año

table(año)
mean(table(año))


#Simulación de las variables por sus respectivas fdps

yajust=rBCPEo(100000,1.073,0.1943,0.1343,0.9788) # box Cox power exponential orig mu es origi
#nalmente negativo, preguntar por qué
mean(yajust)
median(yajust)

zajust=rGG(100000,0.9611,0.292,0.1248) #Generalisaed Gamma " mu es negativo original/
mean(zajust)
median(zajust)

aajust=rBCPEo(100000,1.037,0.3319,0.002123,0.604) #box Cox power exponential orig
mean(aajust)
median(aajust)


#Gráfico de densidades de pérdidas simuladas.

win.graph()
par(mfrow = c(3, 1))
plot(density(yajust),main="Simulación Pérdidas al P95
     Distribución BCPEo",ylab="Densidad")
plot(density(zajust),freq = FALSE,main="Simulación Pérdidas al P99
     Distribución GG",ylab="Densidad")
plot(density(aajust),freq = FALSE,main="Simulación Pérdidas al P99.9
     Distribución BCPEo",ylab="Densidad")

#simulación poisson

freqprom<-rpois(100000,1080)
freqmax<-rpois(100000,1352)

#salidas finales

ensayo1e=quantile(yajust,0.5)*quantile(freqprom,0.5)
ensayo1=quantile(yajust,0.95)*quantile(freqprom,0.95)
ensayo1p=quantile(yajust,0.99)*quantile(freqprom,0.99)
s1=c(ensayo1e,ensayo1,ensayo1p)   #Escenario 1
s1
ensayo2e=quantile(zajust,0.5)*quantile(freqprom,0.5)
ensayo2=quantile(zajust,0.95)*quantile(freqprom,0.95)
ensayo2p=quantile(zajust,0.99)*quantile(freqprom,0.99)
s2=c(ensayo2e,ensayo2,ensayo2p)   #Escenario 2
s2
ensayo3e=quantile(aajust,0.5)*quantile(freqprom,0.5)
ensayo3=quantile(aajust,0.95)*quantile(freqprom,0.95)
ensayo3p=quantile(aajust,0.99)*quantile(freqprom,0.99)
s3=c(ensayo3e,ensayo3,ensayo3p)   #Escenario 3
s3

ensayo11e=quantile(yajust,0.5)*quantile(freqmax,0.5)
ensayo11=quantile(yajust,0.95)*quantile(freqmax,0.95)
ensayo11p=quantile(yajust,0.99)*quantile(freqmax,0.99)
s11=c(ensayo11e,ensayo11,ensayo11p)   #Escenario 4
s11
ensayo22e=quantile(zajust,0.5)*quantile(freqmax,0.55)
ensayo22=quantile(zajust,0.95)*quantile(freqmax,0.95)
ensayo22p=quantile(zajust,0.99)*quantile(freqmax,0.99)
s22=c(ensayo22e,ensayo22,ensayo22p)   #Escenario 5
s22
ensayo33e=quantile(aajust,0.5)*quantile(freqmax,0.5)
ensayo33=quantile(aajust,0.95)*quantile(freqmax,0.95)
ensayo33p=quantile(aajust,0.99)*quantile(freqmax,0.99)
s33=c(ensayo33e,ensayo33,ensayo33p)   #Escenario 6
s33


s1
s2
s3
s11
s22
s33

