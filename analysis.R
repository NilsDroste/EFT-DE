# Panel data analysis of nature conservation area and population density in Germany
# April 2017
# nils.droste@ufz.de
# research question: 
# is there an correlation of population density (pop/sqr(km)) and nature conservation area (in sqr(m) per capita) among German federal states?

setwd("C:/Users/droste/.../EFT-GER/data/")

# start -------------------------------------------------------------------

df <- read.csv("C:/Users/droste/Dropbox/Dokumente/doctorate/fiscal transfers/EFT-GER/data/panel.csv", sep=";", dec=",")
df<-reshape(df, direction='long', varying=list(3:6,7:10,11:14,15:18,19:22,23:26,27:30,31:34,35:38,39:42,43:46,47:50,51:54,55:58), v.names=c("pop.dens","natsch.per","landsch.per","ges.per","natsch.tot","landsch.tot","ges.tot","ausg","BIP.cap","BWSagr", "BWSind","BWSser","pop","area"), idvar="ID", times=c(2012, 2010, 2008, 2006))
df<- df[order(df$No, df$time),] 
df$nat.cap<-df$natsch.tot*1000000/df$pop
df$land.cap<-df$landsch.tot*1000000/df$pop
df$ges.cap<-df$ges.tot*1000000/df$pop
df$ausg.cap<-df$ausg*1000000/df$pop
df<-df[df$time %in% c(2006,2008,2010),]


#summary table
require(stargazer)
stargazer(df[df$time %in% c(2006,2008,2010),c("nat.cap","land.cap","ges.cap","pop.dens","ausg.cap")], type = "html", title="Descriptive statistics", digits=1, out="summary.html",covariate.labels=c("Nature and species conservation area per capita in mÂ² (nat.cap)","Landscape protection area per capita in mÂ² (land.cap)","Total protected area per capita in mÂ² (tot.cap)","Population density in kmÂ²/capita (pop.dens)", "public expenditure environmental protection and nature conservation in â¬ per capita (spend.cap)"), notes = "Sources: authorsâ calculations based on IOER (2015) and Statistisches Bundesamt (personal communciation), monetary values are in constant â¬ 2005 prices.")


# regressions -------------------------------------------------------------
require(plm)
require(lmtest)
require(sandwich)

#multicollinearity
source("C:/Users/droste/.../panel_cor.R")
pairs(~ log(nat.cap) + log(land.cap) + log(ges.cap) + log(pop.dens)+log(BIP.cap)+log(BWSagr)+log(BWSind)+log(BWSser)+log(ausg.cap)+as.integer(time)+ log(area), data=df, upper.panel=panel.cor, lower.panel=panel.smooth, pch=20)
#exclude BIPcap, BWser, area due to correlation coefficients about or greater than 0.7

#regressions

#nat
f1.nat<-plm(log(nat.cap) ~ log(pop.dens)+log(BIP.cap)+as.integer(time), data=df, index=c("ID", "time"), model="within", effect="individual")
summary(f1.nat)
f1.nat.rob<-coeftest(f1.nat, vcov=function(x) vcovSCC(x, type="HC3", maxlag=1)); f1.nat.rob

f2.nat<-plm(log(nat.cap) ~ log(pop.dens)+log(BIP.cap)+log(BWSagr)+log(BWSind)+log(ausg.cap)+as.integer(time), data=df, index=c("ID", "time"), model="within", effect="individual")
summary(f2.nat)
f2.nat.rob<-coeftest(f2.nat, vcov=function(x) vcovSCC(x, type="HC3", maxlag=1)); f2.nat.rob


#land
f1.land<-plm(log(land.cap) ~ log(pop.dens)+log(BIP.cap)+as.integer(time), data=df, index=c("ID", "time"), model="within", effect="individual")
summary(f1.land)
f1.land.rob<-coeftest(f1.land, vcov=function(x) vcovSCC(x, type="HC3", maxlag=1)); f1.land.rob

f2.land<-plm(log(land.cap) ~ log(pop.dens)+log(BIP.cap)+log(BWSagr)+log(BWSind)+log(ausg.cap)+as.integer(time), data=df, index=c("ID", "time"), model="within", effect="individual")
summary(f2.land)
f2.land.rob<-coeftest(f2.land, vcov=function(x) vcovSCC(x, type="HC3", maxlag=1)); f2.land.rob


#total
f1.tot<-plm(log(ges.cap) ~ log(pop.dens)+log(BIP.cap)+as.integer(time), data=df, index=c("ID", "time"), model="within", effect="individual")
summary(f1.tot)
f1.tot.rob<-coeftest(f1.tot, vcov=function(x) vcovSCC(x, type="HC3", maxlag=1)); f1.tot.rob

f2.tot<-plm(log(ges.cap) ~ log(pop.dens)+log(BIP.cap)+log(BWSagr)+log(BWSind)+log(ausg.cap)+as.integer(time), data=df, index=c("ID", "time"), model="within", effect="individual")
summary(f2.tot)
f2.tot.rob<-coeftest(f2.tot, vcov=function(x) vcovSCC(x, type="HC3", maxlag=1)); f2.tot.rob


# analytical plots --------------------------------------------------------
par(mfrow=c(2,2))
plot(density(resid(f1.nat)),main = "residual density")
fitted.f1.nat=(f1.nat$model[[1]] - f1.nat$residuals)
plot(resid(f1.nat)~(fitted.f1.nat), main = "Residuals vs Fitted", ylab="residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted.f1.nat, residuals(f1.nat)), col="red")
qqnorm(f1.nat$resid)
qqline(f1.nat$resid, col="red")
lev = hat(model.matrix(f1.nat))
plot(lev, main="leverage")
#df[lev>0.4,]

par(mfrow=c(2,2))
plot(density(resid(f2.nat)),main = "residual density")
fitted.f2.nat=(f2.nat$model[[1]] - f2.nat$residuals)
plot(resid(f2.nat)~(fitted.f2.nat), main = "Residuals vs Fitted", ylab="residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted.f2.nat, residuals(f2.nat)), col="red")
qqnorm(f2.nat$resid)
qqline(f2.nat$resid, col="red")
lev = hat(model.matrix(f2.nat))
plot(lev, main="leverage")
#df[lev>0.4,]

par(mfrow=c(2,2))
plot(density(resid(f1.land)),main = "residual density")
fitted.f1.land=(f1.land$model[[1]] - f1.land$residuals)
plot(resid(f1.land)~(fitted.f1.land), main = "Residuals vs Fitted", ylab="residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted.f1.land, residuals(f1.land)), col="red")
qqnorm(f1.land$resid)
qqline(f1.land$resid, col="red")
lev = hat(model.matrix(f1.land))
plot(lev, main="leverage")
#df[lev>0.4,]

par(mfrow=c(2,2))
plot(density(resid(f2.land)),main = "residual density")
fitted.f2.land=(f2.land$model[[1]] - f2.land$residuals)
plot(resid(f2.land)~(fitted.f2.land), main = "Residuals vs Fitted", ylab="residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted.f2.land, residuals(f2.land)), col="red")
qqnorm(f2.land$resid)
qqline(f2.land$resid, col="red")
lev = hat(model.matrix(f2.land))
plot(lev, main="leverage")
#df[lev>0.4,]

par(mfrow=c(2,2))
plot(density(resid(f1.tot)),main = "residual density")
fitted.f1.tot=(f1.tot$model[[1]] - f1.tot$residuals)
plot(resid(f1.tot)~(fitted.f1.tot), main = "Residuals vs Fitted", ylab="residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted.f1.tot, residuals(f1.tot)), col="red")
qqnorm(f1.tot$resid)
qqline(f1.tot$resid, col="red")
lev = hat(model.matrix(f1.tot))
plot(lev, main="leverage")
#df[lev>0.4,]

par(mfrow=c(2,2))
plot(density(resid(f2.tot)),main = "residual density")
fitted.f2.tot=(f2.tot$model[[1]] - f2.tot$residuals)
plot(resid(f2.tot)~(fitted.2.tot), main = "Residuals vs Fitted", ylab="residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted.f2.tot, residuals(f2.tot)), col="red")
qqnorm(f2.tot$resid)
qqline(f2.tot$resid, col="red")
lev = hat(model.matrix(f2.tot))
plot(lev, main="leverage")
#df[lev>0.4,]

# outputtable -------------------------------------------------------------
require(stargazer)
rob.se<-list(f1.nat.rob[,2],f2.nat.rob[,2],f1.land.rob[,2],f2.land.rob[,2],f1.tot.rob[,2],f2.tot.rob[,2])
rob.p<-list(f1.nat.rob[,4],f2.nat.rob[,4],f1.land.rob[,4],f2.land.rob[,4],f1.tot.rob[,4],f2.tot.rob[,4])
stargazer(f1.nat.rob, f2.nat.rob,f1.land.rob,f2.land.rob, f1.tot.rob, f2.tot.rob, dep.var.labels=c("ln of 'nature and species conservation' area in m² per capita", "ln of landscape protection area in m² per capita","ln of total protected area in ² per capita"), p=rob.p, se=rob.se, type="html", out="table1.htm", model.numbers = T, p=rob2.p, se=rob2.se, type="html", out="table2.htm", df=F, notes = "The panel data sample is balanced with n=16, T=3, N=nT=48. Robust standard errors are reported in parentheses below the estimated coefficients. Individual coefficients are indicated by a *10%, **5% or ***1% significance level. The models use an individual fixed effects specification.", notes.append = FALSE)


# mapping----
require(sp)
require(maps)                   
require(maptools)
Laender<-readShapeSpatial("Laender.shp", proj4string = CRS("+init=epsg:4326"))

require(latticeExtra)
#control variables
cols<-colorRampPalette(c("gray100","gray75","gray50"))
Laender@data$lnat_cap<-log(Laender@data$nat_cap)
nat.plot<-spplot(Laender, c("lnat_cap"), col.regions=cols(16), main=list(label="A        "))
nat.plot
Laender@data$lland_cap<-log(Laender@data$land_cap)
land.plot<-spplot(Laender, c("land_cap"), col.regions=cols(16), main=list(label="B        "))
land.plot
Laender@data$lges_cap<-log(Laender@data$ges_cap)
tot.plot<-spplot(Laender, c("ges_cap"), col.regions=cols(16), main=list(label="C        "))
tot.plot

#printing plots
require(gridExtra)
gridPlot<-grid.arrange(nat.plot, land.plot, tot.plot, nrow=1, ncol=3)

tiff("Fig1.tif", width = 21, height = 7, units = 'cm', res = 300, compression = "lzw+p")
grid.arrange(nat.plot, land.plot, tot.plot, nrow=1, ncol=3)
dev.off()


