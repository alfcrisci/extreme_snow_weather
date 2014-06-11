##################################################################################
#  
# Authors : Alfonso Crisci,Gianni Messeri,Valerio Capecchi
# IBIMET CNR Institute of Biometeorology Firenze via caproni 8,50145,Italia  
# Consorzio LaMMA CNR - Regione Toscana                             
# mail: a.crisci@ibimet.cnr.it
# mail: messeri@lamma.rete.toscana.it
# file: extrem_snow.r
# github: https://github.com/alfcrisci/extreme_weather 
# language= R
# Licence: https://creativecommons.org/licenses/by/3.0/it/
##################################################################################

#####################################################################################à
# load and install library
if("POT" %in% rownames(installed.packages())==FALSE){install.packages("POT", repos="http://R-Forge.R-project.org")}
library(POT)

#####################################################################################à
# define useful function

value_retlev=function(val,model,npy=1)  {loc <- model$threshold[1]
                                   scale <- model$scale
                                   shape <- model$param["shape"] 
								   p=as.numeric(pgpd(val,loc,scale,shape))
								   return(prob2rp(p,npy=npy)[,"retper"]) 
}
#####################################################################################à
# define working directory
setwd("D:\\")
########################################################################################à
# load data

h_snow_FI=read.csv("h_snow_FI.csv",header=T)
names(h_snow_FI)
x=h_snow_FI$h_snow

##############################################################################
# Analisi con tutti i metodi usando la Pareto Gmeneralizzata per analisi estremi con superamento a soglia. ( 1 cm)
##############################################################################


mom <- fitgpd(x, 1, "moments")$param
mle <- fitgpd(x, 1, "mle")$param
pwmu <- fitgpd(x, 1, "pwmu")$param
pwmb <- fitgpd(x, 1, "pwmb")$param
pickands <- fitgpd(x, 1, "pickands")$param
med <- fitgpd(x, 1, "med", start = list(scale = 2, shape = 0.25))$param
mdpd <- fitgpd(x, 1, "mdpd")$param
mple <- fitgpd(x, 1, "mple")$param
ad2r <- fitgpd(x, 1, "mgf", stat = "AD2R")$param

print( " Tabella valori per i differenti metodi della Librareia POT")
print(rbind(mom, mle, pwmu, pwmb, pickands, med, mdpd, mple,ad2r))


##############################################################################
# Scelgo il metodo MLE soglia 1 cm
##############################################################################

mle_neve_FI <- fitgpd(x, 1, "mle")
h_snow_FI$return_time=value_retlev(h_snow_FI$h_snow,mle_neve_FI)
png("curva_tempi_ritorno.png",width = 800, height = 680, units = "px", pointsize = 15)
retlev(mle_neve_FI, npy = 1,points=T,xlimsup=100,main="Curva di ritorno altezza di neve a Firenze\n R POT Package Metodo MLE",xlab="Tempo di ritorno (Anni)",ylab="Altezza di Neve");
points(h_snow_FI$return_time[51],28,col="red",pch = 23,bg="red")
abline(h = 28, v = 0, col = "gray60")

abline(h = 0, v = h_snow_FI$return_time[51], col = "gray60")
text(25,5, paste0("Tempo di ritorno Neve 2010 Firenze\n Anni Stimati: ",signif(h_snow_FI$return_time[51],3)), col = 2)
dev.off()


write.csv(h_snow_FI,"h_snow_FI_res.csv")

# Reference R POT functions 
# http://www4.stat.ncsu.edu/~mannshardt/st810/Lectures/Lec19Extremes.pdf
# http://www.juergen-grieser.de/extremerain_standalone.pdf
# http://pot.r-forge.r-project.org/docs/guide.pdf
# http://forum.meteonetwork.it/nowcasting-discussioni-climatiche-italia-centrale/149886-massima-altezza-neve-suolo-varie-localita-toscane-3.html

