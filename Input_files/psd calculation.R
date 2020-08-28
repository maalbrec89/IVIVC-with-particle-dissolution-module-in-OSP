rm(list = ls())
library(arules)
library(RColorBrewer)
library(caTools)
require(xlsx)

#############################################
##### FUNCTIONS #############################
#############################################

fitDistribToQ <- function(par,qobs,xobs,distr){
  sum( (qwrapper(qobs,par,distr)-xobs)^2)
}

pwrapper <- function(q,par,distr){
  if(distr=="weibull"){
    res <- pweibull(q = q,shape = par[1],scale = par[2])
  } else if(distr=="lognormal"){
    res <- plnorm(q = q,meanlog = par[1],sdlog = par[2])
  } else if(distr=="normal"){
    res <- pnorm(q = q,mean = par[1],sd = par[2])
  } else if(distr=="gamma"){
    res <- pgamma(q = q,shape = par[1],rate = par[2])
  } else {
    error("Unknown distribution specified.")
  }
  return(res)
}

dwrapper <- function(x,par,distr){
  if(distr=="weibull"){
    res <- dweibull(x = x,shape = par[1],scale = par[2])
  } else if(distr=="lognormal"){
    res <- dlnorm(x = x,meanlog = par[1],sdlog = par[2])
  } else if(distr=="normal"){
    res <- dnorm(x = x,mean = par[1],sd = par[2])
  } else if(distr=="gamma"){
    res <- dgamma(x = x,shape = par[1],rate = par[2])
  } else {
    error("Unknown distribution specified.")
  }
  return(res)
}

qwrapper <- function(p,par,distr){
  if(distr=="weibull"){
    res <- qweibull(p = p,shape = par[1],scale = par[2])
  } else if(distr=="lognormal"){
    res <- qlnorm(p = p,meanlog = par[1],sdlog = par[2])
  } else if(distr=="normal"){
    res <- qnorm(p = p,mean = par[1],sd = par[2])
  } else if(distr=="gamma"){
    res <- qgamma(p = p,shape = par[1],rate = par[2])
  } else {
    error("Unknown distribution specified.")
  }
  return(res)
}

#############################################
##### MEASURED PARTICLE SIZE DATA ###########
#############################################

# set working directory to "<your local path>/IVIVE-with-particle-dissolution-module-in-OSP"
workingDir <- file.path("C:/OSP","IVIVE-with-particle-dissolution-module-in-OSP")
setwd(workingDir)

compound <- "CompoundA"    # compound name
batch <- "BatchX"          # batch name
size_unit <- "µm"          # unit of the particle size
obsDataFile <- "particle_distribution_raw_data.csv"     # file name
obsData <- read.table(file.path("Input_files",obsDataFile),sep=",",header = TRUE, as.is = TRUE)
names(obsData)

obsData$q_obs <- obsData$q_obs/100   # convert percentage to fraction
x_valid <- obsData$q_obs>0 & obsData$q_obs <1
t <- seq(0,obsData$x_obs[length(obsData$x_obs)],by=0.01) # for plotting

# Number of bins
NBins <- 10

# Distributions to test
distributions <- c("weibull","lognormal","gamma","normal")

# Config
color <- brewer.pal(9, "Set1")
colorRaw <- brewer.pal(11, "Set3")
colorBg <- rgb(t(col2rgb(colorRaw)),alpha = 50, maxColorValue = 255)
colorLi <- rgb(t(col2rgb(colorRaw)),alpha = 255, maxColorValue = 255)


#############################################
##### FITTING ###############################
#############################################

models <- list()
res <- data.frame(name=distributions,error=NA)
for(ii in seq_along(distributions)){
  m <- optim(par=c(1,1), fitDistribToQ, qobs=obsData$q_obs[x_valid], xobs=obsData$x_obs[x_valid], distr=distributions[ii])
  res$error[ii] <- m$value
  models <- c(models,list(m))
}
names(models) <- distributions

(res)
cat(as.character(res$name[which.min(res$error)]))

finalModel <- models[[which.min(res$error)]]
finalDistr <- distributions[which.min(res$error)]

png(filename = file.path("Output_files",paste0(paste("FitComp",compound,batch,sep="_"),".png")), width = 12, height = 12, pointsize = 8 ,res = 300, units = "cm")
plot(NA,xlim = c(min(t),max(t)),ylim = c(0,100),las=1,ylab="Cumulative probability [%]",xlab=paste0("Particle size [",size_unit,"]"), paste(compound,batch,sep = " / "))
for(ii in seq_along(distributions)){
  lines(t,100*pwrapper(t,models[[ii]]$par,distr = distributions[ii]),col=color[ii],lwd=2)
}
points(obsData$x_obs,obsData$q_obs*100,pch=16,col="black")
legend("bottomright", legend = c("observed",distributions),
       col = c("black",color),lwd = c(NA,rep(2,length(distributions))),pch=c(16,rep(NA,length(distributions))))
dev.off()

#############################################
##### DISCRETIZATION ########################
#############################################

# equidistant radius bin means
(r_min <- qwrapper(p = 0.001,finalModel$par,distr = finalDistr)) 
(r_max <- qwrapper(p = 0.999,finalModel$par,distr = finalDistr)) 
(r_bin_borders <- r_min+(seq(1,NBins+1,1)-1)*(r_max-r_min)/(NBins))
#(r_bin_means <- diff(r_bin_borders)/2+r_bin_borders[1:(length(r_bin_borders)-1)])

(q_bin_borders <- pwrapper(r_bin_borders,finalModel$par,distr = finalDistr))
(q_bin_means <- diff(q_bin_borders)/2+q_bin_borders[1:(length(q_bin_borders)-1)])
(r_bin_means <- qwrapper(q_bin_means,finalModel$par,distr = finalDistr))

(rel_amountFactor <- dwrapper(r_bin_means,finalModel$par,distr = finalDistr))
trapz(r_bin_means,rel_amountFactor)

dat <- data.frame(bin=seq(1,NBins,1),r_bin_means=r_bin_means,q_bin_means=q_bin_means,rel_amountFactor=rel_amountFactor,
           r_bin_border_low=r_bin_borders[1:(length(r_bin_borders)-1)],r_bin_border_up=r_bin_borders[2:(length(r_bin_borders))],
           distribution=c(finalDistr,rep("",NBins-1)),par1=c(finalModel$par[1],rep("",NBins-1)),par2=c(finalModel$par[2],rep("",NBins-1)))
write.table(x = dat,file = file.path("Output_files",paste0(paste("ParticleDistr",compound,batch,finalDistr,sep="_"),".csv")),
            quote = FALSE,sep = ";",row.names = FALSE,col.names = TRUE)

#############################################
##### PLOT BINS #############################
#############################################

t <- seq(0,r_bin_borders[length(r_bin_borders)],by=0.01) # for plotting

# density:
png(filename = file.path("Output_files",paste0(paste("ProbDensBinned",compound,batch,sep="_"),".png")),
                         width = 12, height = 12, pointsize = 8 ,res = 300, units = "cm")
plot(t,dwrapper(t,finalModel$par,distr = finalDistr),type="l",las=1,
     xlab=paste0("Particle size [",size_unit,"]"),ylab="Density",main=paste(compound,batch,sep = " / "))
#points(r_bin_means,rel_amountFactor)
for(ii in 1:(NBins)){
  tBin <- t[t >= r_bin_borders[ii] & t < r_bin_borders[ii+1]]
  polygon(c(tBin,rev(tBin)),
            c(dwrapper(tBin,finalModel$par,distr = finalDistr),rep(0,length(tBin))),
          col = colorBg[ii])
  lines(rep(r_bin_means[ii],2),c(0,dwrapper(r_bin_means[ii],finalModel$par,distr = finalDistr)),col=colorLi[ii])
  points(r_bin_means[ii],dwrapper(r_bin_means[ii],finalModel$par,distr = finalDistr),col=colorLi[ii],pch=20)
}
dev.off()

# cumulative density:
png(filename = file.path("Output_files",paste0(paste("CumDistrBinned",compound,batch,sep="_"),".png")),
                         width = 12, height = 12, pointsize = 8 ,res = 300, units = "cm")
plot(t,100*pwrapper(t,finalModel$par,distr = finalDistr),type="l",lwd=2,las=1,
     xlab=paste0("Particle size [",size_unit,"]"),ylab="Cumulative probability [%]",col="black",main=paste(compound,batch,sep = " / "))
#points(r_bin_means,rel_amountFactor)
for(ii in 1:(NBins)){
  tBin <- t[t >= r_bin_borders[ii] & t < r_bin_borders[ii+1]]
  polygon(c(tBin,rev(tBin)),
          100*c(pwrapper(tBin,finalModel$par,distr = finalDistr),rep(0,length(tBin))),
          col = colorBg[ii])
  lines(rep(r_bin_means[ii],2),100*c(0,pwrapper(r_bin_means[ii],finalModel$par,distr = finalDistr)),col=colorLi[ii])
  points(r_bin_means[ii],100*pwrapper(r_bin_means[ii],finalModel$par,distr = finalDistr),col=colorLi[ii],pch=20)
}
points(obsData$x_obs,obsData$q_obs*100,pch=16,col="black")
legend("topleft",legend = c("observed","simulated"),pch=c(16,NA),lwd=c(NA,2),col="black")
dev.off()

#############################################
##### EXPORT RESULTS FOR MOBI ###############
#############################################

ContainerPath <- NULL
for(i in 1:NBins){
  ContainerPath <- c(ContainerPath, paste0("Organism|Particle_Size_BIN",i))
}

radius_bin_means <- r_bin_means/2   # particle size divided by 2 to obtain particle radius

PSV <- data.frame("Container Path"=rep(ContainerPath,2), "Parameter Name"=c(rep("radius (at t=0)",NBins),rep("rel_amountFactor",NBins)),
                  Value=c(radius_bin_means,rel_amountFactor), Unit=c(rep("µm",NBins),rep("",NBins)))
write.xlsx(x = PSV,file = file.path("Output_files",paste0(paste("PSV",compound,batch,finalDistr,sep="_"),".xlsx")),sheetName="Sheet1",row.names=F,col.names=T)
