# ADathe 2014-11-12 - adapted for class 2022-09-30
# assume no dependence from O2 in Monod kinetics
# following an example of Karline Soetaert
# with the help of Paul Torfs from Wageningen

# only run if necessary!
#rm(list=ls()) # make work directory clean
# rm(list=setdiff(ls(), "df"))

# set or change working directory
setwd("C:/Users/ad273/Documents/Classes/Cornell_PLSCI_7202")

# load libraries
library(ggplot2) #library for plotting
library(reshape2) # library for reshaping data (tall-narrow <-> short-wide)
library(deSolve) # library for solving differential equations
library(minpack.lm) # library for least squares fit using levenberg-marquart algorithm
library(FME)

#load concentration data
# concentrations are given in mol/flask
df=read.table("flask21n.txt", header=TRUE)

# plot data for O2 and CO2 on two axis
op  <- par(las=1,xaxs="r",mai=c(1,1,1,1), 
           mfrow = c(1,1)) #, mar = c(2, 3, 2, 1))

plot(df[,1],df[,2],
     xlab="time (hours)",
     ylab="O2 (mol flask-1)",
     main="flask 21",
     type='l',
     col=1,
     xlim=c(0,350),
     ylim=c(8.5e-4,1.0e-3)
)

par(new=TRUE)

plot(df[,1],df[,3],xlim=c(0,350),xaxt="n",yaxt="n",xlab="",ylab="",pch=16)
axis(4,at=c(0,20,40,60,80,100))
text(11, 50, "CO2 (mol flask-1)", srt = 270, xpd = TRUE)
par(op)  # reset par

# PG decay rate function
rxnrate=function(t, cinit,pars,positive = TRUE){  
  with (as.list(c(cinit,pars)), {
    # these have fixed values
    # kd  =5.6e-4 # ORCHESTRA kd=1.556e-7 [mol.l-1.s-1]
    yPG =0.32    # ORCHESTRA [molB.molPG-1]
    yO2 =2.102    # 2.575 estimated following Barry et al. 
    yCO2=1.402    # ORCHESTRA [molCO2.molPG-1]
    yBd=0.5
    
    dPG=  muO2*(-1/yPG)     *PG/(kO2PG+PG)*B # dPG/dt
    dS= dPG -    yBd*d*B                     # dead bacteria back in as substrate, with yield factor.
    # negative because dPG is negative and dS later goes in negative
    dB=   yPG   *-dS         -d*B                        # dB/dt
    dO2=  yO2   * dS  -kd*B  -(1-yBd)*5*d*B              # dO2/dt
    dCO2= yCO2  *-dS  +kd*B  +(1-yBd)*5*d*B              # dCO2/dt
    
       return(list(c(dPG,dB,dO2,dCO2)))
  })
}

t=df$time
pars <- list(muO2=5.4e-2, kO2PG=1.0e-5, B_0=2.0e-8, kd =2.0e-3, d=5.5e-4)

bactModel=function(pars,t){
  cinit <- c(PG=2.74e-5, B=pars[[3]], O2=9.59e-4, CO2=2.99e-6)
  P<-pars
  P["muO2"]  <-pars[1]
  P["kO2PG"] <-pars[2]
  P["B"]     <-pars[3]
  P["kd"]    <-pars[4]
  P["d"]     <-pars[5]
  
  ## ode solves the model by integration ...
  return(as.data.frame(ode(y = cinit, times = t, func = rxnrate,
                           parms=P)))
}

out <- bactModel(as.numeric(pars),t)
#write.table(out,"result-f21.txt")
head(out)
tail(out)

# plot the simulation results
old.par = par(no.readonly=TRUE)
par(mfrow = c(2,2), mar = c(4, 5, 1, 1))

plot(out$time, out$B, main = "Bacteria",
     xlab = "time, hour", ylab = "mol/flask", type = "l", lwd = 2)

plot(out$time, out$PG, main = "Propylene Glycol",
     xlab = "time, hour", ylab = "mol/flask", type = "l", lwd = 2)

plot(out$time, out$O2, main = "O2",
     xlab = "time, hour", ylab = "mol/flask", type = "l", lwd = 2)

plot(out$time, out$CO2, main = "CO2",
     xlab = "time, hour", ylab = "mol/flask", type = "l", lwd = 2)

par(old.par)
