install.packages(c("FLCore", "ggplotFL","FLAssess", "FLXSA","FLBRP","FLash","FLa4a"), repos="http://flr-project.org/R")
library(FLCore) 
library(ggplotFL)
library(reshape2)
library(FLXSA)
library(FLAssess)
library(FLBRP)
library(FLash)
library(FLa4a)
library(readxl)

setwd("C:/Users/au732610/OneDrive - Aarhus Universitet/Desktop")

#import landing/discard number data
land.n <- read_excel("Fish.xlsx",sheet="landings.n")
land.n<-melt(land.n,id.vars = "age",variable.name = "year",value.name = "data")
discard.n <- read_excel("Fish.xlsx", sheet = "discards.n")
discard.n<-melt(discard.n,id.vars = "age",variable.name = "year",value.name = "data")
dat.n<-rbind(land.n,discard.n)
dat.n$slot<-c(rep("landings.n",420),rep("discards.n",420))
dat.n$units=1000

#import landing/discard weight data
land.w <- read_excel("Fish.xlsx", sheet = "landings.wt")
land.w<-melt(land.w,id.vars = "age",variable.name = "year",value.name = "data")
discard.w <- read_excel("Fish.xlsx", sheet = "discards.wt")
discard.w<-melt(discard.w,id.vars = "age",variable.name = "year",value.name = "data")
dat.w<-rbind(land.w,discard.w)
dat.w$slot<-c(rep("landings.wt",420),rep("discards.wt",420))
dat.w$units="kg"

#Combine number and weight and convert it to FLS class

data<-rbind(dat.n,dat.w)
cod <- as.FLStock(data)

#import natural mortality and maturity
m(cod) <- c(1.054,0.89,0.284,0.2,0.2,0.2,0.2)
m.spwn(cod) <- harvest.spwn(cod) <- 0
mat(cod) <- c(0.071, 0.185, 0.248,0.587,0.927,1,1)

#Compute caught biomass
landings(cod) <- computeLandings(cod)
discards(cod) <- computeDiscards(cod)
catch(cod) <- computeCatch(cod, slot="all")
stock.wt(cod) <- catch.wt(cod)

#import plusgroup (7+) and f bar (2-4)
range(cod, c("plusgroup","minfbar", "maxfbar")) <- c(1,2, 4)
harvest(cod)[ac(range(cod)["max"]), ]     <- 1
harvest(cod)[, ac(range(cod)["maxyear"])] <- 1

#Build the vpa model and symthesize it with cod data
vpa_model <- VPA(cod, fratio = 0.5, fit.plusgroup = T)
vpa <- cod + vpa_model

#Import EPCU of NW, S, Viking and N3
NW <- read_excel("Fish.xlsx", sheet = "NW")
NW<-melt(NW,id.vars = "age",variable.name = "year",value.name = "data")
NW$slot<-c(rep("landings.wt",280))
 NW <- as.FLStock(NW)
 
S <- read_excel("Fish.xlsx", sheet = "S")
 S<-melt(S,id.vars = "age",variable.name = "year",value.name = "data")
 S$slot<-c(rep("landings.wt",280))
 S <- as.FLStock(S)
 
Viking <- read_excel("Fish.xlsx", sheet = "Viking")
 Viking<-melt(Viking,id.vars = "age",variable.name = "year",value.name = "data")
 Viking$slot<-c(rep("landings.wt",280))
 Viking <- as.FLStock(Viking)
 
N3 <- read_excel("Fish.xlsx", sheet = "N3")
 N3<-melt(N3,id.vars = "age",variable.name = "year",value.name = "data")
 N3$slot<-c(rep("landings.wt",217))
 N3 <- as.FLStock(N3)
 
 #Convert those EPCU to index class. Here I use ple4.index as the template
 data("ple4.index")
 NW.index<-ple4.index
 NW.index@index<-NW.index@index.var<-NW.index@catch.n<-NW.index@catch.wt<-NW.index@effort<-NW.index@sel.pattern<-NW.index@index.q<-NW@landings.wt 
 NW.index@index.var[1:7]<-NW.index@catch.wt[1:7]<-NW.index@sel.pattern[1:7]<-NW.index@index.q[1:7]<-NA
 NW.index@effort<-NW.index@effort[1,]
 NW.index@effort[1,]<-1
 NW.index@range["max"]<-7
 NW.index@range["minyear"]<-1983
 NW.index@range["maxyear"]<-2022
 
 S.index<-NW.index
 S.index@index<-S.index@catch.n<-S@landings.wt
 
 Viking.index<-NW.index
 Viking.index@index<-Viking.index@catch.n<-Viking@landings.wt
 
 N3.index<-NW.index
 N3.index@index<-N3.index@catch.n<-N3@landings.wt
 N3.index@index.var<-N3.index@catch.wt<-N3.index@effort<-N3.index@sel.pattern<-N3.index@index.q<-N3.index@index.var[,10:40]
 N3.index@effort<-N3.index@effort[1,]
 N3.index@effort[1,]<-1
 N3.index@range["minyear"]<-1992
 N3.index@index[7,13]<-0.1
 
 #Fit the EPCU indexes to vpa to tune the model
 fit.NW <- sca(vpa, NW.index)
 fit.S<-sca(vpa,S.index)
 fit.Viking<-sca(vpa,Viking.index)
 fit.N3<-sca(vpa,N3.index)
 
 #Plot tuning resultes in 3-d wireframe, in which F is the z axis
 wireframe(harvest(fit.NW), zlab="F")
 wireframe(harvest(fit.S), zlab="F")
 wireframe(harvest(fit.Viking), zlab="F")
 wireframe(harvest(fit.N3),zlab="F")
 
 #compare the AIC value of different tuning
 print(paste0('NW: ',round(AIC(fit.NW),4),' ',
              'S: ',round(AIC(fit.S),4),' ',
              'Viking: ',round(AIC(fit.Viking),4),' ',
              'N3:',round(AIC(fit.N3),4)))
 
 
# plot ssb-rec relationship 

ggplot(aes(ssb, rec), data=model.frame(FLQuants(vpa, "ssb", "rec"))) +
  geom_point() + geom_smooth(method="loess")

#convert vpa to ssb-rec model
vpasr <- as.FLSR(vpa)
ssb(vpasr)[,1]
rec(vpasr)[,1]

#Predict future recruitment using ricker model
model(vpasr) <- ricker()
#Fit the ricker model through MLE
vpasr<-fmle(vpasr)

vpasr_ri <- vpasr
# change model to bevholt
model(vpasr) <- bevholt()
# fit through MLE
vpasr_bh <- fmle(vpasr)
# change model to geometric mean
model(vpasr) <- geomean()
# fit through MLE
vpasr_geomean <- fmle(vpasr) 
#Compare the 3 models
print(paste0('Ricker: ',round(AIC(vpasr_ri),4),' ',
             'Beverton-Holt: ',round(AIC(vpasr_bh),4),' ',
             'Geomean: ',round(AIC(vpasr_geomean),4)))

plot(vpasr_geomean)

par(mfrow=c(1,3))
profile(vpasr_ri, main="Ricker")
profile(vpasr_bh, main="Beverton-Holt")
profile(vpasr_geomean, main="Geomean")

#short term forecast for 3 years (2023-2025)
maxyr_stk <- range(vpa)[["maxyear"]]
vpa_stf <- stf(vpa,nyears=3,wts.nyears=3, na.rm=TRUE)
maxyr_stf <- range(vpa_stf)[["maxyear"]]

#estimate the geometic mean recruitment for the future
mean_rec <- exp(mean(log(rec(vpa))))
vpa_sr <- as.FLSR(vpa, model="geomean")
params(vpa_sr)['a',] <- mean_rec
params(vpa_sr)

#Control the f for projection
round(fbar(vpa),3)
ggplot(fbar(vpa), aes(x=year,y=data)) + geom_line()

#define the f in the last year of the stock 
fbar_SQ <- mean(fbar(vpa)[,as.character(maxyr_stk)])

# Set the control object - year, quantity and value for the moment
ctrl_target <- data.frame(year = 2023:2025, quantity = "f", val = fbar_SQ)
ctrl_f <- fwdControl(ctrl_target)

vpa_sq <- fwd(vpa_stf, ctrl = ctrl_f, sr = vpa_sr)
#Set up different f values and scenarios
fbar_multiplier <- seq(from = 0, to = 2, by = 0.4)
fbar_scenarios <- cbind(rep(fbar_SQ,length(fbar_multiplier)),
                        fbar_multiplier*fbar_SQ,
                        fbar_multiplier*fbar_SQ)

# Add the msy scenario as a final scenario 
msy <- c(refpts(brp(FLBRP(vpa)))["msy","harvest"])
fbar_scenarios <- rbind(fbar_scenarios, c(fbar_SQ,msy,msy))



colnames(fbar_scenarios) <- c("2023","2024","2025")
rownames(fbar_scenarios) <- c(fbar_multiplier, "msy")

stf_results <- matrix(NA,nrow = nrow(fbar_scenarios),ncol = 8)
# Set some column names
colnames(stf_results) <- c('Fbar',
                           paste0('Catch',maxyr_stk+1), 
                           paste0('Catch',maxyr_stk+2),
                           paste0('Catch',maxyr_stk+3),
                           paste0('SSB',maxyr_stk+2),
                           paste0('SSB',maxyr_stk+3),
                           paste0('SSB_change_',maxyr_stk+2,'-',maxyr_stk+3,'(%)'),
                           paste0('Catch_change_',maxyr_stk,'-',maxyr_stk+2,'(%)'))

stk_stf <- FLStocks()
# Loop over the scenarios (each row in the fbar_scenarios table)
for (scenario in 1:nrow(fbar_scenarios)) {
  cat("Scenario: ", scenario, "\n")
  flush.console()
  # Make a target object with F values for that scenario
  # Set the control object 
  ctrl_target <- data.frame(year = (maxyr_stf-2):maxyr_stf,
                            quantity = "f",
                            val = fbar_scenarios[scenario,])
  # ctrl_target
  ctrl_f <- fwdControl(ctrl_target)
  # Run the forward projection. We could include an additional argument, maxF.
  # By default the value of maxF is 2.0. It could be increased to 10.0, say,
  # so that F is less limited, and the bound is not hit (not a problem here).
  vpa_fwd <- fwd(vpa_stf, ctrl = ctrl_f, sr = vpa_sr)#, maxF = 10.0)
  ## Check it has worked - uncomment out to check scenario by scenario
  # plot(ple4_fwd[,ac(2001:2011)])
  # Store the result - if you want to, comment out if unnecessary
  stk_stf[[as.character(scenario)]] <- vpa_fwd
  
  # Fill results table
  stf_results[scenario,1] <- round(fbar(vpa_fwd)[,ac(2025)],3) # final stf year
  stf_results[scenario,2] <- catch(vpa_fwd)[,ac(maxyr_stk+1)] # 1st stf year
  stf_results[scenario,3] <- catch(vpa_fwd)[,ac(maxyr_stk+2)] # 2nd stf year
  stf_results[scenario,4] <- catch(vpa_fwd)[,ac(maxyr_stk+3)] # final stf year
  stf_results[scenario,5] <- ssb(vpa_fwd)[,ac(maxyr_stk+2)] # 2nd stf year
  stf_results[scenario,6] <- ssb(vpa_fwd)[,ac(maxyr_stk+3)] # final stf year
  
  # change in ssb in last two stf years
  stf_results[scenario,7] <- round((ssb(vpa_fwd)[,ac(maxyr_stk+3)]-ssb(vpa_fwd)[,ac(maxyr_stk+2)])/
                                     ssb(vpa_fwd)[,ac(maxyr_stk+2)]*100,1) 
  
  # change in catch from true year, to 2nd to last stf year
  stf_results[scenario,8] <- round((catch(vpa_fwd)[,ac(maxyr_stk+2)]-catch(vpa_fwd)[,ac(maxyr_stk)])/
                                     catch(vpa_fwd)[,ac(maxyr_stk)]*100,1)  
}
# Give the FLStocks object some names
names(stk_stf) <- rownames(fbar_scenarios)

# Plotting
plot(stk_stf)

#Biological reference points
brp <- FLBRP(vpa)
summary(brp)

#Long-term prognosis: Build stock-recrument model for biological reference points
vpaSR <- transform(as.FLSR(vpa, model=ricker), ssb=ssb/100, rec=rec/100)
vpaSR <- fmle(vpaSR,control=list(silent=T))
params(vpaSR)['b',] <- params(vpaSR)['b',] / 100
vpaSR <- transform(vpaSR, ssb=ssb*100, rec=rec*100)
brpRi <- FLBRP(vpa, sr=vpaSR)
brpRi <- brp(brpRi)

plot(brpRi)

#Adding price
price<-c(8.165,9.260,9.705,10.625,21.34,21.34,21.34)
#converting price to euros
price<-price*0.089
price(brpRi)<-price
price(brpRi)@units <- "keuro/ton"
# variable costs per F 
vcost(brpRi) <- 70000
vcost(brpRi)@units <- "keuro/unit F"
# fixed costs per F 
fcost(brpRi) <- 35000
fcost(brpRi)@units <- "keuro/unit F"
brpEco <- brp(brpRi)
plot(brpEco)