# Species: Cattle
# Ivermectin B1a, with a molecular weight of 875.1 g/mol, is highly lipophilic, exhibits slow absorption depending on dose and formulation, undergoes tissue accumulation, is minimally metabolized, and is primarily eliminated via bile and feces. (https://doi.org/10.1002/psp4.12950)



# Clear current environment -------------------------------------------------------------------------------------------------------
install.packages("pkgbuild")
library(pkgbuild)
pkgbuild::find_rtools()
rm(list=ls(all=TRUE))
getwd()
setwd("C:/Users/jeongminchoi/Desktop/Ivermectin_PK datasets_for calibration and evaluation")

{
  library(readxl)      # Package for reading Excel files, supporting both .xls and .xlsx file formats
  library(mrgsolve)    # Package for Loading mrgsolve code into r via mcode from the 'mrgsolve' package
  library(magrittr)    # Package for the pipe, %>% 
  library(dplyr)       # Package for data manipulation and transformation using a grammar of data manipulation
  library(ggplot2)     # Package for plots
  library(FME)         # Package for MCMC simulation and model fitting
  library(minpack.lm)  # Package for model fitting
  library(reshape)     # Package for melt function to reshape the table
  library(truncnorm)   # Package for the truncated normal distribution function   
  library(EnvStats)    # Package for Environmental Statistics, including the US EPA Guidance
  library(invgamma)    # Package for inverse gamma distribution function
  library(foreach)     # Package for executing loop iterations in parallel, facilitating efficient computation
  library(doParallel)  # Package for parallel computing
  library(bayesplot)   # Package for MCMC traceplot
  library(tidyr)       # Package for tidy messy data
  library(tidyverse)   # Package for tidy messy data
  library(ggpubr)      # Package for plotting the data
}



# Build mrgsolve-based PBPK Model --------------------------------------------------------------------------------------------------
SolvePBPK <- '
$PARAM @annotated

// IV infusion rate constants
Timeiv     : 0.01       : /h,                 IV injection/infusion time

// Oral absorption rate constants
Kst        : 0.06       : /h,                 Gastric emptying rate constant (from stomach to small intestine)
Ka         : 0.035      : /h,                 Intestinal absorption rate constant (from intestine to systemic circulation)
Kint       : 0.028      : /h,                 Intestinal transit rate constant (from intestine to colon)

// SC absorption rate constants
Ksc        : 0.2        : /h,                 Absorption rate constant of SC administration; 0.0458 for sheep value from Lo et al., 1985
Kdiss      : 0.012      : /h,                 Dissolution rate constant of SC administration

// Elimination rate constants
KbileC     : 0.017      : L/h/kg,             Biliary clearance (from liver to intestine)
Kfeces     : 0.012      : /h,                 Fecal elimination rate constant (from colon to feces)
KurineC    : 0.0001     : L/h/kg,             Urinary clearance

// Physiological parameters
// Blood flow rates
QCC        : 5.45       : L/h/kg,             Cardiac output; Lin et al., 2020
QLCa       : 0.44       : Unitless,           Fraction of blood flow to the liver; value from Lin et al., 2020
QKCa       : 0.11 	    : Unitless,           Fraction of blood flow to the kidney; value from Lin et al., 2020
QMCa       : 0.28       : Unitless,           Fraction of blood flow to the muscle; value from Lin et al., 2020
QFCa       : 0.08       : Unitless,           Fraction of blood flow to the fat; value from Li et al., 2019
QLuCa      : 1          : Unitless,           Fraction of blood flow to the lung; value from Lin et al., 2020; Brown et al., 1997
QRestCa    : 0.09       : Unitless,           Fractional rest of body; value calculated from 1-QLC-QKC-QMC-QFC

// Tissue  volumes
VLCa       : 0.0122     : L/kg bw,            Fraction volume of liver; valule from Lin et al., 2020
VKCa       : 0.0021     : L/kg bw,            Fraction volume of kidney; valule from Lin et al., 2020
VMCa       : 0.361      : L/kg bw,            Fraction volume of muscle tissue; valule from Lin et al., 2020
VFCa       : 0.1218     : L/kg bw,            Fraction volume of fat tissue; valule from Lin et al., 2020
VLuCa      : 0.0085     : L/kg bw,            Fraction volume of lung tissue; valule from Lin et al., 2020
VbloodCa   : 0.0399     : L/kg bw,            Fraction volume of blood; valule from Lin et al., 2020
VvenCa     : 0.01737    : L/kg bw,            Fraction volume of venous plasma; 70% of total plasma
VartCa     : 0.00745    : L/kg bw,            Fraction volume of arterial plasma; 30% of total plasma
VRestCa    : 0.4696     : L/kg bw,            Fractional rest of body, total sum to 1 (1-VLC-VKC-VMC-VFC-VLuC-VplasmaC)
Hct        : 0.378      : Unitless,           Hematocrit
BW         : 655.4      : kg                  Body weight of cattle; 655.4 kg (≈ 1445 lbs) for the average market weight for beef cattle (USDA Daily Cattle Reports, Aug 2026); 552 kg for adult cattle (Meyer et al., 2014)
MW         : 875.1      : g/mol,              Ivermectin B1a molecular weight

// Chemical-specific parameters
// Partition coefficients (PC, tissue:plasma)
PL         : 17         : Unitless,           Liver:plasma PC
PK         : 1          : Unitless,           Kidney:plasma PC 
PM         : 0.65       : Unitless,           Muscle:plasma PC
PF         : 8.75       : Unitless,           Fat:plasma PC; 5.75 for human model from Yeo & Wesche, 2023
PLu        : 0.5        : Unitless,           Lung:plasma PC
PRest      : 1.5        : Unitless,           Rest of the body:plasma PC
PB         : 0.93       : Unitless,           Fraction of plasma protein binding; PubChem

$MAIN

// Blood flows to organs
double sumQ = QLCa + QKCa + QMCa + QFCa + QRestCa;                            // Unitless, Sum of cardiac output fractions
double QLC = QLCa/sumQ;                                                       // Unitless, Adjusted blood flow rate fraction to liver
double QKC = QKCa/sumQ;                                                       // Unitless, Adjusted blood flow rate fraction to kidney
double QMC = QMCa/sumQ;                                                       // Unitless, Adjusted blood flow rate fraction to muscle
double QFC = QFCa/sumQ;                                                       // Unitless, Adjusted blood flow rate fraction to fat
double QLuC = QLuCa;                                                          // Unitless, Adjusted blood flow rate fraction to lung
double QRestC = QRestCa/sumQ;                                                 // Unitless, Adjusted blood flow rate fraction to rest of body
double QC = QCC*BW;                                                           // L/h, Cardiac output
double QL = QLC*QC;                                                           // L/h, Blood flow in liver
double QK = QKC*QC;                                                           // L/h, Blood flow in kidney
double QM = QMC*QC;                                                           // L/h, Blood flow in muscle
double QF = QFC*QC;                                                           // L/h, Blood flow in fat
double QLu = QLuC*QC;                                                         // L/h, Blood flow in lung
double QRest = QRestC*QC;                                                     // L/h, Blood flow in in rest of body

//  Tissue volumes
double sumV = VLCa + VKCa + VMCa + VFCa + VLuCa + VvenCa + VartCa + VRestCa;  // L/kg bw, Sum of vascular and tissue volumes
double VLC = VLCa/sumV;                                                       // Unitless, Adjusted fraction of liver volume
double VKC = VKCa/sumV;                                                       // Unitless, Adjusted fraction of kidney volume
double VMC = VMCa/sumV;                                                       // Unitless, Adjusted fraction of muscle volume
double VFC = VFCa/sumV;                                                       // Unitless, Adjusted fraction of fat volume
double VLuC = VLuCa/sumV;                                                     // Unitless, Adjusted fraction of lung volume
double VvenC = VvenCa/sumV;                                                   // Unitless, Adjusted fraction of venous plasma volume
double VartC = VartCa/sumV;                                                   // Unitless, Adjusted fraction of arterial plasma volume
double VRestC = VRestCa/sumV;                                                 // Unitless, Adjusted fraction of rest of body volume
double VL = VLC*BW;                                                           // L, Liver volume
double VK = VKC*BW;                                                           // L, Kidney volume
double VM = VMC*BW;                                                           // L, Muscle volume
double VF = VFC*BW;                                                           // L, Fat volume
double VLu = VLuC*BW;                                                         // L, Lung volume
double VRest = VRestC*BW;                                                     // L, Rest of body volume
double Vven = VvenC*BW;                                                       // L, Venous plasma volume
double Vart = VartC*BW;                                                       // L, Arterial plasma volume
double Free = 1-PB;                                                           // Unitless, Fraction of free drug in plasma
double Kurine = KurineC*BW;                                                   // L/h, Urinary clearance
double Kbile = KbileC*BW;                                                     // L/h, Biliary clearance

$CMT AA AV AIV ADOSE ADOSEsc Absorbsc Asitesc AST AI AAO AL AK AM AF ALu ARest AUCCV AUCCL AUCCK AUCCM AUCCF AUCCLu Aurine Afeces Abile Acolon ABileGut

$ODE

dxdt_AIV=0;
dxdt_ADOSE=0;

// Concentrations in plasma, tissues, and venous outflow from tissues (mg/L)
double CV = AV/Vven;                                                          // Venous plasma concentration
double CA = AA/Vart;                                                          // Arterial plasma concentration
double CL = AL/VL;                                                            // Liver concentration
double CK = AK/VK;                                                            // Kidney concentration
double CM = AM/VM;                                                            // Muscle concentration
double CF = AF/VF;                                                            // Fat concentration
double CLu = ALu/VLu;                                                         // Lung concentration
double CRest = ARest/VRest;                                                   // Rest of body concentration
double CVL = AL/(VL * PL);                                                    // Liver venous outflow concentration
double CVK = AK/(VK * PK);                                                    // Kidney venous outflow concentration
double CVM = AM/(VM * PM);                                                    // Muscle venous outflow concentration
double CVF = AF/(VF * PF);                                                    // Fat venous outflow concentration
double CVLu = ALu/(VLu * PLu);                                                // Lung venous outflow concentration
double CVRest = ARest/(VRest * PRest);                                        // Rest of body venous outflow concentration

// Rate of changes in amount of drug in each compartment (mg/h)    
double Rsc = Ksc * Asitesc;                                                   // Subcutaneous absorption rate
double Rsitesc = - Rsc + Kdiss * ADOSEsc;                                     // Drug amount change at the subcutaneous injection site
double Rbile = Kbile * CVL;                                                   // Biliary excretion rate
double RAST = - Kst * AST;                                                    // Gastric emptying rate
double RAI = Kst * AST - Kint * AI - Ka * AI;                                 // Drug amount change at the intestine
double RAO = Ka * AI;                                                         // Oral absorption rate
double RBileGut = Rbile - Kint * ABileGut;                                    // Bile-derived drug to intestine; no reabsoprtion
double Rcolon = Kint * AI + Kint * ABileGut - Kfeces * Acolon;                // Intestinal transit to colon
double Rfeces = Kfeces * Acolon;                                              // Fecal excretion rate
double Rurine = Kurine * CVK;                                                 // Urinary excretion rate
double RL = QL * (CA - CVL) * Free - Rbile;                                   // Rate of change in amount of the drug in liver
double RK = QK * (CA - CVK) * Free - Rurine;                                  // Rate of change in amount of the drug in kidney
double RM = QM * (CA - CVM) * Free;                                           // Rate of change in amount of the drug in muscle
double RF = QF * (CA - CVF) * Free;                                           // Rate of change in amount of the drug in fat
double RLu = QLu * (CV - CVLu) * Free;                                        // Rate of change in amount of the drug in lung
double RRest = QRest * (CA - CVRest) * Free;                                  // Rate of change in amount of the drug in rest of the body
double RV = (QL * CVL  + QK  * CVK  + QM * CVM + QF * CVF + 
              QRest * CVRest)*Free - QC*CV*Free + RAO + Rsc;                  // Rate of change in amount of the drug in systemic venous plasma
double RA = QC * (CVLu - CA) * Free;                                          // Rate of change in amount of the drug in arterial plasma

// The amount in each compartment (mg)
dxdt_AV = RV;                                                                 // Amount of drug in venous plasma
dxdt_AA = RA;                                                                 // Amount of drug in arterial plasma
dxdt_Absorbsc = Rsc;                                                          // Amount of drug in SC absorption
dxdt_Asitesc = Rsitesc;                                                       // Amount of drug in SC injection site
dxdt_ADOSEsc = - Kdiss * ADOSEsc;                                             // Amount of drug in SC depot
dxdt_AST = RAST;                                                              // Amount of drug in stomach
dxdt_AI = RAI;                                                                // Amount of drug in intestine
dxdt_ABileGut = RBileGut;                                                     // Amount of drug in intestine from biliary excretion; no reabsoprtion
dxdt_Acolon = Rcolon;                                                         // Amount of drug in colon
dxdt_AAO = RAO;                                                               // Amount of absorbed drug
dxdt_Aurine = Rurine;                                                         // Amount of drug in urine
dxdt_Afeces = Rfeces;                                                         // Amount of drug in feces
dxdt_Abile = Rbile;                                                           // Amount of drug in bile
dxdt_AL = RL;                                                                 // Amount of drug in liver
dxdt_AK = RK;                                                                 // Amount of drug in kidney
dxdt_AM = RM;                                                                 // Amount of drug in muscle
dxdt_AF = RF;                                                                 // Amount of drug in fat
dxdt_ALu = RLu;                                                               // Amount of drug in lung
dxdt_ARest = RRest;                                                           // Amount of drug in rest of the body

// The AUC of drug in each compartment (mg*h/L) 
dxdt_AUCCV = CV;                                                              // AUC of drug in venous plasma
dxdt_AUCCL = CL;                                                              // AUC of drug in liver
dxdt_AUCCK = CK;                                                              // AUC of drug in kidney
dxdt_AUCCM = CM;                                                              // AUC of drug in muscle
dxdt_AUCCF = CF;                                                              // AUC of drug in fat
dxdt_AUCCLu = CLu;                                                            // AUC of drug in lung

double Qbal = QC - QL - QK - QM - QF - QRest;
double Tmass = AA + AV + AL + AK + AM + AF + ALu + ARest + Afeces + Aurine + Acolon + ABileGut + Asitesc + ADOSEsc + AI + AST;
double bal   = ADOSE - Tmass;


$TABLE

capture Plasma = CV*1000;                                                     // ng/mL, multiplied by 1000 to convert mg/L to ng/mL, Venous plasma concentration
capture Liver  = CL*1000;                                                     // ng/mL, multiplied by 1000 to convert mg/L to ng/mL, Liver concentration
capture Kidney = CK*1000;                                                     // ng/mL, multiplied by 1000 to convert mg/L to ng/mL, Kidney concentration
capture Muscle = CM*1000;                                                     // ng/mL, multiplied by 1000 to convert mg/L to ng/mL, Muscle concentration
capture Fat    = CF*1000;                                                     // ng/mL, multiplied by 1000 to convert mg/L to ng/mL, Fat concentration
capture Lung   = CLu*1000;                                                    // ng/mL, multiplied by 1000 to convert mg/L to ng/mL, Lung concentration
capture Rest   = CRest*1000;                                                  // ng/mL, multiplied by 1000 to convert mg/L to ng/mL, Rest of body concentration
capture Bal    = bal;
capture tmass  = Tmass;
'



# Load Model ------------------------------------------------------------------------------------------------------------------------
mod <- mcode_cache("pbpk", SolvePBPK) # refer to mcode function in mrgsolve user guide 3.1.2 Inline

pred<- function (pars, tinterval, BW, TDOSE, route, Dose){
  
MW = 875.1                 # g/mol, molecular weight of ivermectin B1a
DOSE = Dose*BW/MW          # mmol
pars <- exp(pars)          # Get out of log domain

  if (route == "oral admin") {
    ev_oral1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval,
                    addl = TDOSE - 1, cmt  = "AST", replicate = FALSE)
    ev_oral2 <- ev (ID   = 1, amt  = DOSE, ii = tinterval,
                    addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
    ex <- ev_oral1 + ev_oral2 }
  
  if (route == "iv") {
    ev_iv1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, tinf = 0.01,
                  addl = TDOSE - 1, cmt  = "AV", replicate = FALSE)
    ev_iv2 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, tinf = 0.01,
                  addl = TDOSE - 1, cmt  = "ADOSE",  replicate = FALSE)
    ex <- ev_iv1 + ev_iv2 }
  
  if (route == "sc") {
    ev_sc1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval,
                  addl = TDOSE - 1, cmt  = "ADOSEsc",  replicate = FALSE)
    ev_sc2 <- ev (ID   = 1, amt  = DOSE, ii = tinterval,
                  addl = TDOSE - 1, cmt  = "ADOSE",  replicate = FALSE)
    ex <- ev_sc1 + ev_sc2 }
  
# Set up the exposure time
tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 24*55, 0.01)

# Simulation
out <- mod %>% param(pars) %>% mrgsim_d (data = ex, tgrid = tsamp)

output <- cbind.data.frame(Time = out$time,
                           CP  = out$Plasma*MW,
                           CL  = out$Liver*MW,
                           CK  = out$Kidney*MW,
                           CM  = out$Muscle*MW,
                           CF  = out$Fat*MW,
                           CLu = out$Lung*MW,
                           Crest = out$Rest*MW,
                           AUC_V = out$AUCCV*MW,
                           AUC_L = out$AUCCL*MW,
                           AUC_K = out$AUCCK*MW,
                           AUC_M = out$AUCCM*MW,
                           AUC_F = out$AUCCF*MW,
                           AUC_Lu = out$AUCCLu*MW,
                           A_AA      = out$AA,      
                           A_AV      = out$AV,      
                           A_AL      = out$AL,       
                           A_AK      = out$AK,      
                           A_AM      = out$AM,       
                           A_AF      = out$AF,      
                           A_ALu     = out$ALu,     
                           A_ARest   = out$ARest,    
                           A_Aurine  = out$Aurine,  
                           A_Afeces  = out$Afeces,  
                           A_Abile   = out$Abile,
                           Bal = out$Bal)
  return(output)
}



# Balance check --------------------------------------------------------------------------------------------------------------------
pars <- c(
  Ka = 0.035,
  Kst = 0.06,
  Kint = 0.028,
  Ksc = 0.2,
  Kdiss = 0.012,
  KbileC = 0.017,
  Kfeces = 0.012,
  KurineC = 0.0001,
  PL = 17,
  PK = 1,
  PM = 0.65,
  PF = 8.75,
  PLu = 0.5,
  PRest = 1.5
)

out <- pred(pars, tinterval = 24, BW = 655.4, TDOSE = 1, route = "iv", Dose = 0.2)
plot(out$Time, out$Bal, type='l', col='red', ylab='Balance (mg)', xlab='Time (h)', ylim=c(-0.001,0.001))

out <- pred(pars, tinterval = 24, BW = 655.4, TDOSE = 1, route = "oral admin", Dose = 0.2)
plot(out$Time, out$Bal, type='l', col='purple', ylab='Balance (mg)', xlab='Time (h)', ylim=c(-0.001,0.001))

out <- pred(pars, tinterval = 24, BW = 655.4, TDOSE = 1, route = "sc", Dose = 0.2)
plot(out$Time, out$Bal, type='l', col='green', ylab='Balance (mg)', xlab='Time (h)', ylim=c(-0.001,0.001))



# Input dataset -------------------------------------------------------------------------------------------------------------------
  ## IV dataset for calibration
  data1 <- read_excel("study1_0.3 mpk (IV single)_plasma.xlsx", col_types = c("numeric", "numeric", "numeric"))
  data1.P <- cbind.data.frame(Time = data1$Time_hours, CP = data1$Concentration_ppb, SD = data1$SD) #BW: NA
  
  ## Oral dataset for calibration
  data2 <- read_excel("study2_0.3 mpk (Oral single)_plasma.xlsx", col_types = c("numeric", "numeric", "numeric"))
  data2.P <- cbind.data.frame(Time = data2$Time_hours, CP = data2$Concentration_ppb, SD = data2$SD) #BW: 240 - 290 kg, Intraruminal
  
  ## SC dataset for calibration
  data3 <- read_excel("study3_0.2 mpk (SC single)_plasma.xlsx", col_types = c("numeric", "numeric", "numeric"))
  data3.P <- cbind.data.frame(Time = data3$Time_hours, CP = data3$Concentration_ppb, SD = data3$SD) #BW: 201 - 289 kg
  
  #data4 <- read_excel("study4_0.2 mpk (SC single)_plasma.xlsx", col_types = c("numeric", "numeric", "numeric"))  # Holstein
  #data4.P <- cbind.data.frame(Time = data4$Time_hours, CP = data4$Concentration_ppb) #BW: 321 - 447 kg
  
  data5 <- read_excel("study5_0.3 mpk (SC single)_liver.xlsx", col_types = c("numeric", "numeric", "numeric"))
  data5.L <- cbind.data.frame(Time = data5$Time_hours, CL  = data5$Concentration_ppb, SD = data5$SD) #BW: NA
  
  data6 <- read_excel("study5_0.3 mpk (SC single)_kidney.xlsx", col_types = c("numeric", "numeric", "numeric"))
  data6.K <- cbind.data.frame(Time = data6$Time_hours, CK  = data6$Concentration_ppb, SD = data6$SD) #BW: NA
  
  data7 <- read_excel("study5_0.3 mpk (SC single)_muscle.xlsx", col_types = c("numeric", "numeric", "numeric"))
  data7.M <- cbind.data.frame(Time = data7$Time_hours, CM  = data7$Concentration_ppb, SD = data7$SD) #BW: NA
  
  data8 <- read_excel("study5_0.3 mpk (SC single)_fat.xlsx", col_types = c("numeric", "numeric", "numeric"))
  data8.F <- cbind.data.frame(Time = data8$Time_hours, CF  = data8$Concentration_ppb, SD = data8$SD) #BW: NA
  
  data9 <- read_excel("study5_0.3 mpk (SC single)_lung.xlsx", col_types = c("numeric", "numeric", "numeric"))
  data9.Lu <- cbind.data.frame(Time = data9$Time_hours, CLu  = data9$Concentration_ppb, SD = data9$SD) #BW: NA
  
  data10 <- read_excel("study6_0.2 mpk (SC single)_liver_NADA.xlsx", col_types = c("numeric", "numeric", "numeric"))
  data10.L <- cbind.data.frame(Time = data10$Time_hours, CL  = data10$Concentration_ppb, SD = data10$SD) #BW: 287 - 401 kg
  
  
# Define the cost function ---------------------------------------------------------------------------------------------------------
  Cost <- function(pars, w) {
  
  #pars <- exp(pars)
  
  out_A <- pred (pars, tinterval = 24, BW = 655.4, TDOSE = 1, route = "iv", Dose = 0.3)
  out_B <- pred (pars, tinterval = 24, BW = 265, TDOSE = 1, route = "oral admin", Dose = 0.3)
  out_C <- pred (pars, tinterval = 24, BW = 245, TDOSE = 1, route = "sc", Dose = 0.2)
  #out_D <- pred (pars, tinterval = 24, BW = 384, TDOSE = 1, route = "sc", Dose = 0.2)  # Holstein excluded
  out_E <- pred (pars, tinterval = 24, BW = 655.4, TDOSE = 1, route = "sc", Dose = 0.3)
  out_F <- pred (pars, tinterval = 24, BW = 344, TDOSE = 1, route = "sc", Dose = 0.2)
  
  cost <- NULL
  
  cost <- modCost (model = out_A, obs = data1.P[, c("Time", "CP")],, x ="Time", cost = cost, weight = w)
  cost <- modCost (model = out_B, obs = data2.P[, c("Time", "CP")],, x ="Time", cost = cost, weight = w)
  cost <- modCost (model = out_C, obs = data3.P[, c("Time", "CP")],, x ="Time", cost = cost, weight = w)
  #cost <- modCost (model = out_D, obs = data4.P, x ="Time", cost = cost, weight = w)  # Holstein excluded
  cost <- modCost (model = out_E, obs = data5.L[, c("Time", "CL")],, x ="Time", cost = cost, weight = w)
  cost <- modCost (model = out_E, obs = data6.K[, c("Time", "CK")],, x ="Time", cost = cost, weight = w)
  cost <- modCost (model = out_E, obs = data7.M[, c("Time", "CM")],, x ="Time", cost = cost, weight = w)
  cost <- modCost (model = out_E, obs = data8.F[, c("Time", "CF")],, x ="Time", cost = cost, weight = w)
  cost <- modCost (model = out_E, obs = data9.Lu[, c("Time", "CLu")],, x ="Time", cost = cost, weight = w)
  cost <- modCost (model = out_F, obs = data10.L[, c("Time", "CL")],, x ="Time", cost = cost, weight = w)
  
  return(cost)
}

# Check initial cost
Cost(log(pars), w="mean")




# Create the plots for calibration --------------------------------------------------------------------------------------------------
  ## Make the cost and dataframe
Cost_cal <- Cost(log(pars), w = "mean")   

PDat <- cbind.data.frame(OBS = Cost_cal$residuals$obs,
                         PRE = Cost_cal$residuals$mod,
                         RES = Cost_cal$residuals$res)
summary(PDat$OBS)

  ## log-log regression (R-squared value)
PDat <- PDat %>% mutate(Log.OBS = log(OBS, 10), Log.PRE = log(PRE, 10))

fit <- lm(Log.OBS ~ Log.PRE, data = PDat)
summary(fit)

  ## MAPE
MAPE=100*mean(abs((PDat$OBS-PDat$PRE)/PDat$OBS))
summary(MAPE)

  ## RMSE
RMSE <- sqrt(mean((PDat$OBS-PDat$PRE)^2))
summary(RMSE)

  ## AAFE (Average Absolute Fold Error)
AAFE <- 10^(mean(abs(log10(PDat$PRE / PDat$OBS))))
summary(AAFE)

  ## Segments (data$Time, data$low,data$Time, data$high)
plot(PDat$PRE,PDat$OBS)
plot(PDat$RES)

  ## Observed-to-Predicted Ratio (O/P Ratio)
plot1 <- PDat %>% mutate(prediction = predict(fit), OPR = PRE/OBS)


  ## OBS vs PRE Plot (p1.cal) ##########################################################################

# R²
fit <- lm(Log.OBS ~ Log.PRE, data = plot1)
R2 <- summary(fit)$r.squared
R2_label <- paste0("R² = ", round(R2, 2))


# p1.cal
p1.cal <- ggplot(plot1, aes(Log.PRE, Log.OBS)) +
  geom_point(colour = "black", fill = "black", shape = 21, size = 4.8, stroke = 0.2) +
  geom_abline(intercept = 0, slope = 1, colour = "grey38", size = 2) + 

  # 2-fold lines
  geom_abline(intercept =  log10(2), slope = 1, linetype = "dashed", size = 0.8, colour = "grey38") +
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", size = 0.8, colour = "grey38") +
  
  # 3-fold lines
  geom_abline(intercept =  log10(3), slope = 1, linetype = "dotted", size = 0.8, colour = "grey38") +
  geom_abline(intercept = -log10(3), slope = 1, linetype = "dotted", size = 0.8, colour = "grey38") +
  
  scale_y_continuous(limits = c(-1.5,3), labels = scales::math_format(10^.x)) +
  scale_x_continuous(limits = c(-1.5,3), labels = scales::math_format(10^.x)) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", size = 1.5),
    text = element_text(family = "Times"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 22, colour = "black", face = "bold"),
    axis.ticks.length = unit(0.4, "cm"),
    axis.ticks = element_line(size = 1.5, colour = "black"),
    axis.title = element_text(size = 22, colour = "black", face = "bold"),
    legend.position = "none"
  ) +
  labs(x = "Prediction (ng/mL or ng/g)",
       y = "Observation (ng/mL or ng/g)") +
  
  annotate("text", x = -1.2, y = 2.5, label = R2_label, size = 10, colour = "black", hjust = 0)

p1.cal

ggsave("p1.cal.tiff", scale = 1,
       plot = p1.cal,
       width = 22, height = 15, units = "cm", dpi=320)


  ## Calibration result plots ##########################################################################

Cal_OBS <- rbind.data.frame(
  data1.P %>% select(Time = Time, Conc = CP, SD = SD) %>% mutate(STUDY = '(A) Plsama (IV, 0.2 mg/kg)'),
  data2.P %>% select(Time = Time, Conc = CP, SD = SD) %>% mutate(STUDY = '(B) Plasma (PO, 0.3 mg/kg)'),
  data3.P %>% select(Time = Time, Conc = CP, SD = SD) %>% mutate(STUDY = '(C) Plasma (SC, 0.2 mg/kg)'),
  #data4.P was excluded as it is Holstein data
  data5.L %>% select(Time = Time, Conc = CL, SD = SD) %>% mutate(STUDY = '(D) Liver (SC, 0.3 mg/kg)'),
  data6.K %>% select(Time = Time, Conc = CK, SD = SD) %>% mutate(STUDY = '(E) Kidney (SC, 0.3 mg/kg)'),
  data7.M %>% select(Time = Time, Conc = CM, SD = SD) %>% mutate(STUDY = '(F) Muscle (SC, 0.3 mg/kg)'),
  data8.F %>% select(Time = Time, Conc = CF, SD = SD) %>% mutate(STUDY = '(G) Fat (SC, 0.3 mg/kg)'),
  data9.Lu %>% select(Time = Time, Conc = CLu, SD = SD) %>% mutate(STUDY = '(H) Lung (SC, 0.3 mg/kg)'),
  data10.L %>% select(Time = Time, Conc = CL, SD = SD) %>% mutate(STUDY = '(I) Liver (SC, 0.2 mg/kg)')
)

Cal_PRE <- rbind.data.frame(
  out_1  <- pred (pars=log(pars), tinterval = 24, BW = 655.4, TDOSE = 1,  route = "iv", Dose = 0.3) %>% select(Time = Time, Conc = CP) %>% mutate(STUDY = '(A) Plsama (IV, 0.2 mg/kg)'),
  out_2  <- pred (pars=log(pars), tinterval = 24, BW = 265, TDOSE = 1,  route = "oral admin", Dose = 0.3) %>% select(Time = Time, Conc = CP) %>% mutate(STUDY = '(B) Plasma (PO, 0.3 mg/kg)'),
  out_3  <- pred (pars=log(pars), tinterval = 24, BW = 245, TDOSE = 1,  route = "sc", Dose = 0.2) %>% select(Time = Time, Conc = CP) %>% mutate(STUDY = '(C) Plasma (SC, 0.2 mg/kg)'),
  #out_4 was excluded as it is Holstein data
  out_5  <- pred (pars=log(pars), tinterval = 24, BW = 655.4, TDOSE = 1,  route = "sc", Dose = 0.3) %>% select(Time = Time, Conc = CL) %>% mutate(STUDY = '(D) Liver (SC, 0.3 mg/kg)'),
  out_6  <- pred (pars=log(pars), tinterval = 24, BW = 655.4, TDOSE = 1,  route = "sc", Dose = 0.3) %>% select(Time = Time, Conc = CK) %>% mutate(STUDY = '(E) Kidney (SC, 0.3 mg/kg)'),
  out_7  <- pred (pars=log(pars), tinterval = 24, BW = 655.4, TDOSE = 1,  route = "sc", Dose = 0.3) %>% select(Time = Time, Conc = CM) %>% mutate(STUDY = '(F) Muscle (SC, 0.3 mg/kg)'),
  out_8  <- pred (pars=log(pars), tinterval = 24, BW = 655.4, TDOSE = 1,  route = "sc", Dose = 0.3) %>% select(Time = Time, Conc = CF) %>% mutate(STUDY = '(G) Fat (SC, 0.3 mg/kg)'),
  out_9  <- pred (pars=log(pars), tinterval = 24, BW = 655.4, TDOSE = 1,  route = "sc", Dose = 0.3) %>% select(Time = Time, Conc = CLu) %>% mutate(STUDY = '(H) Lung (SC, 0.3 mg/kg)'),
  out_10 <- pred (pars=log(pars), tinterval = 24, BW = 344, TDOSE = 1,  route = "sc", Dose = 0.2) %>% select(Time = Time, Conc = CL) %>% mutate(STUDY = '(I) Liver (SC, 0.2 mg/kg)')
)

Levels <- c('(A) Plsama (IV, 0.2 mg/kg)','(B) Plasma (PO, 0.3 mg/kg)','(C) Plasma (SC, 0.2 mg/kg)','(D) Liver (SC, 0.3 mg/kg)','(E) Kidney (SC, 0.3 mg/kg)','(F) Muscle (SC, 0.3 mg/kg)','(G) Fat (SC, 0.3 mg/kg)','(H) Lung (SC, 0.3 mg/kg)','(I) Liver (SC, 0.2 mg/kg)')

Cal_OBS$Study <- factor(Cal_OBS$STUDY, levels = Levels)
Cal_PRE$Study <- factor(Cal_PRE$STUDY, levels = Levels)



p22 <- ggplot(data = Cal_OBS, aes(x = Time, y = Conc)) + 
  
  geom_errorbar(aes(ymin = 10^(log10(Conc) - (log10(Conc + SD) - log10(Conc))),
                    ymax = 10^(log10(Conc) + (log10(Conc + SD) - log10(Conc)))),
                width = 20,
                colour = "#F0625D",
                size = 0.8,
                na.rm = TRUE) +
  
  geom_point(shape = 21, colour = "#F0625D", fill = "#F0625D", size = 4.8, stroke = 0.2) + 
  geom_line(data = data.frame(Cal_PRE), aes(x = Time, y = Conc), 
            size = 2, colour = "#00529B") +
  scale_y_log10 (labels = function(x) x)  + 
  facet_wrap(~Study, scales = "free", ncol = 3) +
  theme_bw() +
  theme (
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", size = 1.5),
    
    strip.text = element_text(size = 25, face = "bold"),
    
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    strip.background        = element_blank(),
    #strip.text              = element_blank(),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(),
    axis.text               = element_text (size = 22, colour = "black", face = "bold"),    # tick labels along axes
    axis.ticks.length = unit(0.4, "cm"),
    axis.ticks = element_line(size = 1.5, colour = "black"),
    axis.title              = element_text (size = 22, colour = "black", face = "bold"),   # label of axes
    legend.position         ='none') +
  labs (x = "Time (hours)",  y = "Concentration (ng/g)")

p22


ggsave("p22.tiff", scale = 1, plot = p22, width = 45, height = 32, units = "cm", dpi=320)
#individual graph: width =20, height = 15, Whole graph: w=45, h=32







# Model evaluation ------------------------------------------------------------------------------------------------------------------
  ## Input IV dataset
data10 <- read_excel("evaluation_IV_0.2.xlsx", col_types = c("numeric", "numeric", "numeric"))
data10.P <- cbind.data.frame(Time = data10$Time_hours, CP = data10$Concentration_ppb, SD = data10$SD) #BW: 224 - 292 kg

## Input Oral dataset
data11 <- read_excel("evaluation_Oral_0.05.xlsx", col_types = c("numeric", "numeric", "numeric"))
data11.P <- cbind.data.frame(Time = data11$Time_hours, CP = data11$Concentration_ppb, SD = data11$SD) #BW: 163 - 209 kg

data12 <- read_excel("evaluation_Oral_0.1.xlsx", col_types = c("numeric", "numeric", "numeric"))
data12.P <- cbind.data.frame(Time = data12$Time_hours, CP = data12$Concentration_ppb, SD = data12$SD) #BW: 163 - 209 kg

data13 <- read_excel("evaluation_Oral_0.2.xlsx", col_types = c("numeric", "numeric", "numeric"))
data13.P <- cbind.data.frame(Time = data13$Time_hours, CP = data13$Concentration_ppb, SD = data13$SD) #BW: 163 - 209 kg

data14 <- read_excel("evaluation_Oral_0.075.xlsx", col_types = c("numeric", "numeric", "numeric"))
data14.P <- cbind.data.frame(Time = data14$Time_hours, CP = data14$Concentration_ppb, SD = data14$SD) #BW: 200 - 232 kg

## Input SC dataset
data15 <- read_excel("evaluation_SC_0.005.xlsx", col_types = c("numeric", "numeric", "numeric"))
data15.P <- cbind.data.frame(Time = data15$Time_hours, CP = data15$Concentration_ppb, SD = data15$SD) #BW: 163 - 209 kg

data16 <- read_excel("evaluation_SC_0.02.xlsx", col_types = c("numeric", "numeric", "numeric"))
data16.P <- cbind.data.frame(Time = data16$Time_hours, CP = data16$Concentration_ppb, SD = data16$SD) #BW: 163 - 209 kg

data17 <- read_excel("evaluation_SC_0.08.xlsx", col_types = c("numeric", "numeric", "numeric"))
data17.P <- cbind.data.frame(Time = data17$Time_hours, CP = data17$Concentration_ppb, SD = data17$SD) #BW: 163 - 209 kg

data18 <- read_excel("evaluation_SC_0.2_liver.xlsx", col_types = c("numeric", "numeric", "numeric"))
data18.L <- cbind.data.frame(Time = data18$Time_hours, CL = data18$Concentration_ppb, SD = data18$SD) #BW: 297 - 401

data19 <- read_excel("evaluation_SC_0.2_kidney.xlsx", col_types = c("numeric", "numeric", "numeric"))
data19.K <- cbind.data.frame(Time = data19$Time_hours, CK = data19$Concentration_ppb, SD = data19$SD) #BW: 297 - 401

data20 <- read_excel("evaluation_SC_0.2_muscle.xlsx", col_types = c("numeric", "numeric", "numeric"))
data20.M <- cbind.data.frame(Time = data20$Time_hours, CM = data20$Concentration_ppb, SD = data20$SD) #BW: 297 - 401

data21 <- read_excel("evaluation_SC_0.2_fat.xlsx", col_types = c("numeric", "numeric", "numeric"))
data21.F <- cbind.data.frame(Time = data21$Time_hours, CF = data21$Concentration_ppb, SD = data21$SD) #BW: 297 - 401

data22 <- read_excel("evaluation_SC_0.3_liver.xlsx", col_types = c("numeric", "numeric", "numeric"))
data22.L <- cbind.data.frame(Time = data22$Time_hours, CL = data22$Concentration_ppb, SD = data22$SD) #BW: 240 - 290

data23 <- read_excel("evaluation_SC_0.3_kidney.xlsx", col_types = c("numeric", "numeric", "numeric"))
data23.K <- cbind.data.frame(Time = data23$Time_hours, CK = data23$Concentration_ppb, SD = data23$SD) #BW: 240 - 290

data24 <- read_excel("evaluation_SC_0.3_muscle.xlsx", col_types = c("numeric", "numeric", "numeric"))
data24.M <- cbind.data.frame(Time = data24$Time_hours, CM = data24$Concentration_ppb, SD = data24$SD) #BW: 240 - 290

data25 <- read_excel("evaluation_SC_0.3_fat.xlsx", col_types = c("numeric", "numeric", "numeric"))
data25.F <- cbind.data.frame(Time = data25$Time_hours, CF = data25$Concentration_ppb, SD = data25$SD) #BW: 240 - 290

data26 <- read_excel("evaluation_SC_0.3_lung.xlsx", col_types = c("numeric", "numeric", "numeric"))
data26.Lu <- cbind.data.frame(Time = data26$Time_hours, CLu = data26$Concentration_ppb, SD = data26$SD) #BW: 240 - 290


  ## Define the cost function
Cost2 <- function(pars, w) {
  
  #pars <- exp(pars)
  
  out_F <- pred (pars,  tinterval = 24, BW = 258, TDOSE = 1, route = "iv", Dose = 0.2)
  out_G_0.05 <- pred (pars,  tinterval = 24, BW = 186, TDOSE = 28, route = "oral admin", Dose = 0.05)  #capsule
  out_G_0.1 <- pred (pars,  tinterval = 24, BW = 186, TDOSE = 28, route = "oral admin", Dose = 0.1)  #capsule
  out_G_0.2 <- pred (pars,  tinterval = 24, BW = 186, TDOSE = 28, route = "oral admin", Dose = 0.2)  #capsule
  out_H <- pred (pars,  tinterval = 72, BW = 216, TDOSE = 8, route = "oral admin", Dose = 0.075)   #capsule
  out_I_0.005 <- pred (pars,  tinterval = 24, BW = 186, TDOSE = 28, route = "sc", Dose = 0.005)
  out_I_0.02 <- pred (pars,  tinterval = 24, BW = 186, TDOSE = 28, route = "sc", Dose = 0.02)
  out_I_0.08 <- pred (pars,  tinterval = 24, BW = 186, TDOSE = 28, route = "sc", Dose = 0.08)
  out_J <- pred (pars,  tinterval = 24, BW = 349, TDOSE = 1, route = "sc", Dose = 0.2)
  out_K <- pred (pars,  tinterval = 24, BW = 265, TDOSE = 1, route = "sc", Dose = 0.3)

  cost <- NULL
  
  cost <- modCost (model = out_F, obs = data10.P[, c("Time", "CP")],, x ="Time", cost = cost, weight = w)
  cost <- modCost (model = out_G_0.05, obs = data11.P[, c("Time", "CP")],, x ="Time", cost = cost, weight = w)
  cost <- modCost (model = out_G_0.1, obs = data12.P[, c("Time", "CP")],, x ="Time", cost = cost, weight = w)  
  cost <- modCost (model = out_G_0.2, obs = data13.P[, c("Time", "CP")],, x ="Time", cost = cost, weight = w)
  cost <- modCost (model = out_H, obs = data14.P[, c("Time", "CP")],, x ="Time", cost = cost, weight = w)
  cost <- modCost (model = out_I_0.005, obs = data15.P[, c("Time", "CP")],, x ="Time", cost = cost, weight = w) 
  cost <- modCost (model = out_I_0.02, obs = data16.P[, c("Time", "CP")],, x ="Time", cost = cost, weight = w)  
  cost <- modCost (model = out_I_0.08, obs = data17.P[, c("Time", "CP")],, x ="Time", cost = cost, weight = w) 
  cost <- modCost (model = out_J, obs = data18.L[, c("Time", "CL")],, x ="Time", cost = cost, weight = w) 
  cost <- modCost (model = out_J, obs = data19.K[, c("Time", "CK")],, x ="Time", cost = cost, weight = w) 
  cost <- modCost (model = out_J, obs = data20.M[, c("Time", "CM")],, x ="Time", cost = cost, weight = w) 
  cost <- modCost (model = out_J, obs = data21.F[, c("Time", "CF")],, x ="Time", cost = cost, weight = w) 
  cost <- modCost (model = out_K, obs = data22.L[, c("Time", "CL")],, x ="Time", cost = cost, weight = w) 
  cost <- modCost (model = out_K, obs = data23.K[, c("Time", "CK")],, x ="Time", cost = cost, weight = w) 
  cost <- modCost (model = out_K, obs = data24.M[, c("Time", "CM")],, x ="Time", cost = cost, weight = w) 
  cost <- modCost (model = out_K, obs = data25.F[, c("Time", "CF")],, x ="Time", cost = cost, weight = w) 
  cost <- modCost (model = out_K, obs = data26.Lu[, c("Time", "CLu")],, x ="Time", cost = cost, weight = w)
  
  return(cost)
}

Cost2(log(pars), w = "mean")


  ## Make the cost and dataframe
Cost_eval <- Cost2(log(pars), w = "mean")

PDat_eval <- cbind.data.frame (OBS1 = Cost_eval$residuals$obs,
                               PRE1 = Cost_eval$residuals$mod,
                               RES1 = Cost_eval$residuals$res)

summary(PDat_eval$OBS)


  ## log-log regression (R-squared value)
PDat_eval <- PDat_eval %>% mutate(Log.OBS1 = log(OBS1, 10), Log.PRE1 = log(PRE1, 10), Species = "Cattle")

fit_eval <- lm(Log.OBS1 ~ Log.PRE1, data = PDat_eval)
summary(fit_eval)

  ## MAPE
MAPE=100*mean(abs((PDat_eval$OBS1-PDat_eval$PRE1)/PDat_eval$OBS1))
summary(MAPE)

  ## RMSE
RMSE <- sqrt(mean((PDat_eval$OBS1-PDat_eval$PRE1)^2))
summary(RMSE)

  ## AAFE (Average Absolute Fold Error)
AAFE <- 10^(mean(abs(log10(PDat_eval$PRE1 / PDat_eval$OBS1))))
summary(AAFE)

  ## 2-, 3-, and 5-fold range
folds <- c(2, 3, 5) 
for(f in folds){
  within_fold  <- sum(PDat_eval$PRE1 / PDat_eval$OBS1 <= f & PDat_eval$PRE1 / PDat_eval$OBS1 >= 1/f)
  total_points <- nrow(PDat_eval)
  percent      <- within_fold / total_points * 100
  cat(sprintf("%d-fold: %d out of %d (%.2f%%)\n", f, within_fold, total_points, percent))
}


  ## Observed-to-Predicted Ratio (O/P Ratio)
plot2 <- PDat_eval %>% mutate(prediction = predict(fit_eval), OPR = PRE1/OBS1)

  ## OBS vs PRE Plot (p1.eval) ##########################################################################

# R²
fit <- lm(Log.OBS1 ~ Log.PRE1, data = plot2)
R2 <- summary(fit_eval)$r.squared
R2_label <- paste0("R² = ", round(R2, 2))  


# p1.eval
p1.eval <- ggplot(plot2, aes(Log.PRE1, Log.OBS1)) +
  geom_point(colour = "black", fill = "black", shape = 21, size = 4.8, stroke = 0.2) +
  geom_abline(intercept = 0, slope = 1, colour = "grey38", size = 2) + 

  # 2-fold lines
  geom_abline(intercept =  log10(2), slope = 1, linetype = "dashed", size = 0.8, colour = "grey38") +
  geom_abline(intercept = -log10(2), slope = 1, linetype = "dashed", size = 0.8, colour = "grey38") +
  
  # 3-fold lines
  geom_abline(intercept =  log10(3), slope = 1, linetype = "dotted", size = 0.8, colour = "grey38") +
  geom_abline(intercept = -log10(3), slope = 1, linetype = "dotted", size = 0.8, colour = "grey38") +
  
  scale_y_continuous(limits = c(-1.5,3), labels = scales::math_format(10^.x)) +
  scale_x_continuous(limits = c(-1.5,3), labels = scales::math_format(10^.x)) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", size = 1.5),
    text = element_text(family = "Times"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 22, colour = "black", face = "bold"),
    axis.ticks.length = unit(0.4, "cm"),
    axis.ticks = element_line(size = 1.5, colour = "black"),
    axis.title = element_text(size = 22, colour = "black", face = "bold"),
    legend.position = "none"
  ) +
  labs(x = "Prediction (ng/mL or ng/g)",
       y = "Observation (ng/mL or ng/g)") +
  
  annotate("text", x = -1.2, y = 2.5, label = R2_label, size = 10, colour = "black", hjust = 0)

p1.eval

ggsave("p1.eval.tiff", scale = 1,
       plot = p1.eval,
       width = 22, height = 15, units = "cm", dpi=320)


  ## Evaluation result plots ##########################################################################

Eval_OBS <- rbind.data.frame(
  data10.P %>% select(Time = Time, Conc = CP, SD = SD) %>% mutate(STUDY = '(A2) Plasma (IV, 0.2 mg/kg)'),
  data11.P %>% select(Time = Time, Conc = CP, SD = SD) %>% mutate(STUDY = 'Oral_0.05_Plasma'),
  data12.P %>% select(Time = Time, Conc = CP, SD = SD) %>% mutate(STUDY = 'Oral_0.1_Plamsa'),
  data13.P %>% select(Time = Time, Conc = CP, SD = SD) %>% mutate(STUDY = 'Oral_0.2_Plamsa'),
  data14.P %>% select(Time = Time, Conc = CP, SD = SD) %>% mutate(STUDY = '(B2) Plasma (PO, 0.075 mg/kg, 7 doses)'),
  data15.P %>% select(Time = Time, Conc = CP, SD = SD) %>% mutate(STUDY = '(C2) Plasma (SC, 0.005 mg/kg, 28 doses)'),
  data16.P %>% select(Time = Time, Conc = CP, SD = SD) %>% mutate(STUDY = 'SC_0.02_Plasma'),
  data17.P %>% select(Time = Time, Conc = CP, SD = SD) %>% mutate(STUDY = 'SC_0.08_Plasma'),
  data18.L %>% select(Time = Time, Conc = CL, SD = SD) %>% mutate(STUDY = 'SC_0.2_Liver'),
  data19.K %>% select(Time = Time, Conc = CK, SD = SD) %>% mutate(STUDY = 'SC_0.2_Kidney'),
  data20.M %>% select(Time = Time, Conc = CM, SD = SD) %>% mutate(STUDY = 'SC_0.2_Muscle'),
  data21.F %>% select(Time = Time, Conc = CF, SD = SD) %>% mutate(STUDY = 'SC_0.2_Fat'),
  data22.L %>% select(Time = Time, Conc = CL, SD = SD) %>% mutate(STUDY = '(D2) Liver (SC, 0.3 mg/kg)'),
  data23.K %>% select(Time = Time, Conc = CK, SD = SD) %>% mutate(STUDY = '(E2) Kidney (SC, 0.3 mg/kg)'),
  data24.M %>% select(Time = Time, Conc = CM, SD = SD) %>% mutate(STUDY = '(F2) Muscle (SC, 0.3 mg/kg)'),
  data25.F %>% select(Time = Time, Conc = CF, SD = SD) %>% mutate(STUDY = '(G2) Fat (SC, 0.3 mg/kg)'),
  data26.Lu %>% select(Time = Time, Conc = CLu, SD = SD) %>% mutate(STUDY = '(H2) Lung (SC, 0.3 mg/kg)')
)

Eval_PRE <- rbind.data.frame(
  out_10  <- pred (pars=log(pars), tinterval = 24, BW = 258, TDOSE = 1,  route = "iv", Dose = 0.2) %>% select(Time = Time, Conc = CP) %>% mutate(STUDY = '(A2) Plasma (IV, 0.2 mg/kg)'),
  out_11  <- pred (pars=log(pars), tinterval = 24, BW = 186, TDOSE = 28,  route = "oral admin", Dose = 0.05) %>% select(Time = Time, Conc = CP) %>% mutate(STUDY = 'Oral_0.05_Plasma'),
  out_12  <- pred (pars=log(pars), tinterval = 24, BW = 186, TDOSE = 28,  route = "oral admin", Dose = 0.1) %>% select(Time = Time, Conc = CP) %>% mutate(STUDY = 'Oral_0.1_Plamsa'),
  out_13  <- pred (pars=log(pars), tinterval = 24, BW = 186, TDOSE = 28,  route = "oral admin", Dose = 0.2) %>% select(Time = Time, Conc = CP) %>% mutate(STUDY = 'Oral_0.2_Plamsa'),
  out_14  <- pred (pars=log(pars), tinterval = 72, BW = 216, TDOSE = 8,  route = "oral admin", Dose = 0.075) %>% select(Time = Time, Conc = CP) %>% mutate(STUDY = '(B2) Plasma (PO, 0.075 mg/kg, 7 doses)'),
  out_15  <- pred (pars=log(pars), tinterval = 24, BW = 186, TDOSE = 28,  route = "sc", Dose = 0.005) %>% select(Time = Time, Conc = CP) %>% mutate(STUDY = '(C2) Plasma (SC, 0.005 mg/kg, 28 doses)'),
  out_16  <- pred (pars=log(pars), tinterval = 24, BW = 186, TDOSE = 28,  route = "sc", Dose = 0.02) %>% select(Time = Time, Conc = CP) %>% mutate(STUDY = 'SC_0.02_Plasma'),  
  out_17  <- pred (pars=log(pars), tinterval = 24, BW = 186, TDOSE = 28,  route = "sc", Dose = 0.08) %>% select(Time = Time, Conc = CP) %>% mutate(STUDY = 'SC_0.08_Plasma'),
  out_18  <- pred (pars=log(pars), tinterval = 24, BW = 349, TDOSE = 1,  route = "sc", Dose = 0.2) %>% select(Time = Time, Conc = CL) %>% mutate(STUDY = 'SC_0.2_Liver'),
  out_19  <- pred (pars=log(pars), tinterval = 24, BW = 349, TDOSE = 1,  route = "sc", Dose = 0.2) %>% select(Time = Time, Conc = CK) %>% mutate(STUDY = 'SC_0.2_Kidney'),
  out_20  <- pred (pars=log(pars), tinterval = 24, BW = 349, TDOSE = 1,  route = "sc", Dose = 0.2) %>% select(Time = Time, Conc = CM) %>% mutate(STUDY = 'SC_0.2_Muscle'),
  out_21  <- pred (pars=log(pars), tinterval = 24, BW = 349, TDOSE = 1,  route = "sc", Dose = 0.2) %>% select(Time = Time, Conc = CF) %>% mutate(STUDY = 'SC_0.2_Fat'),
  out_22  <- pred (pars=log(pars), tinterval = 24, BW = 265, TDOSE = 1,  route = "sc", Dose = 0.3) %>% select(Time = Time, Conc = CL) %>% mutate(STUDY = '(D2) Liver (SC, 0.3 mg/kg)'),
  out_23  <- pred (pars=log(pars), tinterval = 24, BW = 265, TDOSE = 1,  route = "sc", Dose = 0.3) %>% select(Time = Time, Conc = CK) %>% mutate(STUDY = '(E2) Kidney (SC, 0.3 mg/kg)'),
  out_24  <- pred (pars=log(pars), tinterval = 24, BW = 265, TDOSE = 1,  route = "sc", Dose = 0.3) %>% select(Time = Time, Conc = CM) %>% mutate(STUDY = '(F2) Muscle (SC, 0.3 mg/kg)'),
  out_25  <- pred (pars=log(pars), tinterval = 24, BW = 265, TDOSE = 1,  route = "sc", Dose = 0.3) %>% select(Time = Time, Conc = CF) %>% mutate(STUDY = '(G2) Fat (SC, 0.3 mg/kg)'),
  out_26  <- pred (pars=log(pars), tinterval = 24, BW = 265, TDOSE = 1,  route = "sc", Dose = 0.3) %>% select(Time = Time, Conc = CLu) %>% mutate(STUDY = '(H2) Lung (SC, 0.3 mg/kg)')
)

Levels <- c('(A2) Plasma (IV, 0.2 mg/kg)','Oral_0.05_Plasma','Oral_0.1_Plamsa','Oral_0.2_Plamsa','(B2) Plasma (PO, 0.075 mg/kg, 7 doses)','(C2) Plasma (SC, 0.005 mg/kg, 28 doses)','SC_0.02_Plasma','SC_0.08_Plasma','SC_0.2_Plasma',
            'SC_0.45_Plasma','SC_0.63_Plasma','SC_0.2_Liver','SC_0.2_Kidney','SC_0.2_Muscle','SC_0.2_Fat','SC_0.3_Plasma','(D2) Liver (SC, 0.3 mg/kg)','(E2) Kidney (SC, 0.3 mg/kg)','(F2) Muscle (SC, 0.3 mg/kg)',
            '(G2) Fat (SC, 0.3 mg/kg)','(H2) Lung (SC, 0.3 mg/kg)')

Eval_OBS$Study <- factor(Eval_OBS$STUDY, levels = Levels)
Eval_PRE$Study <- factor(Eval_PRE$STUDY, levels = Levels)



p222 <- ggplot(data = Eval_OBS, aes(x = Time, y = Conc)) + 
  
  geom_errorbar(aes(ymin = 10^(log10(Conc) - (log10(Conc + SD) - log10(Conc))),
                    ymax = 10^(log10(Conc) + (log10(Conc + SD) - log10(Conc)))),
                width = 20,
                colour = "#F0625D",
                size = 0.8,
                na.rm = TRUE) +
  
  geom_point(shape = 21, colour = "#F0625D", fill = "#F0625D", size = 4.8, stroke = 0.2) + 
  geom_line(data = data.frame(Eval_PRE), aes(x = Time, y = Conc), 
            size = 2, colour = "#00529B") +
  scale_y_log10 (labels = function(x) x)  + 
  facet_wrap(~Study, scales = "free", ncol = 4) +
  theme_bw() +
  theme (
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", size = 1.5),
    
    strip.text = element_text(size = 25, face = "bold"),  # Title size
    
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    strip.background        = element_blank(),
    #strip.text              = element_blank(),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(),
    axis.text               = element_text (size = 22, colour = "black", face = "bold"),    # tick labels along axes
    axis.ticks.length = unit(0.4, "cm"),
    axis.ticks = element_line(size = 1.5, colour = "black"),
    axis.title              = element_text (size = 22, colour = "black", face = "bold"),   # label of axes
    legend.position         ='none') +
  labs (x = "Time (hours)",  y = "Concentration (ng/g)")

p222

ggsave("p222.tiff", scale = 1, plot = p222, width = 85, height = 70, units = "cm", dpi=320)
#individual graph: width =20, height = 15;   Whole graph: w=45, h=32




# Save the fitting results -----------------------------------------------------------------------------------------------------------
saveRDS(Cost_eval, file = 'Cattle_IVE_evaluation.rds')







