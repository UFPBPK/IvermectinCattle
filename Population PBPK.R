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
  library(knitr)
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
PL         : 17         : Unitless,           Liver:plasma PC // 15.69
PK         : 1          : Unitless,           Kidney:plasma PC 
PM         : 0.65       : Unitless,           Muscle:plasma PC
PF         : 8.75       : Unitless,           Fat:plasma PC; 5.75 for human model from Yeo & Wesche, 2023
PLu        : 0.5        : Unitless,           Lung:plasma PC
PRest      : 1.5        : Unitless,           Rest of the body:plasma PC // 1
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

$CMT AA AV AIV ADOSE ADOSEsc Aabsorbsc Alocalsc AST AI AAO AL AK AM AF ALu ARest AUCCV AUCCL AUCCK AUCCM AUCCF AUCCLu Aurine Afeces Abile Acolon ABileGut

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
double Rsc = Ksc * Alocalsc;                                                  // Subcutaneous absorption rate
double Rdepotsc = - Kdiss * ADOSEsc;                                          // Drug amount change at the injection-site depot
double Rlocalsc = Kdiss * ADOSEsc - Rsc;                                      // Drug amount change at the local SC tissue
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
dxdt_ADOSEsc = Rdepotsc;                                                      // Amount of drug in SC injection-site depot
dxdt_Alocalsc = Rlocalsc;                                                     // Amount of drug in local SC tissue
dxdt_Aabsorbsc = Rsc;                                                         // Amount of drug in SC absorption
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
double Tmass = AA + AV + AL + AK + AM + AF + ALu + ARest + Afeces + Aurine + Acolon + ABileGut + Alocalsc + ADOSEsc + AI + AST;
double bal   = ADOSE - Tmass;
$TABLE

capture Plasma = CV;                                                          // mmol/L, Venous plasma concentration
capture Liver  = CL;                                                          // mmol/L, Liver concentration
capture Kidney = CK;                                                          // mmol/L, Kidney concentration
capture Muscle = CM;                                                          // mmol/L, Muscle concentration
capture Fat    = CF;                                                          // mmol/L, Fat concentration
capture Lung   = CLu;                                                         // mmol/L, Lung concentration
capture Rest   = CRest;                                                       // mmol/L, Rest of body concentration
capture Bal    = bal;
capture tmass  = Tmass;
'



# Load Model ------------------------------------------------------------------------------------------------------------------------
mod <- mcode_cache("pbpk", SolvePBPK) # refer to mcode function in mrgsolve user guide 3.1.2 Inline



## Set Up Simulation Subjects
## There will be N simulation subjects, with their own sampled parameters that will override the default values.  The PBPK model does not use the Omega block to simulate invidual values, instead we fill the idata dataframe with these values.  Here we only change body weight, but other parameters could also be varied.

N = 1000
ndoses = 1
doseinterval = 24
MW = 875.1            # g/mol molecular weight of ivermectin B1a
End_time = 55         # days
PDOSEsc = 0.2         # mg/kg bw  # Labeled = 0.2 mg/kg; Extralabel = 0.3 mg/kg
PDOSEoral = 0.2       # mg/kg bw  # Labeled = 0.2 mg/kg; Extralabel = 0.3 mg/kg

## Chemical-specific parameters  *0.2
sd.log.PL = sqrt(log(1 + 3.4 ^ 2 / 17 ^ 2))
m.log.PL = log(17) - 0.5 * sd.log.PL ^ 2

sd.log.PK = sqrt(log(1 + 0.2 ^ 2 / 1 ^ 2))
m.log.PK = log(1) - 0.5 * sd.log.PK ^ 2

sd.log.PM = sqrt(log(1 + 0.13 ^ 2 / 0.65 ^ 2))
m.log.PM = log(0.65) - 0.5 * sd.log.PM ^ 2

sd.log.PF = sqrt(log(1 + 1.75 ^ 2 / 8.75 ^ 2))
m.log.PF = log(8.75) - 0.5 * sd.log.PF ^ 2

sd.log.PLu = sqrt(log(1 + 0.1 ^ 2 / 0.5 ^ 2))
m.log.PLu = log(0.5) - 0.5 * sd.log.PLu ^ 2

sd.log.PRest = sqrt(log(1 + 0.3 ^ 2 / 1.5 ^ 2))
m.log.PRest = log(1.5) - 0.5 * sd.log.PRest ^ 2

sd.log.Ka = sqrt(log(1 + 0.0105 ^ 2 / 0.035 ^ 2)) # *0.3
m.log.Ka = log(0.035) - 0.5 * sd.log.Ka ^ 2

sd.log.Kst = sqrt(log(1 + 0.018 ^ 2 / 0.06 ^ 2))
m.log.Kst = log(0.06) - 0.5 * sd.log.Kst ^ 2

sd.log.Kint = sqrt(log(1 + 0.0084 ^ 2 / 0.028 ^ 2))
m.log.Kint = log(0.028) - 0.5 * sd.log.Kint ^ 2

sd.log.Ksc = sqrt(log(1 + 0.06 ^ 2 / 0.2 ^ 2))
m.log.Ksc = log(0.2) - 0.5 * sd.log.Ksc ^ 2

sd.log.Kdiss = sqrt(log(1 + 0.0036 ^ 2 / 0.012 ^ 2))
m.log.Kdiss = log(0.012) - 0.5 * sd.log.Kdiss ^ 2

sd.log.KbileC = sqrt(log(1 + 0.0051 ^ 2 / 0.017 ^ 2))
m.log.KbileC = log(0.017) - 0.5 * sd.log.KbileC ^ 2

sd.log.Kfeces = sqrt(log(1 + 0.0036 ^ 2 / 0.012 ^ 2))
m.log.Kfeces = log(0.012) - 0.5 * sd.log.Kfeces ^ 2

sd.log.KurineC = sqrt(log(1 + 0.00003 ^ 2 / 0.0001 ^ 2))
m.log.KurineC = log(0.0001) - 0.5 * sd.log.KurineC ^ 2

sd.log.PB = sqrt(log(1 + 0.093 ^ 2 / 0.93 ^ 2))
m.log.PB = log(0.93) - 0.5 * sd.log.PB ^ 2



## Physiological parameters
set.seed(42)
idata <- 
  data_frame(ID=1:N) %>% 
  mutate(
    BW = rtruncnorm(n = N, 
                    a = qnorm(0.025, mean = 655.4, sd = 0.157*655.4), 
                    b = qnorm(0.975, mean = 655.4, sd = 0.157*655.4), 
                    mean = 655.4,
                    sd = 0.157*655.4),  # CV = 15.7 (%) (from Table 6, Lin et al., 2020)
    
    # Fractional blood flow to the tissue  
    QCC = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 5.45, sd = 1.47),
      b = qnorm(0.975, mean = 5.45, sd = 1.47),
      mean = 5.45,
      sd = 1.47),  # SD = 1.47 (L/h/kg bw) (from Table 13, Lin et al., 2020)
    
    QLCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.44, sd = 0.44*0.3),
      b = qnorm(0.975, mean = 0.44, sd = 0.44*0.3),
      mean = 0.44,
      sd = 0.44*0.3),  # SD = 25 (%) (from Table 21, Lin et al., 2020), however, setting it this way may produce negative values, CV of 30% was applied instead.
    
    QKCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.11, sd = 0.11*0.3),
      b = qnorm(0.975, mean = 0.11, sd = 0.11*0.3),
      mean = 0.11,
      sd = 0.11*0.3),  # SD = 8 (%) (from Table 21, Lin et al., 2020), however, setting it this way may produce negative values, CV of 30% was applied instead.
    
    QMCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.28, sd = 0.09),
      b = qnorm(0.975, mean = 0.28, sd = 0.09),
      mean = 0.28,
      sd = 0.09),  # SD = 9 (%) (from Table 21, Lin et al., 2020)
    
    QFCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.08, sd = 0.024),
      b = qnorm(0.975, mean = 0.08, sd = 0.024),
      mean = 0.08,
      sd = 0.024),  # SD = 2.400E-02 (unitless) (from Table S2, Li et al., 2019)
    
    QLuCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 1, sd = 1*0.3),
      b = qnorm(0.975, mean = 1, sd = 1*0.3),
      mean = 1,
      sd = 1*0.3),  # CV = 30%
    
    #QRestCa = pmax(1 - QLCa - QKCa - QMCa - QFCa, 0.01),  
        QRestCa = rtruncnorm(
          n = N,
          a = qnorm(0.025, mean = 0.09, sd = 0.09*0.3),
          b = qnorm(0.975, mean = 0.09, sd = 0.09*0.3),
          mean = 0.09,
          sd = 0.09*0.3),  # CV = 30%
    
    # Fractional organ tissue volumes
    VLCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.0122, sd = 0.0018),
      b = qnorm(0.975, mean = 0.0122, sd = 0.0018) ,
      mean = 0.0122,
      sd = 0.0018),  # SD = 0.18 (% BW) (from Table 2, Lin et al., 2020)
    
    VKCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.0021, sd = 0.0005),
      b = qnorm(0.975, mean = 0.0021, sd = 0.0005),
      mean = 0.0021,
      sd = 0.0005),  # SD = 0.05 (% BW) (from Table 2, Lin et al., 2020)
    
    VMCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.361, sd = 0.1173),
      b = qnorm(0.975, mean = 0.361, sd = 0.1173),
      mean = 0.361,
      sd = 0.1173),  # SD = 11.73 (% BW) (from Table 2, Lin et al., 2020)
    
    VFCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.1218, sd = 0.0506),
      b = qnorm(0.975, mean = 0.1218, sd = 0.0506),
      mean = 0.1218,
      sd = 0.0506),  # SD = 5.06 (% BW) (from Table 2, Lin et al., 2020)
    
    VLuCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.0085, sd = 0.0025),
      b = qnorm(0.975, mean = 0.0085, sd = 0.0025),
      mean = 0.0085,
      sd = 0.0025),  # SD = 0.25 (% BW) (from Table 2, Lin et al., 2020)
    
    VvenCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.01737, sd = 0.002961),
      b = qnorm(0.975, mean = 0.01737, sd = 0.002961),
      mean = 0.01737,
      sd = 0.002961),  # Hct=0.378, 70% of plasma blood (from Table 2, Lin et al., 2020)
    
    VartCa = rtruncnorm(
      n = N,
      a = qnorm(0.025, mean = 0.00745, sd = 0.001269),
      b = qnorm(0.975, mean = 0.00745, sd = 0.001269),
      mean = 0.00745,
      sd = 0.001269),  # Hct=0.378, 30% of plasma blood (from Table 2, Lin et al., 2020)
    
    #VRestCa = pmax(1 - VLCa - VKCa - VMCa - VFCa - VLuCa - VvenCa - VartCa, 0.01),
        VRestCa = rtruncnorm(
          n = N,
          a = qnorm(0.025, mean = 0.4696, sd = 0.4696*0.3),
          b = qnorm(0.975, mean = 0.4696, sd = 0.4696*0.3),
          mean = 0.4696,
          sd = 0.4696*0.3),
    
    
    
    ## Generate individual parameters for Monte Carlo simulation
    PL = rlnormTrunc(N, meanlog = m.log.PL, sdlog = sd.log.PL, 
                      min = qlnorm(0.025, meanlog = m.log.PL, sdlog = sd.log.PL), 
                      max = qlnorm(0.975, meanlog = m.log.PL, sdlog = sd.log.PL)),
    
    PK = rlnormTrunc(N, meanlog = m.log.PK, sdlog = sd.log.PK, 
                      min = qlnorm(0.025, meanlog = m.log.PK, sdlog = sd.log.PK), 
                      max = qlnorm(0.975, meanlog = m.log.PK, sdlog = sd.log.PK)),
    
    PM = rlnormTrunc(N, meanlog = m.log.PM, sdlog = sd.log.PM, 
                      min = qlnorm(0.025, meanlog = m.log.PM, sdlog = sd.log.PM), 
                      max = qlnorm(0.975, meanlog = m.log.PM, sdlog = sd.log.PM)),
    
    PF = rlnormTrunc(N, meanlog = m.log.PF, sdlog = sd.log.PF, 
                      min = qlnorm(0.025, meanlog = m.log.PF, sdlog = sd.log.PF), 
                      max = qlnorm(0.975, meanlog = m.log.PF, sdlog = sd.log.PF)),
    
    PLu = rlnormTrunc(N, meanlog = m.log.PLu, sdlog = sd.log.PLu, 
                       min = qlnorm(0.025, meanlog = m.log.PLu, sdlog = sd.log.PLu), 
                       max = qlnorm(0.975, meanlog = m.log.PLu, sdlog = sd.log.PLu)),
    
    PRest = rlnormTrunc(N, meanlog = m.log.PRest, sdlog = sd.log.PRest, 
                         min = qlnorm(0.025, meanlog = m.log.PRest, sdlog = sd.log.PRest), 
                         max = qlnorm(0.975, meanlog = m.log.PRest, sdlog = sd.log.PRest)),
    
    Ka = rlnormTrunc(N, meanlog = m.log.Ka, sdlog = sd.log.Ka, 
                      min = qlnorm(0.025, meanlog = m.log.Ka, sdlog = sd.log.Ka), 
                      max = qlnorm(0.975, meanlog = m.log.Ka, sdlog = sd.log.Ka)),
    
    Kst = rlnormTrunc(N, meanlog = m.log.Kst, sdlog = sd.log.Kst, 
                       min = qlnorm(0.025, meanlog = m.log.Kst, sdlog = sd.log.Kst), 
                       max = qlnorm(0.975, meanlog = m.log.Kst, sdlog = sd.log.Kst)),
    
    Kint = rlnormTrunc(N, meanlog = m.log.Kint, sdlog = sd.log.Kint, 
                        min = qlnorm(0.025, meanlog = m.log.Kint, sdlog = sd.log.Kint), 
                        max = qlnorm(0.975, meanlog = m.log.Kint, sdlog = sd.log.Kint)),
    
    Ksc = rlnormTrunc(N, meanlog = m.log.Ksc, sdlog = sd.log.Ksc, 
                       min = qlnorm(0.025, meanlog = m.log.Ksc, sdlog = sd.log.Ksc), 
                       max = qlnorm(0.975, meanlog = m.log.Ksc, sdlog = sd.log.Ksc)),
    
    Kdiss = rlnormTrunc(N, meanlog = m.log.Kdiss, sdlog = sd.log.Kdiss, 
                         min = qlnorm(0.025, meanlog = m.log.Kdiss, sdlog = sd.log.Kdiss), 
                         max = qlnorm(0.975, meanlog = m.log.Kdiss, sdlog = sd.log.Kdiss)),
    
    KbileC = rlnormTrunc(N, meanlog = m.log.KbileC, sdlog = sd.log.KbileC,
                         min = qlnorm(0.025, meanlog = m.log.KbileC, sdlog = sd.log.KbileC),
                         max = qlnorm(0.975, meanlog = m.log.KbileC, sdlog = sd.log.KbileC)),
    
    Kfeces = rlnormTrunc(N, meanlog = m.log.Kfeces, sdlog = sd.log.Kfeces,
                          min = qlnorm(0.025, meanlog = m.log.Kfeces, sdlog = sd.log.Kfeces),
                          max = qlnorm(0.975, meanlog = m.log.Kfeces, sdlog = sd.log.Kfeces)),
    
    KurineC = rlnormTrunc(N, meanlog = m.log.KurineC, sdlog = sd.log.KurineC,
                          min = qlnorm(0.025, meanlog = m.log.KurineC, sdlog = sd.log.KurineC),
                          max = qlnorm(0.975, meanlog = m.log.KurineC, sdlog = sd.log.KurineC)),
    
    PB = rlnormTrunc(N, meanlog = m.log.PB, sdlog = sd.log.PB,
                     min = 0.5,
                     max = 1),
    
    
    
    ## Set up dosing
    DOSEsc = PDOSEsc*BW, DOSEoral = PDOSEoral*BW)
idata %>% head %>% kable

ev_sc1 <- ev(ID = 1:N, amt = idata$DOSEsc, ii = doseinterval, addl = ndoses-1, cmt = "ADOSEsc", replicate = FALSE)
ev_sc2 <- ev(ID = 1:N, amt = idata$DOSEsc, ii = doseinterval, addl = ndoses-1, cmt = "ADOSE", replicate = FALSE)
ex_sc <- ev_sc1 + ev_sc2

dim(as.data.frame(ex_sc))

ev_oral1 <- ev(ID = 1:N, amt = idata$DOSEoral, ii = doseinterval, addl = ndoses-1, cmt = "AST", replicate = FALSE)
ev_oral2 <- ev(ID = 1:N, amt = idata$DOSEoral, ii = doseinterval, addl = ndoses-1, cmt = "ADOSE", replicate = FALSE)
ex_oral <- ev_oral1 + ev_oral2

dim(as.data.frame(ex_oral))



## set up time grid
tsamp = tgrid(0, doseinterval*ndoses + 55*24, 1)  # Sample every hour for 55 days after the last dose



## Combine data and run the simulation (SC) --------------------------------------------------------------------------------------------------
set.seed(11)
system.time({
  out_sc <- mod %>% data_set(ex_sc) %>% update(atol= 1E-10, rtol=1E-5, maxsteps=50000, hmax = 0.01) %>% idata_set(idata) %>%
    mrgsim(obsonly=TRUE, tgrid=tsamp) 
}) -> run.time_sc

out.sum_sc <- out_sc %>% as.data.frame %>%
  mutate(Time1 = time/24 - (ndoses-1)*doseinterval/24) %>%
  dplyr::select(ID, Time1, Plasma:Fat) %>%
  arrange(Time1)  %>%
  group_by(Time1) %>%
  dplyr::summarise(CV1 = quantile(Plasma, 0.01), CV50 = median(Plasma), CV99 = quantile(Plasma, 0.99), 
                   CL1 = quantile(Liver, 0.01), CL50 = median(Liver), CL95 = quantile(Liver, 0.95), CL99 = quantile(Liver, 0.99), # US = 99th, EU = 95th
                   CK1 = quantile(Kidney, 0.01), CK50 = median(Kidney), CK95 = quantile(Kidney, 0.95),
                   CM1 = quantile(Muscle, 0.01), CM50 = median(Muscle), CM95 = quantile(Muscle, 0.95), CM99 = quantile(Muscle, 0.99),  # US = 99th, EU = 95th
                   CF1 = quantile(Fat, 0.01), CF50 = median(Fat), CF95 = quantile(Fat, 0.95))

data_tissue_sc <- as.data.frame(out.sum_sc)

end_time = 55

## Tissue concentration plots for withdrawal time estimation (SC)
##############################################################################################################################################################
  ## Liver - FDA tolerance (99th percentile) #################################################################################################################
withdrawal_Liver_sc   <- min(subset(data_tissue_sc, CL99 <= 0.1 & Time1 > 2)$Time1)  # (ug/g tissue), Historical FDA tolerance = 0.1 ug/g, Current FDA tolerance = 1.6 ug/g
has_cross_Liver_sc    <- any(data_tissue_sc$CL99 > 0.1)

p1_sc <- ggplot(data_tissue_sc, aes(Time1)) +
  geom_ribbon(aes(ymin = CL1, ymax = CL99), fill = '#B0C4DE', alpha = 0.3, size = 0.2, show.legend = TRUE) +
  geom_line(aes(y = CL50, color = '50th percentile'), size = 1.5, show.legend = TRUE) +  
  geom_line(aes(y = CL99, color = '99th percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) +  
  geom_line(aes(y = CL1, color = '1st percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) + 
  scale_color_manual(values = c('1st percentile' = '#616161', '50th percentile' = '#424242', '99th percentile' = '#616161')) +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dashed"), linewidth = 0.3))) +
  scale_x_continuous(name = 'Time (days)', breaks = seq(0, end_time, by = 5)) + 
  scale_y_log10(breaks = c(0.00001, 0.001, 0.1, 1, 10), labels = scales::label_number(drop0trailing = TRUE), expand = expansion(mult = c(0, 0.15))) + 
  ylab("Concentration (µg/mL or µg/g)") + 
  
  geom_line(aes(y = 0.1), color = '#E69F00', size = 1.5, linetype = 'dashed', show.legend = F) +  # add tolerance line / yellow (historical FDA tolerance): #E69F00, purple (current FDA tolerance): #9C27B0, blue (Codex MRL): #0057B8
  
  list(
    if (has_cross_Liver_sc) geom_vline(aes(xintercept = withdrawal_Liver_sc), size = 1.5, color = "#F0625D", linetype = 2, show.legend = F),
    if (has_cross_Liver_sc) annotate("text", x = withdrawal_Liver_sc, y = 0.1, label = paste("WDI:", round(withdrawal_Liver_sc, 2), "days"), hjust = -0.1, vjust=-0.5, angle = 0, size = 5)
  ) +

  theme_bw() +
  theme(
    axis.text = element_text(size = 23, face = "bold", colour = "black"), 
    axis.title = element_text(size = 22, face = "bold", colour = "black"), 
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.title.position = "panel"
  ) +
  labs(title = "Liver (SC, labeled) — Historical FDA tolerance")

p1_sc

ggsave("Lableled_SC_Liver_FDA1995.tiff", plot = p1_sc, width=25, height=13, units='cm', dpi=320)

  ## Liver - Codex MRL (95th percentile) ###################################################################################################################
withdrawal_Liver_sc   <- min(subset(data_tissue_sc, CL95 <= 0.8 & Time1 > 3)$Time1)  # (ug/g tissue), JECFA MRL = 0.8 ug/g
has_cross_Liver_sc    <- any(data_tissue_sc$CL95 > 0.8)

p11_sc <- ggplot(data_tissue_sc, aes(Time1)) +
  geom_ribbon(aes(ymin = CL1, ymax = CL95), fill = '#B0C4DE', alpha = 0.3, size = 0.2, show.legend = TRUE) +
  geom_line(aes(y = CL50, color = '50th percentile'), size = 1.5, show.legend = TRUE) +  
  geom_line(aes(y = CL95, color = '95th percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) +  
  geom_line(aes(y = CL1, color = '1st percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) + 
  scale_color_manual(values = c('1st percentile' = '#616161', '50th percentile' = '#424242', '95th percentile' = '#616161')) +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dashed"), linewidth = 0.3))) +
  scale_x_continuous(name = 'Time (days)', breaks = seq(0, end_time, by = 5)) + 
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) + 
  ylab("Concentration (µg/mL or µg/g)") + 
  
  geom_line(aes(y = 0.8), color = '#0057B8', size = 1.5, linetype = 'dashed', show.legend = F) +  # add tolerance line / yellow (historical FDA tolerance): #E69F00, purple (current FDA tolerance): #9C27B0, blue (Codex MRL): #0057B8
  
  list(
    if (has_cross_Liver_sc) geom_vline(aes(xintercept = withdrawal_Liver_sc), size = 1.5, color = "#F0625D", linetype = 2, show.legend = F),
    if (has_cross_Liver_sc) annotate("text", x = withdrawal_Liver_sc, y = 0.8, label = paste("WDI:", round(withdrawal_Liver_sc, 2), "days"), hjust = -0.1, vjust=-0.5, angle = 0, size = 5)
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 23, face = "bold", colour = "black"),
    axis.title = element_text(size = 22, face = "bold", colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.title.position = "panel"
  ) +
  labs(title = "Liver (SC, labeled) — Codex MRL")

p11_sc

ggsave("Lableled_SC_Liver_Codex.tiff", plot = p11_sc, width=25, height=13, units='cm', dpi=320)

##############################################################################################################################################################
  ## Kidney - Not established by the US FDA ##################################################################################################################

  ## Kidney - Codex MRL (95th percentile) ####################################################################################################################
withdrawal_Kidney_sc <- min(subset(data_tissue_sc, CK95 <= 0.1 & Time1 > 2)$Time1)  # (ug/g tissue), JECFA MRL = 0.1 ug/g
has_cross_Kidney_sc  <- any(data_tissue_sc$CK95 > 0.1)

p22_sc <- ggplot(data_tissue_sc, aes(Time1)) +
  geom_ribbon(aes(ymin = CK1, ymax = CK95), fill = '#B0C4DE', alpha = 0.3, size = 0.2, show.legend = TRUE) +
  geom_line(aes(y = CK50, color = '50th percentile'), size = 1.5, show.legend = TRUE) +  
  geom_line(aes(y = CK95, color = '95th percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) +  
  geom_line(aes(y = CK1, color = '1st percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) + 
  scale_color_manual(values = c('1st percentile' = '#616161', '50th percentile' = '#424242', '95th percentile' = '#616161')) +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dashed"), linewidth = 0.3))) +
  scale_x_continuous(name = 'Time (days)', breaks = seq(0, end_time, by = 5)) + 
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) + 
  ylab("Concentration (µg/mL or µg/g)") + 
  
  geom_line(aes(y = 0.1), color = '#0057B8', size = 1.5, linetype = 'dashed', show.legend = F) +  # add tolerance line / yellow (historical FDA tolerance): #E69F00, purple (current FDA tolerance): #9C27B0, blue (Codex MRL): #0057B8
  
  list(
    if (has_cross_Kidney_sc) geom_vline(aes(xintercept = withdrawal_Kidney_sc), size = 1.5, color = "#F0625D", linetype = 2, show.legend = F),
    if (has_cross_Kidney_sc) annotate("text", x = withdrawal_Kidney_sc, y = 0.1, label = paste("WDI:", round(withdrawal_Kidney_sc, 2), "days"), hjust = -0.1, vjust=-0.5, angle = 0, size = 5)
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 23, face = "bold", colour = "black"), 
    axis.title = element_text(size = 22, face = "bold", colour = "black"), 
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.title.position = "panel"
  ) +
  labs(title = "Kidney (SC, labeled) — Codex MRL")

p22_sc

ggsave("Lableled_SC_Kidney_Codex.tiff", plot = p22_sc, width=25, height=13, units='cm', dpi=320)

##############################################################################################################################################################
  ## Muscle - FDA tolerance (99th percentile) ################################################################################################################
withdrawal_Muscle_sc <- min(subset(data_tissue_sc, CM99 <= 0.65 & Time1 > 5)$Time1)  # (ug/g tissue), Current FDA tolerance = 0.65 ug/g
has_cross_Muscle_sc  <- any(data_tissue_sc$CM99 > 0.65)

p3_sc <- ggplot(data_tissue_sc, aes(Time1)) +
  geom_ribbon(aes(ymin = CM1, ymax = CM99), fill = '#B0C4DE', alpha = 0.3, size = 0.2, show.legend = TRUE) +
  geom_line(aes(y = CM50, color = '50th percentile'), size = 1.5, show.legend = TRUE) +  
  geom_line(aes(y = CM99, color = '99th percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) +  
  geom_line(aes(y = CM1, color = '1st percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) + 
  scale_color_manual(values = c('1st percentile' = '#616161', '50th percentile' = '#424242', '99th percentile' = '#616161')) +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dashed"), linewidth = 0.3))) +
  scale_x_continuous(name = 'Time (days)', breaks = seq(0, end_time, by = 5)) + 
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) + 
  ylab("Concentration (µg/mL or µg/g)") + 
  
  geom_line(aes(y = 0.65), color = '#9C27B0', size = 1.5, linetype = 'dashed', show.legend = F) +  # add tolerance line / yellow (historical FDA tolerance): #E69F00, purple (current FDA tolerance): #9C27B0, blue (Codex MRL): #0057B8
  
  list(
    if (has_cross_Muscle_sc) geom_vline(aes(xintercept = withdrawal_Muscle_sc), size = 1.5, color = "#F0625D", linetype = 2, show.legend = F),
    if (has_cross_Muscle_sc) annotate("text", x = withdrawal_Muscle_sc, y = 0.65, label = paste("WDI:", round(withdrawal_Muscle_sc, 2), "days"), hjust = -0.1, vjust=-0.5, angle = 0, size = 5)
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 23, face = "bold", colour = "black"), 
    axis.title = element_text(size = 22, face = "bold", colour = "black"), 
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.title.position = "panel"
  ) +
  labs(title = "Muscle (SC, labeled) — Current FDA tolerance")

p3_sc

ggsave("Lableled_SC_Muscle_FDA2020.tiff", plot = p3_sc, width=25, height=13, units='cm', dpi=320)

  ## Muscle - Codex MRL (95th percentile) ###################################################################################################################
withdrawal_Muscle_sc <- min(subset(data_tissue_sc, CM95 <= 0.03 & Time1 > 5)$Time1)  # (ug/g tissue), JECFA MRL = 0.03 ug/g
has_cross_Muscle_sc  <- any(data_tissue_sc$CM95 > 0.03)

p33_sc <- ggplot(data_tissue_sc, aes(Time1)) +
  geom_ribbon(aes(ymin = CM1, ymax = CM95), fill = '#B0C4DE', alpha = 0.3, size = 0.2, show.legend = TRUE) +
  geom_line(aes(y = CM50, color = '50th percentile'), size = 1.5, show.legend = TRUE) +  
  geom_line(aes(y = CM95, color = '95th percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) +  
  geom_line(aes(y = CM1, color = '1st percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) + 
  scale_color_manual(values = c('1st percentile' = '#616161', '50th percentile' = '#424242', '95th percentile' = '#616161')) +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dashed"), linewidth = 0.3))) +
  scale_x_continuous(name = 'Time (days)', breaks = seq(0, end_time, by = 5)) + 
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) + 
  ylab("Concentration (µg/mL or µg/g)") + 
  
  geom_line(aes(y = 0.03), color = '#0057B8', size = 1.5, linetype = 'dashed', show.legend = F) +  # add tolerance line / yellow (historical FDA tolerance): #E69F00, purple (current FDA tolerance): #9C27B0, blue (Codex MRL): #0057B8
  
  list(
    if (has_cross_Muscle_sc) geom_vline(aes(xintercept = withdrawal_Muscle_sc), size = 1.5, color = "#F0625D", linetype = 2, show.legend = F),
    if (has_cross_Muscle_sc) annotate("text", x = withdrawal_Muscle_sc, y = 0.03, label = paste("WDI:", round(withdrawal_Muscle_sc, 2), "days"), hjust = -0.1, vjust=-0.5, angle = 0, size = 5)
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 23, face = "bold", colour = "black"), 
    axis.title = element_text(size = 22, face = "bold", colour = "black"), 
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.title.position = "panel"
  ) +
  labs(title = "Muscle (SC, labeled) — Codex MRL")

p33_sc

ggsave("Lableled_SC_Muscle_Codex.tiff", plot = p33_sc, width=25, height=13, units='cm', dpi=320)

##############################################################################################################################################################
  ## Fat - Not established by the US FDA #####################################################################################################################

  ## Fat - Codex MRL (95th percentile) #######################################################################################################################
withdrawal_Fat_sc <- min(subset(data_tissue_sc, CF95 <= 0.4 & Time1 > 5)$Time1)  # (ug/g tissue), JECFA MRL = 0.4 ug/g
has_cross_Fat_sc  <- any(data_tissue_sc$CF95 > 0.4)

p44_sc <- ggplot(data_tissue_sc, aes(Time1)) +
  geom_ribbon(aes(ymin = CF1, ymax = CF95), fill = '#B0C4DE', alpha = 0.3, size = 0.2, show.legend = TRUE) +
  geom_line(aes(y = CF50, color = '50th percentile'), size = 1.5, show.legend = TRUE) +  
  geom_line(aes(y = CF95, color = '95th percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) +  
  geom_line(aes(y = CF1, color = '1st percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) + 
  scale_color_manual(values = c('1st percentile' = '#616161', '50th percentile' = '#424242', '95th percentile' = '#616161')) +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dashed"), linewidth = 0.3))) +
  scale_x_continuous(name = 'Time (days)', breaks = seq(0, end_time, by = 5)) + 
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) + 
  ylab("Concentration (µg/mL or µg/g)") + 
  
  geom_line(aes(y = 0.4), color = '#0057B8', size = 1.5, linetype = 'dashed', show.legend = F) +  # add tolerance line / yellow (historical FDA tolerance): #E69F00, purple (current FDA tolerance): #9C27B0, blue (Codex MRL): #0057B8
  
  list(
    if (has_cross_Fat_sc) geom_vline(aes(xintercept = withdrawal_Fat_sc), size = 1.5, color = "#F0625D", linetype = 2, show.legend = F),
    if (has_cross_Fat_sc) annotate("text", x = withdrawal_Fat_sc, y = 0.4, label = paste("WDI:", round(withdrawal_Fat_sc, 2), "days"), hjust = -0.1, vjust=-0.5, angle = 0, size = 5)
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 23, face = "bold", colour = "black"),
    axis.title = element_text(size = 22, face = "bold", colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.title.position = "panel"
  ) +
  labs(title = "Fat (SC, labeled) — Codex MRL")

p44_sc

ggsave("Lableled_SC_Fat_Codex.tiff", plot = p44_sc, width=25, height=13, units='cm', dpi=320)




## Combine data and run the simulation (Oral) --------------------------------------------------------------------------------------------------
set.seed(11)
system.time({
  out_oral <- mod %>% data_set(ex_oral) %>% update(atol= 1E-10, rtol=1E-5, maxsteps=50000, hmax = 0.01) %>% idata_set(idata) %>%
    mrgsim(obsonly=TRUE, tgrid=tsamp)
}) -> run.time_oral

out.sum_oral <- out_oral %>% as.data.frame %>%
  mutate(Time1 = time/24 - (ndoses-1)*doseinterval/24) %>%
  dplyr::select(ID, Time1, Plasma:Fat) %>%
  arrange(Time1)  %>%
  group_by(Time1) %>%
  dplyr::summarise(CV1 = quantile(Plasma, 0.01), CV50 = median(Plasma), CV99 = quantile(Plasma, 0.99), 
                   CL1 = quantile(Liver, 0.01), CL50 = median(Liver), CL95 = quantile(Liver, 0.95), CL99 = quantile(Liver, 0.99), # US = 99th, EU/Codex = 95th
                   CK1 = quantile(Kidney, 0.01), CK50 = median(Kidney), CK95 = quantile(Kidney, 0.95),
                   CM1 = quantile(Muscle, 0.01), CM50 = median(Muscle), CM95 = quantile(Muscle, 0.95), CM99 = quantile(Muscle, 0.99), # US = 99th, EU/Codex = 95th
                   CF1 = quantile(Fat, 0.01), CF50 = median(Fat), CF95 = quantile(Fat, 0.95))

data_tissue_oral <- as.data.frame(out.sum_oral)


## Tissue concentration plots for withdrawal time estimation (Oral)
##############################################################################################################################################################
  ## Liver - FDA tolerance (99th percentile) #################################################################################################################
withdrawal_Liver_oral   <- min(subset(data_tissue_oral, CL99 <= 0.1 & Time1 > 2)$Time1)  # (ug/g tissue), Historical FDA tolerance = 0.1 ug/g, Current FDA tolerance = 1.6 ug/g
has_cross_Liver_oral  <- any(data_tissue_oral$CL99 > 0.1)

p1_oral <- ggplot(data_tissue_oral, aes(Time1)) +
  geom_ribbon(aes(ymin = CL1, ymax = CL99), fill = '#B0C4DE', alpha = 0.3, size = 0.2, show.legend = TRUE) +
  geom_line(aes(y = CL50, color = '50th percentile'), size = 1.5, show.legend = TRUE) +  
  geom_line(aes(y = CL99, color = '99th percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) +  
  geom_line(aes(y = CL1, color = '1st percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) + 
  scale_color_manual(values = c('1st percentile' = '#616161', '50th percentile' = '#424242', '99th percentile' = '#616161')) +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dashed"), linewidth = 0.3))) +
  scale_x_continuous(name = 'Time (days)', breaks = seq(0, end_time, by = 5)) + 
  scale_y_log10(breaks = c(0.00001, 0.001, 0.1, 1, 10), labels = scales::label_number(drop0trailing = TRUE), expand = expansion(mult = c(0, 0.15))) + 
  ylab("Concentration (µg/mL or µg/g)") + 
  
  geom_line(aes(y = 0.1), color = '#E69F00', size = 1.5, linetype = 'dashed', show.legend = F) +  # add tolerance line / yellow (historical FDA tolerance): #E69F00, purple (current FDA tolerance): #9C27B0, blue (Codex MRL): #0057B8
  
  list(
    if (has_cross_Liver_oral) geom_vline(aes(xintercept = withdrawal_Liver_oral), size = 1.5, color = "#F0625D", linetype = 2, show.legend = F),
    if (has_cross_Liver_oral) annotate("text", x = withdrawal_Liver_oral, y = 0.1, label = paste("WDI:", round(withdrawal_Liver_oral, 2), "days"), hjust = -0.1, vjust=-0.5, angle = 0, size = 5)
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 23, face = "bold", colour = "black"), 
    axis.title = element_text(size = 22, face = "bold", colour = "black"), 
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.title.position = "panel"
  ) +
  labs(title = "Liver (PO, labeled) — Historical & Current FDA tolerance")

p1_oral

ggsave("Lableled_PO_Liver_FDA1995.tiff", plot = p1_oral, width=25, height=13, units='cm', dpi=320)

  ## Liver - Codex MRL (95th percentile) #####################################################################################################################
withdrawal_Liver_oral   <- min(subset(data_tissue_oral, CL95 <= 0.8 & Time1 > 1)$Time1)  # (ug/g tissue), Codex MRL = 0.8 ug/g
has_cross_Liver_oral  <- any(data_tissue_oral$CL95 > 0.8)

p11_oral <- ggplot(data_tissue_oral, aes(Time1)) +
  geom_ribbon(aes(ymin = CL1, ymax = CL95), fill = '#B0C4DE', alpha = 0.3, size = 0.2, show.legend = TRUE) +
  geom_line(aes(y = CL50, color = '50th percentile'), size = 1.5, show.legend = TRUE) +  
  geom_line(aes(y = CL95, color = '95th percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) +  
  geom_line(aes(y = CL1, color = '1st percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) + 
  scale_color_manual(values = c('1st percentile' = '#616161', '50th percentile' = '#424242', '95th percentile' = '#616161')) +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dashed"), linewidth = 0.3))) +
  scale_x_continuous(name = 'Time (days)', breaks = seq(0, end_time, by = 5)) + 
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) + 
  ylab("Concentration (µg/mL or µg/g)") + 
  
  geom_line(aes(y = 0.8), color = '#0057B8', size = 1.5, linetype = 'dashed', show.legend = F) +  # add tolerance line / yellow (historical FDA tolerance): #E69F00, purple (current FDA tolerance): #9C27B0, blue (Codex MRL): #0057B8
  
  list(
    if (has_cross_Liver_oral) geom_vline(aes(xintercept = withdrawal_Liver_oral), size = 1.5, color = "#F0625D", linetype = 2, show.legend = F),
    if (has_cross_Liver_oral) annotate("text", x = withdrawal_Liver_oral, y = 0.8, label = paste("WDI:", round(withdrawal_Liver_oral, 2), "days"), hjust = -0.1, vjust=-0.5, angle = 0, size = 5)
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 23, face = "bold", colour = "black"), 
    axis.title = element_text(size = 22, face = "bold", colour = "black"), 
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.title.position = "panel"
  ) +
  labs(title = "Liver (PO, labeled) — Codex MRL")

p11_oral

ggsave("Lableled_PO_Liver_Codex.tiff", plot = p11_oral, width=25, height=13, units='cm', dpi=320)

##############################################################################################################################################################
## Kidney - Not established by the US FDA ####################################################################################################################

## Kidney - Codex MRL (95th percentile) ######################################################################################################################
withdrawal_Kidney_oral   <- min(subset(data_tissue_oral, CK95 <= 0.1 & Time1 > 2)$Time1)  # (ug/g tissue), JECFA MRL = 100 ug/kg 
has_cross_Kidney_oral  <- any(data_tissue_oral$CK95 > 0.1)

p22_oral <- ggplot(data_tissue_oral, aes(Time1)) +
  geom_ribbon(aes(ymin = CK1, ymax = CK95), fill = '#B0C4DE', alpha = 0.3, size = 0.2, show.legend = TRUE) +
  geom_line(aes(y = CK50, color = '50th percentile'), size = 1.5, show.legend = TRUE) +  
  geom_line(aes(y = CK95, color = '95th percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) +  
  geom_line(aes(y = CK1, color = '1st percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) + 
  scale_color_manual(values = c('1st percentile' = '#616161', '50th percentile' = '#424242', '95th percentile' = '#616161')) +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dashed"), linewidth = 0.3))) +
  scale_x_continuous(name = 'Time (days)', breaks = seq(0, end_time, by = 5)) + 
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) + 
  ylab("Concentration (µg/mL or µg/g)") + 
  
  geom_line(aes(y = 0.1), color = '#0057B8', size = 1.5, linetype = 'dashed', show.legend = F) +  # add tolerance line / yellow (historical FDA tolerance): #E69F00, purple (current FDA tolerance): #9C27B0, blue (Codex MRL): #0057B8
  
  list(
    if (has_cross_Kidney_oral) geom_vline(aes(xintercept = withdrawal_Kidney_oral), size = 1.5, color = "#F0625D", linetype = 2, show.legend = F),
    if (has_cross_Kidney_oral) annotate("text", x = withdrawal_Kidney_oral, y = 0.1, label = paste("WDI:", round(withdrawal_Kidney_oral, 2), "days"), hjust = -0.1, vjust=-0.5, angle = 0, size = 5)
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 23, face = "bold", colour = "black"),
    axis.title = element_text(size = 22, face = "bold", colour = "black"), 
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.title.position = "panel"
  ) +
  labs(title = "Kidney (PO, labeled) — Codex MRL")

p22_oral

ggsave("Lableled_PO_Kidney_Codex.tiff", plot = p22_oral, width=25, height=13, units='cm', dpi=320)

##############################################################################################################################################################
  ## Muscle - FDA tolerance (99th percentile) ################################################################################################################
withdrawal_Muscle_oral   <- min(subset(data_tissue_oral, CM99 <= 0.65 & Time1 > 5)$Time1)  # (ug/g tissue), Current FDA tolerance = 0.65 ug/g
has_cross_Muscle_oral  <- any(data_tissue_oral$CM99 > 0.65)

p3_oral <- ggplot(data_tissue_oral, aes(Time1)) +
  geom_ribbon(aes(ymin = CM1, ymax = CM99), fill = '#B0C4DE', alpha = 0.3, size = 0.2, show.legend = TRUE) +
  geom_line(aes(y = CM50, color = '50th percentile'), size = 1.5, show.legend = TRUE) +  
  geom_line(aes(y = CM99, color = '99th percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) +  
  geom_line(aes(y = CM1, color = '1st percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) + 
  scale_color_manual(values = c('1st percentile' = '#616161', '50th percentile' = '#424242', '99th percentile' = '#616161')) +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dashed"), linewidth = 0.3))) +
  scale_x_continuous(name = 'Time (days)', breaks = seq(0, end_time, by = 5)) + 
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) + 
  ylab("Concentration (µg/mL or µg/g)") + 
  
  geom_line(aes(y = 0.65), color = '#9C27B0', size = 1.5, linetype = 'dashed', show.legend = F) +  # add tolerance line / yellow (historical FDA tolerance): #E69F00, purple (current FDA tolerance): #9C27B0, blue (Codex MRL): #0057B8
  
  list(
    if (has_cross_Muscle_oral) geom_vline(aes(xintercept = withdrawal_Muscle_oral), size = 1.5, color = "#F0625D", linetype = 2, show.legend = F),
    if (has_cross_Muscle_oral) annotate("text", x = withdrawal_Muscle_oral, y = 0.65, label = paste("WDI:", round(withdrawal_Muscle_oral, 2), "days"), hjust = -0.1, vjust=-0.5, angle = 0, size = 5)
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 23, face = "bold", colour = "black"), 
    axis.title = element_text(size = 22, face = "bold", colour = "black"), 
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.title.position = "panel"
  ) +
  labs(title = "Muscle (PO, labeled) — Current FDA tolerance")

p3_oral

ggsave("Lableled_PO_Muscle_FDA2020.tiff", plot = p3_oral, width=25, height=13, units='cm', dpi=320)

  ## Muscle - Codex MRL (95th percentile) ###################################################################################################################
withdrawal_Muscle_oral   <- min(subset(data_tissue_oral, CM95 <= 0.03 & Time1 > 2)$Time1)  # (ug/g tissue), JECFA MRL = 0.03 ug/g
has_cross_Muscle_oral  <- any(data_tissue_oral$CM95 > 0.03)

p33_oral <- ggplot(data_tissue_oral, aes(Time1)) +
  geom_ribbon(aes(ymin = CM1, ymax = CM95), fill = '#B0C4DE', alpha = 0.3, size = 0.2, show.legend = TRUE) +
  geom_line(aes(y = CM50, color = '50th percentile'), size = 1.5, show.legend = TRUE) +  
  geom_line(aes(y = CM95, color = '95th percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) +  
  geom_line(aes(y = CM1, color = '1st percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) + 
  scale_color_manual(values = c('1st percentile' = '#616161', '50th percentile' = '#424242', '95th percentile' = '#616161')) +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dashed"), linewidth = 0.3))) +
  scale_x_continuous(name = 'Time (days)', breaks = seq(0, end_time, by = 5)) + 
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) + 
  ylab("Concentration (µg/mL or µg/g)") + 
  
  geom_line(aes(y = 0.03), color = '#0057B8', size = 1.5, linetype = 'dashed', show.legend = F) +  # add tolerance line / yellow (historical FDA tolerance): #E69F00, purple (current FDA tolerance): #9C27B0, blue (Codex MRL): #0057B8
  
  list(
    if (has_cross_Muscle_oral) geom_vline(aes(xintercept = withdrawal_Muscle_oral), size = 1.5, color = "#F0625D", linetype = 2, show.legend = F),
    if (has_cross_Muscle_oral) annotate("text", x = withdrawal_Muscle_oral, y = 0.03, label = paste("WDI:", round(withdrawal_Muscle_oral, 2), "days"), hjust = -0.1, vjust=-0.5, angle = 0, size = 5)
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 23, face = "bold", colour = "black"), 
    axis.title = element_text(size = 22, face = "bold", colour = "black"), 
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.title.position = "panel"
  ) +
  labs(title = "Muscle (PO, labeled) — Codex MRL")

p33_oral

ggsave("Lableled_PO_Muscle_Codex.tiff", plot = p33_oral, width=25, height=13, units='cm', dpi=320)

##############################################################################################################################################################
## Fat - Not established by the US FDA #######################################################################################################################

## Fat - Codex MRL (95th percentile) #########################################################################################################################
withdrawal_Fat_oral   <- min(subset(data_tissue_oral, CF95 <= 0.4 & Time1 > 2.5)$Time1)  # (ug/g tissue), JECFA MRL = 0.4 ug/g
has_cross_Fat_oral  <- any(data_tissue_oral$CF95 > 0.4)

p44_oral <- ggplot(data_tissue_oral, aes(Time1)) +
  geom_ribbon(aes(ymin = CF1, ymax = CF95), fill = '#B0C4DE', alpha = 0.3, size = 0.2, show.legend = TRUE) +
  geom_line(aes(y = CF50, color = '50th percentile'), size = 1.5, show.legend = TRUE) +  
  geom_line(aes(y = CF95, color = '95th percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) +  
  geom_line(aes(y = CF1, color = '1st percentile'), size = 1.3, linetype = "dashed", show.legend = TRUE) + 
  scale_color_manual(values = c('1st percentile' = '#616161', '50th percentile' = '#424242', '95th percentile' = '#616161')) +
  guides(color = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dashed"), linewidth = 0.3))) +
  scale_x_continuous(name = 'Time (days)', breaks = seq(0, end_time, by = 5)) + 
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) + 
  ylab("Concentration (µg/mL or µg/g)") + 
  
  geom_line(aes(y = 0.4), color = '#0057B8', size = 1.5, linetype = 'dashed', show.legend = F) +  # add tolerance line / yellow (historical FDA tolerance): #E69F00, purple (current FDA tolerance): #9C27B0, blue (Codex MRL): #0057B8 
  
  list(
    if (has_cross_Fat_oral) geom_vline(aes(xintercept = withdrawal_Fat_oral), size = 1.5, color = "#F0625D", linetype = 2, show.legend = F),
    if (has_cross_Fat_oral) annotate("text", x = withdrawal_Fat_oral, y = 0.4, label = paste("WDI:", round(withdrawal_Fat_oral, 2), "days"), hjust = -0.1, vjust=-0.5, angle = 0, size = 5)
  ) +
  
  theme_bw() +
  theme(
    axis.text = element_text(size = 23, face = "bold", colour = "black"), 
    axis.title = element_text(size = 22, face = "bold", colour = "black"), 
    axis.ticks = element_line(colour = "black", linewidth = 1.2),
    axis.ticks.length = unit(0.25, "cm"),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.title.position = "panel"
  ) +
  labs(title = "Fat (PO, labeled) — Codex MRL")

p44_oral

ggsave("Lableled_PO_Fat_Codex.tiff", plot = p44_oral, width=25, height=13, units='cm', dpi=320)



