## ======================================================================================= ##
## PCPF-1R model Version 2.2 (Last update: 2022/01/04)                                     ##
## Simulating fate and transport processes of pesticide and its metabolite in rice paddy   ##
## Coded by K. Kondo (The Institute of Environmental Toxicology)                           ##
## Contact adresss; kondoh@iet.or.jp                                                       ##
## ======================================================================================= ##

####################
#       Log        #
####################

# 2016/02/20: Main program was developed (except nursery box application)
# 2016/02/29: Mass balance error due to spray and multiple application was improved.
# 2016/03/09: Metabolite formation in paddy soil was added.
# 2016/03/24: Consideration of temperature dependency.
# 2016/04/21: Change the assumption of flowable application (if feq=0 -> no partition)
# 2016/09/29: Improve feq multiply to both soild phase and liquid phase in soil
# 2017/02/08: Correct water balance bag
# 2017/05/23: Kinetic sorption process was substituted by linear diffusion exchange
# 2017/06/13: Pseudo equilibrium fraction "frac" was installed at dissolution term
# 2017/06/18: Integrate pore water diffusion and water sediment kinetic sorption
# 2017/07/26: Modify sorption term to use original assumption (kdiff = 0 and fdes=0)
# 2019/03/27: Add fast sorption term for flowable or emulsion application 
#             (during fast adsorption, reversible sorption term is off).
# 2019/04/03: cwDES and csDES term are simplified as 
#             (kdiff/24 + dDLi * bulk * Kd * (ksp/24)) -> (om / 24)
#             when fdes = 2. 
# 2019/04/15: kbiow -> kbiow (hydrolysis + biodegradation) + Kphot(photolysis)
# 2019/07/23: Correct unit conversion of diffusion term
# 2020/12/16: kbiow -> khyd + kbiow + kphot (need to add kbulk option)
#             pH dependency for khyd 
#             temperature dependency for khyd and kbiow
#             irradiation dependency for kphot  
# 2021/01/29: bulk degradation option is revived (set fdegw = 1 as default)
# 2022/01/04: User specified initial condition "ini_state" was introduced
#             Logical parameter "rBC" was installed for coupled simulation

####################
#   Main program   #
####################

PCPF1R <- function(

# Subsets of input parameters
   parms,
   parms2,

# User specified initial conditions (this can be ommitted)
   ini_state,

# Water balance data as a data frame
   water,

# Maximum simulation period
   DMAX,

# Switching parameter for kinetic sorption
   fdes = 1,      # = 0;  diffusion + constant desorption (original assumption)
                  # = 1;  diffusion + reversible sorption with sigle rate constant
                  # = 2;  overall mass trasnport constant (diffusion + reversible sorption)

# Switching parameter for degradation in paddy water
   fdegw = 1,     # = 1;  bulk degradation (specify kbiow only)
                  # = 2;  process based degradation (specify khyd, kbiow, kphot)

# Switching parameter for degradation in paddy soil -> under construction
#   fdegs = 1,     # = 1;  degradation in soil particle (original assumption; specify kbios only)
#                  # = 2;  include processes in pore water degradation (specify khyd, kbiow)

# Application method of pesticide 
   type,          # = 1;  granule
                  # = 2;  flowable or emulsion
                  # = 3;  nursery box application -> under construction
                  # = 4;  spray application

# Number of metabolite (user can ignore this parameter if user doesn't use this option)
   meta = 0,      # default is 0; if > 1, execute the metabolite computation -> under improvement

# Number of multiple application (user can ignore this parameter if user doesn't use this option)
   mapp = 0,      # default is 0
   mapplist,      # application schedule
   SMassC = 0,    # cumulative mass for multiple applicatons

# Environmnetal covariates dependences of degradation parameters (user can ignore this parameter if user doesn't use this option)
   temp = FALSE,  # temperature dependency frag (default is FALSE)
   pH   = FALSE,  # pH dependency frag (default is FALSE)
   uv   = FALSE,  # uv dependency frag (default is FALSE)
   Tdef = 20,     # reference temperature
   meteolist,     # all environmnetal covariates considerd should be included in this list

# Other optional parameters (user can ignore this parameter if user doesn't use this option)
   rBC = FALSE    # If TRUE, boundary condition is returned(for coupled modeling, default is FALSE)
){

# -------------------------------- Initialization  -------------------------------------------
if (missing(ini_state)) {
   if (meta > 0) {
     state <- c( Cpwi         = 0,
                 Csdli        = 0,
                 Mpwdegi      = 0,
                 Msdldegi     = 0,
                 CpwH         = 0,
                 CsdlH        = 0,
                 Cmpwi        = 0,
                 Cmsdli       = 0,
                 Mmpwdegi     = 0,
                 Mmsdldegi    = 0,
                 CmpwH        = 0,
                 CmsdlH       = 0,
                 PMdissi      = 0,
                 TPMdissi     = 0,
                 CRFloss      = 0,
                 CPERCloss    = 0,
                 CLSEEPloss   = 0,
                 CMpwdeg      = 0,
                 CMsdldeg     = 0,
                 CMCdeg       = 0,
                 CRFmloss     = 0,
                 CPERCmloss   = 0,
                 CLSEEPmloss  = 0,
                 CMmpwdeg     = 0,
                 CMmsdldeg    = 0,
                 CMmCdeg      = 0,
                 Dperci       = 0,
                 dDLi         = 0,
                 Fdiss        = 1,
                 dfdDLi       = 0,
                 DriMass      = 1
                 )
   } else {  
     state <- c( Cpwi         = 0,
                 Csdli        = 0,
                 Mpwdegi      = 0,
                 Msdldegi     = 0,
                 CpwH         = 0,
                 CsdlH        = 0,
                 PMdissi      = 0,
                 TPMdissi     = 0,
                 CRFloss      = 0,
                 CPERCloss    = 0,
                 CLSEEPloss   = 0,
                 CMpwdeg      = 0,
                 CMsdldeg     = 0,
                 CMCdeg       = 0,
                 Dperci       = 0,
                 dDLi         = 0,
                 Fdiss        = 1,
                 dfdDLi       = 0,
                 DriMass      = 1
                 ) 
   } 
} else {
   state <-  ini_state
}

#
# Define time step
DT      <- 1

# Assume cyclic irrigation is not practiced
cwIRR   <- 0

# Simulation period (unit: hour)
DMAX    <- DMAX * 24

# Calculate application amount (unit: ug)
APPMass <- with(as.list(c(parms, parms2)), AppR * Area * 1000000)
CUMMass <- 0

# 
# -------------------- Define sub-function solving mass balance equation -----------------------------
DISS <- function(parms,
                 wbdata1, 
                 wbdata2, 
                 state, 
                 DT = 1, 
                 APPMass,
                 CUMMass,
                 DriMass,
                 SMass = 0,
                 mapp
  ){
  with(as.list(c(parms, wbdata1, wbdata2, state)), {

# Initial mass partitioning for flowerble, emulsion and spray
      if ((type == 2) && (Fdiss == 1) || (type == 4) && (Fdiss == 1)) {
          if (type==2) {
            if ((mapp > 0) || (mapp < 0)) {
              PMass    <- AppRd * Area * 1000000 # flowable and emulsion
            } else {
              PMass    <- AppR * Area * 1000000 # flowable and emulsion
            }
          } else if (type==4) {
            if ((mapp > 0) || (mapp < 0)) {
              PMass    <- AppRd * (1 - DRIFT) * (1 - COVER) * Area * 1000000 # spray
            } else {
              PMass    <- AppR * (1 - DRIFT) * (1 - COVER) * Area * 1000000 # spray
            }
            SMass    <- PMass
          }
          if (f == 0) { 
            frac     <- 0
            PmassDL  <- 0
          } else {
 #           frac     <- HpwH / (f *(Kd * bulk + SatWC))
 #           PmassDL  <- PMass / (1 + frac)
            frac     <- 0
            PmassDL  <- 0
          }
          PmassPW  <- PMass - PmassDL
          Cpwi     <- Cpwi + PmassPW / (10000 * Area * Hpw2)
          Csdli    <- Csdli + PmassDL / (10000 * Area * dDL * (bulk + SatWC / Kd))
          Fdiss    <- 0
          dfdDLi   <- 0
          dfdDLH   <- 0
          dfdDL    <- 0
          ksp      <- 0
      } else {
          PMass    <- APPMass
      }
# Loop for dissolution terminator 
      for (k in 1:2) {
# Forth-order Runge-Kutta computation
# k1
        cwDISS   <- ((kdiss / 24) * (CSLB - Cpwi) + Cpwi / Hpw1 * dHpw1 / DT) * Fdiss
        if (fdes == 2) {
          cwDES    <- 1 / Hpw1 * (om / 24) *(Csdli / Kd - Cpwi)
        } else {
          cwDES    <- 1 / Hpw1 * ((100*kdiff)/24 + dDLi * bulk * Kd * (ksorp / 24)) *(Csdli / Kd - fdes * Cpwi)
        }
        if ((type == 2) || (type == 4)) {
          cwPERT   <- 1 / Hpw1 * dDLi * bulk * (alp / 24) *(Csdli - f * Kd * Cpwi)
        } else {
          cwPERT   <- 0
        }        
        cwIRR    <- 1 / Hpw1 * irr * cwIRR
        cwDPL    <- 1 / Hpw1 * Cpwi * (drain + perc + seep)
        cwVOL    <- 1 / Hpw1 * 100 * (kvol / 24) * Cpwi
        if (fdegw == 2) {
          cwDEG    <- ((khyd +kbiow + kphot) / 24) * Cpwi
        } else {
          cwDEG    <- (kbiow / 24) * Cpwi
        }
        cwdhdt   <- 1 / Hpw1 * dHpw1 / DT * Cpwi
        K1w      <- cwDISS + cwDES + cwPERT + cwIRR - cwDPL - cwVOL - cwDEG - cwdhdt 
        K1wdeg   <- -(cwVOL + cwDEG) * 10000 * Area * Hpw2
#
        if (dDLi == 0) {
          csDISS    <- (f * Kd * (kdiss / 24) * (CSLB - Cpwi)) * Fdiss
          csPERC    <- 0
          csDES     <- 0
          csdDLdt   <- 0
          csPERT    <- 0
        } else if (dDLi < 1) {
          csDISS    <- (f * Kd * (kdiss / 24) * (CSLB - Cpwi) + f * Kd / dDLi * Cpwi * dfdDLi / DT) * Fdiss
          csPERC    <- (Kd / ((SatWC + bulk * Kd) * dDLi)) * perc * (Cpwi)
          if (fdes == 2) {
            csDES    <- (Kd / ((SatWC + bulk * Kd) * dDLi)) * (om / 24) * (Csdli / Kd - Cpwi)
          } else { 
            csDES    <- (Kd / ((SatWC + bulk * Kd) * dDLi)) * ((100*kdiff)/24 + dDLi * bulk * Kd * (ksorp / 24)) * (Csdli / Kd - fdes * Cpwi)
          }
          if ((type == 2) || (type == 4)) {
            csPERT   <- (Kd / ((SatWC + bulk * Kd) * dDLi)) * (dDLi * bulk * (alp / 24)) * (Csdli - f * Kd * Cpwi)
          } else {
            csPERT   <- 0
          } 
          csdDLdt  <- Csdli * dfdDLi / (dDLi * DT)
        } else {
          csDISS    <- (f *Kd * (kdiss / 24) * (CSLB - Cpwi) + f * Kd / dDLi * Cpwi * dfdDLi / DT) * Fdiss
          csPERC   <- (Kd / ((SatWC + bulk * Kd) * dDLi)) * perc * (Cpwi - Csdli / Kd)
          if (fdes == 2) {
            csDES    <- (Kd / ((SatWC + bulk * Kd) * dDLi)) * (om / 24) * (Csdli / Kd - Cpwi)
          } else {
            csDES    <- (Kd / ((SatWC + bulk * Kd) * dDLi)) * ((100*kdiff)/24 + dDLi * bulk * Kd * (ksorp / 24)) * (Csdli / Kd - fdes * Cpwi)
          }
          if ((type == 2) || (type == 4)) {
            csPERT   <- (Kd / ((SatWC + bulk * Kd) * dDLi)) * (dDLi * bulk * (alp / 24)) * (Csdli - f * Kd * Cpwi)
          } else {
            csPERT   <- 0
          } 
          csdDLdt  <- Csdli * dfdDLi / (dDLi * DT)
        } 
        csBIO    <- Kd / (SatWC + bulk * Kd) * bulk * (kbios / 24) * Csdli
        K1s      <- csDISS + csPERC - csBIO - csDES - csPERT - csdDLdt
        K1sdeg   <- -csBIO / (Kd / (SatWC + bulk * Kd))  * 10000 * Area * dDL
#
# k2
        cwDISS   <- ((kdiss / 24) * (CSLB - (Cpwi + 0.5 * K1w)) + (Cpwi + 0.5 * K1w) / HpwH * dHpwH / DT) * Fdiss
        if (fdes == 2) {
          cwDES    <- 1 / HpwH * (om / 24) *((Csdli + 0.5 * K1s) / Kd - (Cpwi + 0.5 * K1w))
        } else {
          cwDES    <- 1 / HpwH * ((100*kdiff)/24 + dDLH * bulk * Kd *  (ksorp / 24)) *((Csdli + 0.5 * K1s) / Kd - fdes * (Cpwi + 0.5 * K1w))
        }
        if ((type == 2) || (type == 4)) {
          cwPERT   <- 1 / Hpw1 * dDLH * bulk * (alp / 24) *((Csdli + 0.5 * K1s) - f * Kd * (Cpwi + 0.5 * K1w))
        } else {
          cwPERT   <- 0
        }   
        cwIRR    <- 1 / HpwH * irr * cwIRR
        cwDPL    <- 1 / HpwH * (Cpwi + 0.5 * K1w) * (drain + perc + seep)
        cwVOL    <- 1 / HpwH * 100 * (kvol / 24) * (Cpwi + 0.5 * K1w)
        cwDEG    <- ((khyd +kbiow + kphot) / 24) * (Cpwi + 0.5 * K1w)
        if (fdegw == 2) {
          cwDEG    <- ((khyd +kbiow + kphot) / 24) * (Cpwi + 0.5 * K1w)
        } else {
          cwDEG    <- (kbiow / 24) * (Cpwi + 0.5 * K1w)
        }
        cwdhdt   <- 1 / HpwH * dHpwH / DT * (Cpwi + 0.5 * K1w)
        K2w      <- cwDISS + cwDES + cwPERT + cwIRR - cwDPL - cwVOL - cwDEG - cwdhdt 
        K2wdeg   <- -(cwVOL + cwDEG) * 10000 * Area * Hpw2
#
        csDISS   <- (f * Kd * (kdiss / 24) * (CSLB - (Cpwi + 0.5 * K1w)) + f * Kd / dDLH * (Cpwi + 0.5 * K1w) * dfdDLH / DT) * Fdiss
        if (dDLH != 1) {
          csPERC   <- (Kd / ((SatWC + bulk * Kd) * dDLH)) * perc * (Cpwi + 0.5 * K1w)
          csdDLdt  <- (Csdli + 0.5 * K1s) * dfdDLH / (dDLH * DT)
        } else {
          csPERC   <- (Kd / ((SatWC + bulk * Kd) * dDLH)) * perc * ((Cpwi + 0.5 * K1w) - (Csdli + 0.5 * K1s) / Kd)
          csdDLdt  <- (Csdli + 0.5 * K1s) * dfdDLH / (dDLH * DT)
        }          
        csBIO    <- Kd / (SatWC + bulk * Kd) * bulk * (kbios / 24) * (Csdli + 0.5 * K1s)
        if (fdes == 2) {
          csDES    <- (Kd / ((SatWC + bulk * Kd) * dDLH)) * (om / 24) * ((Csdli + 0.5 * K1s) / Kd - (Cpwi + 0.5 * K1w))
        } else {
          csDES    <- (Kd / ((SatWC + bulk * Kd) * dDLH)) * ((100*kdiff)/24 + dDLH * bulk * Kd *  (ksorp / 24)) * ((Csdli + 0.5 * K1s) / Kd - fdes * (Cpwi + 0.5 * K1w))
        }
        if ((type == 2) || (type == 4)) {
          csPERT   <- (Kd / ((SatWC + bulk * Kd) * dDLH)) * (dDLH * bulk * (alp / 24)) * ((Csdli + 0.5 * K1s) - f * Kd * (Cpwi + 0.5 * K1w))
        } else {
          csPERT   <- 0
        } 
        K2s      <- csDISS + csPERC - csBIO - csDES - csPERT - csdDLdt
        K2sdeg   <- -csBIO / (Kd / (SatWC + bulk * Kd)) * 10000 * Area * dDL
#
# k3
        cwDISS   <- ((kdiss / 24) * (CSLB - (Cpwi + 0.5 * K2w)) + (Cpwi + 0.5 * K2w) / HpwH * dHpwH / DT) * Fdiss
        if (fdes == 2) {
          cwDES    <- 1 / HpwH * (om / 24) * ((Csdli + 0.5 * K2s) / Kd - (Cpwi + 0.5 * K2w))
        } else {
          cwDES    <- 1 / HpwH * ((100*kdiff)/24 + dDLH * bulk * Kd *  (ksorp/24)) * ((Csdli + 0.5 * K2s) / Kd - fdes * (Cpwi + 0.5 * K2w))
        }
        if ((type == 2) || (type == 4)) {
          cwPERT   <- 1 / Hpw1 * dDLH * bulk * (alp / 24) *((Csdli + 0.5 * K2s) - f * Kd * (Cpwi + 0.5 * K2w))
        } else {
          cwPERT   <- 0
        }   
        cwIRR    <- 1 / HpwH * irr * cwIRR
        cwDPL    <- 1 / HpwH * (Cpwi + 0.5 * K2w) * (drain + perc + seep)
        cwVOL    <- 1 / HpwH * 100 * (kvol/24) * (Cpwi + 0.5 * K2w)
        if (fdegw == 2) {
          cwDEG    <- ((khyd + kbiow + kphot) / 24) * (Cpwi + 0.5 * K2w)
        } else {
          cwDEG    <- (kbiow / 24) * (Cpwi + 0.5 * K2w)
        }
        cwdhdt   <- 1 / HpwH * dHpwH / DT * (Cpwi + 0.5 * K2w)
        K3w      <- cwDISS + cwDES + cwPERT + cwIRR - cwDPL - cwVOL - cwDEG - cwdhdt 
        K3wdeg   <- -(cwVOL + cwDEG) * 10000 * Area * Hpw2
#
        csDISS   <- (f * Kd * (kdiss / 24) * (CSLB - (Cpwi + 0.5 * K2w)) + f * Kd / dDLH * (Cpwi + 0.5 * K2w) * dfdDLH / DT) * Fdiss
        if (dDLH != 1) {
          csPERC   <- (Kd / ((SatWC + bulk * Kd) * dDLH)) * perc * (Cpwi + 0.5 * K2w)
          csdDLdt  <- (Csdli + 0.5 * K2s) * dfdDLH / (dDLH * DT)
        } else {
          csPERC   <- (Kd / ((SatWC + bulk * Kd) * dDLH)) * perc * ((Cpwi + 0.5 * K2w) - (Csdli + 0.5 * K2s) / Kd)
          csdDLdt  <- (Csdli + 0.5 * K2s) * dfdDLH / (dDLH * DT)
        }   
        csBIO    <- Kd / (SatWC + bulk * Kd) * bulk * (kbios / 24) * (Csdli + 0.5 * K2s)
        if (fdes == 2) {
          csDES    <- (Kd / ((SatWC + bulk * Kd) * dDLH)) * (om / 24) * ((Csdli + 0.5 * K2s) / Kd - (Cpwi + 0.5 * K2w))
        } else {
          csDES    <- (Kd / ((SatWC + bulk * Kd) * dDLH)) * ((100*kdiff)/24 + dDLH * bulk * Kd *  (ksorp/24)) * ((Csdli + 0.5 * K2s) / Kd - fdes * (Cpwi + 0.5 * K2w))
        }
        if ((type == 2) || (type == 4)) {
          csPERT   <- (Kd / ((SatWC + bulk * Kd) * dDLH)) * (dDLH * bulk * (alp / 24)) * ((Csdli + 0.5 * K2s) - f * Kd * (Cpwi + 0.5 * K2w))
        } else {
          csPERT   <- 0
        } 
        K3s      <- csDISS + csPERC - csBIO - csDES - csPERT - csdDLdt
        K3sdeg   <- -csBIO / (Kd / (SatWC + bulk * Kd)) * 10000 * Area * dDL
#
# k4
        cwDISS   <- ((kdiss / 24) * (CSLB - (Cpwi + K3w)) + (Cpwi + K3w) / HpwH * dHpwH / DT) * Fdiss
        if (fdes == 2) {
          cwDES    <- 1 / Hpw2 * (om / 24) *((Csdli + K3s) / Kd - (Cpwi + K3w))
        } else {
          cwDES    <- 1 / Hpw2 * ((100*kdiff) / 24 + dDL * bulk * Kd *  (ksorp/24)) *((Csdli + K3s) / Kd - fdes * (Cpwi + K3w))
        }
        if ((type == 2) || (type == 4)) {
          cwPERT   <- 1 / Hpw1 * dDL * bulk * (alp / 24) *((Csdli + K3s) - f * Kd * (Cpwi + K3w))
        } else {
          cwPERT   <- 0
        }  
        cwIRR    <- 1 / Hpw2 * irr * cwIRR
        cwDPL    <- 1 / Hpw2 * (Cpwi + K3w) * (drain + perc + seep)
        cwVOL    <- 1 / Hpw2 * 100 * (kvol / 24) * (Cpwi + K3w)
        if (fdegw == 2) {
          cwDEG    <- ((khyd +kbiow + kphot) / 24) * (Cpwi + K3w)
        } else {
          cwDEG    <- (kbiow / 24) * (Cpwi + K3w)
        }
        cwdhdt   <- 1 / Hpw2 * dHpw2 / DT * (Cpwi + K3w)
        K4w      <- cwDISS + cwDES + cwPERT + cwIRR - cwDPL - cwVOL - cwDEG - cwdhdt 
        K4wdeg   <- -(cwVOL + cwDEG) * 10000 * Area * Hpw2
#
        csDISS   <- (f * Kd * (kdiss / 24) * (CSLB - (Cpwi + K3w)) + f * Kd / dDL * (Cpwi + K3w) * dfdDL / DT) * Fdiss
        if (dDL != 1) {
          csPERC   <- (Kd / ((SatWC + bulk * Kd) * dDL)) * perc * (Cpwi + K3w)
          csdDLdt  <- (Csdli + K3s) * dfdDL / (dDL * DT)
        } else {
          csPERC   <- (Kd / ((SatWC + bulk * Kd) * dDL)) * perc * ((Cpwi + K3w) - (Csdli + K3s) / Kd)
          csdDLdt  <- (Csdli + K3s) * dfdDL / (dDL * DT)
        }   
        csBIO    <- Kd / (SatWC + bulk * Kd) * bulk * (kbios/24) * (Csdli + K3s)
        if (fdes == 2) {
          csDES    <- (Kd / ((SatWC + bulk * Kd) * dDL)) * (om / 24) * ((Csdli + K3s) / Kd - (Cpwi + K3w))
        } else {
          csDES    <- (Kd / ((SatWC + bulk * Kd) * dDL)) * ((100*kdiff)/24 + dDL * bulk * Kd *  (ksorp/24)) * ((Csdli + K3s) / Kd - fdes * (Cpwi + K3w))
        }
        if ((type == 2) || (type == 4)) {
          csPERT   <- (Kd / ((SatWC + bulk * Kd) * dDL)) * (dDL * bulk * (alp / 24)) * ((Csdli + K3s) - f * Kd * (Cpwi + K3w))
        } else {
          csPERT   <- 0
        } 
        K4s      <- csDISS + csPERC - csBIO - csDES - csPERT - csdDLdt
        K4sdeg   <- -csBIO / (Kd / (SatWC + bulk * Kd)) * 10000 * Area * dDL
  #
  # For metabolite
        if (meta > 0) {
          # k1
          cwDES     <- 1 / Hpw1 * (kmdiff/24) * (Cmsdli / Kmd - Cmpwi)
          cwIRR     <- 1 / Hpw1 * irr * cwIRR
          cwDPL     <- 1 / Hpw1 * Cmpwi * (drain + perc + seep)
          cwVOL     <- 1 / Hpw1 * 100 * (kmvol / 24) * Cmpwi
          cwPHOTO   <- 0
#          cwPHOTO  <- Kphot * dUVBH / DT * Cpwi
          cwBIO     <- (kmbiow / 24) * Cmpwi
          cwFORM    <- ffw * (kbiow / 24) * Cpwi * Mf 
          cwdhdt    <- 1 / Hpw1 * dHpw1 / DT * Cmpwi
          K1mw      <- cwFORM + cwDES + cwIRR - cwDPL - cwVOL - cwPHOTO - cwBIO - cwdhdt 
          K1mwdeg   <- -(cwVOL + cwPHOTO + cwBIO) * 10000 * Area * Hpw2
          #
          if (dDLi == 0) {
            csPERC   <- 0
            csdDLdt  <- 0
          } else if (dDLi < 1) {
            csPERC   <- (Kmd / ((SatWC + bulk * Kmd) * dDLi)) * perc * (Cmpwi)
            csdDLdt  <- Cmsdli * dfdDLi / (dDLi * DT)
          } else {
            csPERC   <- (Kmd / ((SatWC + bulk * Kmd) * dDLi)) * perc * (Cmpwi - Cmsdli / Kmd)
            csdDLdt  <- Cmsdli * dfdDLi / (dDLi * DT)
          }
          csFORM    <- ffs * (Kmd / (SatWC + bulk * Kmd) * bulk * (kbios / 24) * Csdli) * Mf
          csBIO     <- Kmd / (SatWC + bulk * Kmd) * bulk * (kmbios / 24) * Cmsdli
          csDES     <- (Kmd / ((SatWC + bulk * Kmd) * dDLi)) * (kmdiff / 24) * (Cmsdli / Kmd - Cmpwi)
          K1ms      <- csFORM + csPERC - csBIO - csDES - csdDLdt
          K1msdeg   <- -csBIO / (Kmd / (SatWC + bulk * Kmd)) * 10000 * Area * dDL
          #
          # k2
          cwDES    <- 1 / HpwH * (kmdiff/24) *( (Cmsdli + 0.5 * K1ms) / Kmd - (Cmpwi + 0.5 * K1mw))
          cwIRR    <- 1 / HpwH * irr * cwIRR
          cwDPL    <- 1 / HpwH * (Cmpwi + 0.5 * K1mw) * (drain + perc + seep)
          cwVOL    <- 1 / HpwH * 100 * (kmvol / 24) * (Cmpwi + 0.5 * K1mw)
          cwPHOTO  <- 0
#          cwPHOTO <- Kphot * dUVBH / DT * (Cpwi + 0.5 * K1w)
          cwBIO    <- (kmbiow / 24) * (Cmpwi + 0.5 * K1mw)
          cwFORM   <- ffw * (kbiow / 24) * (Cpwi + 0.5 * K1w) * Mf
          cwdhdt   <- 1 / HpwH * dHpwH / DT * (Cmpwi + 0.5 * K1mw)
          K2mw      <- cwFORM + cwDES + cwIRR - cwDPL - cwVOL - cwPHOTO - cwBIO - cwdhdt 
          K2mwdeg   <- -(cwVOL + cwPHOTO + cwBIO) * 10000 * Area * Hpw2
#  
          if (dDLH != 1) {
            csPERC <- (Kmd / ((SatWC + bulk * Kmd) * dDLH)) * perc * (Cmpwi + 0.5 * K1mw)
            csdDLdt  <- (Cmsdli + 0.5 * K1ms) * dfdDLH / (dDLH * DT)
          } else {
            csPERC <- (Kmd / ((SatWC + bulk * Kmd) * dDLH)) * perc * ((Cmpwi + 0.5 * K1mw) - (Cmsdli + 0.5 * K1ms) / Kmd)
            csdDLdt  <- (Cmsdli + 0.5 * K1ms) * dfdDLH / (dDLH * DT)
          }
          csFORM    <- ffs * (Kmd / (SatWC + bulk * Kmd) * bulk * (kbios / 24) * (Csdli + 0.5 * K1s)) * Mf            
          csBIO     <- Kmd / (SatWC + bulk * Kmd) * bulk * (kmbios / 24) * (Cmsdli + 0.5 * K1ms)
          csDES     <- (Kmd / ((SatWC + bulk * Kmd) * dDLH)) * (kmdiff / 24) * ((Cmsdli + 0.5 * K1ms) / Kmd - (Cmpwi + 0.5 * K1mw))
          K2ms      <- csFORM + csPERC - csBIO - csDES - csdDLdt
          K2msdeg   <- -csBIO / (Kmd / (SatWC + bulk * Kmd)) * 10000 * Area * dDL
#
# k3
          cwDES    <- 1 / HpwH * (kmdiff/24) *( (Cmsdli + 0.5 * K2ms) / Kmd - (Cmpwi + 0.5 * K2mw))
          cwIRR    <- 1 / HpwH * irr * cwIRR
          cwDPL    <- 1 / HpwH * (Cmpwi + 0.5 * K2mw) * (drain + perc + seep)
          cwVOL    <- 1 / HpwH * 100 * (kmvol/24) * (Cmpwi + 0.5 * K2mw)
          cwPHOTO  <- 0
#          cwPHOTO <- Kphot * dUVBH / DT * (Cpwi + 0.5 * K2w)
          cwBIO    <- (kmbiow / 24) * (Cmpwi + 0.5 * K2mw)
          cwFORM   <- ffw * (kbiow / 24) * (Cpwi + 0.5 * K2w) * Mf
          cwdhdt   <- 1 / HpwH * dHpwH / DT * (Cmpwi + 0.5 * K2mw)
          K3mw      <- cwFORM + cwDES + cwIRR - cwDPL - cwVOL - cwPHOTO - cwBIO - cwdhdt 
          K3mwdeg   <- -(cwVOL + cwPHOTO + cwBIO) * 10000 * Area * Hpw2
          #
          if (dDLH != 1) {
            csPERC <- (Kmd / ((SatWC + bulk * Kmd) * dDLH)) * perc * (Cmpwi + 0.5 * K2mw)
            csdDLdt  <- (Cmsdli + 0.5 * K2ms) * dfdDLH / (dDLH * DT)
          } else {
            csPERC <- (Kmd / ((SatWC + bulk * Kmd) * dDLH)) * perc * ((Cmpwi + 0.5 * K2mw) - (Cmsdli + 0.5 * K2ms) / Kmd)
            csdDLdt  <- (Cmsdli + 0.5 * K2ms) * dfdDLH / (dDLH * DT)
          }
          csFORM    <- ffs * (Kmd / (SatWC + bulk * Kmd) * bulk * (kbios / 24) * (Csdli + 0.5 * K2s)) * Mf  
          csBIO     <- Kmd / (SatWC + bulk * Kmd) * bulk * (kmbios / 24) * (Cmsdli + 0.5 * K2ms)
          csDES     <- (Kmd / ((SatWC + bulk * Kmd) * dDLH)) * (kmdiff / 24) * ((Cmsdli + 0.5 * K2ms) / Kmd - (Cmpwi + 0.5 * K2mw))
          K3ms      <- csFORM + csPERC - csBIO - csDES - csdDLdt
          K3msdeg   <- -csBIO / (Kmd / (SatWC + bulk * Kmd)) * 10000 * Area * dDL
#
# k4
          cwDES    <- 1 / Hpw2 * (kmdiff / 24) *((Cmsdli + K3ms) / Kmd - (Cmpwi + K3mw))
          cwIRR    <- 1 / Hpw2 * irr * cwIRR
          cwDPL    <- 1 / Hpw2 * (Cmpwi + K3mw) * (drain + perc + seep)
          cwVOL    <- 1 / Hpw2 * 100 * (kmvol / 24) * (Cmpwi + K3mw)
          cwPHOTO  <- 0
#          cwPHOTO <- Kphot * dUVBH / DT * (Cpwi + K3w)
          cwBIO    <- (kmbiow / 24) * (Cmpwi + 0.5 * K3mw)
          cwFORM   <- ffw * (kbiow / 24) * (Cpwi + 0.5 * K3w) * Mf 
          cwdhdt   <- 1 / Hpw2 * dHpw2 / DT * (Cmpwi + K3mw)
          K4mw      <- cwFORM + cwDES + cwIRR - cwDPL - cwVOL - cwBIO - cwdhdt 
          K4mwdeg   <- -(cwVOL + cwBIO) * 10000 * Area * Hpw2
          #
          if (dDL != 1) {
            csPERC <- (Kmd / ((SatWC + bulk * Kmd) * dDL)) * perc * (Cmpwi + K3mw)
            csdDLdt  <- (Cmsdli + K3ms) * dfdDL / (dDL * DT)
          } else {
            csPERC <- (Kmd / ((SatWC + bulk * Kmd) * dDL)) * perc * ((Cmpwi + K3mw) - (Cmsdli + K3ms) / Kmd)
            csdDLdt  <- (Cmsdli + K3ms) * dfdDL / (dDL * DT)
          }
          csFORM    <- ffs * (Kmd / (SatWC + bulk * Kmd) * bulk * (kbios/24) * (Csdli + K3s)) * Mf   
          csBIO     <- Kmd / (SatWC + bulk * Kmd) * bulk * (kmbios/24) * (Cmsdli + K3ms)
          csDES     <- (Kmd / ((SatWC + bulk * Kmd) * dDL)) * (kmdiff/24) * ( (Cmsdli + K3ms) / Kmd - (Cmpwi + K3mw))
          K4ms      <- csFORM + csPERC - csBIO - csDES - csdDLdt
          K4msdeg   <- -csBIO / (Kmd / (SatWC + bulk * Kmd)) * 10000 * Area * dDL
  # 
  # For metabolite concentration
          Cmpw     <- Cmpwi + (K1mw + 2 * K2mw + 2 * K3mw + K4mw) / 6
          Cmsdl    <- Cmsdli + (K1ms + 2 * K2ms + 2 * K3ms + K4ms) / 6
          Mmpwdeg  <- (Mmpwdegi + (K1mwdeg + 2 * K2mwdeg + 2 * K3mwdeg + K4mwdeg) / 6) 
          Mmsdldeg <- (Mmsdldegi + (K1msdeg + 2 * K2msdeg + 2 * K3msdeg + K4msdeg) / 6) 
          if (Cmsdl < 0) {Cmsdl  <- 0}
          CmpwH    <- 0.5 * (Cmpwi + Cmpw)
          CmsdlH   <- 0.5 * (Cmsdli + Cmsdl)
        } # End of metabolite calculation
#
# Caluculate the pesticide concentration
        Cpw     <- Cpwi + (K1w + 2 * K2w + 2 * K3w + K4w) / 6
        if (Cpw < 0) {Cpw  <- 0}
        Csdl    <- Csdli + (K1s + 2 * K2s + 2 * K3s + K4s) / 6
        if (Csdl < 0) {Csdl  <- 0}
        Mpwdeg  <- (Mpwdegi + (K1wdeg + 2 * K2wdeg + 2 * K3wdeg + K4wdeg) / 6) 
        Msdldeg <- (Msdldegi + (K1sdeg + 2 * K2sdeg + 2 * K3sdeg + K4sdeg) / 6) 
        CpwH    <- 0.5 * (Cpwi + Cpw)
        CsdlH   <- 0.5 * (Csdli + Csdl)
#
# Mass balance calculation
        CMpwdeg    <- Mpwdeg 
        CMsdldeg   <- Msdldeg
        CMCdeg     <- Mpwdeg + Msdldeg
        PmassPW    <- 10000 * Area * Cpw * Hpw2
        PmassDL    <- 10000 * Area * Csdl * dDL * (bulk + SatWC / Kd)
        Pmass      <- PmassPW + PmassDL
        RFloss     <- 10000 * Area * drain * CpwH
        CRFloss    <- CRFloss + RFloss
        LSEEPloss  <- 10000 * Area * seep * CpwH
        CLSEEPloss <- CLSEEPloss + LSEEPloss
        if (dDL != 1) {
          PERCloss <- 0
        } else {
          PERCloss <- 10000 * Area * perc * CsdlH / Kd
        }
        CPERCloss  <- CPERCloss + PERCloss
        TPMdiss    <- Pmass + CRFloss + CPERCloss + CLSEEPloss - CMCdeg
#
# Dissolution terminator
        if (TPMdiss > PMass) {
          if (Fdiss == 0) break
          Fdiss      <- (PMass - TPMdissi) / (TPMdiss - TPMdissi)
          CRFloss    <- CRFloss - RFloss
          CLSEEPloss <- CLSEEPloss - LSEEPloss
          CPERCloss  <- CPERCloss - PERCloss   
        } else {
          break
        }
      } # End of loop for dissolution terminator 
#
# Mass balance calculation for metabolite
      if (meta > 0) { 
        CMmpwdeg      <- Mmpwdeg / Mf 
        CMmsdldeg     <- Mmsdldeg / Mf 
        CMmCdeg       <- Mmpwdeg + Mmsdldeg
        PmmassPW      <- 10000 * Area * Cmpw  / Mf * Hpw2
        PmmassDL      <- 10000 * Area * Cmsdl / Mf * dDL * (bulk + SatWC / Kmd)
        Pmmass        <- PmmassPW + PmmassDL
        RFmloss       <- 10000 * Area * drain * CmpwH / Mf 
        CRFmloss      <- CRFmloss + RFmloss
        LSEEPmloss    <- 10000 * Area * seep * CmpwH / Mf 
        CLSEEPmloss   <- CLSEEPmloss + LSEEPmloss
        if (dDL != 1) {
          PERCmloss   <- 0 
        } else {
          PERCmloss   <- 10000 * Area * perc * CmsdlH / Kmd / Mf 
        }
        CPERCmloss    <- CPERCmloss + PERCmloss
        TPMmdiss      <- Pmmass + CRFmloss + CPERCmloss + CLSEEPmloss - CMmCdeg
      }
      PmassDLD   <- PmassDL + PMass - TPMdiss
#
      CsdlD      <- PmassDLD / (10000 * Area * dDL * (bulk + SatWC / Kd))
      ACsdl      <- Csdl * dDL
      ACsdlD     <- CsdlD * dDL
      if (Fdiss < 1) {
        PMdiss <- PMass
        if (CUMMass > 0) PMdiss <- CUMMass + SMass
        if ((type==4) && (SMass == 0)) PMdiss <- DriMass * APPMass
        Fdiss  <- 0
      } else {
        PMdiss <- TPMdiss
      }
      MBE        <- (PMdiss- TPMdiss)  / PMdiss * 100
      if (meta > 0) {
        as.list(c(Fdiss=Fdiss, PMdiss=PMdiss, TPMdiss=TPMdiss, MBE=MBE, PmassDL=PmassDL, PmassPW=PmassPW,
                  Cpwi=Cpwi, Csdli=Csdli, Cpw=Cpw, Csdl=Csdl, CsdlD=CsdlD, Mpwdeg=Mpwdeg, Msdldeg=Msdldeg,
                  CMpwdeg= CMpwdeg,CMsdldeg=CMsdldeg,CMCdeg=CMCdeg,APPMass=APPMass,ACsdl=ACsdl,ACsdlD=ACsdlD,
                  PmassPW=PmassPW, PmassDL= PmassDL,PmassDLD=PmassDLD,Pmass=Pmass,CRFloss=CRFloss,
                  CLSEEPloss=CLSEEPloss,CPERCloss=CPERCloss,Cmpwi=Cmpwi, Cmsdli=Cmsdli,Cmpw=Cmpw, Cmsdl=Cmsdl,
                  Mmpwdeg=Mmpwdeg, Mmsdldeg=Mmsdldeg,CMmCdeg=CMmCdeg,PmmassPW=PmmassPW, PmmassDL= PmmassDL,
                  Pmmass=Pmmass,CRFmloss=CRFmloss,CLSEEPmloss=CLSEEPmloss,CPERCmloss=CPERCmloss,TPMmdiss=TPMmdiss,
                  CMmpwdeg= CMmpwdeg,CMmsdldeg=CMmsdldeg,SMass=SMass, DriMass = DriMass
                ))
      } else {       
        as.list(c(Fdiss=Fdiss, PMdiss=PMdiss,TPMdiss=TPMdiss,MBE=MBE, PmassDL=PmassDL,PmassPW=PmassPW,
                  Cpwi=Cpwi, Csdli=Csdli,Cpw=Cpw, Csdl=Csdl,CsdlD=CsdlD,Mpwdeg=Mpwdeg, Msdldeg=Msdldeg,
                  CMpwdeg= CMpwdeg,CMsdldeg=CMsdldeg,CMCdeg=CMCdeg,APPMass=APPMass,ACsdl=ACsdl,ACsdlD=ACsdlD,
                  PmassPW=PmassPW, PmassDL= PmassDL,PmassDLD=PmassDLD,Pmass=Pmass,CRFloss=CRFloss,
                  CLSEEPloss=CLSEEPloss,CPERCloss=CPERCloss,SMass=SMass, DriMass = DriMass
                ))
      }      
  })  # End of the statement "with"         
} # End of sub-function "DISS"

# --------------------------------------- End of sub-functions ------------------------------------------------
# Combine parameter set
para <- c(parms, parms2)


# Start time loop for the computation of PCPF-1
for (i in 1:DMAX){
   if (i == 1) {
        WB  <- with(as.data.frame(water), {
                 Hpw1   <-  water[i   ,"h"]
                 perc   <-  water[i+1 ,"perc"] / 24
                 rain   <-  water[i+1 ,"rain"] / 24
                 irr    <-  water[i+1 ,"irr"] / 24
                 drain  <-  water[i+1 ,"drain"] / 24
                 seep   <-  water[i+1 ,"seep"] / 24
                 et     <-  water[i+1 ,"et"] / 24
                 Cj     <-  rain + irr - perc - seep - et - drain
                 as.list(c(Hpw1=Hpw1, perc = perc, rain = rain, irr = irr, drain = drain, 
                           seep = seep, et = et, Cj = Cj))      
        })
        Hpw1 <- with(as.list(WB), Hpw1)          
        BC <- state
   } else if (i %% 24 == 1) {
        WB  <- with(as.data.frame(water), {
                 perc   <-  water[(i %/% 24)+2 ,"perc"] / 24
                 rain   <-  water[(i %/% 24)+2 ,"rain"] / 24
                 irr    <-  water[(i %/% 24)+2 ,"irr"] / 24
                 drain  <-  water[(i %/% 24)+2 ,"drain"] / 24
                 seep   <-  water[(i %/% 24)+2 ,"seep"] / 24
                 et     <-  water[(i %/% 24)+2 ,"et"] / 24
                 Cj     <-  rain + irr - perc - seep - et - drain
                 as.list(c(perc = perc, rain = rain, irr = irr, drain = drain, 
                           seep = seep, et = et, Cj = Cj))      
        })
        Hpw1 <- with(as.list(BC), Hpw1)
   } else {   
        Hpw1 <- with(as.list(BC), Hpw1)          
   }
   Cj   <- with(as.list(WB), Cj)
   Hpw2 <- Hpw1 + Cj

# Forward difference:water balance
   HpwH  <- (Hpw1 + Hpw2) * 0.5
   dHpw1 <- Cj
   dHpw2 <- Cj
   dHpwH <- (dHpw1 + dHpw2) * 0.5

# Advance solute front via percolation
   Dperc <-with(as.list(c(WB,BC)), Dperci + perc)
   Cfront <-with(as.list(c(WB,para)), Dperc / SatWC)
   if (Cfront > 1) {
     dDL <- 1
   } else {
     dDL <- Cfront
   }
   dDLH  <- with(as.list(BC), 0.5 * (dDL + dDLi))
   dfdDL <- with(as.list(BC), dDL - dDLi)
   dfdDLH <- with(as.list(BC), dDL - dDLi)
   wbdata <- cbind(Hpw1,Hpw2,HpwH,dHpw1,dHpw2,dHpwH,dDL,dDLH,dfdDL,dfdDLH)

# Multiple application option
   if (mapp > 0) {
     if (i == 1) {
       mapp          <-  mapp - 1
       AppRd         <-  with(as.data.frame(mapplist), mapplist[1, "appmass"])
       APPMass       <-  with(as.list(para), AppRd * Area * 1000000)
       AppRd         <-  as.list(AppRd)
       names(AppRd)  <-  "AppRd"
       para          <-  c(para, AppRd)
       APPTime       <-  with(as.data.frame(mapplist), mapplist[2, "apptime"])
       m             <-  2
     } else if (i/24 == APPTime) {
       mapp          <-  mapp - 1 
       if (mapp == 0) mapp <- -1
       para          <-  para[-length(para)] 
       AppRd         <-  with(as.data.frame(mapplist), mapplist[m, "appmass"])
       APPMassd      <-  with(as.list(para), AppRd * Area * 1000000)
       CUMMass       <-  with(as.list(BC), APPMass * DriMass)
       AppRd         <-  as.list(AppRd)
       names(AppRd)  <-  "AppRd"
       para          <-  c(para, AppRd)  
       if (mapp > 0) {
         APPTime   <-  with(as.data.frame(mapplist), mapplist[m + 1, "apptime"])
         m         <-  m + 1
         if (mapp == 0) CUMMass  <-  with(as.data.frame(mapplist), sum(mapplist$apptime))
       }
       Fdiss         <- with(as.list(BC), Fdiss)
       if (Fdiss > 0) {
         APPMass   <-  APPMass + APPMassd
       } else {
         APPMass   <-  APPMass + APPMassd
         BC[4]     <-  1
       }
     }
   }
#browser()
# pH dependence of hydrolysis parameter
   if (pH == TRUE) {      
     pH1 <- with(as.data.frame(meteolist), meteolist[i+1,"pH"])
     khyd_ref   <-  with(as.list(para), ka_ref * 10^-pH1 + kb_ref * 10^pH1 + kn_ref)
     if (exists('khyd_ref', where=para) == TRUE) {
       para["khyd_ref"] <- khyd_ref
     } else {
       para_d <- list(khyd_ref)
       names(para_d) <- "khyd_ref"
       para <- c(para, para_d)
     }
   }

# Temperature dependence of degradation parameter
   if (temp == TRUE) {
     T1 <- with(as.data.frame(meteolist), meteolist[i+1,"temp"])
     kbiow   <-  with(as.list(para), kbiow_ref * exp(Ea/0.008314*(1/(273.14+Tdef)-1/(273.14+T1))))
     if (exists('kbiow', where=para) == TRUE) {
       para["kbiow"] <- kbiow
     } else {
       para_d <- list(kbiow)
       names(para_d) <- "kbiow"
       para <- c(para, para_d)
     }
     khyd   <-  with(as.list(para), khyd_ref * exp(Ea/0.008314*(1/(273.14+Tdef)-1/(273.14+T1))))
     if (exists('khyd', where=para) == TRUE) {
       para["khyd"] <- khyd
     } else {
       para_d <- list(khyd)
       names(para_d) <- "khyd"
       para <- c(para, para_d)
     }
     parms["khyd"] <- khyd
   }          

# UV dependence of degradation parameter
   if (uv == TRUE) {
     uv1 <- with(as.data.frame(meteolist), meteolist[i+1,"uv"])
     kphot   <-  with(as.list(c(parms, parms2)), fcut * kphot_ref * uv1/I_ref)
     if (exists('kphot', where=para) == TRUE) {
       para["kphot"] <- kphot
     } else {
       para_d <- list(kphot)
       names(para_d) <- "kphot"
       para <- c(para, para_d)
     }
   }     

# solve mass balane
   x <- DISS(parms=para, wbdata1=WB,wbdata2=wbdata,state=BC,APPMass=APPMass,
             CUMMass=CUMMass,DriMass=DriMass, mapp=mapp)

# Convert sorption parameters 
   if ((type == 2) || (type == 4)) {
     alp <- with(as.list(c(parms, parms2)), alp)
     if (alp != 0) {
       Cpw <-  with(as.list(c(parms, parms2, x)), f * Kd * Cpw)
       Cs  <-  with(as.list(c(parms, parms2, x)), Csdl)
       if (Cs >= Cpw) {
         if ("alp" %in% names(parms) == TRUE) {
           parms["alp"] <- 0
         } else {
           parms2["alp"] <- 0
         }
       }
     }
   }

# Model output
   if (i == 1) {
     if (meta > 0) {
       x0 <- with(as.list(c(x,para,WB,BC)), {
            Time         <-  0  
            Hpw          <-  Hpw1
            dDL          <-  dDLi
            Drain        <-  0
            Perc         <-  0
            Seep         <-  0
            Cpw          <-  0
            Csdl         <-  0
            Cmpw         <-  0
            Cmsdl        <-  0
            ACsdl        <-  0
            CsdlD        <-  0
            ACsdlD       <-  0
            Mformdiss    <-  0
            TPMmdiss     <-  0
            PmmassPW     <-  0
            PmmassDL     <-  0
            CMmCdeg      <-  0
            CRFmloss     <-  0
            CPERCmloss   <-  0
            CLSEEPmloss  <-  0
            if ((type == 4) && (SMass > 0)) {
              SMassC       <-  0
              DriMass      <-  0 
            }
            TPMdiss      <-  0
            PmassPW      <-  0
            PmassDL      <-  0
            CMCdeg       <-  0
            CRFloss      <-  0
            CPERCloss    <-  0
            CLSEEPloss   <-  0
            MBE          <-  0
            if (fdegw == 2) {
              khyd         <-  khyd
              kbiow        <-  kbiow
              kphot        <-  kphot
              kbulk        <-  khyd + kbiow + kphot
            } 
            if (type == 4) {
              as.list(c(Time        =  Time, 
                        Hpw         =  Hpw, 
                        dDL         =  dDL,
                        Drain       =  drain,
                        Perc        =  perc,
                        Seep        =  seep, 
                        Cpw         =  Cpw, 
                        Csdl        =  Csdl,
                        ACsdl       =  ACsdl, 
                        TPMdiss     =  TPMdiss,
                        PmassPW     =  PmassPW,
                        PmassDL     =  PmassDL, 
                        CMCdeg      =  CMCdeg, 
                        CRFloss     =  CRFloss,
                        CPERCloss   =  CPERCloss, 
                        DriMass     =  DriMass,
                        MBE         =  MBE, 
                        Fdiss       =  Fdiss,
                        Cmpw        =  Cmpw,
                        Cmsdl       =  Cmsdl,
                        TPMmdiss    =  TPMmdiss,
                        PmmassPW    =  PmmassPW, 
                        PmmassDL    =  PmmassDL,
                        CMmCdeg     =  CMmCdeg, 
                        CRFmloss    =  CRFmloss,
                        CPERCmloss  =  CPERCmloss,  
                        Mformdiss   =  Mformdiss
                     ))
            } else {
              as.list(c(Time        =  Time, 
                        Hpw         =  Hpw, 
                        dDL         =  dDL,
                        Drain       =  drain,
                        Perc        =  perc,
                        Seep        =  seep, 
                        Cpw         =  Cpw, 
                        Csdl        =  Csdl,
                        ACsdl       =  ACsdl, 
                        TPMdiss     =  TPMdiss,
                        PmassPW     =  PmassPW,
                        PmassDL     =  PmassDL, 
                        CMCdeg      =  CMCdeg, 
                        CRFloss     =  CRFloss,
                        CPERCloss   =  CPERCloss, 
                        MBE         =  MBE, 
                        Fdiss       =  Fdiss,
                        Cmpw        =  Cmpw,
                        Cmsdl       =  Cmsdl,
                        TPMmdiss    =  TPMmdiss,
                        PmmassPW    =  PmmassPW, 
                        PmmassDL    =  PmmassDL,
                        CMmCdeg     =  CMmCdeg, 
                        CRFmloss    =  CRFmloss,
                        CPERCmloss  =  CPERCmloss,  
                        Mformdiss   =  Mformdiss
                     ))
            }   
       }) 
     } else {
       x0 <- with(as.list(c(x,para,WB,BC)), {
            Time         <-  0  
            Hpw          <-  Hpw1
            dDL          <-  dDLi
            Drain        <-  0
            Perc         <-  0
            Seep         <-  0 
            Cpw          <-  0
            Csdl         <-  0
            ACsdl        <-  0
            CsdlD        <-  0
            ACsdlD       <-  0
            if ((type == 4) && (SMass > 0)) {
              SMassC       <-  0
              DriMass      <-  0
            }
            TPMdiss      <-  0
            PmassPW      <-  0
            PmassDL      <-  0
            CMCdeg       <-  0
            CRFloss      <-  0
            CPERCloss    <-  0
            CLSEEPloss   <-  0
            MBE          <-  0
            if (fdegw == 2) {
              khyd         <-  khyd
              kbiow        <-  kbiow
              kphot        <-  kphot
              kbulk        <-  khyd + kbiow + kphot
            }
            if (type == 4) {
              as.list(c(Time      =  Time,
                        Hpw       =  Hpw, 
                        dDL       =  dDL,
                        Drain     =  drain,
                        Perc      =  perc,
                        Seep      =  seep,  
                        Cpw       =  Cpw,
                        Csdl      =  Csdl,
                        ACsdl     =  ACsdl,
                        TPMdiss   =  TPMdiss,
                        PmassPW   =  PmassPW,
                        PmassDL   =  PmassDL,
                        CMCdeg    =  CMCdeg, 
                        CRFloss   =  CRFloss,
                        CPERCloss =  CPERCloss,
                        DriMass   =  DriMass * 100, 
                        MBE       =  MBE,
                        Fdiss     =  Fdiss
                     ))
            } else {
              as.list(c(Time      =  Time,
                        Hpw       =  Hpw,
                        dDL       =  dDL,
                        Drain     =  drain,
                        Perc      =  perc,
                        Seep      =  seep,  
                        Cpw       =  Cpw,
                        Csdl      =  Csdl,
                        ACsdl     =  ACsdl,
                        TPMdiss   =  TPMdiss,
                        PmassPW   =  PmassPW,
                        PmassDL   =  PmassDL,
                        CMCdeg    =  CMCdeg, 
                        CRFloss   =  CRFloss,
                        CPERCloss =  CPERCloss, 
                        MBE       =  MBE,
                        Fdiss     =  Fdiss,
                        khyd      =  khyd,
                        kbiow     =  kbiow,
                        kphot     =  kphot,
                        kbulk     =  khyd + kbiow + kphot 
                     ))
            }        
       }) 
     }
   }

   if (meta > 0) {
     xx <- with(as.list(c(x,para,WB)), {
          Time         <-  i / 24  
          Hpw          <-  Hpw2
          dDL          <-  dDL
          Drain        <-  drain
          Perc         <-  perc
          Seep         <-  seep 
          Cpw          <-  Cpw
          Csdl         <-  Csdl
          Cmpw         <-  Cmpw
          Cmsdl        <-  Cmsdl
          ACsdl        <-  ACsdl
          CsdlD        <-  CsdlD
          ACsdlD       <-  ACsdlD
          Mformdiss    <-  TPMmdiss / APPMass * 100
          TPMmdiss     <-  -TPMmdiss / CMCdeg * 100
          PmmassPW     <-  -PmmassPW / CMCdeg * 100
          PmmassDL     <-  -PmmassDL / CMCdeg * 100
          CMmCdeg      <-  CMmCdeg / CMCdeg * 100
          CRFmloss     <-  -CRFmloss / CMCdeg * 100
          CPERCmloss   <-  -CPERCmloss / CMCdeg * 100
          CLSEEPmloss  <-  -CLSEEPmloss / CMCdeg * 100
          if ((type == 4) && (SMass > 0)) {
            SMassC       <-  PMdiss
            DriMass      <-  SMassC / APPMass *100 
          }
          TPMdiss      <-  TPMdiss / APPMass * 100
          PmassPW      <-  PmassPW / APPMass * 100
          PmassDL      <-  PmassDL / APPMass * 100
          CMCdeg       <-  CMCdeg / APPMass * 100
          CRFloss      <-  CRFloss / APPMass * 100
          CPERCloss    <-  CPERCloss / APPMass * 100
          CLSEEPloss   <-  CLSEEPloss / APPMass * 100
          MBE          <-  MBE
          if (fdegw == 2) {
            khyd         <-  khyd
            kbiow        <-  kbiow
            kphot        <-  kphot
            kbulk        <-  khyd + kbiow + kphot
          }
          if (type == 4) {
            as.list(c(Time        =  Time, 
                      Hpw         =  Hpw, 
                      dDL         =  dDL,
                      Drain       =  drain,
                      Perc        =  perc,
                      Seep        =  seep, 
                      Cpw         =  Cpw, 
                      Csdl        =  Csdl,
                      ACsdl       =  ACsdl, 
                      TPMdiss     =  TPMdiss,
                      PmassPW     =  PmassPW,
                      PmassDL     =  PmassDL, 
                      CMCdeg      =  CMCdeg, 
                      CRFloss     =  CRFloss,
                      CPERCloss   =  CPERCloss, 
                      DriMass     =  DriMass,
                      MBE         =  MBE, 
                      Fdiss       =  Fdiss,
                      Cmpw        =  Cmpw,
                      Cmsdl       =  Cmsdl,
                      TPMmdiss    =  TPMmdiss,
                      PmmassPW    =  PmmassPW, 
                      PmmassDL    =  PmmassDL,
                      CMmCdeg     =  CMmCdeg, 
                      CRFmloss    =  CRFmloss,
                      CPERCmloss  =  CPERCmloss,  
                      Mformdiss   =  Mformdiss
                   ))
          } else {
            as.list(c(Time        =  Time, 
                      Hpw         =  Hpw, 
                      dDL         =  dDL,
                      Drain       =  drain,
                      Perc        =  perc,
                      Seep        =  seep, 
                      Cpw         =  Cpw, 
                      Csdl        =  Csdl,
                      ACsdl       =  ACsdl, 
                      TPMdiss     =  TPMdiss,
                      PmassPW     =  PmassPW,
                      PmassDL     =  PmassDL, 
                      CMCdeg      =  CMCdeg, 
                      CRFloss     =  CRFloss,
                      CPERCloss   =  CPERCloss, 
                      MBE         =  MBE, 
                      Fdiss       =  Fdiss,
                      Cmpw        =  Cmpw,
                      Cmsdl       =  Cmsdl,
                      TPMmdiss    =  TPMmdiss,
                      PmmassPW    =  PmmassPW, 
                      PmmassDL    =  PmmassDL,
                      CMmCdeg     =  CMmCdeg, 
                      CRFmloss    =  CRFmloss,
                      CPERCmloss  =  CPERCmloss,  
                      Mformdiss   =  Mformdiss
                   ))
          }   
     }) 
   } else {
     xx <- with(as.list(c(x,para,WB)), {
          Time         <-  i / 24  
          Hpw          <-  Hpw2
          dDL          <-  dDL
          Drain        <-  drain
          Perc         <-  perc
          Seep         <-  seep 
          Cpw          <-  Cpw
          Csdl         <-  Csdl
          ACsdl        <-  ACsdl
          CsdlD        <-  CsdlD
          ACsdlD       <-  ACsdlD
          if ((type == 4) && (SMass > 0)) {
            SMassC       <-  PMdiss
            DriMass      <-  SMassC / APPMass
          }
          TPMdiss      <-  (TPMdiss / APPMass * 100) / DriMass
          PmassPW      <-  (PmassPW / APPMass * 100) / DriMass
          PmassDL      <-  (PmassDL / APPMass * 100) / DriMass
          CMCdeg       <-  (CMCdeg / APPMass * 100) / DriMass
          CRFloss      <-  (CRFloss / APPMass * 100) / DriMass
          CPERCloss    <-  (CPERCloss / APPMass * 100) / DriMass
          CLSEEPloss   <-  (CLSEEPloss / APPMass * 100) / DriMass
          MBE          <-  MBE
          if (fdegw == 2) {
            khyd         <-  khyd
            kbiow        <-  kbiow
            kphot        <-  kphot
            kbulk        <-  khyd + kbiow + kphot
          }
          if (type == 4) {
            as.list(c(Time      =  Time,
                      Hpw       =  Hpw, 
                      dDL       =  dDL,
                      Drain     =  drain,
                      Perc      =  perc,
                      Seep      =  seep,  
                      Cpw       =  Cpw,
                      Csdl      =  Csdl,
                      ACsdl     =  ACsdl,
                      TPMdiss   =  TPMdiss,
                      PmassPW   =  PmassPW,
                      PmassDL   =  PmassDL,
                      CMCdeg    =  CMCdeg, 
                      CRFloss   =  CRFloss,
                      CPERCloss =  CPERCloss,
                      DriMass   =  DriMass * 100, 
                      MBE       =  MBE,
                      Fdiss     =  Fdiss
                   ))
          } else {
            as.list(c(Time      =  Time,
                      Hpw       =  Hpw,
                      dDL       =  dDL,
                      Drain     =  drain,
                      Perc      =  perc,
                      Seep      =  seep,  
                      Cpw       =  Cpw,
                      Csdl      =  Csdl,
                      ACsdl     =  ACsdl,
                      TPMdiss   =  TPMdiss,
                      PmassPW   =  PmassPW,
                      PmassDL   =  PmassDL,
                      CMCdeg    =  CMCdeg, 
                      CRFloss   =  CRFloss,
                      CPERCloss =  CPERCloss, 
                      MBE       =  MBE,
                      Fdiss     =  Fdiss,
                      khyd      =  khyd,
                      kbiow     =  kbiow,
                      kphot     =  kphot,
                      kbulk     =  khyd + kbiow + kphot 
                   ))
          }        
     }) 
   }

# Update the boundary condition
   if (meta > 0) {
     BC <- with(as.list(x), {
          CpwH         <-  0.5 * (Cpwi + Cpw)
          CsdlH        <-  0.5 * (Csdli + Csdl)        
          Hpw1         <-  Hpw2
          Cpwi         <-  Cpw
          Csdli        <-  Csdl
          Cmpwi        <-  Cmpw
          Cmsdli       <-  Cmsdl
          Mpwdegi      <-  Mpwdeg
          Msdldegi     <-  Msdldeg
          Mmpwdegi     <-  Mmpwdeg
          Mmsdldegi    <-  Mmsdldeg
          PMdissi      <-  PMdiss
          TPMdissi     <-  TPMdiss
          CRFloss      <-  CRFloss
          CPERCloss    <-  CPERCloss
          CLSEEPloss   <-  CLSEEPloss
          CMpwdeg      <-  CMpwdeg
          CMsdldeg     <-  CMsdldeg
          CMCdeg       <-  CMCdeg
          CRFmloss     <-  CRFmloss
          CPERCmloss   <-  CPERCmloss
          CLSEEPmloss  <-  CLSEEPmloss
          CMmpwdeg     <-  CMmpwdeg
          CMmsdldeg    <-  CMmsdldeg
          CMmCdeg      <-  CMmCdeg
          Dperci       <-  Dperc
          dDLi         <-  dDL
          Fdiss        <-  Fdiss
          dfdDLi       <-  dfdDL
          if ((type == 4) && (SMass > 0)) {
            SMassC       <-  PMdiss
            DriMass      <-  SMassC / APPMass
          }
          as.list(c(Hpw1         =  Hpw1,
                    Cpwi         =  Cpwi,
                    Csdli        =  Csdli,
                    Fdiss        =  Fdiss,
                    Mpwdegi      =  Mpwdegi,
                    Msdldegi     =  Msdldegi,
                    PMdissi      =  PMdissi,
                    TPMdissi     =  TPMdissi,
                    CRFloss      =  CRFloss,
                    Cmpwi        =  Cmpwi,
                    Cmsdli       =  Cmsdli,
                    Mmpwdegi     =  Mmpwdegi,
                    Mmsdldegi    =  Mmsdldegi,
                    CPERCloss    =  CPERCloss, 
                    CLSEEPloss   =  CLSEEPloss,
                    CMpwdeg      =  CMpwdeg,
                    CMsdldeg     =  CMsdldeg,
                    CMCdeg       =  CMCdeg,
                    Dperci       =  Dperci,
                    dDLi         =  dDLi,
                    dfdDLi       =  dfdDLi,
                    CPERCmloss   =  CPERCmloss, 
                    CLSEEPmloss  =  CLSEEPmloss,
                    CMmpwdeg     =  CMmpwdeg,
                    CMmsdldeg    =  CMmsdldeg,
                    CMmCdeg      =  CMmCdeg,
                    CRFmloss     =  CRFmloss,
                    DriMass      =  DriMass
                 ))      
     })
   } else {           
     BC <- with(as.list(x), {
          CpwH       <-  0.5 * (Cpwi + Cpw)
          CsdlH      <-  0.5 * (Csdli + Csdl)        
          Hpw1       <-  Hpw2
          Cpwi       <-  Cpw
          Csdli      <-  Csdl
          Mpwdegi    <-  Mpwdeg
          Msdldegi   <-  Msdldeg
          PMdissi    <-  PMdiss
          TPMdissi   <-  TPMdiss
          CRFloss    <-  CRFloss
          CPERCloss  <-  CPERCloss
          CLSEEPloss <-  CLSEEPloss
          CMpwdeg    <-  CMpwdeg
          CMsdldeg   <-  CMsdldeg
          CMCdeg     <-  CMCdeg
          Dperci     <-  Dperc
          dDLi       <-  dDL
          Fdiss      <-  Fdiss
          dfdDLi     <-  dfdDL
          if ((type == 4) && (SMass > 0)) {
            SMassC       <-  PMdiss
            DriMass      <-  SMassC / APPMass
          }
          as.list(c(Hpw1        =  Hpw1,
                    Cpwi        =  Cpwi,
                    Csdli       =  Csdli,
                    Fdiss       =  Fdiss,
                    Mpwdegi     =  Mpwdegi,
                    Msdldegi    =  Msdldegi,
                    PMdissi     =  PMdissi,
                    TPMdissi    =  TPMdissi,
                    CRFloss     =  CRFloss,
                    CPERCloss   =  CPERCloss,
                    CLSEEPloss  =  CLSEEPloss,
                    CMpwdeg     =  CMpwdeg,
                    CMsdldeg    =  CMsdldeg,
                    CMCdeg      =  CMCdeg,
                    Dperci      =  Dperci,
                    dDLi        =  dDLi,
                    dfdDLi      =  dfdDLi,
                    DriMass     =  DriMass
                 ))      
     })
   } 
   if (rBC == TRUE) {
     out <- BC
   } else {       
     if (i == 1) {
        out <- as.data.frame(x0)
        out <- rbind(out,xx)
     } else {
        out <- rbind(out,xx)
     }
   } 
 }
 return(out)
}

# define data splitting function for PCPF-1R output
XX <- function(X){
      a <-X[X[,"Time"] ==  0.125, c("Time","Cpw")]
      out <- a
      b <-X[X[,"Time"] ==  1, c("Time","Cpw")]
      out <-rbind(out,b)
      c <-X[X[,"Time"] ==  2, c("Time","Cpw")]
      out <-rbind(out,c)
      d <-X[X[,"Time"] ==  3, c("Time","Cpw")]
      out <-rbind(out,d)
      e <-X[X[,"Time"] ==  5, c("Time","Cpw")]
      out <-rbind(out,e)
      f <-X[X[,"Time"] ==  7, c("Time","Cpw")]
      out <-rbind(out,f)
      g <-X[X[,"Time"] ==  8, c("Time","Cpw")]
      out <-rbind(out,g)
      h <-X[X[,"Time"] ==  10, c("Time","Cpw")]
      out <-rbind(out,h)
      i <-X[X[,"Time"] ==  14, c("Time","Cpw")]
      out <-rbind(out,i)
      j <-X[X[,"Time"] ==  21, c("Time","Cpw")]
      out <-rbind(out,j)
      return(out)
      } 

# define data splitting function for sensRange output
XX2 <- function(X){
       a <-X[X[,"x"] ==  0.125, c("x","q05")]
       out <- a
       b <-X[X[,"x"] ==  1, c("x","q05")]
       out <-rbind(out,b)
       c <-X[X[,"x"] ==  2, c("x","q05")]
       out <-rbind(out,c)
       d <-X[X[,"x"] ==  3, c("x","q05")]
       out <-rbind(out,d)
       e <-X[X[,"x"] ==  5, c("x","q05")]
       out <-rbind(out,e)
       f <-X[X[,"x"] ==  7, c("x","q05")]
       out <-rbind(out,f)
       g <-X[X[,"x"] ==  8, c("x","q05")]
       out <-rbind(out,g)
       h <-X[X[,"x"] ==  10, c("x","q05")]
       out <-rbind(out,h)
       i <-X[X[,"x"] ==  14, c("x","q05")]
       out <-rbind(out,i)
       j <-X[X[,"x"] ==  21, c("x","q05")]
       out <-rbind(out,j)
return(out)
} 

XXL <- function(X,Y,Z){
      for(i in 1:nrow(Y)){
        if (i == 1){
          a <-X[X[,Z[1]] ==  Y[i,1], Z]
          out <- a
        } else {
          b <-X[X[,Z[1]] ==  Y[i,1], Z]
          out <-rbind(out,b)
        }
      }
      return(out)
}
