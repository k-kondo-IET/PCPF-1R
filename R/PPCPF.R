# Sub-function for pre-processing of raw data
# Setting up the basic parameter set for PCPF-1R modeling

PPCPF <- function(input_EXP, input_PHYSCHEM, input_LAB){
  with(as.list(c(input_EXP, input_PHYSCHEM, input_LAB)), {
    
    # Initialization
    kxfer <- 0.00000001 # empirical coefficient of overall water-sediment mass transfer
    
    # Parameterization
    AppR  <-  APPLICATION*CONTENT/100/1000      # Application rate (g/m2)
    Cmax  <-  AppR*1000/50                      # Theoretical max conc. with 5cm water depth (mg/L)
    Henry <-  16.04*MW[1]*VP[1]*(1/133.3)/(CSLB[1]*298)  # Dimensionless Henry's law constant (-)
    kvol  <-  (1/(4.752*(44/MW[1])^(1/2))+1/(Henry*(720*(18/MW[1])^(1/2))))^(-1)  # Mass transfer rate (m/day) 
    Kd    <-  median(KF*Cmax^(FN-1))            # Equilibrium partitioning coefficient (L/kg)
    kdiff <-  69.35/365*SAT_WC*MW[1]^(-2/3)     # Diffusion rate constant (m/day)
    ksorp <-  kxfer*(SAT_WC+BULK_D*Kd)*86400    # First-order sorption rate constant (1/day)
    om    <-  kdiff+1*BULK_D*Kd*ksorp           # Overall mass transfer constant (m/day)
    khyd  <-  log(2)/DT50_HYD[1]
    kphot <-  log(2)/DT50_PHT[1]
    kbiow <-  log(2)/DT50_BIOW[1]
    kbulk <-  khyd+kphot+kbiow                  # First-order bulk degradation rate constant in PW (1/day)
    kbios <-  log(2)/DT50_BIOS[1]               # First-order degradation rate constant in PSL (1/day)
    
    # Create parameter list
    out   <-  list(
                name   =  NAME[1],
                plot   =  PLOT,
                Area   =  AREA,
                bulk   =  BULK_D,
                SatWC  =  SAT_WC,
                AppR   =  AppR,
                kdiss  =  0.01,                 # First-order dissolution rate constant (1/day)
                alp    =  1,                    # First-order mixing rate constant (1/day)
                CSLB   =  CSLB[1],
                Kd     =  Kd,
                f      =  1,                    # Fraction for apparent equilibrium partitioning 
                kdiff  =  kdiff,
                ksorp  =  ksorp,
                om     =  om,
                kvol   =  kvol,
                khyd   =  khyd,
                kphot  =  kphot,
                kbiow  =  kbiow,
                kbulk  =  kbulk,
                kbios  =  kbios
    )
    
    return(out)
  })
}