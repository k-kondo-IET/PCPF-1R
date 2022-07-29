PCPFBR <- function(pcpfin, waterin, blockin, flow, meteo, DMAX, type){
  # start main computation
  # Determine pesticide application schedule
  with(as.list(c(pcpfin, waterin, blockin)), {
    INI   <-  flow[flow$date == Dcal, 1]
    Dapp  <-  as.Date(Dapp)
    Dcal  <-  as.Date(Dcal)
    z     <-  as.numeric(Dapp - Dcal)
    zz    <-  z + AppP - 1
    AppA  <-  NULL
    DOHA  <-  NULL
    APFA  <-  TArea * pfr 
    for (j in z:zz){
      AppA  <- append(AppA, APFA*(1/(AppStdev*sqrt(2*3.14)))*exp(-0.5*(j-z-AppMean)*(j-z-AppMean)/(AppStdev*AppStdev)))
      DOHA  <- append(DOHA, j)
    }
    coeff <-  APFA / sum(AppA)
    AppA  <-  AppA * coeff
    pas   <-  cbind(DOHA, AppA)
    #
    # Number of PCPF-1R run 
    AppEvent   <-  length(pas[,1])
    # Start block computation
    for (i in 1:AppEvent){
      DST   <-  as.numeric(pas[i,1])
      # water balance computation
      meteoin   <-  meteo[DST:(DST+DMAX-1),]
      DSIM      <-  DMAX-DST+1
      wb        <-  WBcalc (pars=waterin, DMAX=DSIM, meteo=meteoin)
      # Plot scale computation
      pcpfin$Area    <-  as.numeric(pas[i,2])*usage
      out            <-  PCPF1R(pcpfin, pcpfin, water=wb, DMAX=DSIM, fdes=1, type=type)
      Y              <-  as.data.frame(seq(1:DSIM))
      Z1             <-  c("Time","Cpw","CRFloss","CSEEPloss")
      outd           <-  XXL(out, Y, Z1)
      DASS           <-  as.vector(Y)+DST-1
      Conc           <-  as.data.frame(outd$Cpw)
      AppArea        <-  as.data.frame(c(pcpfin$Area, rep(0,nrow(Y)-1)))
      AppMass        <-  pcpfin$AppR*pcpfin$Area
      AppAmount      <-  as.data.frame(c(AppMass, rep(0,nrow(Y)-1)))
      DRFloss        <-  as.data.frame(c(outd$CRFloss[1]*AppMass/100, diff(outd$CRFloss)*AppMass/100))
      DLSEEPloss     <-  as.data.frame(c(outd$CSEEPloss[1]*AppMass/100/Klevee, diff(outd$CSEEPloss)*AppMass/100/Klevee))
      outb           <-  cbind(DASS,Conc,AppArea,AppAmount,DRFloss,DLSEEPloss)
      colnames(outb) <-  c("DASS","Conc","AppArea","AppAmount","DRFloss","DLSEEPloss")
      if (i == 1){
        x <- outb
      } else {
        x <- rbind(x,outb)
      }
    } # end of for loop of i
    for (i in 1:DMAX){    # create final data
      if (length(x[x$DASS == i, 1]) == 0) {
        xx <- as.list(c(DASS = i, AppArea = 0, AppAmount=0, DRFLoss = 0, DLSEEPLoss = 0))
      } else {
        x2 <-  apply(x[x$DASS == i,], 2, sum)
        xx <- as.list(c(DASS = i, 
                        AppArea = as.numeric(x2[3]), 
                        AppAmount = as.numeric(x2[4]), 
                        DRFLoss = as.numeric(x2[5]), 
                        DLSEEPLoss = as.numeric(x2[6])))
      }
      if (i == 1) { 
        list <- as.data.frame(xx)
      } else {
        list <- rbind(list,xx)
      }
    } # end of for loop of i
    # extract flow data
    END  <- INI+DMAX-1
    if (blockin$hflow == 0){
      FD   <- flow[INI:END, 2:4]
    } else {
      FD   <- cbind(flow[INI:END, 2:3], flow=rep(blockin$hflow*86400, DMAX))
    }
    list <- cbind(list,FD)
    # PEC computation
    preoutput  <- transform(list, PEC=(list$DRFLoss+list$DLSEEPLoss)/list$flow*1000)
    #preoutput   <- transform(list, PEC=list$DRFLoss/list$flow*1000)
        output  <- preoutput[,c(1,6,7,8,2,3,4,5,9)]
    return(output)
  }) # end of with statement
} # end of program PCPFB
