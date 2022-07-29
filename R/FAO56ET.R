FAO56ET <- function(pars, DMAX, meteo){
   # par includes; lat, z, SN, Kcini, Kcmid, Kcend, D1:D4
 with(as.list(c(pars,as.data.frame(meteo))), {
   
   # Initialization
   ap   <- 101.3 * (((293 - 0.0065 * z) / 293) ^ 5.26)
   gama <- 0.000655 * ap
   if (SN == "N") {
     Phi <- 3.14 * lat/180
   } else {
     Phi = -3.14 * lat/180
   }

   # Start loop
   for (i in 1:DMAX){
     Date  <-  as.Date(meteo[i ,"Date"])
     R     <-  meteo[i+1 ,"Rain"]
     Tmax  <-  meteo[i+1 ,"Tmax"]
     Tmin  <-  meteo[i+1 ,"Tmin"]
     Tave0 <-  meteo[i ,"Tave"]
     Tave  <-  meteo[i+1 ,"Tave"]
     U2    <-  meteo[i+1 ,"Wind"]
     Hum   <-  meteo[i+1 ,"Hum"]
     Solar <-  meteo[i+1 ,"Solar"]
     if (i <= D1) {
       Kcrop   <-  Kcini
     } else if (i <= D2) {
       Kcrop   <-  Kcini + (i - D1) * (Kcmid - Kcini) / D2
     } else if (i <= (D2+D3)) {
       Kcrop   <-  Kcmid
     } else {
       Kcrop   <-  Kcmid - (i - D2 - D3) * (Kcmid - Kcend) / D4
     }
     jj    <-  as.integer(275 * month(Date) / 9 - 30 + day(Date)) - 2
     if (month(Date) < 3) jj + 2
     if (month(Date) > 2) jj + 1
     Denta   <- 4098 * (0.6108 * exp(17.27 * Tave/(Tave + 237.3))) / ((Tave + 237.3) ^ 2)
     Estmax  <- 0.6108 * exp(17.27 * Tmax/(237.3 + Tmax))
     Estmin  <- 0.6108 * exp(17.27 * Tmin/(237.3 + Tmin))
     Es      <- 0.5 * (Estmax + Estmin)
     ea      <- 0.5 * (Estmax + Estmin) * Hum * 0.01
     Ge      <- 0.38 * (Tave - Tave0)
     sig     <- 0.409 * sin(0.0172 * jj - 1.39)
     dr      <- 1 + 0.033 * cos(0.0172 * jj)
     XX      <- 1 - tan(Phi) * tan(Phi) * tan(sig) * tan(sig)
     omes    <- 1.57 - atan((-tan(Phi) * tan(sig))/(XX ^ 0.5))
     Ra      <- 37.6 * dr * (omes * sin(Phi) * sin(sig) + cos(Phi) * cos(sig) * sin(omes))
     Rso     <- 0.75 * Ra
     Rnl     <- (2.45 * 10 ^ -9) * (((Tmax + 273.3) ^ 4) + ((Tmin + 273.3) ^ 4)) * ((1.35 * Solar / Rso) - 0.35) * (0.34 - 0.14 * (ea ^ 0.5))
     Rns     <- 0.77 * Solar
     Rn      <- Rns - Rnl
     Eto     <- (0.408 * Denta * (Rn - Ge) + (900 * gama * U2 * (Es - ea)/(273 + Tave)))/(Denta + gama * (1 + 0.34 * U2))
     Etc     <- Kcrop * Eto 
   
     if (i == 1) {
       out   <-  data.frame(Date=Date, Day=i, ET0=Eto, ETC=Etc)
     } else {
       outd  <-  data.frame(Date=Date, Day=i, ET0=Eto, ETC=Etc)
       rownames(outd)  <- i + 1
       out <- rbind(out, outd)       
     }
   }
 return(out)
 })
}