WBcalc <- function(pars, DMAX, meteo){
   # par includes; Hmax, Hmin, Hini, ddrain, WHP, dperc, dseep
 with(as.list(c(pars,as.data.frame(meteo))), {
   # Initialization
   out <- data.frame(day=0,rain=NA,irr=NA,drain=NA,perc=NA,seep=NA,et=NA,h=Hini)
   rownames(out)  <- 0
   # Start loop
   for (i in 1:DMAX){
     Hpw    <-  out[i   ,"h"]
     rain   <-  as.numeric(meteo[i ,"rain"])
     et     <-  as.numeric(meteo[i ,"et"])
     def    <-  rain - et - dperc - dseep
     Hpw1   <-  Hpw + def
     if (Hpw1 > Hmax) {
       drain   <-  Hpw1 - Hmax
       irr     <-  0
       Hpw2    <-  Hmax
     } else if (Hmin > Hpw1) {
       drain   <-  0
       irr     <-  Hmin-Hpw1
       Hpw2    <-  Hmin
     } else if (i <= WHP){
       drain   <-  0
       irr     <-  0
       Hpw2    <-  Hpw1
     } else {
       drain   <-  ddrain
       irr     <-  Hmax-Hpw1+drain
       Hpw2    <-  Hmax
     }
     outd <- data.frame(day=i,rain=rain,irr=irr,drain=drain,perc=dperc,seep=dseep,et=et,h=Hpw2)
     rownames(outd)  <- i
     out <- rbind(out, outd)       
   }
   
 return(out)
 })
}