# Sub-function for pre-processing of raw data
# handling the concentration data below LOQ

PCONC <- function(data_conc){
  with(as.list(data_conc), {
    for (i in 1:(nrow(data_conc))) {
      if (data_conc[i,2]=="="){
        Timed = data_conc[i,1]
        Cpwd  = data_conc[i,3]
      } else if(data_conc[i,2]=="<") {
        if (data_conc[i-1,2]=="=") {
          Timed = data_conc[i,1]
          Cpwd  = 0.5*data_conc[i,3]
        } else {
          for (j in i:(nrow(data_conc))) {
            if (data_conc[j,2]=="=") {
              Timed = data_conc[i,1]
              Cpwd  = 0.5*data_conc[i,3]
              break
            }
          }
          break
        }
      } else {
        stop("Qualifier is incorrect")
      }
     if (i==1) {
       Time = Timed
       Cpw  = as.numeric(Cpwd)
     } else {
       Time = cbind(Time, Timed)
       Cpw  = cbind(Cpw, Cpwd)
     }
    }
    out = data.frame(Time=as.numeric(Time), Cpw=as.numeric(Cpw))
    return(out)
  })
}