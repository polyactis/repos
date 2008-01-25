############################################################################

getSEQ2TR<-function(outRoot,N,fileNamesMULT){
  
  setwd(outRoot)                       
  for(tr in 1:nrow(N)){
    if (file.exists(fileNamesMULT[tr])) {  
      hello1 <- system(sprintf("seq2tr.exe tree%d tree%d.tr",tr,tr),intern=TRUE)
    } else {
      return(FALSE)
    }
  }
  return(TRUE)
}

getTREEPIC<-function(outRoot,N,fileNamesTR){

  setwd(outRoot)                       
  for(tr in 1:nrow(N)){
    if (file.exists(fileNamesTR[tr])) {
     hello2 <- system(sprintf("treepic.exe tree%d.tr tree%d.ps",tr,tr),intern=TRUE)
    } else {
      return(FALSE)
    }
  }
  return(TRUE)
}

plotTree<-function(outRoot,tr, R, M){

  setwd(outRoot)
  myfile <- paste(outRoot,"pickedtr", sep = "")                                     
  sink(myfile)
  for(i in 1:nrow(R[[tr]])){
    cat(paste(M[[tr]][i],":"))
    cat(paste(" ",R[[tr]][i,]))
    cat("\n")
  }
  sink() 
  hello1 <- system("seq2tr.exe pickedtr pickedtr.tr",intern=TRUE)
  hello2 <- system("treepic.exe pickedtr.tr pickedtr.ps",intern=TRUE)
  unlink(paste(outRoot,"pickedtr", sep = ""))   
  unlink(paste(outRoot,"pickedtr.tr", sep = ""))   
                                       
}

