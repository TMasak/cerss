get_last_no <- function(aa){
  ind <- m <- nchar(aa)
  for(j in rev(1:m)){
    if(substring(aa,j,j) != ' '){
      ind <- ind-1
    }else{
      x <- substring(aa,j+1,m)
      break
    }
  }
  return(x)
}
read_one_file <- function(X){
  Res <- array(0,c(64,256))
  i <- j <- 1
  pom <- FALSE
  for(lineno in 1:length(X$V1)){
    aa <- X$V1[lineno]
    if(substring(aa,1,1) != "#"){
      pom <- TRUE
      Res[i,j] <- as.double(get_last_no(aa))
      j <- j+1
    }else if(pom){
      i <- i+1
      j <- 1
    }
  }
  return(Res)
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
### unwrap all patients
fnames <- list.files("./data/tar")
for(j in 1:length(fnames)){
  untar(paste0("./data/tar/",fnames[j]), exdir="./data/untar")
}

### unwrap all trials for all patients
library(R.utils)
setwd("./data/untar")
fnames <- list.files()
for(j in 1:length(fnames)){
  fsubnames <- list.files(fnames[j])
  for(k in 1:length(fsubnames)){
    gunzip(paste0("./", fnames[j], "/", fsubnames[k]))
  }
}

### read data from the unwrapped files

setwd("./data/untar")
pnames <- list.files()
Data <- array(0,c(122,120,64,256)) # patients, trials, channels, time
Data2 <- array(0,c(122,120)) # patients, trials - indicating 0 for missing trial,
                             # 1 for S1 obj, 2 for S2 match, and 3 for S2 nomatch

for(i in 109:length(pnames)){
  tnames <- list.files(paste0("./",pnames[i]))
  for(j in 1:length(tnames)){
    print(c(i,j))
    X <- read.delim(paste0("./", pnames[i],"/", tnames[j]), header=F)
    tno <- as.double(substring(tnames[j],nchar(tnames[j])-2,nchar(tnames[j])))
    Data[i,tno+1,,] <- read_one_file(X)
    cname <- substring(X$V1[4], 3, 8)
    if(is.na(cname)) cname <- "bla"
    if(cname=="S1 obj"){
      Data2[i,tno+1] <- 1
    }else if(cname=="S2 mat"){
      Data2[i,tno+1] <- 2
    }else if(cname=="S2 nom"){
      Data2[i,tno+1] <- 3
    }
  }
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./data")
save(Data,file="surfaces.RData")
save(Data2,file="trial_info.RData")





