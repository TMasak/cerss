maxage <- 99
minage <- 60
maxyear <- 2013
data <- array(0,c(32,maxyear-1963,maxage-minage+1))
for(i in 1:32){
  print(i)
  country <- read.table(paste("death_rates/",i,".txt",sep=""),header=T,na.strings = ".")
  country <- country[country$Year >= 1964,]
  country <- country[country$Year <= maxyear,]
  country <- country[country$Age <= maxage,]
  country <- country[country$Age >= minage,]
  data[i,,] <- matrix(country$Total,ncol=maxage-minage+1,byrow=T) # rows = year , columns = age
}
save(data,file="mortality_data.RData")