p = read.table("/Users/shuyuzhang/OneDrive/000_Desktop_documents/Documents/Russia Meteorological Stations_Bulygina and Razuvaev/Russia_stations/wr225330/wr225330.txt", fill = T)
colnames(p) <- c("station","year","month","day","Tmin","Ta","Tmax","P")
## reorganize the dataframe 

syear <- aggregate(year ~ station, data = p, FUN = range)

syear$length = syear$year[,2]-syear$year[,1]
## check if there is any records have longer than 30days of NA
station = unique(p$station)
P= list()
for(i in 1:600){
  data = p[p$station == station[i],]
  P[[i]] = data
}
goodstation = array()
n = 1
for(i in 1:length(station)){
  x = p[p$station == station[i],]
  if(min(x$year) <= 1961 & max(x$year >= 2022)){
    x = x[x$year>=1961 & x$year <= 2022,]
    y = sum(is.na(x$P))/nrow(x)
    if( y <= 0.05 ){
      goodstation[n] = x$station[1]
      n = n+1
    }
  }
}
#### fill the na value
P= list()
for(i in 1:length(goodstation)){
  data = p[p$station == goodstation[i],]
  P[[i]] = data[data$year>=1961 & data$year <= 2022,]
}

for(i in 1:length(goodstation)){
  data = P[[i]]
  pmean = aggregate(P~month+day, data = data , FUN = mean)
  tmmean = aggregate(Tmin~month+day, data = data , FUN = mean)
  tmean = aggregate(Ta~month+day, data = data , FUN = mean)
  txmean = aggregate(Tmax~month+day, data = data , FUN = mean)
  
  Avg = data.frame(month = pmean$month,
                    day = pmean$day,
                    Pm = pmean$P, 
                    Tm = tmmean$Tmin,
                    T = tmean$Ta,
                    Tx = txmean$Tmax)
  
  x2 = merge(x = data, y = Avg, by = c("month","day"), all.x = TRUE)
  na.id = which(is.na(x2$P))
  x2[na.id, "P"] = x2[na.id, "Pm"] 
  
  na.id = which(is.na(x2$Ta))
  x2[na.id, "Ta"] = x2[na.id, "T"] 
  
  na.id = which(is.na(x2$Tmax))
  x2[na.id, "Tmax"] = x2[na.id, "Tx"] 
  
  na.id = which(is.na(x2$Tmin))
  x2[na.id, "Tmin"] = x2[na.id, "Tm"] 
  
  P[[i]] = x2
}

save(P, file = "/Users/shuyuzhang/OneDrive/Russia Meteorological Stations_Bulygina and Razuvaev/00-Russia STDF/408station_fillNA.RData")


setwd("/Users/shuyuzhang/OneDrive/Russia Meteorological Stations_Bulygina and Razuvaev/00-Russia STDF")
load("408station_fillNA.RData")

library(standaRdized)
library(xts)
dates <- seq(from = as.Date('1961-01-01'), to = as.Date('2022-12-31'), by = 1)
# for( i in 1:length(P)){
#   data = P[[i]]
#   data = data[with(data,order(year,month,day)),]
#   P[[i]] = data
# }

SPI = list()
for( i in 1:length(P)){
  data1 = P[[i]]$P
  SPI_1 = spi(data1, 30)
  SPI[[i]] = SPI_1$fitted
}
##3 remove 37. 230. 304 stations 
goodstation = goodstation[-c(37,230,304)]

p = P[-c(37,230,304)]

## fill other stations with missing data ### 
# 4   79  81  88 108 112 249 252 300 307 321 403
miss.id = c(4,79,81,88,108,112,249,252,300,307,321,403)

corr.date = P[[1]][,c("year","month","day")]

for(i in miss.id){
  data = P[[i]][,c(1:8)]
  pmean = aggregate(P~month+day, data = data , FUN = mean)
  tmmean = aggregate(Tmin~month+day, data = data , FUN = mean)
  tmean = aggregate(Ta~month+day, data = data , FUN = mean)
  txmean = aggregate(Tmax~month+day, data = data , FUN = mean)
  
  Avg = data.frame(month = pmean$month,
                   day = pmean$day,
                   Pm = pmean$P, 
                   Tm = tmmean$Tmin,
                   T = tmean$Ta,
                   Tx = txmean$Tmax)
  
  x2 = merge(x = data, y = corr.date, by = c("year","month","day"), all.y = TRUE)
  x2 = merge(x = x2, y = Avg, by = c("month","day"), all = TRUE)
  na.id = which(is.na(x2$P))
  x2[na.id, "P"] = x2[na.id, "Pm"] 
  
  na.id = which(is.na(x2$Ta))
  x2[na.id, "Ta"] = x2[na.id, "T"] 
  
  na.id = which(is.na(x2$Tmax))
  x2[na.id, "Tmax"] = x2[na.id, "Tx"] 
  
  na.id = which(is.na(x2$Tmin))
  x2[na.id, "Tmin"] = x2[na.id, "Tm"] 
  
  x2 = x2[with(x2,order(year,month,day)),]

  P[[i]] = x2
}

for( i in miss.id){
  data1 = P[[i]]$P
  SPI_1 = spi(data1, 30)
  SPI[[i]] = SPI_1$fitted
}

##### calculate the 3-days cumulative precipitation 
library(zoo)
rollP = list()
for(i in 1:408){
  data1 = P[[i]]$P
  rollP[[i]] = rollsum(data1, 3)
}

#### build the stdf data.frame ####
DATA = list()
for(i in 1:408){
  DATA[[i]] = data.frame(year =  P[[i]]$year, month = P[[i]]$month, 
                         day = P[[i]]$day, station = P[[i]]$station, 
                         pre = P[[i]]$P, spi = SPI[[i]], 
                         rollP = c(NA,NA,rollP[[i]]))
}

goodstation = unlist(lapply(P, FUN = function(x) {x$station[1]}))
### analysis the distribution of spi and p 
Pinfo = read.csv("/Users/shuyuzhang/OneDrive/Russia Meteorological Stations_Bulygina and Razuvaev/00-Russia STDF/WMO_stations.csv", header = T)
Pinfo = Pinfo[Pinfo$station %in% goodstation, ]

## the quantile value and QR of P at first
## the 90th and 10th of P for each station 
Pinfo$p90 = Pinfo$p10 = array()
DATA = list()
for ( i in 1:408){
  data = DATA[[i]]
  Pinfo$p90[i] = quantile(data$pre, prob = 0.9, na.rm = T)
  Pinfo$p10[i] = quantile(data$pre, prob = 0.1, na.rm = T)
}

library(raster)
library(ggplot2)


data("wrld_simpl", package = "maptools")  
wm <- crop(wrld_simpl, extent(-180, 180, 30, 90)) 
rus <- wrld_simpl[wrld_simpl$NAME == "Russia",]
p = ggplot()+
  geom_polygon(data = wm,aes(x = long, y = lat,group = group),
               fill = NA, colour = "grey50", size =0.15)+
  geom_polygon(data = rus,aes(x = long, y = lat,group = group),
               fill = "grey95", colour = "grey10", size =0.5)+
  coord_map("orthographic",orientation = c(45,100,0),
            xlim = c(-15,180),ylim = c(40,90))+
  theme_bw()+
  theme(
    panel.border=element_blank(),
    axis.text=element_blank(),
    axis.title=element_blank(),
    panel.grid.minor = element_line(size = 0.2,
                                    linetype = 'solid',
                                    colour = "black"))

p + geom_point(data = Pinfo, aes(x = longitude, y = latitude, 
                                 col = p90, size = p90, alpha = 0.2))+
  scale_colour_gradient2(low = "orange", mid = "cyan",
                         high = "blue", space = "Lab",
                         guide = "colourbar", aesthetics = "colour")+
  ggtitle("90th 3-days precipitation")

## filter the stdf events 
## larger than 90th quantile levels 
## first check the continuous drought, and within 10 days
STDF = drought = flood = STDF2= drought2 = STDF3= drought3 = list()

for ( i in 1:408){
  data = DATA[[i]]
  p90 = quantile(data$rollP, prob = 0.9, na.rm = T)
   
  flood[[i]] = data[data$rollP >= p90, c(1:4)] ## flooding events 
  
  sp = rle(data$spi <= -0.5)
  sqc = cumsum(sp$lengths)
  dr.id = sp$lengths >= 40 
  drought3[[i]] = data.frame(data[sqc[dr.id],c(1:4)], 
                            length = sp$lengths[sp$lengths >= 40]) ## drought events
  
  stdf.id = array()
  k = 1
  for( n in sqc[dr.id]){
    if(sum(data$rollP[n:(n+10)] >= p90, na.rm = T) > 0){
      stdf.id[k] = n ## also, the ending point of the drought events
      k = k+1
    }
  }
 STDF3[[i]] = data[stdf.id,c(1:4)]
}


### analysis the STDF evetns for each station ####
colfunc<-colorRampPalette(c("red","magenta","blue","cyan","green","yellow","orangered"))

Pinfo$total.stdf = unlist(lapply(STDF, nrow))
Pinfo$total.stdf2 = unlist(lapply(STDF2, nrow))
Pinfo$total.stdf3 = unlist(lapply(STDF3, nrow))

## total STDF events 
p + geom_point(data = Pinfo, aes(x = longitude, y = latitude, 
                                 col = total.stdf3/61),alpha = 0.8)+
  scale_colour_gradientn(colours=c("#fbeaff","#845EC2","#ff6f91","#ffc75f","#f9f871"),
                         space = "Lab")+
  ggtitle("annual mean STDF frequency (spi<-0.5 for 40days)")
## calculate the percentage of STDF/drought
pct.stdf2 = list()
year = data.frame(year = 1962:2022)
for( i in 1:408){
  stdf = aggregate(station~year, STDF2[[i]], length)
  dry =  aggregate(station~year, drought2[[i]], length)
  x = merge(stdf, dry, by = "year", all =T)
  x$pct = x$station.x/x$station.y
  x = merge(x, year, by.y = "year", all =T)
  pct.stdf2[[i]] = x
}

Pinfo$pct.stdf2 = unlist(lapply(pct.stdf2, FUN = function(x){
  mean(x$pct, na.rm  = T)}))


p + geom_point(data = Pinfo, aes(x = longitude, y = latitude, 
                                 col = pct.stdf2),alpha = 0.8)+
  scale_colour_gradientn(colours=c("#fbeaff","#845EC2","#ff6f91","#ffc75f","#f9f871"),
                         space = "Lab")+
  ggtitle("percentage of STDF over Drought (spi<-1 for 20days)")


### calculate the trend of stdf
Pinfo$sp.co2 = Pinfo$sp.p2 = array()
for (i in 1:408){
  x = STDF2[[i]]
  x.ts =  aggregate(station~year, x, length)
  a = cor.test(x.ts$station + rnorm(nrow(x.ts))/10000, x.ts$year, method = "spearman")
  Pinfo$sp.co2[i] = a$estimate
  Pinfo$sp.p2[i] = a$p.value
}

Pinfo$lm.co2 = Pinfo$lm.p2 = array()
for (i in 1:408){
  x = STDF2[[i]]
  x.ts =  aggregate(station~year, x, length)
  a = lm(station ~ year, data = x.ts)
  a = summary(a)
  Pinfo$lm.co2[i] = a$coefficients[2,1]
  Pinfo$lm.p2[i] = a$coefficients[2,4]
}

p + geom_point(data = Pinfo[Pinfo$sp.p2>0.1,], aes(x = longitude, y = latitude, 
                                                  col = sp.co2),size = 1.5, alpha = 0.5)+
  geom_point(data = Pinfo[Pinfo$sp.p2<= 0.1,], aes(x = longitude, y = latitude, 
                                 col = sp.co2),size = 4, alpha = 0.9)+
  scale_colour_gradient2(low = "royalblue1",mid = "white",high = "orangered",
                                   midpoint = 0,  space = "Lab")+
  ggtitle("spearman cor of STDF frequency (spi<= -1 for 20days)")

## calculate the trend of pct.stdf

Pinfo$lm.co = Pinfo$lm.p = array()
for (i in 1:408){
  x = STDF2[[i]]
  x.ts =  aggregate(station~year, x, length)
  x.ts = merge(x.ts, year, by = "year", all = T)
  a = lm(station ~ year, data = x)
  a = summary(a)
  Pinfo$lm.co[i] = a$coefficients[2,1]
  Pinfo$lm.p[i] = a$coefficients[2,4]
}


p + geom_point(data = Pinfo[Pinfo$lm.p<=0.1,], aes(x = longitude, y = latitude, 
                                                   col = lm.co*10),size = 4, alpha = 0.8)+
  geom_point(data = Pinfo[Pinfo$lm.p>0.1,], aes(x = longitude, y = latitude, 
                                                col = lm.co*10),size = 1, alpha = 0.5)+
  scale_colour_gradient2(low = "blue",mid = "white",high = "red",
                         midpoint = 0,  space = "Lab")+
  ggtitle("trend of STDF frequency")

## compare the frequency of stdf before and after 1990

Pinfo$total.2.early = unlist(lapply(STDF2, FUN = function(x) {
  nrow(x[x$year <= 1990,])}))
Pinfo$total.2.late = unlist(lapply(STDF2, FUN = function(x) {
  nrow(x[x$year > 1990,])}))

p + geom_point(data = Pinfo, aes(x = longitude, y = latitude, 
                                 col = total.2.early),size = 2, alpha = 0.8)+
  scale_colour_gradientn(colours=colfunc(12),
                         space = "Lab")+
  ggtitle("total number of STDF for 1961-1990")


p + geom_point(data = Pinfo, aes(x = longitude, y = latitude, 
                                 col = total.2.late),size = 2, alpha = 0.8)+
  scale_colour_gradientn(colours=colfunc(12),
                         space = "Lab")+
  ggtitle("total number of STDF for 1991-2022")


p + geom_point(data = Pinfo, 
               aes(x = longitude, y = latitude, 
                   col = total.2.late - total.2.early),
               size = 2)+
  scale_colour_gradientn(colours= c("orangered","orange","white","blue","navy"),
                         space = "Lab")+
  ggtitle("total number changes of STDF")




# meaningless to calculate the trend of individual stations, 
# separate them into clusters maybe meaningful 
# and also analysis the average timing of the STDF

#### change the dates into Julian days 
STDF2.avg.timing = list()
for(i in 1:408){
  data = STDF2[[i]]
  data$dates = as.Date(paste(data$year, data$month, data$day, sep = "-"))
  for(k in 1:nrow(data)){
    data$Julian[k] = as.numeric(julian(data$dates[k], 
                            origin = as.Date(paste(data$year[k], "1-1", sep = "-"))))
  }
  STDF2.avg.timing[[i]] = aggregate(Julian~year, data, mean)
}

Pinfo$avg.timing = unlist(lapply(STDF2.avg.timing, FUN= function(x){
  mean(x$Julian)
}))

Pinfo$avg.timing.early = unlist(lapply(STDF2.avg.timing, FUN= function(x){
  mean(x[x$year<=1990,"Julian"])
}))

Pinfo$avg.timing.late = unlist(lapply(STDF2.avg.timing, FUN= function(x){
  mean(x[x$year>1990,"Julian"])
}))

p + geom_point(data = Pinfo, aes(x = longitude, y = latitude, 
                                col = avg.timing.early),size = 2, alpha = 0.8)+
  scale_colour_gradientn(colours=colfunc(12),
                         space = "Lab")+
  ggtitle("average timing of STDF for 1961-1990")


p + geom_point(data = Pinfo, aes(x = longitude, y = latitude, 
                                 col = avg.timing.late),size = 2, alpha = 0.8)+
  scale_colour_gradientn(colours=colfunc(12),
                         space = "Lab")+
  ggtitle("average timing of STDF for 1991-2022")


p + geom_point(data = Pinfo, 
               aes(x = longitude, y = latitude, 
                  col = avg.timing.late - avg.timing.early),
               size = 2)+
  scale_colour_gradientn(colours= c("orangered","orange","white","blue","navy"),
                         space = "Lab")+
  ggtitle("average timing changes of STDF")


## cluster them into groups ####


som_grid <- somgrid(xdim = 3, ydim=3, topo="rectangular")
traj_som <- som(scale(som_mx),  grid=som_grid, 
                rlen=2000, 
                alpha=c(0.05,0.01),
                keep.data = TRUE )





