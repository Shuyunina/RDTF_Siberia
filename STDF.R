p = read.table("/", fill = T)
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

save(P, file = "")


setwd("")
load("")

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
rollP = list()
for(i in 1:408){
  data1 = P[[i]]$P
  n = length(data1)
  data2 = cbind(data1[1:(n-2)],data1[2:(n-1)],data1[3:n])
  rollP[[i]] = rowSums(data2)
}

#### build the stdf data.frame ####
DATA = list()
DATA[[i]] = data.frame(year =  P[[i]]$year, month = P[[i]]$month, 
                       day = P[[i]]$day, station = P[[i]]$station, 
                       pre = P[[i]]$P, spi = SPI[[i]], 
                       Ta = P[[i]]$Ta, Tmin = P[[i]]$Tmin,Tmax = P[[i]]$Tmax,
                       rollP = c(NA,NA,rollP[[i]]))


goodstation = unlist(lapply(P, FUN = function(x) {x$station[1]}))
### analysis the distribution of spi and p 
Pinfo = read.csv("", header = T)
Pinfo = Pinfo[Pinfo$station %in% goodstation, ]


Pinfo$p90 = Pinfo$p10 = array()

library(raster)
library(ggplot2)


data("wrld_simpl", package = "maptools")  
wm <- crop(wrld_simpl, extent(-180, 180, 30, 90)) 
rus <- wrld_simpl[wrld_simpl$NAME == "Russia",]
p = ggplot()+
  geom_polygon(data = wm,aes(x = long, y = lat,group = group),
               fill = NA, colour = "grey50", size =0.15)+
  geom_polygon(data = rus,aes(x = long, y = lat,group = group),
               fill = "grey95", colour = "grey20", size =0.25)+
  coord_map("orthographic",orientation = c(45,100,0),
            xlim = c(-15,180),ylim = c(40,90))+
  theme_bw()+
  theme(
    panel.border=element_blank(),
    axis.text=element_blank(),
    axis.title=element_blank(),
    panel.grid.major = element_line(size = 0.1,
                                    linetype = 'solid',
                                    colour = "black"))

p1 = ggplot()+
  geom_polygon(data = rus,aes(x = long, y = lat,group = group),
               fill = "grey95", colour = "grey10", size =0.25)+
  coord_map("orthographic",orientation = c(45,100,0),
            xlim = c(-15,180),ylim = c(40,90))+
  theme_bw()+
  theme(
    panel.border=element_blank(),
    axis.text=element_blank(),
    axis.title=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank())




drought = flood = STDF2 = list()

for ( i in 1:408){
  data = DATA[[i]]
  p90 = quantile(data$rollP, prob = 0.9, na.rm = T)
  
  flood[[i]] = data[data$rollP >= p90, c(1:4,7)] ## flooding events 
  
  sp = rle(data$spi <= -1)
  sp = data.frame(values = sp[[2]], lengths = sp[[1]])
  sp$values[is.na(sp$values)] = FALSE
  sp$dr.id = (sp$values == T) & (sp$lengths >=20)
  dr = sp$lengths[sp$dr.id]
  loc = cumsum(sp$lengths) [sp$dr.id] ## the end of a drought evnets
  
  drought2[[i]] = dry = data.frame(data[loc,], ## the end of a drought evnets
                                   duration = dr) ## drought events
  
  stdf.id = itv = rp = dr1 = dry.in = dry.T =  dry.Tmax =  dry.Tmin = array()
  k = 1
  num = c(1:11)
  for( n in 1:length(loc)){
    a = data$rollP[loc[n]:(loc[n]+10)] 
    m = a >= p90
    
    if(sum(m, na.rm = T) > 0){
      stdf.id[k] = loc[n] ## also, the ending point of the drought events
      dr1[k] = dr[n]
      itv[k] = min(num[m])
      rp[k] = data$rollP[loc[n]:(loc[n]+10)][m][1]
      dry.in[k] = sum(data$spi[loc[n]:(loc[n]-dr[n]+1)]+1, na.rm = T)
      
      dry.T[k] = mean(data$Ta[loc[n]:(loc[n]-dr[n]+1)], na.rm = T)
      dry.Tmax[k] = max(data$Ta[loc[n]:(loc[n]-dr[n]+1)], na.rm = T)
      dry.Tmin[k] = min(data$Ta[loc[n]:(loc[n]-dr[n]+1)], na.rm = T)
      
      k = k+1
    }
  }
  STDF2[[i]] =  data.frame(data[stdf.id, c(1:4,6)], dry.in,
                           dr1, dry.T, dry.Tmax, dry.Tmin,
                           interval = itv, fld = rp)
}
## 别删
flood = lapply(flood, function(x){
  x[3:nrow(x),]
})

### analysis the STDF events for each station ####
colfunc<-colorRampPalette(c("red","magenta","blue","cyan","green","yellow","orangered"))

Pinfo$total.stdf2 = unlist(lapply(STDF2, nrow)) ### THE BEST CRITERION OF STDF

## TOTAL STDF events  and interval days between the drought and the extreme P 
p + geom_point(data = Pinfo, aes(x = longitude, y = latitude, 
                                 col = total.stdf2, size = -itv),alpha = 0.7)+
  scale_colour_gradientn(colours=c("#fbeaff","#845EC2","#ff6f91","#ffc75f","#f9f871","green","cyan","blue"),
                         space = "Lab")+
  ggtitle("total number of STDF (spi<-1 for 20days)")

## total drought event 
Pinfo$total.drought2 = unlist(lapply(drought2, nrow)) ### total number of drought and the duration
Pinfo$stdf.dr = unlist(lapply(STDF2, function(x){
  mean(x$dr1)
}))

Pinfo$stdf.Ta = unlist(lapply(STDF2, function(x){
  mean(x$dry.T)
}))

Pinfo$stdf.fld = unlist(lapply(STDF2, function(x){
  mean(x$fld)
}))

Pinfo$stdf.Ta.dis = cut(Pinfo$stdf.Ta, breaks = c(-Inf, -20,-15,-10,-5,0,5,10,15,20,Inf), 
                        labels = c("<-20","-20","-15","-10","-5","5","10","15","20",">20"), 
                        right = FALSE)



p1 + geom_point(data = Pinfo, aes(x = longitude, y = latitude, 
                                  col = stdf.Ta.dis, size = stdf.dr),alpha = 0.8)+
  scale_colour_manual(values= c("#3F51B5","#5C6BC0","#7986CB","#9FA8DA","#E8EAF6",
                                "#FFF8E1","#FFE082","#FFA000","#FF6F00","#dc4730"))+
  ggtitle("mean T and duration of stdf (spi<-1 for 20days)")

## percentage of STDF/drought
pct.stdf2 = list()
year = data.frame(year = 1961:2022)
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
                                 col = pct.stdf2), size = 4,alpha = 0.8)+
  scale_colour_gradientn(colours=c("#fbeaff","#845EC2","#ff6f91","#ffc75f","#f9f871","green","cyan","blue"),
                         space = "Lab")+
  ggtitle("percentage of STDF over Drought (spi<-1 for 20days)")

## the mean extreme P volume
Pinfo$avg.3d.P = unlist(lapply(flood, FUN = function(x){
  mean(x$rollP)
})
)

## frequency per year
Pinfo$f.3d.P = unlist(lapply(flood, nrow))
p + geom_point(data = Pinfo, aes(x = longitude, y = latitude, 
                                 col = avg.3d.P, size = f.3d.P/62),alpha = 0.7)+
  scale_colour_gradient2(low = "white",high = "blue")+
  ggtitle("frequency and volume of 3-day P >= 90th")


large.mx = do.call(rbind.data.frame, STDF2)

large = merge(large.mx, Pinfo[,c(2:4)], by = "station")

large$lon.cl = as.factor(cut(large$longitude, 
                             breaks = c(-Inf, 60,120,Inf), 
                             labels = c("-60E","60-120E",
                                        "120E-"), 
                             right = FALSE))


ggplot(data = large, aes(x = dry.in/dr1, y = fld))+
  geom_point(aes(col = lon.cl, size = fld), alpha = 0.25)+
  scale_colour_manual(values= c("#479003", "#eb9911","#845EC2"))+
  geom_smooth(method = "lm", col = "orangered")+
  facet_grid(lon.cl ~ month)+
  ggtitle("scatter of intensity of Drought and Pre volume, with lon")+
  theme_bw()+
  theme(
    panel.border = element_rect(fill = NULL, linetype = 'solid', colour = "black"),
    panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "grey"),
    panel.grid.minor = element_blank())


subtotal = aggregate(duration~month, large, mean)
subtotal$freq = aggregate(duration~month, large, length)$duration /408

ggplot(subtotal)  +  
  geom_bar(aes(x=month, y=duration),stat="identity", fill="#479003",colour="#006000")+ 
  geom_line(aes(x=month, y=4*freq),stat="identity",color="red",size=1)+ 
  labs(title= "Duration vs frequency for months", 
       x="Month",y="Duration of Drougth")+ 
  scale_y_continuous(sec.axis=sec_axis(~.*0.25,name="Frequency")) 



#### the timing of stdf  
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
Pinfo$avg.timing.month = as.factor(cut(Pinfo$avg.timing, 
                                       breaks = c(-Inf, 60, 91, 121, 152, 182, 213, 244, 275), 
                                       labels = c("FEB", "MAM", "APR", "MAY","JUN","JUL","AUG","SEP"), 
                                       right = FALSE))
p1 + geom_point(data = Pinfo, 
                aes(x = longitude, y = latitude, col  = avg.timing.month, shape = as.factor(cluster.4)),
                alpha = 0.9, size = 2,  stroke = 0.25)+
  scale_colour_manual(values= c("#479003","#a4d13a","#eee296","#eb9911","#fc5e1d","#680e03"))+
  scale_shape_manual(values = c(15:18))+
  ggtitle("average timing of STDF")


  
########################################################################################
### -------------- TREND of the STDF drought and extreme P and the pct--------------####
########################################################################################



LM = Pinfo[, c("station", "latitude", "longitude")]

## Trend of stdf frequency
LM$stdf.co = LM$stdf.p = array()
for (i in 1:408){
  x = STDF2[[i]]
  x.ts =  aggregate(station~year, x, length)
  x.ts = merge(x.ts, year, by = "year", all = T)
  x.ts$station[is.na(x.ts$station)] = 0
  if(nrow(x.ts) >= 20){
    a = lm(station ~ year, data = x.ts)
    a = summary(a)
    LM$stdf.co[i] = a$coefficients[2,1]
    LM$stdf.p[i] = a$coefficients[2,4]
    
  }else{
    LM$stdf.co[i] = LM$stdf.p[i] = NA
  }
}

p + geom_point(data = LM[LM$stdf.p>0.1,], aes(x = longitude, y = latitude, 
                                              col = stdf.co),size = 1.5, alpha = 0.7)+
  geom_point(data = LM[LM$stdf.p<= 0.1,], aes(x = longitude, y = latitude, 
                                              col = stdf.co),size = 4, alpha = 0.95)+
  scale_colour_gradient2(low = "royalblue1",mid = "white",high = "orangered",
                         midpoint = 0,  space = "Lab")+
  ggtitle("lm of stdf frequency (spi<= -1 for 20days)")



## Trend of 
# drought, frequency and duration, 
LM$f.co = LM$f.p = LM$d.co = LM$d.p = array()
for (i in 1:408){
  x = drought2[[i]]
  x.ts =  aggregate(station~year, x, length)
  x.ts = merge(x.ts, year, by = "year", all = T)
  a = lm(station ~ year, data = x.ts)
  a = summary(a)
  LM$f.co[i] = a$coefficients[2,1]
  LM$f.p[i] = a$coefficients[2,4]
  
  x.ts =  aggregate(duration~year, x, mean)
  x.ts = merge(x.ts, year, by = "year", all = T)
  a = lm(duration ~ year, data = x.ts)
  a = summary(a)
  LM$d.co[i] = a$coefficients[2,1]
  LM$d.p[i] = a$coefficients[2,4]
}

p + geom_point(data = LM[LM$f.p>0.1,], aes(x = longitude, y = latitude, 
                                           col = f.co),size = 1.5, alpha = 0.7)+
  geom_point(data = LM[LM$f.p<= 0.1,], aes(x = longitude, y = latitude, 
                                           col = f.co),size = 4, alpha = 0.95)+
  scale_colour_gradient2(low = "royalblue1",mid = "white",high = "orangered",
                         midpoint = 0,  space = "Lab")+
  ggtitle("lm of drought frequency (spi<= -1 for 20days)")


p + geom_point(data = LM[LM$d.p>0.1,], aes(x = longitude, y = latitude, 
                                           col = d.co),size = 1.5, alpha = 0.7)+
  geom_point(data = LM[LM$d.p<= 0.1,], aes(x = longitude, y = latitude, 
                                           col = d.co),size = 4, alpha = 0.95)+
  scale_colour_gradient2(low = "royalblue1",mid = "white",high = "orangered",
                         midpoint = 0,  space = "Lab")+
  ggtitle("lm of drought duration (spi<= -1 for 20days)")

# P frequency and volume, intervals
LM$f.flood.co = LM$f.flood.p = LM$p.flood.co = LM$p.flood.p = array()
for (i in 1:408){
  x = flood[[i]]
  x.ts =  aggregate(station~year, x, length)
  a = lm(station ~ year, data = x.ts)
  a = summary(a)
  LM$f.flood.co[i] = a$coefficients[2,1]
  LM$f.flood.p[i] = a$coefficients[2,4]
  
  x.ts =  aggregate(rollP~year, x, mean)
  a = lm(rollP ~ year, data = x.ts)
  a = summary(a)
  LM$p.flood.co[i] = a$coefficients[2,1]
  LM$p.flood.p[i] = a$coefficients[2,4]
}

p + geom_point(data = LM[LM$f.flood.p>0.1,], aes(x = longitude, y = latitude, 
                                                 col = f.flood.co),size = 1.5, alpha = 0.7)+
  geom_point(data = LM[LM$f.flood.p<= 0.1,], aes(x = longitude, y = latitude, 
                                                 col = f.flood.co),size = 4, alpha = 0.95)+
  scale_colour_gradient2(low = "royalblue1",mid = "white",high = "orangered",
                         midpoint = 0,  space = "Lab")+
  ggtitle("lm of extreme P frequency (spi<= -1 for 20days)")


p + geom_point(data = LM[LM$p.flood.p>0.1,], aes(x = longitude, y = latitude, 
                                                 col = p.flood.co),size = 1.5, alpha = 0.7)+
  geom_point(data = LM[LM$p.flood.p<= 0.1,], aes(x = longitude, y = latitude, 
                                                 col = p.flood.co),size = 4, alpha = 0.95)+
  scale_colour_gradient2(low = "royalblue1",mid = "white",high = "orangered",
                         midpoint = 0,  space = "Lab")+
  ggtitle("lm of  extreme P volume(spi<= -1 for 20days)")



qr.p = list()
library(stats)
for ( i in 1:408){
  data = flood[[i]]
  ecdf_p = ecdf(data$rollP)
  
  x = STDF2[[i]]
  x.ts = aggregate(fld ~ year, x, mean)
  x.ts$pos = ecdf_p(x.ts$fld)
  qr.p[[i]] = x.ts
}

plot(unlist(lapply(qr.p, function(x) median(x$pos))))

## meaningless to draw the linear regression of quantile levels


# timing 
LM$timing.co = LM$timing.p = array()
for (i in 1:408){
  x = STDF2.avg.timing[[i]]
  if(nrow(x) >= 20){
    a = lm(Julian ~ year, data = x)
    a = summary(a)
    LM$timing.co[i] = a$coefficients[2,1]
    LM$timing.p[i] = a$coefficients[2,4]
  }else{
    LM$timing.co[i] = LM$timing.p[i] = NA
  }
}

p + geom_point(data = LM[LM$timing.p>0.1,], aes(x = longitude, y = latitude, 
                                                col = timing.co),size = 1.5, alpha = 0.7)+
  geom_point(data = LM[LM$timing.p<= 0.1,], aes(x = longitude, y = latitude, 
                                                col = timing.co),size = 4, alpha = 0.95)+
  scale_colour_gradient2(low = "royalblue1",mid = "white",high = "orangered",
                         midpoint = 0,  space = "Lab")+
  ggtitle("lm of stdf timing (spi<= -1 for 20days)")



drought2.avg.timing =  list()
for(i in 1:408){
  data = drought2[[i]]
  data$dates = as.Date(paste(data$year, data$month, data$day, sep = "-"))
  for(k in 1:nrow(data)){
    data$Julian[k] = as.numeric(julian(data$dates[k], 
                                       origin = as.Date(paste(data$year[k], "1-1", sep = "-"))))
  }
  drought2.avg.timing[[i]] = aggregate(Julian~year, data, mean)
}

LM$dry.timing = unlist(lapply(drought2.avg.timing, FUN= function(x){
  mean(x$Julian)
}))
LM$avg.timing.month = as.factor(cut(LM$dry.timing, 
                                    breaks = c(-Inf, 60, 91, 121, 152, 182, 213, 244, 275), 
                                    labels = c("FEB", "MAM", "APR", "MAY","JUN","JUL","AUG","SEP"), 
                                    right = FALSE))
Pinfo$dry.timing = LM$dry.timing
p + geom_point(data = Pinfo, 
               aes(x = longitude, y = latitude, 
                   col = avg.timing - dry.timing),size = 3, alpha = 0.9)+
  scale_colour_gradient2(low = "royalblue1",mid = "white",high = "orangered",
                         midpoint = 0,  space = "Lab")+
  ggtitle("lag time of STDF to drought")


LM$timing.co = LM$timing.p = array()
for (i in 1:408){
  x = STDF2.avg.timing[[i]]
  if(nrow(x) >= 20){
    a = lm(Julian ~ year, data = x)
    a = summary(a)
    LM$timing.co[i] = a$coefficients[2,1]
    LM$timing.p[i] = a$coefficients[2,4]
  }else{
    LM$timing.co[i] = LM$timing.p[i] = NA
  }
}

p + geom_point(data = LM[LM$timing.p>0.1,], aes(x = longitude, y = latitude, 
                                                col = timing.co),size = 1.5, alpha = 0.7)+
  geom_point(data = LM[LM$timing.p<= 0.1,], aes(x = longitude, y = latitude, 
                                                col = timing.co),size = 4, alpha = 0.95)+
  scale_colour_gradient2(low = "royalblue1",mid = "white",high = "orangered",
                         midpoint = 0,  space = "Lab")+
  ggtitle("lm of stdf timing (spi<= -1 for 20days)")
########################################################################################
### ----------------------   BEFORE & AFTER      ---------------------------------####
########################################################################################


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








########################################################################################
### -----------------------------   CLUSTERING     ---------------------------------####
########################################################################################

#### cluster them into groups ####
library(kohonen)
library(aweSOM)

data =   scale(Pinfo[,c("latitude","longitude","total.stdf2","avg.timing")])
k <- list()
for ( i in 1:15){
  som_grid <- somgrid(xdim = 1, ydim=i, topo="rectangular")
  som.results <- som(data,  grid=som_grid, 
                     rlen=1000, 
                     alpha=c(0.05,0.01),
                     keep.data = TRUE )
  k[[i]] <- somQuality(som.results, data)
  print(i)
}

qe = unlist(lapply(k, function(x){
  return(x$err.quant)
}))

plot(qe)


som_grid <- somgrid(xdim = 1, ydim=4, topo="rectangular")
som.results <- som(data,  grid=som_grid, 
                   rlen=1000, 
                   alpha=c(0.05,0.01),
                   keep.data = TRUE )

Pinfo$cluster.4 = som.results$unit.classif

p + geom_point(data = Pinfo, aes(x = longitude, y = latitude,
                                 col = as.factor(cluster.4))) +
  # geom_label(data = Pinfo,aes(x = longitude, y = latitude,label = station))+
  scale_color_discrete()+
  ggtitle("4 clusters of stations")


### try which cluster would have better results 
ggplot()+
  geom_point(data = Pinfo, aes(x = longitude, y = latitude,
                               col = total.stdf2)) +
  # geom_label(data = Pinfo,aes(x = longitude, y = latitude,label = station))+
  scale_colour_gradientn(colours= c("orangered","orange","white","blue","navy"),
                         space = "Lab")+
  facet_grid(. ~ cluster.4)+
  ggtitle("4 clusters of stations")

#### re-organize
Pinfo$cluster.4[Pinfo$cluster.4 == 1 & Pinfo$longitude<50 ] = 4
Pinfo$cluster.4[Pinfo$cluster.4 == 1 & Pinfo$longitude<80 ] = 3

write.csv(Pinfo, file = "Pinfo with 4cluster.csv")


p + geom_point(data = Pinfo, aes(x = longitude, y = latitude, 
                                 col  = total.stdf2, shape = as.factor(cluster.4)),
               alpha = 0.9, size = 2)+
  scale_color_gradientn(colours=c("#FF5722","#FF9800","#ffc75f","#FFEB3B","#CDDC39","#4CAF50","#2196F3","#3F51B5","#9C27B0"),
                        space = "Lab")+
  scale_shape_manual(values = c(15:18)) + 
  ggtitle("total number of STDF (spi<-1 for 20days)")






### classified to analysis the temporal changes
large = merge(large, Pinfo[,c("station","cluster.4")], by = "station")




ggplot(data = large, aes(x = dr1, y = fld))+
  geom_point(aes(col = as.factor(cluster.4), size = dry.T), alpha = 0.25)+
  scale_colour_manual(values= c("#479003", "#eb9911","#845EC2","#9bddf9"))+
  geom_smooth(method = "lm", col = "orangered")+
  facet_grid(cluster.4 ~ month)+
  ggtitle("Duration-Drought VS Pre, with cluster")+
  theme_bw()+
  theme(
    panel.border = element_rect(fill = NULL, linetype = 'solid', colour = "black"),
    panel.grid.major = element_line(size = 0.2, linetype = 'solid', colour = "grey"),
    panel.grid.minor = element_blank())


subtotal <- aggregate(station ~ month+cluster.4, large, function(x){
  length(x)/length(unique((x)))})
subtotal$dr1 <- aggregate(dr1 ~ month+cluster.4, large, mean)$dr1
subtotal$dry.in <- aggregate(dry.in ~ month+cluster.4, large, mean)$dry.in

library(dplyr)
library(ggplot2)
library(stringr)




ggplot(subtotal) +
  # geom_hline(aes(yintercept = station), color = "grey40", size = 0.25) + 
  geom_col(aes(x = month, y = station, fill = as.factor(month)),
           position = "dodge2",
           show.legend = TRUE,
           alpha = .9) +
  scale_fill_manual(values =  c("#1976D2","#B2EBF2", "#AED581","#81C784","#CDDC39","#FDD835",
                                "#FFA000","#F4511E", "#BF360C","#D81B60","#AB47BC","#1565C0"))+
  labs(title= "frequency of STDF per month") + 
  facet_grid(.~ cluster.4)+
  # geom_vline(aes(xintercept = month), color = "lightgrey") + 
  coord_polar()+
  theme_bw()+
  theme(
    panel.border=element_blank(),
    # axis.text=element_blank(),
    axis.title=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.major.y = element_line(size = 0.2,
                                      linetype = 'dashed',
                                      colour = "black"))


ggplot(subtotal)  +  
  geom_bar(aes(x=month, y= dr1, 
               fill = as.factor(cluster.4)), alpha = 0.5,
           stat="identity")+ 
  xlim(1,12)+
  scale_fill_manual(values= c("#479003", "#eb9911","#845EC2","#9bddf9"))+
  geom_line(aes(x=month, y= -dry.in),stat="identity",color="#EF9A9A",size= 0.75)+ 
  geom_point(aes(x=month, y= -dry.in), col = "#C62828")+
  facet_grid(.~ cluster.4)+
  labs(title= "frequency, duration & severity for months", 
       x="Month",y="Duration") + 
  scale_y_continuous(sec.axis=sec_axis(~.*1.5,name="Severity")) +
  scale_x_continuous(breaks = c(1:12)) +
  theme_bw()+
  theme(
    panel.border=element_blank(),
    # axis.text=element_blank(),
    # axis.title=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.major.y = element_line(size = 0.2,
                                      linetype = 'dashed',
                                      colour = "black"))





# plot the time series #
ts.stdf = ts.drought = ts.pct = list()
year = data.frame(year = 1961:2022)
som.id = Pinfo$cluster.4

for(n in 1:4){
  ts.stdf[[n]] = ts.drought[[n]] = ts.pct[[n]] =year
  
  for( i in which(som.id == n)){
    pre = aggregate(station~year, STDF2[[i]], length)
    dry =  aggregate(station~year, drought2[[i]], length)
    x1 = merge(stdf, year, by = "year", all =T)
    x2 = merge(dry, year, by.y = "year", all =T)
    x3 = data.frame(year = x1$year, station = x1$station/x2$station)
    
    x1$station[is.na(x1$station)] = 0
    x2$station[is.na(x2$station)] = 0
    x3$station[is.na(x3$station)] = 0
    
    ts.pre[[n]] = cbind(ts.pre[[n]],x1$station)
    ts.drought[[n]] = cbind(ts.drought[[n]],x2$station)
    ts.pct[[n]] = cbind(ts.pct[[n]],x3$station)
  }
}
