### analyze the SST, wind field, divergence, omega, GPH500/850, 
# latent/sensible net flux(shtfl)
## 4-WEEKS before and 2 week after the turning point
library(abind)
library(ncdf4)
library(raster)
load("~")
data("wrld_simpl", package = "maptools")  
wm <- crop(wrld_simpl, extent(0, 180, 0, 90)) 
rus <- wrld_simpl[wrld_simpl$NAME == "Russia",]


setwd("")

rotate <- function(x) t(apply(x, 2, rev))

var = c("hgt", "air.sig","pr_wtr", "pr", 
        "omega","lhtfl","shtfl","div",
        "uwnd","vwnd")## wind [78:152]

### extract anomalies ####

## SST-air.gis
file = list.files(pattern = "air.sig")
file = file[-c(1:13)]
file = file[-c(63,64)]

ncin = nc_open(file[1])
div = ncvar_get(ncin, "div")[1:72,1:36,] # 0:180 E, 0:90N
nc_close(ncin)

for(i in 2:length(file)){
  ncin = nc_open(file[i])
  a =  ncvar_get(ncin, "air")[1:72,1:36,]
  nc_close(ncin)
  sst = abind(sst, a, along= 3)
}



### extract those days 4 weeks before and 2 weeks after the STDF
### extract the days of STDF
days = c(rep(c(365,365,365,366),15),365,365)
days = c(0,cumsum(days))
days = days[-63]
year = 1961:2022

info1 = info[info$seg %in% c(1,2),]
info2 = info[info$seg %in% c(3,4),]
dm = c(72,36,1)

Z1 = Z2 =  list()
for ( i in 1:6 ){
  Z1[[i]] = Z2[[i]] =  array(NA, dim = dm)
}
rg = list()

### dealing with SST at first
z = sst
for(i in 1:62){
  y0 = year[i]
  id1 = which(info1$year == y0)
  id2 = which(info2$year == y0)
  
  jd1 = unique(info1$julian[id1])
  jd2 = unique(info2$julian[id2])
  for( k in jd1){
    rg[[1]] = (k+7):(k+13)
    rg[[2]] = k:(k+6)
    rg[[3]] = (k-6):k
    rg[[4]] = (k-13):(k-7)
    rg[[5]] = (k-20):(k-14)
    rg[[6]] = (k-27):(k-21)
    
    for ( w in 1:6){
      a = rg[[w]]
      z.1 = apply(z[,, (a+ days[i])],c(1,2), mean) 
      rg.1 = a+days[1]
      for ( n in 1:62){
        rg.1 = c(rg.1, a+days[n])
      }
      rg.1 = rg.1[rg.1>0]
      rg.1 = rg.1[rg.1<= dim(z)[3]]
      
      z.clim = apply(z[,, rg.1],c(1,2), mean) 
      z.1 = z.1-z.clim
      dim(z.1) <- dm
      Z1[[w]] = abind(Z1[[w]], z.1, along = 3)
    }
    }
  
  for( k in jd2){
    rg[[1]] = (k+7):(k+13)
    rg[[2]] = k:(k+6)
    rg[[3]] = (k-6):k
    rg[[4]] = (k-13):(k-7)
    rg[[5]] = (k-20):(k-14)
    rg[[6]] = (k-27):(k-21)
    
    for ( w in 1:6){
      a = rg[[w]]
      z.1 = apply(z[,, (a+ days[i])],c(1,2), mean) 
      rg.1 = a+days[1]
      
      for ( n in 1:62){
        rg.1 = c(rg.1, a+days[n])
      }
      rg.1 = rg.1[rg.1>0]
      rg.1 = rg.1[rg.1<= dim(z)[3]]
      
      z.clim = apply(z[,, rg.1],c(1,2), mean) 
      z.1 = z.1-z.clim
      dim(z.1) <- dm
      Z2[[w]] = abind(Z2[[w]], z.1, along = 3)  
      }    
  }
  print(i)
}


SST1 = Z1
SST2 = Z2
setwd("")
save(SST1, SST2, file = "SST.RData")




## dealing with multiple levels ####
file = list.files(pattern = "uwnd")
file = file[-c(1:77)]
file = file[14:75]

ncin = nc_open(file[1])
uwnd.850 = ncvar_get(ncin, "uwnd")[1:72,1:24,3, ] # 850(3), 500(6), 200mb(10) GPH at N.H 30-90N
uwnd.500 = ncvar_get(ncin, "uwnd")[1:72,1:24,6, ] # 850(3), 500(6), 200mb(10) GPH at N.H 30-90N
uwnd.200 = ncvar_get(ncin, "uwnd")[1:72,1:24,10, ] # 850(3), 500(6), 200mb(10) GPH at N.H 30-90N
nc_close(ncin)

for(i in 2:length(file)){
  ncin = nc_open(file[i])
  a = ncvar_get(ncin, "uwnd")
  nc_close(ncin)
  uwnd.850 = abind(uwnd.850, a[1:72,1:24,3, ], along= 3)
  uwnd.500 = abind(uwnd.500, a[1:72,1:24,6, ], along= 3)
  uwnd.200 = abind(uwnd.200, a[1:72,1:24,10, ], along= 3)
}


### extract those days 4 weeks before and 2 weeks after the STDF
### extract the days of STDF
days = c(rep(c(365,365,365,366),15),365,365)
days = c(0,cumsum(days))
days = days[-63]
year = 1961:2022

info1 = info[info$seg %in% c(1,2),]
info2 = info[info$seg %in% c(3,4),]
dm = c(72,24,1)

Z1.1 = Z2.1 =  Z1.2 = Z2.2 =  Z1.3 = Z2.3 =  list()
for ( i in 1:6 ){
  Z1.1[[i]] = Z2.1[[i]] = Z1.2[[i]] = Z2.2[[i]] = Z1.3[[i]] = Z2.3[[i]] =  array(NA, dim = dm)
}
rg = list()

### dealing with 
z = uwnd.850
z1 = uwnd.500
z2 = uwnd.200

for(i in 1:62){
  y0 = year[i]
  id1 = which(info1$year == y0)
  id2 = which(info2$year == y0)
  
  jd1 = unique(info1$julian[id1])
  jd2 = unique(info2$julian[id2])
  for( k in jd1){
    rg[[1]] = (k+7):(k+13)
    rg[[2]] = k:(k+6)
    rg[[3]] = (k-6):k
    rg[[4]] = (k-13):(k-7)
    rg[[5]] = (k-20):(k-14)
    rg[[6]] = (k-27):(k-21)
    
    for ( w in 1:6){
      a = rg[[w]]
      rg.1 = a+days[1]
      for ( n in 1:62){
        rg.1 = c(rg.1, a+days[n])
      }
      rg.1 = rg.1[rg.1>0]
      rg.1 = rg.1[rg.1<= dim(z)[3]]
      
      z.1 = apply(z[,, (a+ days[i])],c(1,2), mean) 
      z.clim = apply(z[,, rg.1],c(1,2), mean) 
      z.1 = z.1-z.clim
      dim(z.1) <- dm
      Z1.1[[w]] = abind(Z1.1[[w]], z.1, along = 3)
      
      
      z.1 = apply(z1[,, (a+ days[i])],c(1,2), mean) 
      z.clim = apply(z1[,, rg.1],c(1,2), mean) 
      z.1 = z.1-z.clim
      dim(z.1) <- dm
      Z1.2[[w]] = abind(Z1.2[[w]], z.1, along = 3)
      
      ## 
      z.1 = apply(z1[,, (a+ days[i])],c(1,2), mean) 
      z.clim = apply(z1[,, rg.1],c(1,2), mean) 
      z.1 = z.1-z.clim
      dim(z.1) <- dm
      Z1.2[[w]] = abind(Z1.2[[w]], z.1, along = 3)
    }
  }
  
  for( k in jd2){
    rg[[1]] = (k+7):(k+13)
    rg[[2]] = k:(k+6)
    rg[[3]] = (k-6):k
    rg[[4]] = (k-13):(k-7)
    rg[[5]] = (k-20):(k-14)
    rg[[6]] = (k-27):(k-21)
    
    for ( w in 1:6){
      a = rg[[w]]
      rg.1 = a+days[1]
      
      for ( n in 1:62){
        rg.1 = c(rg.1, a+days[n])
      }
      rg.1 = rg.1[rg.1>0]
      rg.1 = rg.1[rg.1<= dim(z)[3]]
      
      z.1 = apply(z[,, (a+ days[i])],c(1,2), mean) 
      z.clim = apply(z[,, rg.1],c(1,2), mean) 
      z.1 = z.1-z.clim
      dim(z.1) <- dm
      Z2.1[[w]] = abind(Z2.1[[w]], z.1, along = 3)  
      
      z.1 = apply(z1[,, (a+ days[i])],c(1,2), mean) 
      z.clim = apply(z1[,, rg.1],c(1,2), mean) 
      z.1 = z.1-z.clim
      dim(z.1) <- dm
      Z2.2[[w]] = abind(Z2.2[[w]], z.1, along = 3)  
      
      z.1 = apply(z2[,, (a+ days[i])],c(1,2), mean) 
      z.clim = apply(z2[,, rg.1],c(1,2), mean) 
      z.1 = z.1-z.clim
      dim(z.1) <- dm
      Z2.3[[w]] = abind(Z2.3[[w]], z.1, along = 3)  
    }    
  }
  print(i)
}


U1 = list(U850 = Z1.1, U500 = Z1.2, U200 = Z1.3)
U2 = list(U850 = Z2.1, U500 = Z2.2, U200 = Z2.3)

save(U1, U2, file = "UWND.RData")



### separate into each months and calculate the anomalies for each month ####
mycrs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
newcrs <- "+proj=laea +lat_0=90 +lon_0=100 +x_0=0 +y_0=0 +a=6371228 +b=6371228 +units=m +no_defs "
library(raster)


## plot sst ####
load("SST.RData")

sst1 = lapply(SST1, FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

sst2 = lapply(SST2, FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

mycol=colorRampPalette(c("#0D47A1","#64B5F6","#E1F5FE","#FFFDE7", "#FF9800","#BF360C")) 

par(mfrow = c(2,3),mar= c(.5,.5,.5,.5),  bty = 'n')
for ( i in 1:6){
  x = sst1[[i]]
  x = rotate(x) 
  x = t(apply(x, 1, rev))
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 0, ymx = 90, crs = mycrs)
  plot (x, main = paste("group 1/2 week-",i), breaks = seq(-0.6,0.6,0.1),
        col = mycol(12))
  plot(wm, lwd = 0.25, add = T, col = "white")
  plot(rus, lwd = 1, add = T)
}

for ( i in 1:6){
  x = sst2[[i]]
  x = rotate(x) 
  x = t(apply(x, 1, rev))
  
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 0, ymx = 90, crs = mycrs)
  plot (x, main = paste("group 3/4 week-",i), breaks = seq(-0.9,0.9,0.1),
        col = mycol(18))
  plot(wm, lwd = 0.25, add = T, col = "white")
  plot(rus, lwd = 1, add = T)
  
}

## plot Pr_wtr ####

load("Pr_wtr.RData")

pw1 = lapply(Pr_wtr1, FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

pw2 = lapply(Pr_wtr2, FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

mycol1=colorRampPalette(c("#E65100","#FBC02D","#FFF59D","#C8E6C9", "#4CAF50","#2E7D32")) 

par(mfrow = c(2,3),mar= c(.5,.5,.5,.5),  bty = 'n')
for ( i in 1:6){
  x = pw1[[i]]
  x = rotate(x) 
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 30, ymx = 90, crs = mycrs)
  plot (x, main = paste("Pr_wtr group 1/2 week-",i), breaks = seq(-0.5,0.5,0.1),
        col = mycol1(10))
  plot(wm, lwd = 0.25, add = T, col = NA)
  plot(rus, lwd = 1, add = T, col = NA)
  
  contour(x, col = "black", lwd = 1, drawlabels = T,
          add = T,levels = seq(0,0.5,0.1),
          line = list(smoothing = 2))
  contour(x, col = "black", lwd = 1, drawlabels = T,
          add = T,levels = seq(-0.6,0,0.1),lty = "dotted",
          line = list(smoothing = 2))
}

for ( i in 1:6){
  x = pw2[[i]]
  x = rotate(x) 
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 30, ymx = 90, crs = mycrs)
  plot (x, main = paste("group 3/4 week-",i), breaks = seq(-0.6,0.6,0.1),
        col = mycol1(12))
  plot(wm, lwd = 0.25, add = T,col = NA)
  plot(rus, lwd = 1, add = T, col = NA)
  
}

#### sensbible/latent heat flex and latent heat flux ####
load("lhtfl.RData")

lh1 = lapply(Lhtfl1, FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

lh2 = lapply(Lhtfl2, FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})


mycol=colorRampPalette(c("#1976D2","#42A5F5","#E3F2FD","#F3E5F5", "#7E57C2","#4527A0")) 

par(mfrow = c(2,3),mar= c(.5,.5,.5,.5),  bty = 'n')
for ( i in 1:6){
  x = lh1[[i]]
  x = rotate(x) 
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 30, ymx = 90, crs = mycrs)
  plot (x, main = paste("Lensible Heat group 1/2 week-",i), breaks = seq(-8,8,1),
        col = mycol(16))
  plot(wm, lwd = 0.25, add = T, col = NA)
  plot(rus, lwd = 1, add = T, col = NA)
}

for ( i in 1:6){
  x = lh2[[i]]
  x = rotate(x) 
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 30, ymx = 90, crs = mycrs)
  plot (x, main = paste("Lensible Heat group 3/4 week-",i), breaks = seq(-8,8,1),
        col = mycol(16))
  plot(wm, lwd = 0.25, add = T,col = NA)
  plot(rus, lwd = 1, add = T, col = NA)
  
}



#### 
load("lhtfl.RData")

lh1 = lapply(Lhtfl1, FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

lh2 = lapply(Lhtfl2, FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})


mycol=colorRampPalette(c("#1976D2","#42A5F5","#E3F2FD","#F3E5F5", "#7E57C2","#4527A0")) 

par(mfrow = c(2,3),mar= c(.5,.5,.5,.5),  bty = 'n')
for ( i in 1:6){
  x = lh1[[i]]
  x = rotate(x) 
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 30, ymx = 90, crs = mycrs)
  plot (x, main = paste("Lensible Heat group 1/2 week-",i), breaks = seq(-8,8,1),
        col = mycol(16))
  plot(wm, lwd = 0.25, add = T, col = NA)
  plot(rus, lwd = 1, add = T, col = NA)
}

for ( i in 1:6){
  x = lh2[[i]]
  x = rotate(x) 
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 30, ymx = 90, crs = mycrs)
  plot (x, main = paste("Lensible Heat group 3/4 week-",i), breaks = seq(-8,8,1),
        col = mycol(16))
  plot(wm, lwd = 0.25, add = T,col = NA)
  plot(rus, lwd = 1, add = T, col = NA)
  
}


#### sensbible/latent heat flex and latent heat flux ####
load("lhtfl.RData")

lh1 = lapply(Lhtfl1, FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

lh2 = lapply(Lhtfl2, FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})


mycol=colorRampPalette(c("#1976D2","#42A5F5","#E3F2FD","#F3E5F5", "#7E57C2","#4527A0")) 

par(mfrow = c(2,3),mar= c(.5,.5,.5,.5),  bty = 'n')
for ( i in 1:6){
  x = lh1[[i]]
  x = rotate(x) 
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 30, ymx = 90, crs = mycrs)
  plot (x, main = paste("Lensible Heat group 1/2 week-",i), breaks = seq(-8,8,1),
        col = mycol(16))
  plot(wm, lwd = 0.25, add = T, col = NA)
  plot(rus, lwd = 1, add = T, col = NA)
}

for ( i in 1:6){
  x = lh2[[i]]
  x = rotate(x) 
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 30, ymx = 90, crs = mycrs)
  plot (x, main = paste("Lensible Heat group 3/4 week-",i), breaks = seq(-8,8,1),
        col = mycol(16))
  plot(wm, lwd = 0.25, add = T,col = NA)
  plot(rus, lwd = 1, add = T, col = NA)
  
}

#### temperature && gph ####
load("Ta.RData")
load("hgt.RData")

ta1 = lapply(Ta1[[3]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

ta2 = lapply(Ta2[[3]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})


hgt1 = lapply(HGT1[[3]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

hgt2 = lapply(HGT2[[3]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})


mycol=colorRampPalette(c("#3949AB","#7986CB","#E8EAF6","#FFCCBC", "#FF8A65","#FF5722")) 

par(mfrow = c(2,3),mar= c(.5,.5,.5,.5),  bty = 'n')
for ( i in 1:6){
  x = ta1[[i]]
  x = rotate(x) 
  x = t(apply(x, 1, rev))
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 30, ymx = 90, crs = mycrs)
  
  x1 = hgt1[[i]]
  x1 = rotate(x1) 
  x1 = t(apply(x1, 1, rev))
  x1 = raster(x1, xmn = 0, xmx = 180, 
             ymn = 30, ymx = 90, crs = mycrs)
  plot (x, main = paste("Ta 200 mb group 1/2 week-",i), breaks = seq(-1,1,0.1),
        col = mycol(20))
  plot(wm, lwd = 0.25, add = T, col = NA)
  plot(rus, lwd = 1, add = T, col = NA)
  contour(x1, col = "black", lwd = 1, drawlabels = T,
          add = T,levels = seq(2,20,2),
          line = list(smoothing = 2))
  contour(x1, col = "black", lwd = 1, drawlabels = T,
          add = T,levels = seq(-20,0,2),lty = "dotted",
          line = list(smoothing = 2))
}

for ( i in 1:6){
  x = ta2[[i]]
  x = rotate(x) 
  x = t(apply(x, 1, rev))
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 30, ymx = 90, crs = mycrs)
  
  x1 = hgt2[[i]]
  x1 = rotate(x1) 
  x1 = t(apply(x1, 1, rev))
  x1 = raster(x1, xmn = 0, xmx = 180, 
              ymn = 30, ymx = 90, crs = mycrs)
  plot (x, main = paste("Ta 200 mb group 3/4 week-",i), breaks = seq(-1,1,0.1),
        col = mycol(20))
  plot(wm, lwd = 0.25, add = T, col = NA)
  plot(rus, lwd = 1, add = T, col = NA)
  contour(x1, col = "black", lwd = 1, drawlabels = T,
          add = T,levels = seq(4,32,4),
          line = list(smoothing = 2))
  contour(x1, col = "black", lwd = 1, drawlabels = T,
          add = T,levels = seq(-28,0,4),lty = "dotted",
          line = list(smoothing = 2))
}
#### wind field ####
load("UWND.RData")
load("VWND.RData")
load("Wind200.RData")
### 200mb N.H
u1 = lapply(U200[[1]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

u2 = lapply(U200[[2]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

v1 = lapply(V200[[1]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

v2 = lapply(V200[[2]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})


## 850mb
u1.1 = lapply(U1[[1]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

u2.1 = lapply(U2[[1]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

v1.1 = lapply(V1[[1]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

v2.1 = lapply(V2[[1]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

## 500mb
u1.2 = lapply(U1[[2]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

u2.2 = lapply(U2[[2]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

v1.2 = lapply(V1[[2]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

v2.2 = lapply(V2[[2]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

## 200mb
u1.3 = lapply(U1[[3]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

u2.3 = lapply(U2[[3]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

v1.3 = lapply(V1[[3]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

v2.3 = lapply(V2[[3]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

library(rasterVis)

for ( i in 1:6){
  x = u1.2[[i]]
  x = rotate(x) 
  x = t(apply(x, 1, rev))
  
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 30, ymx = 90, crs = mycrs)
  
  y = v1.2[[i]]
  y = rotate(y) 
  y = t(apply(y, 1, rev))
  
  y = raster(y, xmn = 0, xmx = 180, 
             ymn = 30, ymx = 90, crs = mycrs)
  
  # x <- projectRaster(x,crs=newcrs)
  # y <- projectRaster(y,crs=newcrs)
  
  
  slope <- sqrt(x^2 +y^2)
  # slope[slope < 0.1] <- NA
  aspect <- atan2(x, y)
  # aspect[slope < 0.1] <- NA
 p = vectorplot(stack(slope,aspect), isField = TRUE,lwd.arrows=0.5,
                 region = F , #par.settings = rasterTheme(region = mycol(20)),  
                 margin = FALSE, narrows = 1000, 
                 axes =FALSE,aspX=3,aspY=3,main =paste("wind 500 group 1/2 week-",i)) 
  
  
 
  # pdf(paste("wind 200 group 3/4 week-",i,".pdf"))
  print(p)
  # dev.off() #turn off develop
  }


#### temperature at 2m ####
load("T2m.RData")

ta1 = lapply(T2m[[1]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})

ta2 = lapply(T2m[[2]], FUN = function(x){
  apply(x, c(1,2), mean, na.rm = T)
})


mycol=colorRampPalette(c("#3949AB","#7986CB","#E8EAF6","#FFCCBC", "#FF8A65","#FF5722")) 

par(mfrow = c(2,3),mar= c(.5,.5,.5,.5),  bty = 'n')
for ( i in 1:6){
  x = ta1[[i]][]
  x = rotate(x) 
  x = t(apply(x, 1, rev))
  x = x[1:32,]
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 0, ymx = 90, crs = mycrs)
  
  plot (x, main = paste("T2m group 1/2 week-",i), breaks = seq(-1,1,0.1),
        col = mycol(20))
  plot(wm, lwd = 0.25, add = T, col = NA)
  plot(rus, lwd = 1, add = T, col = NA)
}

for ( i in 1:6){
  x = ta2[[i]]
  x = rotate(x) 
  x = t(apply(x, 1, rev))
  x = x[1:32,]
  
  x = raster(x, xmn = 0, xmx = 180, 
             ymn = 30, ymx = 90, crs = mycrs)
  
  plot (x, main = paste("T2m group 3/4 week-",i), breaks = seq(-1,1,0.1),
        col = mycol(20))
  plot(wm, lwd = 0.25, add = T, col = NA)
  plot(rus, lwd = 1, add = T, col = NA)
}