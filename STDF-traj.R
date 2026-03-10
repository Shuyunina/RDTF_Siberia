### let's check the trajectories where the water comes from ###

data = large[,c("year","month","day","longitude","latitude","interval","julian")]
data$date = data$julian+data$interval-1

for(i in 1:nrow(data)){
  data$date2[i] =  as.Date(data$date[i], origin = as.Date(paste(data$year[i], "1-1", sep = "-")))
}

write.csv(data, file = "4hysplit.csv")
library(splitr)
for(i in nrow(data)){
  traj[[k]] <-  try(hysplit_trajectory(lat = data$latitude[i],
                                       lon = data$longitude[i],
                                       days = c(data$date2[i]-1,data$date2[i], data$date2[i]+1),
                                       height = 5,
                                       daily_hours =  c(0,3, 6,9, 12, 15,18,21),
                                       model_height = 5000,
                                       duration = 192,
                                       direction = "backward",
                                       met_type = "reanalysis",
                                       vert_motion = 0,
                                       extended_met = T,
                                       met_dir = "",
                                       exec_dir = "",
                                       config = NULL,
                                       ascdata = NULL))
  if(i %% 100 == 0 ){
    print(paste(Sys.time(),i,"I'm---GOOD!"))
  }
  save(traj, file = paste("traj_",i,".RData"))
}


###

