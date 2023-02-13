#### d18O Seawater ####

library(tidyverse)
library(ncdf4)
library(fuzzyjoin) # Allows non-exact joining

# Read in location file

locs <- read.csv("Data/03_DIC_Isotopes/Unique_Locations.csv")

# Set depth NA to 0

locs$Depth_Original <- locs$Depth

locs$Depth[is.na(locs$Depth)] <- 0

#### Load the NC file ####

d18O_NC <- nc_open("Data/04_SW_Oxygen_Isotopes/calculated_d18O_v1_1.nc")

# Get NC variables

Longitude <- ncvar_get(d18O_NC, 'lon')
Latitude <- ncvar_get(d18O_NC, 'lat')
Depth <- ncvar_get(d18O_NC, 'depth')
d18O_SW <- ncvar_get(d18O_NC, 'd18o')

#### Turn into data frame ####

# Turn into vectors

long_v <- as.vector(Longitude)
lat_v <- as.vector(Latitude)
Depth_v <- as.vector(Depth)
d18O_SW_v<- as.vector(d18O_SW)

## Get dimensions

long_d <- dim(Longitude)
lat_d <- dim(Latitude)
Depth_d <- dim(Depth)

# Combine vectors and dimensions

Depth_v <- rep(Depth_v, times = lat_d[1]*long_d[1])
Depth <- sort(Depth_v, decreasing = FALSE)

Long <- rep(long_v, times = Depth_d[1]*lat_d[1])

Lat <- rep(lat_v, times = long_d[1])
Lat <- sort(Lat, decreasing = FALSE)
Lat <- rep(Lat, times = Depth_d[1])

whole <- cbind(d18O_SW_v, Depth, Lat, Long)
whole <- as.data.frame(whole)

d18O_SW <- whole
colnames(d18O_SW) <- c("d18O_SW", "Depth", "lat_dec", "long_dec")

d18O_SW <- filter(d18O_SW, !is.na(d18O_SW))

#### Fuzzy join ####

result <- data.frame()

for(i in 1:nrow(locs)){
  Stn <- as.data.frame(slice(locs, i))
  r <- difference_left_join(Stn, d18O_SW, by = c("lat_dec",
                                                 "long_dec",
                                                 "Depth"),
                            max_dist = c(1, 1, 100))
  d18O_SW_mean <- mean(r$d18O_SW, na.rm = TRUE)
  d18O_SW_sd <- sd(r$d18O_SW, na.rm = TRUE)
  dat <- data.frame(lat_dec = Stn$lat_dec,
                    long_dec = Stn$long_dec,
                    Depth = Stn$Depth_Original,
                    Year = Stn$Year,
                    d18O_SW_mean = d18O_SW_mean,
                    d18O_SW_sd = d18O_SW_sd)
  result <- rbind(result, dat)
  print(dat)
}

#### Get missing data ####

res_missing <- filter(result, is.na(d18O_SW_mean))

result_2 <- data.frame()

for(i in 1:nrow(res_missing)){
  Stn <- as.data.frame(slice(res_missing, i))
  r <- difference_left_join(Stn, d18O_SW, by = c("lat_dec",
                                                 "long_dec",
                                                 "Depth"),
                            max_dist = c(1, 1, 500))
  d18O_SW_mean <- mean(r$d18O_SW, na.rm = TRUE)
  d18O_SW_sd <- sd(r$d18O_SW, na.rm = TRUE)
  dat <- data.frame(lat_dec = Stn$lat_dec,
                    long_dec = Stn$long_dec,
                    Depth = Stn$Depth,
                    Year = Stn$Year,
                    d18O_SW_mean = d18O_SW_mean,
                    d18O_SW_sd = d18O_SW_sd)
  result_2 <- rbind(result_2, dat)
  print(result_2)
}

#### Get missing data ####

res_missing <- filter(result_2, is.na(d18O_SW_mean))

result_3 <- data.frame()

for(i in 1:nrow(res_missing)){
  Stn <- as.data.frame(slice(res_missing, i))
  r <- difference_left_join(Stn, d18O_SW, by = c("lat_dec",
                                                 "long_dec",
                                                 "Depth"),
                            max_dist = c(1, 1, 1000))
  d18O_SW_mean <- mean(r$d18O_SW, na.rm = TRUE)
  d18O_SW_sd <- sd(r$d18O_SW, na.rm = TRUE)
  dat <- data.frame(lat_dec = Stn$lat_dec,
                    long_dec = Stn$long_dec,
                    Depth = Stn$Depth,
                    Year = Stn$Year,
                    d18O_SW_mean = d18O_SW_mean,
                    d18O_SW_sd = d18O_SW_sd)
  result_3 <- rbind(result_3, dat)
  print(result_3)
}

#### Combine all ####

SW_1 <- filter(locs, !is.na(d18O_SW_mean))
SW_2 <- filter(result, !is.na(d18O_SW_mean))
SW_3 <- filter(result_2, !is.na(d18O_SW_mean))
SW_4 <- result_3

SW_combined <- rbind(SW_1, SW_2, SW_3, SW_4)
glimpse(SW_combined)

write.csv(SW_combined, "Outputs/06_Temp_Calculations/d18O_SW_combined.csv", row.names = F)
