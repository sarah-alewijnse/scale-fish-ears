#### Myctophid joining ####

library(tidyverse)

# Load compiled data

comp <- read.csv("Outputs/04_C_resp_Components/Compiled/All_Data.csv")

#### Join to myctophids estimates ####

# Load in myctophids estimates

myct <- read.csv("Data/04_SW_Oxygen_Isotopes/Myct_d18O.csv")

# Join to new column

comp_myct <- left_join(comp, myct, by = c("My_Number", "Spp_Code"))

#### Replace for myctophids ####

comp_rep <- comp_myct

# ELN

comp_rep$mean_d18O_SW[comp_rep$Spp_Code == "ELN"] <- comp_rep$d18O_SW[comp_rep$Spp_Code == "ELN"]
comp_rep$sd_d18O_SW.x[comp_rep$Spp_Code == "ELN"] <- comp_rep$sd_d18O_SW.y[comp_rep$Spp_Code == "ELN"]

# ELC

comp_rep$mean_d18O_SW[comp_rep$Spp_Code == "ELC"] <- comp_rep$d18O_SW[comp_rep$Spp_Code == "ELC"]
comp_rep$sd_d18O_SW.x[comp_rep$Spp_Code == "ELC"] <- comp_rep$sd_d18O_SW.y[comp_rep$Spp_Code == "ELC"]

# GYR

comp_rep$mean_d18O_SW[comp_rep$Spp_Code == "GYR"] <- comp_rep$d18O_SW[comp_rep$Spp_Code == "GYR"]
comp_rep$sd_d18O_SW.x[comp_rep$Spp_Code == "GYR"] <- comp_rep$sd_d18O_SW.y[comp_rep$Spp_Code == "GYR"]

# GYN

comp_rep$mean_d18O_SW[comp_rep$Spp_Code == "GYN"] <- comp_rep$d18O_SW[comp_rep$Spp_Code == "GYN"]
comp_rep$sd_d18O_SW.x[comp_rep$Spp_Code == "GYN"] <- comp_rep$sd_d18O_SW.y[comp_rep$Spp_Code == "GYN"]

# KRA

comp_rep$mean_d18O_SW[comp_rep$Spp_Code == "KRA"] <- comp_rep$d18O_SW[comp_rep$Spp_Code == "KRA"]
comp_rep$sd_d18O_SW.x[comp_rep$Spp_Code == "KRA"] <- comp_rep$sd_d18O_SW.y[comp_rep$Spp_Code == "KRA"]

# PRM

comp_rep$mean_d18O_SW[comp_rep$Spp_Code == "PRM"] <- comp_rep$d18O_SW[comp_rep$Spp_Code == "PRM"]
comp_rep$sd_d18O_SW.x[comp_rep$Spp_Code == "PRM"] <- comp_rep$sd_d18O_SW.y[comp_rep$Spp_Code == "PRM"]

#### Tidy ####

comp_rep <- select(comp_rep, -c(d18O_SW, sd_d18O_SW.y))
colnames(comp_rep)[14] <- "sd_d18O_SW"
glimpse(comp_rep)

# Write into file

write.csv(comp_rep, "Outputs/04_C_resp_Components/Compiled/All_Data_Combined.csv", row.names = F)
