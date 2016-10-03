##############################################################################
##
##
##  Knudsen Plume Modeling
##
##
##############################################################################

# Substrate Definition
Substrates <- read.csv("~/R/Flux Modeling/data/Substrates.csv", row.names=1) # this is a set of rectangular substrates
View(Substrates)

# Load Reference Tables
  Periodic_Table <- read.csv("~/R/Flux Modeling/data/Periodic_Table.csv", row.names=1)
  # 1 degree =0.0174532925 radians
  degrees_to_radians <- 0.0174532925


#Choose the substrate you want to work with
Substrate <- as.list(Substrates[4,])  # 4 is the index for the G6H

# Define the resultution of the model: the length and width of the matrix that will be assigned to the plume shape
Resolution <- 100.0

# Define the scale of each point
Substrate_Scale <- as.double(Substrate$Diagonal)/100

# Define System Parameters
  # Offset of source from center of substrate in cm
  Offset <- c(26,90,90,75)                                                 
  # Vertical distance from the source to the substrate in cm
  Height <- c(60,60,60,60) 
  # Source Angle - Degrees then convert to Radians
  Phi <- c(0,0,0,0)
  Phi <- Phi*degrees_to_radians
  # area of the opening of the source in square cm
  Source_Area <- c(3^2*pi,3^2*pi,3^2*pi,3^2*pi)                 
  # Atomic Mass
  Molecular_Weight <- c(Periodic_Table["Silver","Atomic.Mass"],
                        Periodic_Table["Silver","Atomic.Mass"],
                        Periodic_Table["Silver","Atomic.Mass"],
                        Periodic_Table["Magnesium","Atomic.Mass"])  
  # Equilibrium Pressure at the source in Torr
  Pressure <- c(9*10^-4, 4*10^- 3, 4*10^-3, 2.5*10^-4)                                    
  # Temperature of source material in Kelvin
  #TODO: generate this using Vapor Pressure Curves) https://en.wikipedia.org/wiki/Vapor_pressures_of_the_elements_(data_page)
  Temperature <- c(800+273.15, 900+273.15, 900+273.15, 350+273.15)                                      

  ####################################################
  # Update this whenever a parameter above is updated
  ####################################################
  System <- data.frame(Offset, Height, Phi, Source_Area, Molecular_Weight, Pressure, Temperature)
                              
# Calculate R_A, R_B, Phi, Theta, I_prime_B
  # Calculate the Flux (Atoms/cm^2/sec) Assumes Vertical orientation
  I_A <- 1.118*10^22*System$Pressure*System$Source_Area/
    (System$Height^2*(System$Molecular_Weight*System$Temperature)^0.5)
  
  # Generate the Substrate Matrix
  X <- c((Resolution/2-Resolution):(Resolution/2))
  Y <- c((Resolution/2-Resolution):(Resolution/2))
  
  X <- as.double(X*Substrate_Scale)
  Y <- as.double(Y*Substrate_Scale)

  R_A <- (System$Height^2 + (System$Height*sin(System$Phi))^2)^0.5
  
  
  