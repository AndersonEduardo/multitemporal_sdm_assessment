######################################################################
####SCRIPT FOR THE STUDY 'Assessing multitemporal calibration for ####
####              species distribution models'                    ####
######################################################################


##necessary packages
library(virtualspecies)
library(maptools)
library(dismo)
library(raster)
library(phyloclim)
# library(SDMTools)


#########################################################
#### BEFORE BEGINING: setting up project file system ####
#########################################################

# Hint: make sure your R-session is set to a convenient 
# folder for the project

folders = c('climate_data', 'results', 'models', 'samples', 
            'virtual_sps_niche', 'virtual_sps_range')

for (folder in folders){
  
  if(!file.exists(folder)){
    
    dir.create(
      folder,
      recursive=TRUE
    )
    
    cat('[STATUS] Directory created:', 
        folder,
        '\n\n')
    
  }
}


##############################################
#### FIRST PART: creating virtual species ####
##############################################


# some necessary Parameters
envVarFolder = "climate_data" #environemntal variables
real_niche_folder = "virtual_sps_niche"
AmSulShape = rgdal::readOGR(dsn="utils/Am_Sul", layer="borders") #shapefile Neotropics
elev = raster('./utils/DEM/DEM.tif') #DEM

# loading the procedure
source('./create_virtual_species.R')

# running the procedure
procedure_create_virtual_species(

  envVarFolder = envVarFolder,
  AmSulShape   = AmSulShape,
  elev         = elev

)


#################################################################
#### SECOND PART: sampling occurrences and building datasets ####
#################################################################


# some necessary Parameters
envVarFolder = "climate_data" #folder with environmental gridfiles
spsNames = c('spHW', 'spCD') #species names
sampleSizes = c(10, 50, 100) #scenarios for sample size
NumRep = 3 #number of replicates (for each scenario implemented)
bgPoints = 10000 #number of background points

# # FOR A MINIMAL EXAMPLE, UNCOMMENT AND USE THESE PARAMETERS #
# sampleSizes = c(10)
# NumRep = 2
# #

## Multitemporal datasets ##

# loading the procedure
source('./multitemporal_sampling.R')

# running the procedure
procedure_multitemporal_sampling(
  
  spsNames       = spsNames,
  sampleSizes    = sampleSizes, 
  NumRep         = NumRep,
  envVarFolder   = envVarFolder,
  bgPoints       = bgPoints
  
)

## Monotemporal datasets ##

# loading the procedure
source('./monotemporal_sampling.R')

# running the procedure
procedure_monotemporal_sampling(
  
  spsNames      = spsNames,
  sampleSizes   = sampleSizes, 
  NumRep        = NumRep,
  envVarFolder  = envVarFolder,
  bgPoints      = bgPoints
  
)


#################################################################
###### THIRD PART: implementing species distribution models #####
#################################################################


###Some necessary Parameters###
envVarFolder = "climate_data" #folder with environmental gridfiles
maxentFolder = 'maxent' #maxent folder
spsNames = c('spHW','spCD') #sps names
sdmTypes = c('multitemporal','monotemporal') #sampling strategy
sampleSizes = c(10,50,100) #sample size
NumRep = 3 #replicates number
Tmax = 22 #max age (here, 22 kyr BP)
mainProjectFolder = getwd() #registering the project folder

# # FOR A MINIMAL EXAMPLE, UNCOMMENT AND USE THESE PARAMETERS #
# sampleSizes = c(10)
# NumRep = 2
# Tmax = 1
# #

# loading the procedure
source('./sdm.R')

# running the procedure
setwd(mainProjectFolder)
procedure_sdm(

  sdmTypes     = sdmTypes, 
  spsNames     = spsNames, 
  envVarFolder = envVarFolder, 
  maxentFolder = maxentFolder

)
setwd(mainProjectFolder)


#################################################################
###### FOURTH PART: comparing SDM projection and the        #####
######      actual spatial distribution of sps              #####
#################################################################


###Some necessary Parameters###
envVarFolder = "climate_data" #folder with environmental gridfiles
AmSulShape = rgdal::readOGR("utils/Am_Sul/borders.shp") #shape for South America
maxentFolder = 'maxent' #maxent folder
spsNames = c('spHW', 'spCD') #sps names
sdmTypes = c('multitemporal', 'monotemporal') #sampling strategy
sampleSizes = c(10,50,100) #sample sizes
NumRep = 3 #replicates number

# # FOR A MINIMAL EXAMPLE, UNCOMMENT AND USE THESE PARAMETERS #
# sampleSizes = c(10)
# NumRep = 2
# #

# loading the procedure
source('./comparing_models.R')

# running the procedure
procedure_for_comparing_models(

  sdmTypes     = sdmTypes, 
  spsNames     = spsNames, 
  sampleSizes  = sampleSizes, 
  NumRep       = NumRep,
  envVarFolder = envVarFolder,
  AmSulShape   = AmSulShape

)
