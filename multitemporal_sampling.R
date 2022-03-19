######################################################################
####SCRIPT FOR THE STUDY 'Assessing multitemporal calibration for ####
####              species distribution models'                    ####
######################################################################


multitemporal_sampling = function(){
  
  ##dataframes to store occurrence data
  sampleData = data.frame()
  sampleDataBG = data.frame()
  
  for (i in 1:length(spsTypes)){ #loop on species
    
    ##creating the directory structure
    if(!file.exists(paste(projectFolder,'/samples','/multitemporal/',spsTypes[i],sep=''))){
      dir.create(paste(projectFolder,'/samples','/multitemporal/',spsTypes[i],sep=''),recursive=TRUE)}
    
    for (sSize in sampleSizes){ #loop on sample sizes
      
      sampledAges = vector()
      sampledAges = round(runif(sSize,0,Tmax)) #selecting 'n' time layers
      rangeRealFolder = paste(projectFolder,'/RangeReal/',spsTypes[i],sep='') #folder with real sps distribution
      rangeRealPath = list.files(path=rangeRealFolder, full.names=TRUE, pattern='.asc') #address list
      
      for (j in 1:NumRep){ #loop on replicates for sampling scenarios
        
        for (sAge in unique(sampledAges)){ #sampling at each time layer in the sample
          
          ## occ pts
          
          sampleData_i = dismo::randomPoints(mask=raster(rangeRealPath[sAge+1])>0.2,prob=TRUE, n=sum(sAge==sampledAges)) #sampling point
          scenarioName = basename(rangeRealPath[1:24][sAge+1]) #time linked to the scenario for environmental variables
          scenarioName = gsub('.asc','',scenarioName) #removing '.asc' from the name
          layers_i = extract(
            x=stack(list.files(path=paste(envVarFolder,'/',scenarioName,sep=''), pattern='bioclim', full.names=TRUE)),
            y=sampleData_i) #extracting environmental variables from the point in its respective time layer
          sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sAge)) #together with the data from the other time layers sampled
          
          ## background points
          
          envVarPath = list.files(path=envVarFolder, full.names=TRUE)[sAge+1] #list of the environmental variables at the time corresponding to the interaction
          envData = list.files(envVarPath, full.names=TRUE, pattern='bioclim')
          sampleDataBG_i = dismo::randomPoints(mask=raster(envData[1], crs='+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'),
                                               n=sum(sAge==sampledAges)*(round(bgPoints/length(sampledAges)))) #sampling points
          scenarioName = list.files(path=paste(envVarFolder))[sAge+1] #scenario name
          layersBG_i = extract(
            x=stack(list.files(path=paste(envVarFolder,'/',scenarioName,sep=''), pattern='asc', full.names=TRUE)),
            y=sampleDataBG_i) #extracting environmental variables from the point in its respective time layer
          sampleDataBG = rbind(sampleDataBG, data.frame(lon=sampleDataBG_i[,1],lat=sampleDataBG_i[,2],layersBG_i,kyrBP=sAge)) #gathering data from the time layers sampled
          
        }
        
        ## occ pts
        names(sampleData) = c('lon', 'lat', names(as.data.frame(layers_i)), 'kyrBP') #adjusting the names
        write.csv(sampleData,paste(projectFolder,'/samples/multitemporal/',spsTypes[i],'/occ_',sSize,'pts_multitemporal_', j ,'rep.csv',sep=''),row.names=FALSE) #saving
        sampleData = data.frame() #returning empty data.frame for the next iteration
        
        ## background pts
        names(sampleDataBG) = c('lon','lat',names(as.data.frame(layersBG_i)),'kyrBP') #adjusting the names
        write.csv(sampleDataBG,paste(projectFolder,'/samples/multitemporal/',spsTypes[i],'/bg_',sSize,'pts_multitemporal_', j ,'rep.csv',sep=''),row.names=FALSE) #saving
        sampleDataBG = data.frame() #returning empty data.frame for the next iteration
        
      }
    }
  }
}
