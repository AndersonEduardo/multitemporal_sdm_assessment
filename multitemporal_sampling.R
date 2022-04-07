######################################################################
####SCRIPT FOR THE STUDY 'Assessing multitemporal calibration for ####
####              species distribution models'                    ####
######################################################################


procedure_multitemporal_sampling = function(spsNames, sampleSizes, NumRep, 
                                            envVarFolder, bgPoints){
  
  
  cat('[STATUS] Running `procedure_multitemporal_sampling`\n\n')
  
  # capturing start time
  timeStart = Sys.time()
  
  ##dataframes for storing occurrence data
  sampleData = data.frame()
  sampleDataBG = data.frame()
  
  for (i in 1:length(spsNames)){ #loop on species
    
    ##creating the directory structure
    folderPath = file.path('samples/multitemporal', spsNames[i])
    if(!file.exists(folderPath)){
      
      dir.create(
        folderPath,
        recursive=TRUE
      )
      
      cat('[STATUS] Directory created:', 
          folderPath,
          '\n\n')
      
    }
    
    for (sSize in sampleSizes){ #loop over sample sizes
      
      rangeRealFolder = file.path('virtual_sps_range', spsNames[i]) #folder with real sps distribution
      rangeRealPath = list.files(path=rangeRealFolder, full.names=TRUE, pattern='.asc') #address list
      
      sampledAges = vector()
      sampledAges = round(runif(sSize, 1, length(rangeRealPath))) #selecting 'n' time layers
      
      for (j in 1:NumRep){ #loop over replicates for sampling scenarios
        
        for (sampledAge in unique(sampledAges)){ #sampling at each time layer in the sample
          
          ## occ pts ##
          
          sampleData_i = dismo::randomPoints(
                                  mask=raster(rangeRealPath[sampledAge]) > 0.2, 
                                  prob=TRUE, 
                                  n=sum(sampledAge==sampledAges)
                         ) #sampling point
          scenarioName = gsub('.asc','', basename(rangeRealPath[sampledAge])) #time linked to the scenario for environmental variables
          projectionAge = as.numeric(scenarioName)
          
          cat('[STATUS] Sampling virtual species occurrences `', spsNames[i], 
              '`, sample size `', sSize, 
              '`, replicate `', j, 
              '`, and `', projectionAge,'` kyr BP...\n\n', sep='')
          
          
          layers_i = extract(
            x=stack(
                list.files(
                  path=file.path(envVarFolder, scenarioName),
                  pattern='bioclim',
                  full.names=TRUE
                )
            ),
            y=sampleData_i
          ) #extracting environmental variables from the point in its respective time layer
          sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,projectionAge)) #together with the data from the other time layers sampled
          
          ## background points ##
          
          envVarPath = list.files(
                        path=envVarFolder, 
                        full.names=TRUE
                       )[sampledAge] #list of the environmental variables at the time corresponding to the interaction
          envData = list.files(
                      envVarPath, 
                      full.names=TRUE, 
                      pattern='bioclim'
                    )
          sampleDataBG_i = dismo::randomPoints(
                                mask=raster(
                                      envData[1], 
                                      crs='+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
                                ),
                                n=sum(sampledAge==sampledAges)*(round(bgPoints/length(sampledAges)))
                            ) #sampling points

          layersBG_i = extract(
            x=stack(
              list.files(
                path=file.path(envVarFolder, scenarioName),
                pattern='asc', 
                full.names=TRUE
              )
            ),
            y=sampleDataBG_i
          ) #extracting environmental variables from the point in its respective time layer

          sampleDataBG = rbind(
                          sampleDataBG, 
                          data.frame(
                            lon=sampleDataBG_i[,1],
                            lat=sampleDataBG_i[,2],
                            layersBG_i,
                            kyrBP=projectionAge
                          )
                        ) #gathering data from the time layers sampled

          cat('[STATUS] ...done.\n\n' )
          
        }
        
        ## occ pts ##
        names(sampleData) = c('lon', 'lat', names(as.data.frame(layers_i)), 'kyrBP') #adjusting the names
        write.csv(
          sampleData, 
          file.path('samples/multitemporal', spsNames[i], paste('occ_', sSize, 'pts_multitemporal_rep', j, '.csv', sep='')),
          row.names=FALSE
        ) #saving
        sampleData = data.frame() #returning empty data.frame for the next iteration
        
        ## background pts ##
        names(sampleDataBG) = c('lon', 'lat', names(as.data.frame(layersBG_i)), 'kyrBP') #adjusting the names
        write.csv(
          sampleDataBG,
          file.path('samples/multitemporal', spsNames[i], paste('bg_', sSize, 'pts_multitemporal_rep', j, '.csv', sep='')),
          row.names=FALSE
        ) #saving
        sampleDataBG = data.frame() #returning empty data.frame for the next iteration
      }
    }
  }
  
  timeEnd = Sys.time()
  
  cat('[STATUS] Multitemporal sampling virtual species occurrences finished. (latency', 
       as.numeric(difftime(timeEnd, timeStart, unit='secs')), 
      'seconds)\n\n')

}
