######################################################################
####SCRIPT FOR THE STUDY 'Assessing multitemporal calibration for ####
####              species distribution models'                    ####
######################################################################


procedure_monotemporal_sampling = function(spsNames, sampleSizes, NumRep, 
                                           envVarFolder, bgPoints){
  
  cat('[STATUS] Running `procedure_monotemporal_sampling`\n\n')
  
  # capturing start time
  timeStart = Sys.time()
  
  ##dataframes to store occurrence data
  sampleData = data.frame()
  sampleDataBg = data.frame()
  
  for (i in 1:length(spsNames)){ #loop on sps
    
    ##creating the directory structure
    
    ##creating the directory structure
    folderPath = file.path('samples/monotemporal', spsNames[i])
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
      
      for (j in 1:NumRep){ #loop over replicates for sampling scenarios

        cat('[STATUS] Sampling virtual species occurrences `', spsNames[i], 
            '`, sample size `', sSize, 
            '` and replicate `', j, '`...\n\n', sep='')
        
        #folder with the real niche maps of sp
        nicheRealFolder = file.path('virtual_sps_niche', spsNames[i]) 

        #list with the addresses of the distribution maps
        nicheRealPath = list.files(
                          path=nicheRealFolder,
                          full.names=TRUE,
                          pattern='.asc'
                        ) 

        #selecting the time layer randomly
        sampledAge = round(runif(1, 1, length(nicheRealPath)))
        scenarioName = gsub('.asc', '', basename(nicheRealPath[sampledAge]))
        projectionAge = as.numeric(scenarioName)
        
        cat('[STATUS] sampled age:', projectionAge, 'kyrBP.')

        #sampling point
        sampleData_i = dismo::randomPoints(
                          mask=raster(nicheRealPath[sampledAge]) > 0.2, 
                          prob=TRUE, 
                          n=sSize
                        )
        
        #extracting environmental variables from the point in its respective time layer
        layers_i = extract(
          x=stack(
              list.files(
                path=file.path(envVarFolder, scenarioName), 
                pattern='asc', 
                full.names=TRUE
              )
          ),
          y=sampleData_i
        )
        
        #gathering data from the time layers sampled
        sampleData = rbind(
                      sampleData, 
                      cbind(sampleData_i,layers_i,projectionAge)
                     )
        
        #adjusting the names
        names(sampleData) = c('lon', 'lat', names(as.data.frame(layers_i)), 'kyrBP')
        
        #saving
        write.csv(
            sampleData,
            file.path(
              'samples/monotemporal', 
              spsNames[i], 
              paste('occ_', sSize, 'pts_monotemporal_', 'rep', j, '.csv', sep='')
            ),
            row.names=FALSE
        )
        
        #returning empty data.frame for the next iteration
        sampleData = data.frame()
        
        ##background points##
        sampleDataBg_i = dismo::randomPoints(
                                mask=raster(nicheRealPath[sampledAge]),
                                n=bgPoints
                          ) #sampling points
        layersBg_i = extract(
          x=stack(
            list.files(
              path=file.path(envVarFolder, scenarioName), 
              pattern='asc', 
              full.names=TRUE
              )
          ),
          y=sampleDataBg_i
        ) #extracting environmental variables from the point in its respective time layer
        sampleDataBg = rbind(
                        sampleDataBg,
                        data.frame(
                          lon=sampleDataBg_i[,1],
                          lat=sampleDataBg_i[,2],
                          layersBg_i,
                          kyrBP=projectionAge
                        )
                      ) #gathering data from the time layers sampled
        names(sampleDataBg) = c('lon', 'lat', names(as.data.frame(layersBg_i)), 'kyrBP') #ajusting names
        write.csv(
            sampleDataBg,
            file.path(
              'samples/monotemporal', 
              spsNames[i], 
              paste('bg_', sSize, 'pts_monotemporal_', 'rep', j, '.csv', sep='')
            ),
            row.names=FALSE
        ) #saving
        sampleDataBg = data.frame() #returning empty data.frame for the next iteration
        
        cat('[STATUS] ...done.\n\n')
        
      }
    }
  }
  
  timeEnd = Sys.time()
  
  cat('[STATUS] Monotemporal sampling virtual species occurrences finished. (latency', 
      as.numeric(difftime(timeEnd, timeStart, unit='secs')), 
      'seconds)\n\n')
  
  
}