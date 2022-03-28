######################################################################
####SCRIPT FOR THE STUDY 'Assessing multitemporal calibration for ####
####              species distribution models'                    ####
######################################################################

##packages
library(biomod2)

options(java.parameters = "-Xmx7g") ###set available memmory to java


procedure_sdm = function(sdmTypes, spsNames, envVarFolder, Tmax, maxentFolder){
  
  cat('[STATUS] Running `procedure_monotemporal_sampling`\n\n')
  
  # capturing start time
  timeStart = Sys.time()
  
  ##adjusting directory
  mainfolderPath = file.path(getwd(), 'models')
  if(!file.exists(mainfolderPath)){
    
    dir.create(mainfolderPath, recursive=TRUE)
    
    cat('[STATUS] Directory created:', 
        mainfolderPath,
        '\n\n')
  }
  
  # local variables
  envVarPaths = list.files(
                  path=file.path(envVarFolder),
                  full.names=TRUE
                ) #lista com os caminhos das camadas no sistema (comp.)
  statResults = data.frame()
  
  
  for (h in 1:length(sdmTypes)){ #loop on SDM types (SDMmono or SDMmulti)
    
    for (i in 1:length(spsNames)){
      
      statResults = data.frame() #table of statistical outputs
      
      for (j in 1:length(sampleSizes)){

        for (k in 1:NumRep){ #loop on replicates
          tryCatch({
            
            ##adjusting directory
            setwd(mainfolderPath)
            folderPath = file.path(sdmTypes[h], spsNames[i])
            if(!file.exists(folderPath)){

              dir.create(folderPath, recursive=TRUE)

              cat('[STATUS] Directory created:', 
                  folderPath,
                  '\n\n')
            }
            setwd(folderPath)

            ##defining local parameters and local variables
            occPoints = read.csv(
              file.path(
                '..','..','..',
                'samples', 
                sdmTypes[h], 
                spsNames[i], 
                paste('occ_', sampleSizes[j], 'pts_', sdmTypes[h], '_', 'rep', k, '.csv', sep='')
              ),
              # paste(mainSampleFolder,'/',sdmTypes[h],'/',spsNames[i],'/occ_',sampleSizes[j],'pts_',sdmTypes[h],'_',k,'rep.csv',sep=''),
              header=TRUE
            ) #occurrence points
            backgroundPoints = read.csv(
              file.path(
                '..','..','..',
                'samples', 
                sdmTypes[h], 
                spsNames[i], 
                paste('bg_', sampleSizes[j], 'pts_', sdmTypes[h], '_', 'rep', k, '.csv', sep='')
              ),
              # paste(mainSampleFolder,'/',sdmTypes[h],'/',spsNames[i],'/bg_',sampleSizes[j],'pts_',sdmTypes[h],'_',k,'rep.csv',sep=''),
              header=TRUE
            ) #background points

            ##consolidating occurrence and background points
            names(backgroundPoints) = names(occPoints) #cerfying column names
            dataSet = data.frame(
                        cbind(
                          rbind(occPoints, backgroundPoints),
                          pres=c(rep(1, nrow(occPoints)), rep(0, nrow(backgroundPoints)))
                        )
                      ) #maxent input data in the format SWD (Sample With Data)
            
            ##parameters for biomod2
            myRespName = paste(spsNames[i], '_sample', sampleSizes[j],'_replicate', k, sep='') #scenario name
            myResp = dataSet[,c('pres')] #response variable
            myRespXY = dataSet[,c('lon','lat')] #coordenates linked to response variable
            myExpl = dataSet[,c('bioclim_01','bioclim_12')]  #predictor variables
            
            ##adjusting data for biomod2
            myBiomodData = BIOMOD_FormatingData(
                                  resp.var = myResp,
                                  expl.var = myExpl,
                                  resp.xy = myRespXY,
                                  resp.name = myRespName
                          )
            
            ## ##inspecionando o objeto gerado pela funcao do biomod2
            ## myBiomodData
            ## plot(myBiomodData)
            
            ##parametrizing model
            myBiomodOption = BIOMOD_ModelingOptions(
              MAXENT.Phillips=list(
                path_to_maxent.jar   = file.path('..','..','..',maxentFolder),
                maximumiterations    = 1000,
                linear               = TRUE,
                quadratic            = TRUE,
                product              = FALSE,
                threshold            = FALSE,
                hinge                = FALSE,
                maximumiterations    = 1000,
                convergencethreshold = 1.0E-5,
                threads              = 2
              )
            )

            ##runing SDM algorithm
            myBiomodModelOut = BIOMOD_Modeling(
              myBiomodData,
              models            = c('MAXENT.Phillips'),
              models.options    = myBiomodOption,
              NbRunEval         = 5, #100,
              DataSplit         = 75,
              VarImport         = 5,
              models.eval.meth  = c('TSS','ROC'),
              SaveObj           = FALSE,
              rescal.all.models = TRUE,
              do.full.models    = FALSE,
              modeling.id       = paste(myRespName)
            )
            
            ##output data
            evaluationScores = get_evaluations(myBiomodModelOut)
            
            ##recording statistical outputs
            statResults = rbind(
                  statResults,
                  cbind(
                      modelType        = sdmTypes[h],
                      sp               = spsNames[i],
                      sampleSize       = sampleSizes[j],
                      replicate        = k,
                      AUC              = mean(evaluationScores['ROC','Testing.data',,,]),
                      TSS              = mean(evaluationScores['TSS','Testing.data',,,]),
                      numbOfTimeLayers = length(unique(occPoints$kyrBP)),
                      medianKyr        =  median(occPoints$kyrBP),
                      minAge           = min(occPoints$kyrBP),
                      maxAge           = max(occPoints$kyrBP)
                  )
            )
            
            write.csv(
              statResults,
              file.path(
                # sdmTypes[h], 
                # spsNames[i], 
                paste('StatisticalResults-', spsNames[i], '.csv', sep='')
              ),
              # file=paste(projectFolder,'/maxent/',sdmTypes[h],'/',spsNames[i],'/StatisticalResults-',spsNames[i],'.csv',sep=''),
              row.names=FALSE
            )

            ##implementing model projections
            for (l in 1:length(envVarPaths[1:Tmax])){
              
              ##local parameters and local variables
              predictors = stack(
                              list.files(
                                path=file.path('..', '..', '..', envVarPaths[l]),
                                full.names=TRUE, 
                                pattern='.asc'
                              )
                            ) #predictors
              predictors = predictors[[c('bioclim_01','bioclim_12')]]
              crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #adjusting CRS
              
              ##selecting best model for projections
              whichModel = names(
                evaluationScores['TSS','Testing.data',,,][which(
                  evaluationScores['TSS','Testing.data',,,]==max(evaluationScores['TSS','Testing.data',,,])
                )]
              )
              modelName = grep(
                            pattern=whichModel, 
                            myBiomodModelOut@models.computed, 
                            value=TRUE
                          )

              ##runing projection algorithm
              myBiomodProj = BIOMOD_Projection(
                modeling.output     = myBiomodModelOut,
                new.env             = predictors,
                proj.name           = paste(l-1,'kyr',sep=''),
                selected.models     = modelName,
                binary.meth         = 'TSS',
                compress            = 'TRUE',
                build.clamping.mask = 'TRUE',
                output.format       = '.grd'
              )
              
            }
            
            # back to models directory
            setwd(mainfolderPath)

          }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        }
      }
    }
  }
  
  # back to main directory
  setwd('../')
  
  # end time
  timeEnd = Sys.time()
  
  cat('[STATUS] Species Distribution Model implementations finished. (latency', 
      as.numeric(difftime(timeEnd, timeStart, unit='secs')), 
      'seconds)\n\n')

}
