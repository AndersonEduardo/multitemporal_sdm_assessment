######################################################################
####SCRIPT FOR THE STUDY 'Assessing multitemporal calibration for ####
####              species distribution models'                    ####
######################################################################

##necessary packages
library(raster)
library(ecospat)

options(java.parameters = "-Xmx7g") ###set available memmory to java

procedure_for_comparing_models = function(sdmTypes, spsNames, sampleSizes, 
                                          NumRep, envVarFolder, AmSulShape, 
                                          Tmax){
  
  cat('[STATUS] Running `procedure_for_comparing_models`\n\n')
  
  # capturing start time
  timeStart = Sys.time()
  
  
  ##creating the directory structure
  folderPath = 'results'
  if(!file.exists(folderPath)){
    
    dir.create(
      folderPath,
      recursive=TRUE
    )
    
    cat('[STATUS] Directory created:', 
        folderPath,
        '\n\n')

  }
  
  # local variables
  envVarPaths = list.files(path=envVarFolder, full.names=TRUE)
  outputData = data.frame()
  
  for (h in 1:length(sdmTypes)){ #loop over SDMs

    for (i in 1:length(spsNames)){ #loop over sps
      
      ##local parameters
      nicheRealFolder = file.path('virtual_sps_niche', spsNames[i])
      nicheRealPath = list.files(
                        path=nicheRealFolder,
                        pattern='.asc',
                        full.names=TRUE
                      )
      
      for (l in 1:length(nicheRealPath[1:Tmax])){ #loop on time layers
        
        #real distribution
        realNiche = nicheRealPath[l] 
        
        ##sampling points of the real distribution for the purchase of SDMs
        binMap = raster(realNiche) > 0.2 #binary map of the real
        realNicheDataOccCoord = dismo::randomPoints(binMap, 1000) #sampling 1000 points of the real sps distribution
        realNicheDataOccPres = extract(
                                    binMap, 
                                    realNicheDataOccCoord, 
                                    na.rm=TRUE
                                ) #extracting occurrences and absences from points
        realNicheDataOcc = data.frame(
                              longitude=realNicheDataOccCoord[,1], 
                              latitude=realNicheDataOccCoord[,2], 
                              pres=realNicheDataOccPres
                           ) #table lon, lat e pres
        predictors = stack(
                      list.files(
                          path=envVarPaths[l],
                          full.names=TRUE, 
                          pattern='.asc'
                      )
                    ) #predictors
        predictors = predictors[[c('bioclim_01','bioclim_12')]] #selecting the variables
        predictors = mask(predictors, AmSulShape) #masking environmental variables
        projection(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #adjusting CRS
        realNicheDataPred = extract(
                              x=predictors,
                              y=realNicheDataOcc[,c('longitude','latitude')],
                              na.rm=TRUE
                            ) #extracting environmental variables from the point in its respective time layer
        realNicheData = data.frame(
                          realNicheDataOcc, 
                          realNicheDataPred
                        ) #gathering data from the time layers sampled
        
        for (m in sampleSizes){ ## loop on sample size scenarios
          
          for(n in 1:NumRep){ ##loop over replicates of each combination of time and sample size
            tryCatch({
              
              # status 
              cat('[STATUS] Running for SDM ', 
                  sdmTypes[h], 
                  ', sps ', spsNames[i], 
                  ', sample size ', m,
                  ', and replicate ', n,
                  '...\n\n',
                  sep='')
              
              sdmNichePath = file.path(
                'models',
                sdmTypes[h],
                spsNames[i],
                paste(spsNames[i], '.sample', m, '.replicate', n, sep=''),
                paste('proj_', l-1, 'kyr', sep=''),
                paste('proj_', l-1, 'kyr_', spsNames[i], '.sample', m, '.replicate', n, '_TSSbin.grd', sep='')
              )
              
              # sdmNichePath = paste(
              #   'models/',
              #   sdmTypes[h],'/',
              #   spsNames[i],'/',
              #   spsNames[i],'.sample',m,'.replica',n,
              #   '/proj_',l-1,'kyr/',
              #   'proj_',l-1,'kyr_',spsNames[i],'.sample',m,'.replica',n,'_TSSbin.grd',sep='') #path to suitability map from SDM

              sdmNicheStack = stack(sdmNichePath) #reading suitability map
              binMapSDM = sdmNicheStack #adjusting object name
              
              SDMDataOccCoord = dismo::randomPoints(binMapSDM, 1000) #sampling points for PCA
              SDMDataOccPres = extract(binMapSDM, SDMDataOccCoord, na.rm=TRUE) #extracting environmental data
              SDMDataOcc = data.frame(
                            longitude=SDMDataOccCoord[,1],
                            latitude=SDMDataOccCoord[,2],
                            pres=as.numeric(SDMDataOccPres)
              ) #building table
              SDMDataPred = extract(
                              x=predictors,
                              y=SDMDataOcc[,c('longitude','latitude')],
                              na.rm=TRUE
              ) #extracting environmental variables from the point in its respective time layer
              SDMData = data.frame(SDMDataOcc, SDMDataPred) #gathering data from the time layers sampled
              SDMData = SDMData[complete.cases(SDMData), ] #some data cleaning
              
              ##The PCA is calibrated on all the sites of the study area
              pca.env = dudi.pca(
                rbind(realNicheData,SDMData)[,c('bioclim_01','bioclim_12')],
                scannf=F,
                nf=2
              )
              #ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) #grafico
              
              ##PCA scores for the whole study area
              scores.globclim = pca.env$li
              ##PCA scores for the species native distribution
              scores.sp.realNiche = suprow(
                pca.env,
                realNicheData[
                  which(realNicheData[,'pres']==1),
                  c('bioclim_01','bioclim_12')
                  ]
              )$li
              
              ##PCA scores for the species invasive distribution
              scores.sp.SDMniche = suprow(
                pca.env,
                SDMData[
                  which(SDMData[,'pres']==1),
                  c('bioclim_01','bioclim_12')
                  ]
              )$li
              
              ##PCA scores for the whole native study area
              scores.clim.realNiche = suprow(
                pca.env,
                realNicheData[,c('bioclim_01','bioclim_12')]
              )$li
              
              ##PCA scores for the whole invaded study area
              scores.clim.SDMniche = suprow(
                pca.env,
                SDMData[,c('bioclim_01','bioclim_12')]
              )$li
              
              ##gridding the native niche
              grid.clim.realNiche = ecospat.grid.clim.dyn(
                glob=scores.globclim,
                glob1=scores.clim.realNiche,
                sp=scores.sp.realNiche, 
                R=100,
                th.sp=0
              )
              
              ##gridding the invasive niche
              grid.clim.SDMniche = ecospat.grid.clim.dyn(
                glob=scores.globclim,
                glob1=scores.clim.SDMniche,
                sp=scores.sp.SDMniche, 
                R=100,
                th.sp=0
              )
              
              ##Niche equivalency
              ## OBS: Compares the observed niche overlap between z1 and z2 to 
              ## overlaps between random niches z1.sim and z2.sim, which are 
              ## built from random reallocations of occurences of z1 and z2.
              ## 'alternative' argument specifies if you want to test for niche 
              ## conservatism (alternative = "greater", i.e.  the niche overlap 
              ## is more equivalent/similar than random) or for niche divergence 
              ## (alternative = "lower", i.e. the niche overlap is less 
              ## equivalent/similar than random).
              eq.test = ecospat.niche.equivalency.test(
                grid.clim.realNiche, 
                grid.clim.SDMniche,
                rep=100,
                #alternative="greater"
              )
              
              ##Niche similarity
              ## OBS: Compares the observed niche overlap between z1 and z2 to 
              ## overlaps between z1 and random niches (z2.sim) as available in 
              ## the range of z2 (z2$Z). z2.sim has the same pattern as z2 but 
              ## the center is randomly translatated in the availabe z2$Z space 
              ## and weighted by z2$Z densities. If rand.type = 1, both z1 and 
              ## z2 are randomly shifted, if rand.type =2, only z2 is randomly 
              ## shifted. 'alternative' specifies if you want to test for niche 
              ## conservatism (alternative = "greater", i.e. the niche overlap is 
              ## more equivalent/similar than random) or for niche divergence 
              ## (alternative = "lower", i.e. the niche overlap is less 
              ## equivalent/similar than random)
              sim.test = ecospat.niche.similarity.test(
                grid.clim.realNiche, 
                grid.clim.SDMniche, 
                rep=100, 
                # alternative="greater"
              )
              
              Dobs_equiv = eq.test$obs$D #index D observed in the niche equivalence test
              Iobs_equiv = eq.test$obs$I #index I observed in the niche equivalence test
              DpValue_equiv = eq.test$p.D #p-valor for D in niche equivalence test
              IpValue_equiv = eq.test$p.I #p-valor for I in niche equivalence test
              ##
              Dobs_simi = sim.test$obs$D #index D observed in the niche similarity test
              Iobs_simi = sim.test$obs$I #indice I observado no teste de similaridade de nicho
              DpValue_simi = sim.test$p.D #p-valor indice D no teste de similaridade de nicho
              IpValue_simi = sim.test$p.I #p-valor indice I no teste de similaridade de nicho
              
              ##opening occ data to extract information of current scenario
              occPoints = read.csv(
                file.path(
                  'samples',
                  sdmTypes[h],
                  spsNames[i], 
                  paste('occ_', m, 'pts_', sdmTypes[h], '_', 'rep', n, '.csv', sep='')
                ),
                # paste(mainSampleFolder,'/',sdmTypes[h],'/',spsNames[i],'/occ_',m,'pts_',sdmTypes[h],'_', n ,'rep.csv',sep=''),
                header=TRUE
              )
              # occPoints[occPoints==0] = NA
              # occPoints = occPoints[complete.cases(occPoints),]
              occPoints = round(occPoints, digits=2)
              occPoints = occPoints[!duplicated(occPoints),]                 
              
              ## output data
              outputData = rbind(
                outputData,
                data.frame(
                  sdmType              = sdmTypes[h],
                  sp                   = spsNames[i],
                  kyrBP                = l-1,
                  sampleSize           = m,
                  replicate            = n,
                  numbOfTimeLayers     = length(unique(occPoints$kyrBP)),
                  medianKyr            = median(occPoints$kyrBP),
                  minAge               = min(occPoints$kyrBP),
                  maxAge               = max(occPoints$kyrBP),
                  Schoeners_D_equiv    = Dobs_equiv,
                  p_value_equiv        = DpValue_equiv,
                  Hellinger_I_equiv    = Iobs_equiv,
                  p_value_equiv        = IpValue_equiv,
                  Schoeners_D_simi     = Dobs_simi,
                  p_value_simi         = DpValue_simi,
                  Hellinger_I_simi     = Iobs_simi,
                  p_value_simi         = IpValue_simi
                )
            )

              # saving
              write.csv(
                outputData, 
                file = file.path('results/output.csv'),
                row.names = FALSE
              ) #saving output data
              
              # status
              cat('[STATUS] ...done.\n\n')

            }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
          }
        }
      }
    }
  }
  
  # end time
  timeEnd = Sys.time()
  
  cat('[STATUS] Model comparison procedure finished. (latency', 
      as.numeric(difftime(timeEnd, timeStart, unit='secs')), 
      'seconds)\n\n')
  
}