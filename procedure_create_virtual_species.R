######################################################################
####SCRIPT FOR THE STUDY 'Assessing multitemporal calibration for ####
####              species distribution models'                    ####
######################################################################

source("./rangeByAC.R")


procedure_create_virtual_species = function(envVarFolder, AmSulShape, elev, 
                                            spsTypes, projectFolder){

  cat('[STATUS] Running `procedure_create_virtual_species`\n\n')
  
  timeStart = Sys.time()
  
  ##resistance to sps moviment
  caminhosCamadasTemp = list.files(
    path = envVarFolder, 
    full.names = TRUE
  )
  elev = mask(x=elev, mask=AmSulShape) #mask for South America
  roug = terrain(x=elev, opt='roughness', unit='degrees') #rater layer for 'roughness'
  resfun = function(ro=roug){ 1/(1 + exp(+0.01*(ro - 200))) }
  resdata = calc(x = roug, fun = resfun)
  
  ##ecological niche and its geographical projection
  for (i in 1:3){ #length(caminhosCamadasTemp)){
    
    cat('[STATUS] Implementing virtual species for',
        basename(caminhosCamadasTemp[i]), 'kyr BP\n\n')
    
    ##predictor variables
    predictors = stack(file.path(caminhosCamadasTemp[i], 'bioclim_01.asc'),
                       file.path(caminhosCamadasTemp[i], 'bioclim_12.asc')) #loading
    predictors = mask(predictors, AmSulShape) #cropping

    ## scenario name
    nameScenario = basename(caminhosCamadasTemp[i])
    
    ##function for hot and wet species (HD species)
    parametersHW = formatFunctions(bioclim_01=c(fun='betaFun',
                                                p1=200,
                                                p2=295,
                                                alpha=1,
                                                gamma=1),
                                   bioclim_12=c(fun='betaFun',
                                                p1=2000,
                                                p2=3500,
                                                alpha=1,
                                                gamma=1)) #sps responses
    
    ##functions for cold and dry species (CD species)
    parametersCD = formatFunctions(bioclim_01=c(fun='betaFun',
                                                p1=50,
                                                p2=220,
                                                alpha=1,
                                                gamma=1),
                                   bioclim_12=c(fun='betaFun',
                                                p1=50,
                                                p2=1800,
                                                alpha=1,
                                                gamma=1)) #sps responses
    
    ##geographical distribution of virtual sps
    spHW = generateSpFromFun(predictors, parametersHW) 
    spCD = generateSpFromFun(predictors, parametersCD)
    
    ##stacking maps
    auxVector = stack(c(spHW$suitab.raster, spCD$suitab.raster))
    names(auxVector) = c('spHW', 'spCD')
    
    ##some adjustments and saving files
    for(j in 1:dim(auxVector)[3]){
      
      #CRS
      projection(auxVector[[j]]) = 
        CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
      
      #saving the sps distribution map
      folderPath = file.path('./virtual_sps_niche', names(auxVector[[j]]))
      if(!file.exists(folderPath)){
        
        dir.create(
          folderPath,
          recursive=TRUE
        )
        
        cat('[STATUS] Directory created:', 
            folderPath,
            '\n\n')
        
      }
      
      writeRaster(
        auxVector[[j]], 
        filename = file.path(
                    './virtual_sps_niche',
                    names(auxVector[[j]]),
                    paste(nameScenario, '.asc', sep='')
                    ),
        overwrite = TRUE,
        prj = TRUE
      )

    }
    
    cat('[STATUS] Virtual species for',
        basename(caminhosCamadasTemp[i]), 
        'kyr BP created.\n\n')
    
  }
  
  ##sps occupancy of suitable areas
  # setwd(projectFolder)
  
  for (spsType in spsTypes){
  
    cat('[STATUS] Running for:', spsType, '\n\n')
    
    folderPath = file.path('./virtual_sps_range', spsType)
    if(!file.exists(folderPath)){

      dir.create(
        folderPath,
        recursive=TRUE
        )

      cat('[STATUS] Directory created:', 
          folderPath,
          '\n\n')

    }
    
    real_niche_filepaths = list.files(
                            file.path(real_niche_folder, spsType),
                            full.names=TRUE, 
                            pattern='.asc'
                          )

    sps_range = rangeByAC(
      envAreas = raster(real_niche_filepaths[grep('002.asc', real_niche_filepaths)]),
      movRes = resdata,
      sps_range = NULL
      ) # assuming start at 120 kyrBP
    
    #updating sps range throughout the years
    for (real_niche_filepath in real_niche_filepaths){

      kyrbp = gsub('.asc', '', basename(real_niche_filepath))
      cat('[STATUS] Running for sps', spsType, 'at', kyrbp, 'KyrBP...')
      
      niche_proj = raster(real_niche_filepath)
  
      range = rangeByAC(envAreas=niche_proj,
                        movRes=resdata,
                        sps_range=sps_range,
                        iter=50)
      crs(range) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')

      writeRaster(
        range, 
        filename = file.path(
                      projectFolder,
                      'virtual_sps_range', 
                      spsType, 
                      paste(kyrbp, '.asc', sep='')
                    ),
        overwrite = TRUE,
        prj = TRUE
      )
      
      cat('done.\n\n')

    }
  }
  
  timeEnd = Sys.time()

  cat('[STATUS] Simulations finished after', 
      as.numeric(difftime(timeEnd, timeStart, unit='secs')),
      'seconds', 
      '\n\n')
  cat('[STATUS] `procedure_create_virtual_species` concluded.\n\n')
  
}
