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
library(SDMTools)


##############################################
#### FIRST PART: creating virtual species ####
##############################################


# some necessary Parameters
envVarFolder = "climate_data" #environemntal variables
projectFolder = "./" #main project folder
real_niche_folder = "virtual_sps_niche"
AmSulShape = rgdal::readOGR(dsn="utils/Am_Sul", layer="borders") #shapefile Neotropics
spsTypes = c('spHW', 'spCD') #species names
elev = raster('./utils/DEM/DEM.tif') #DEM

# loading the procedure
source('./procedure_create_virtual_species.R')

# running the procedure
procedure_create_virtual_species(
  
  envVarFolder=envVarFolder,
  AmSulShape=AmSulShape,
  elev=elev,
  spsTypes=spsTypes,
  projectFolder=projectFolder

)


#################################################################
#### SECOND PART: sampling occurrences and building datasets ####
#################################################################


# some necessary Parameters
envVarFolder = "/home/anderson/Projetos/dados_ambientais/climate_data" #folder with environmental gridfiles
projectFolder = "/home/anderson/Projetos/especies_artificiais/desenvolvimento_2021" #working folder
mainSampleFolder = '/home/anderson/Projetos/especies_artificiais/desenvolvimento_2021/samples' #folde for occurrence/background datasets
AmSulShape = rgdal::readOGR("/home/anderson/Projetos/dados_ambientais/Am_Sul/borders.shp") #shape of America do Sul
crs(AmSulShape) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #adjusting shapefile
biomodFolder = '/home/anderson/Projetos/especies_artificiais/desenvolvimento_2021/biomod/' #folder for maxent outputs
spsTypes = c('spHW', 'spCD') #species names
sampleSizes = c(10, 50, 100) #sample size scenarios
NumRep = 5 #number of replicates (for each scenario implemented)
Tmax = 22 #maximum age (for the past)
bgPoints = 10000 #number of background points


## Multitemporal datasets ##

# loading the procedure
source('./multitemporal_sampling.R')

# running the procedure
multitemporal_sampling(
  
  ###############################
  ####### CONTINUAR DAQUI #######
  ###############################
  
)






## Monotemporal datasets ##

##dataframes to store occurrence data
sampleData = data.frame()
sampleDataBg = data.frame()

for (i in 1:length(spsTypes)){ #loop on sps

    ##creating the directory structure
    if(!file.exists(paste(projectFolder,'/samples','/monotemporal/',spsTypes[i],sep=''))){
        dir.create(paste(projectFolder,'/samples','/monotemporal/',spsTypes[i],sep=''), recursive=TRUE)}
    
    for (sSize in sampleSizes){ #number of points
        
        for (j in 1:NumRep){ #loop on replicates for sampling scenarios

            sampledAge = round(runif(1,0,Tmax)) #selecting the time layer randomly
            nicheRealFolder = paste(projectFolder,'/NichoReal/',spsTypes[i],sep='') #folder with the real niche maps of sp
            nicheRealPath = list.files(path=nicheRealFolder, full.names=TRUE, pattern='.asc') #list with the addresses of the distribution maps
            sampleData_i = dismo::randomPoints(mask=raster(nicheRealPath[sampledAge])>0.2, prob=TRUE, n=sSize) #sampling point
            scenarioName = basename(nicheRealPath[1:24][sampledAge]) #time for environmental variables
            scenarioName = gsub('.asc','',scenarioName) #removing '.asc' from the name
            layers_i = extract(
                x=stack(list.files(path=paste(envVarFolder,'/',scenarioName,sep=''), pattern='asc', full.names=TRUE)),
                y=sampleData_i) #extracting environmental variables from the point in its respective time layer
            sampleData = rbind(sampleData, cbind(sampleData_i,layers_i,sampledAge)) #gathering data from the time layers sampled
            names(sampleData) = c('lon','lat',names(as.data.frame(layers_i)),'kyrBP') #adjusting the names
            write.csv(sampleData,paste(projectFolder,'/samples/monotemporal/',spsTypes[i],'/occ_',sSize,'pts_monotemporal_',j,'rep','.csv',sep=''),row.names=FALSE) #saving
            
            sampleData = data.frame() #returning empty data.frame for the next iteration

            ##background points##
            sampleDataBg_i = dismo::randomPoints(mask=raster(nicheRealPath[sampledAge]),
                                                 n=bgPoints) #sampling points
            layersBg_i = extract(
                x=stack(list.files(path=paste(envVarFolder,'/',scenarioName,sep=''), pattern='asc', full.names=TRUE)),
                y=sampleDataBg_i) #extracting environmental variables from the point in its respective time layer
            sampleDataBg = rbind(sampleDataBg, data.frame(lon=sampleDataBg_i[,1],lat=sampleDataBg_i[,2],layersBg_i,kyrBP=sampledAge)) #gathering data from the time layers sampled
            names(sampleDataBg) = c('lon','lat',names(as.data.frame(layersBg_i)),'kyrBP') #ajusting names
            write.csv(sampleDataBg,paste(projectFolder,'/samples/monotemporal/',spsTypes[i],'/bg_',sSize,'pts_monotemporal_',j,'rep','.csv',sep=''),row.names=FALSE) #saving
            sampleDataBg = data.frame() #returning empty data.frame for the next iteration

        }
    }
}


#### THIRD PART: implementing species distribution models ####


##packages
library(biomod2)

###Some necessary Parameters###
##Workstation
options(java.parameters = "-Xmx7g") ###set available memmory to java
projectFolder =  "/home/anderson/Projetos/especies_artificiais/desenvolvimento_2021" #pasta do projeto
envVarFolder = "/home/anderson/Projetos/dados_ambientais/climate_data" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = rgdal::readOGR("/home/anderson/Projetos/dados_ambientais/Am_Sul/borders.shp") #shape da America do Sul
mainSampleFolder = "/home/anderson/Projetos/especies_artificiais/desenvolvimento_2021/samples" #caminho para pasta onde a planilha
maxentFolder = '/home/anderson/Projetos/maxent' #pasta para resultados do maxent
spsTypes = c('spHW','spCD') #c('spHW', 'spHD', 'spCD') #nomes das especies
sdmTypes = c('monotemporal') #c('multitemporal','monotemporal')
sampleSizes = c(10,50,100) #tamanhos das amostras
NumRep = 5 #numero de replicas (de cada cenario amostral)
statResults = data.frame() #tabela de estatisticas basicas do modelo


for (h in 1:length(sdmTypes)){ #loop on SDM types (SDMmono or SDMmulti)
    for (i in 1:length(spsTypes)){
        
        statResults = data.frame() #table of statistical outputs
        
        for (j in 1:length(sampleSizes)){
            for (k in 1:NumRep){ #loop on replicates
                tryCatch({
                    
                    ##adjusting directory
                    if(!file.exists(file.path(projectFolder,'maxent',sdmTypes[h], spsTypes[i],sep=''))){
                        dir.create(file.path(projectFolder,'maxent',sdmTypes[h],spsTypes[i],sep=''),recursive=TRUE)
                    }
                    setwd(file.path(projectFolder,'maxent',sdmTypes[h],spsTypes[i]))
                    
                    ##defining local parameters and local variables
                    occPoints = read.csv(paste(mainSampleFolder,'/',sdmTypes[h],'/',spsTypes[i],'/occ_',sampleSizes[j],'pts_',sdmTypes[h],'_',k,'rep.csv',sep=''),header=TRUE) #occurrence points
                    backgroundPoints = read.csv(paste(mainSampleFolder,'/',sdmTypes[h],'/',spsTypes[i],'/bg_',sampleSizes[j],'pts_',sdmTypes[h],'_',k,'rep.csv',sep=''),header=TRUE) #background points
                    
                    
                    ##consolidating occurrence and background points
                    names(backgroundPoints) = names(occPoints) #cerfifying column names
                    dataSet = data.frame(cbind(rbind(occPoints,backgroundPoints),pres=c(rep(1,nrow(occPoints)),rep(0,nrow(backgroundPoints))))) #maxent input data in the format SWD (Sample With Data)
                    
                    ##parameters for biomod2
                    myRespName <- paste(spsTypes[i],'_sample',sampleSizes[j],'_replica',k,sep='') #scenario name
                    myResp <- dataSet[,c('pres')] #response variable
                    myRespXY <- dataSet[,c('lon','lat')] #coordenates linked to response variable
                    myExpl = dataSet[,c('bioclim_01','bioclim_12')]  #predictor variables
                    
                    ##adjusting data for biomod2
                    myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                                         expl.var = myExpl,
                                                         resp.xy = myRespXY,
                                                         resp.name = myRespName)
                    
                    ## ##inspecionando o objeto gerado pela funcao do biomod2
                    ## myBiomodData
                    ## plot(myBiomodData)
                    
                    ##parametrizing model
                    myBiomodOption <- BIOMOD_ModelingOptions(
                        MAXENT.Phillips=list(
                            path_to_maxent.jar=maxentFolder,
                            maximumiterations=1000,
                            linear=TRUE,
                            quadratic=TRUE,
                            product=FALSE,
                            threshold=FALSE,
                            hinge=FALSE,
                            maximumiterations=1000,
                            convergencethreshold=1.0E-5,
                            threads=2))
                    
                    ##runing SDM algorithm
                    myBiomodModelOut <- BIOMOD_Modeling(
                        myBiomodData,
                        models = c('MAXENT.Phillips'),
                        models.options = myBiomodOption,
                        NbRunEval = 3, #100,
                        DataSplit = 75,
                        VarImport = 5,
                        models.eval.meth = c('TSS','ROC'),
                        SaveObj = FALSE,
                        rescal.all.models = TRUE,
                        do.full.models = FALSE,
                        modeling.id = paste(myRespName))
                    
                    ##output data
                    evaluationScores = get_evaluations(myBiomodModelOut)
                    
                    ##recording statistical outputs
                    statResults = rbind(statResults,cbind(
                                                        modelType = sdmTypes[h],
                                                        sp = spsTypes[i],
                                                        sampleSize = sampleSizes[j],
                                                        replicate = k,
                                                        AUC = mean(evaluationScores['ROC','Testing.data',,,]),
                                                        TSS = mean(evaluationScores['TSS','Testing.data',,,]),
                                                        numbOfTimeLayers = length(unique(occPoints$kyrBP)),
                                                        medianKyr = median(occPoints$kyrBP),
                                                        minAge = min(occPoints$kyrBP),
                                                        maxAge = max(occPoints$kyrBP)))
                    
                    write.csv(statResults,file=paste(projectFolder,'/maxent/',sdmTypes[h],'/',spsTypes[i],'/StatisticalResults-',spsTypes[i],'.csv',sep=''),row.names=FALSE)
                    
                    ##implementing model projections
                    for (l in 1:length(envVarPaths[1:24])){
                        
                        ##local parameters and local variables
                        predictors = stack(list.files(path=envVarPaths[l],full.names=TRUE, pattern='.asc')) #predictors
                        predictors = predictors[[c('bioclim_01','bioclim_12')]]
                        crs(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #adjusting CRS

                        ##selecting best model for projections
                        whichModel = names(evaluationScores['TSS','Testing.data',,,][which(evaluationScores['TSS','Testing.data',,,]==max(evaluationScores['TSS','Testing.data',,,]) )])
                        modelName = grep(pattern=whichModel, myBiomodModelOut@models.computed, value=TRUE)
                        
                        ##runing projection algorithm
                        myBiomodProj <- BIOMOD_Projection(
                            modeling.output = myBiomodModelOut,
                            new.env = predictors,
                            proj.name = paste(l-1,'kyr',sep=''),
                            selected.models = modelName,
                            binary.meth = 'TSS',
                            compress = 'TRUE',
                            build.clamping.mask = 'TRUE',
                            output.format = '.grd')
                        
                    }
                }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
            }
        }
    }
}

##tempo gasto
print(Sys.time() - timeStart)




#### FOURTH PART: comparing SDM projection and the actual spatial distribution of sps ####


##necessary packages
library(raster)
library(ecospat)

###Some necessary Parameters###
##Workstation
options(java.parameters = "-Xmx7g") ###set available memmory to java
projectFolder =  "J:/Anderson_Eduardo/spsArtificiais" #pasta do projeto
envVarFolder = "J:/Anderson_Eduardo/dados_projeto" #pasta com as variaveis ambientais
envVarPaths = list.files(path=envVarFolder, full.names=TRUE) #lista com os caminhos das camadas no sistema (comp.)
AmSulShape = rgdal::readOGR("J:/Anderson_Eduardo/shapefiles/Am_Sul/borders.shp") #shape da America do Sul
mainSampleFolder = "J:/Anderson_Eduardo/spsArtificiais/Amostras" #caminho para pasta onde a planilha
maxentFolder = 'C:/Users/WS/Documents/R/win-library/3.4/dismo/java' #pasta para resultados do maxent
spsTypes = c('spHW','spCD') #c('spHW', 'spHD', 'spCD') #nomes das especies
sdmTypes = c('multitemporal','monotemporal')
sampleSizes = c(10,100) #c(5,10,20,40,80,160) #tamanhos das amostras
NumRep = 5 #numero de replicas (de cada cenario amostral)
outputData = data.frame()


##algoritmo da analise do projeto
for (h in 1:length(sdmTypes)){ #loop on SDMs
    for (i in 1:length(spsTypes)){ #loop on sps

        ##local parameters
        nicheRealFolder = paste(projectFolder,'NichoReal/',spsTypes[i],sep='') #pasta com os mapas de nicho real da sp
        nicheRealPath = list.files(path=nicheRealFolder,pattern='.asc',full.names=TRUE) #lista com os enderecos dos mapas de distribuicao da sp

        for (l in 1:length(nicheRealPath[1:24])){ #loop on time layers

            realNiche = nicheRealPath[l] #real distribution

            ##sampling points of the real distribution for the purchase of SDMs
            binMap = raster(realNiche)>0.2 #binary map of the real
            realNicheDataOccCoord = dismo::randomPoints(binMap,1000) #sampling 1000 points of the real sps distribution
            realNicheDataOccPres = extract(binMap,realNicheDataOccCoord,na.rm=TRUE) #extracting occurrences and absences from points
            realNicheDataOcc = data.frame(longitude=realNicheDataOccCoord[,1], latitude=realNicheDataOccCoord[,2], pres=realNicheDataOccPres) #table lon, lat e pres
            predictors = stack(list.files(path=envVarPaths[l],full.names=TRUE, pattern='.asc')) #predictors
            predictors = predictors[[c('bioclim_01','bioclim_12')]] #selecting the variables
            predictors = mask(predictors,AmSulShape) #masking environmental variables
            projection(predictors) = CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0') #adjusting CRS
            realNicheDataPred = extract(x=predictors,y=realNicheDataOcc[,c('longitude','latitude')],na.rm=TRUE) #extracting environmental variables from the point in its respective time layer
            realNicheData = data.frame(realNicheDataOcc, realNicheDataPred) #gathering data from the time layers sampled
            

            for (m in sampleSizes){ ## loop on sample size scenarios
            
                for(n in 1:NumRep){ ##loop over replicates of each combination of time and sample size
                  tryCatch({
                                      
                    sdmNichePath = paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',m,'.replica',n,'/proj_',l-1,'kyr/','proj_',l-1,'kyr_',spsTypes[i],'.sample',m,'.replica',n,'_TSSbin.grd',sep='') #path to suitability map from SDM
                    sdmNicheStack = stack(sdmNichePath) #reading suitability map
                    binMapSDM = sdmNicheStack #adjusting object name

                    SDMDataOccCoord = dismo::randomPoints(binMapSDM, 1000) #sampling points for PCA
                    SDMDataOccPres = extract(binMapSDM, SDMDataOccCoord, na.rm=TRUE) #extracting environmental data
                    SDMDataOcc = data.frame(longitude=SDMDataOccCoord[,1],latitude=SDMDataOccCoord[,2],pres=as.numeric(SDMDataOccPres)) #building table
                    SDMDataPred = extract(x=predictors,y=SDMDataOcc[,c('longitude','latitude')],na.rm=TRUE) #extracting environmental variables from the point in its respective time layer
                    SDMData = data.frame(SDMDataOcc, SDMDataPred) #gathering data from the time layers sampled
                    SDMData = SDMData[complete.cases(SDMData),] #some data cleaning
                    
                    ##The PCA is calibrated on all the sites of the study area
                    pca.env <- dudi.pca(rbind(realNicheData,SDMData)[,c('bioclim_01','bioclim_12')],scannf=F,nf=2)
                    #ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig) #grafico

                    ##PCA scores for the whole study area
                    scores.globclim <- pca.env$li
                    ##PCA scores for the species native distribution
                    scores.sp.realNiche <- suprow(pca.env,realNicheData[which(realNicheData[,'pres']==1),c('bioclim_01','bioclim_12')])$li

                    ##PCA scores for the species invasive distribution
                    scores.sp.SDMniche <- suprow(pca.env,SDMData[which(SDMData[,'pres']==1),c('bioclim_01','bioclim_12')])$li

                    ##PCA scores for the whole native study area
                    scores.clim.realNiche <-suprow(pca.env,realNicheData[,c('bioclim_01','bioclim_12')])$li

                    ##PCA scores for the whole invaded study area
                    scores.clim.SDMniche <- suprow(pca.env,SDMData[,c('bioclim_01','bioclim_12')])$li

                    ##gridding the native niche
                    grid.clim.realNiche <-ecospat.grid.clim.dyn(glob=scores.globclim,glob1=scores.clim.realNiche,sp=scores.sp.realNiche, R=100,th.sp=0)

                    ##gridding the invasive niche
                    grid.clim.SDMniche <- ecospat.grid.clim.dyn(glob=scores.globclim,glob1=scores.clim.SDMniche,sp=scores.sp.SDMniche, R=100,th.sp=0)

                    ##Niche equivalency
                    ##OBS: Compares the observed niche overlap between z1 and z2 to overlaps between random niches z1.sim
                    ## and z2.sim, which are built from random reallocations of occurences of z1 and z2.
                    ##'alternative' argument specifies if you want to test for niche conservatism (alternative = "greater", i.e.  the
                    ## niche overlap is more equivalent/similar than random) or for niche divergence (alternative = "lower",
                    ## i.e. the niche overlap is less equivalent/similar than random).
                    eq.test <- ecospat.niche.equivalency.test(grid.clim.realNiche, grid.clim.SDMniche,rep=100, alternative = "greater")

                    ##Niche similarity
                    ##OBS: Compares the observed niche overlap between z1 and z2 to overlaps between z1 and random niches
                    ## (z2.sim) as available in the range of z2 (z2$Z). z2.sim has the same pattern as z2 but the center is
                    ## randomly translatated in the availabe z2$Z space and weighted by z2$Z densities. If rand.type = 1,
                    ## both z1 and z2 are randomly shifted, if rand.type =2, only z2 is randomly shifted.
                    ## 'alternative' specifies if you want to test for niche conservatism (alternative = "greater", i.e.  the
                    ## niche overlap is more equivalent/similar than random) or for niche divergence (alternative = "lower",
                    ## i.e. the niche overlap is less equivalent/similar than random)
                    sim.test <- ecospat.niche.similarity.test(grid.clim.realNiche, grid.clim.SDMniche, rep=100, alternative = "greater")
                    
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
                    occPoints = read.csv(paste(mainSampleFolder,'/',sdmTypes[h],'/',spsTypes[i],'/occ_',m,'pts_',sdmTypes[h],'_', n ,'rep.csv',sep=''),header=TRUE) 
                    occPoints[occPoints==0] = NA
                    occPoints = occPoints[complete.cases(occPoints),]
                    occPoints = round(occPoints, digits=2)
                    occPoints = occPoints[!duplicated(occPoints),]                 
                    
		    ## output data
                    outputData = rbind(outputData,data.frame(sdmType = sdmTypes[h],
                                                        sp = spsTypes[i],
                                                        kyrBP = l-1,
                                                        sampleSize = m,
                                                        replicate = n,
                                                        numbOfTimeLayers = length(unique(occPoints$kyrBP)),
                                                        medianKyr = median(occPoints$kyrBP),
                                                        minAge = min(occPoints$kyrBP),
                                                        maxAge = max(occPoints$kyrBP),
                                                        Schoeners_D_equiv = Dobs_equiv,
                                                        p_value_equiv = DpValue_equiv,
                                                        Hellinger_I_equiv = Iobs_equiv,
                                                        p_value_equiv = IpValue_equiv,
                                                        Schoeners_D_simi = Dobs_simi,
                                                        p_value_simi = DpValue_simi,
                                                        Hellinger_I_simi = Iobs_simi,
                                                        p_value_simi = IpValue_simi))
                    
                    write.csv(outputData, file=paste(projectFolder,'/maxent/output.csv',sep=''),row.names=FALSE) #saving output data
                    
                  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                }
            }
        }
    }
}

##time record
print(Sys.time() - timeStart)




#### FIFTH PART: building results charts ####


###Some necessary Parameters###
##Workstation
spsTypes = c('spHW', 'spCD') #c('spHW', 'spHD', 'spCD') #nomes das especies
outputData = list() #tabela de dados de saida
vetor.nomes = vector()
projectFolder = "J:/Anderson_Eduardo/spsArtificiais" #pasta do projeto

##My notebook###############
##spsTypes = c('spHW', 'spCD') #c('spHW', 'spHD', 'spCD') #nomes das especies
##outputData = list() #tabela de dados de saida
##vetor.nomes = vector()
##projectFolder = "/home/anderson/Projetos/Sps artificiais" #pasta do projeto
#projectFolder = '/media/anderson/PIBi/ANDERSON EDUARDO/Sps artificiais'
############################


### AUC and TSS of the models

##multitemporal, spHW

spHWmulti = read.csv(paste(projectFolder,'/maxent/multitemporal/spHW/StatisticalResults-spHW.csv', sep=''), header=TRUE)

##multitemporal, spCD
spCDmulti = read.csv(paste(projectFolder,'/maxent/multitemporal/spCD/StatisticalResults-spCD.csv', sep=''), header=TRUE)

##monotemporal, spHW
spHWmono = read.csv(paste(projectFolder,'/maxent/monotemporal/spHW/StatisticalResults-spHW.csv', sep=''), header=TRUE)

##monotemporal, spCD
spCDmono = read.csv(paste(projectFolder,'/maxent/monotemporal/spCD/StatisticalResults-spCD.csv', sep=''), header=TRUE)


## boxplots models X AUC and TSS, full dataset

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotModelos&Acuracia_dadosTotais.jpeg', width=800)
par(mfrow=c(1,2), las=2, mar=c(8,5,1,1), cex=1.3)
boxplot(rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$AUC ~  rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$modelType, ylim=c(0,1), ylab='AUC')
boxplot(rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$TSS ~  rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$modelType, ylim=c(0,1), ylab='TSS')
dev.off()

##tests (both statistically significative)

kruskal.test( rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$AUC ~  rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$modelType )
kruskal.test( rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$TSS ~  rbind(spHWmulti,spCDmulti,spHWmono,spCDmono)$modelType )



## boxplots models X AUC and TSS, species

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotModelos&Acuracia_sps.jpeg', height=600)
par(mfrow=c(2,2), las=2, mar=c(8,5,2,1), cex=1.1)
boxplot(rbind(spHWmulti,spHWmono)$AUC ~  rbind(spHWmulti,spHWmono)$modelType, ylim=c(0,1), ylab=c('AUC'), main='HW species')
boxplot(rbind(spCDmulti,spCDmono)$AUC ~  rbind(spCDmulti,spCDmono)$modelType, ylim=c(0,1), ylab=c('AUC'), main='CD species')
boxplot(rbind(spHWmulti,spHWmono)$TSS ~  rbind(spHWmulti,spHWmono)$modelType, ylim=c(0,1), ylab=c('TSS'), main='HW species')
boxplot(rbind(spCDmulti,spCDmono)$TSS ~  rbind(spCDmulti,spCDmono)$modelType, ylim=c(0,1), ylab=c('TSS'), main='CD species')
dev.off()

##tests (all statistically significative)

kruskal.test(rbind(spHWmulti,spHWmono)$AUC ~  rbind(spHWmulti,spHWmono)$modelType)
kruskal.test(rbind(spCDmulti,spCDmono)$AUC ~  rbind(spCDmulti,spCDmono)$modelType)
kruskal.test(rbind(spHWmulti,spHWmono)$TSS ~  rbind(spHWmulti,spHWmono)$modelType)
kruskal.test(rbind(spCDmulti,spCDmono)$TSS ~  rbind(spCDmulti,spCDmono)$modelType)


## boxplots sp X AUC and TSS, full dataset

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotSps&Acuracia.jpeg')
par(mfrow=c(2,2), mar=c(3,4,5,1),cex=1.1)
boxplot(rbind(spHWmulti,spCDmulti)$AUC ~  rbind(spHWmulti,spCDmulti)$sp, ylim=c(0,1), ylab=c('AUC'), main='Multitemporal')
boxplot(rbind(spHWmono,spCDmono)$AUC ~  rbind(spHWmono,spCDmono)$sp, ylim=c(0,1), ylab=c('AUC'), main='Monotemporal')
boxplot(rbind(spHWmulti,spCDmulti)$TSS ~  rbind(spHWmulti,spCDmulti)$sp, ylim=c(0,1), ylab=c('TSS'), main='Multitemporal')
boxplot(rbind(spHWmono,spCDmono)$TSS ~  rbind(spHWmono,spCDmono)$sp, ylim=c(0,1), ylab=c('TSS'), main='Monotemporal')
dev.off()


##tests (both AUC and TSS statistically significative only for SDMmulti)

kruskal.test(rbind(spHWmulti,spCDmulti)$AUC ~  rbind(spHWmulti,spCDmulti)$sp) 
kruskal.test(rbind(spHWmono,spCDmono)$AUC ~  rbind(spHWmono,spCDmono)$sp)
kruskal.test(rbind(spHWmulti,spCDmulti)$TSS ~  rbind(spHWmulti,spCDmulti)$sp)
kruskal.test(rbind(spHWmono,spCDmono)$TSS ~  rbind(spHWmono,spCDmono)$sp)


## boxplots sampleSize X AUC aand TSS, full dataset

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotSampleSize&Acuracia_dadosTotais.jpeg')
par(mfrow=c(2,2))
boxplot(rbind(spHWmulti,spCDmulti)$AUC ~  rbind(spHWmulti,spCDmulti)$sampleSize, ylim=c(0,1), ylab='AUC', main='Multitemporal')
boxplot(rbind(spHWmono,spCDmono)$AUC ~  rbind(spHWmono,spCDmono)$sampleSize, ylim=c(0,1), ylab='AUC', main='Monotemporal')
boxplot(rbind(spHWmulti,spCDmulti)$TSS ~  rbind(spHWmulti,spCDmulti)$sampleSize, ylim=c(0,1), ylab='TSS', main='Multitemporal')
boxplot(rbind(spHWmono,spCDmono)$TSS ~  rbind(spHWmono,spCDmono)$sampleSize, ylim=c(0,1), ylab='TSS', main='Monotemporal')
dev.off()

## boxplots sampleSize X AUC, species

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotSampleSize&Acuracia_spHW.jpeg')
par(mfrow=c(2,2), mar=c(3,4,5,1))
boxplot(rbind(spHWmulti)$AUC ~  rbind(spHWmulti)$sampleSize, ylim=c(0,1), ylab='AUC', main='Multitemporal')
boxplot(rbind(spHWmono)$AUC ~  rbind(spHWmono)$sampleSize, ylim=c(0,1), ylab='AUC', main='Monotemporal')
boxplot(rbind(spHWmulti)$TSS ~  rbind(spHWmulti)$sampleSize, ylim=c(0,1), ylab='TSS', main='Multitemporal')
boxplot(rbind(spHWmono)$TSS ~  rbind(spHWmono)$sampleSize, ylim=c(0,1), ylab='TSS', main='Monotemporal')
title('spHW', outer=TRUE, line=-1)
dev.off()

## boxplots modelo X AUC, species

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotSampleSize&Acuracia_spCD.jpeg')
par(mfrow=c(2,2), mar=c(3,4,5,1))
boxplot(rbind(spCDmulti)$AUC ~  rbind(spCDmulti)$sampleSize, ylim=c(0,1), ylab='AUC', main='Multitemporal')
boxplot(rbind(spCDmono)$AUC ~  rbind(spCDmono)$sampleSize, ylim=c(0,1), ylab='AUC', main='Monotemporal')
boxplot(rbind(spCDmulti)$TSS ~  rbind(spCDmulti)$sampleSize, ylim=c(0,1), ylab='TSS', main='Multitemporal')
boxplot(rbind(spCDmono)$TSS ~  rbind(spCDmono)$sampleSize, ylim=c(0,1), ylab='TSS', main='Monotemporal')
title('spCD', outer=TRUE, line=-1)
dev.off()


### niche overlap

outputData = read.csv(file=paste(projectFolder,'/maxent/output.csv',sep=''), header=TRUE)
##outputData = read.csv(file=paste(projectFolder,'/Resultados Lorien/output.csv',sep=''),header=TRUE)
##outputData = read.csv(file=paste(projectFolder,'/maxent/output_uni.csv',sep=''), header=TRUE)
#vetor.nomes = append(vetor.nomes,paste(spsTypes[i],sep=''))


## Schoener e Hellinger for full dataset
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotDadosTotais.jpeg', width=600)
par(mfrow=c(1,2), mar=c(8,3,3,1), cex=1.4, las=2)
boxplot(outputData$Schoeners_D_simi~ outputData$sdmType, ylim=c(0,1), main="Schoeners' D")
boxplot(outputData$Hellinger_I_simi~ outputData$sdmType, ylim=c(0,1), main='Hellinger')
dev.off()

##tests (both statistically non-significative)
kruskal.test(Schoeners_D_simi ~ sdmType, data = outputData)
kruskal.test(Hellinger_I_simi ~ sdmType, data = outputData)


## bosplots for sps
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplotSps.jpeg', height=650)
par(mfrow=c(2,2), mar=c(7,4.5,6,1), cex=1.1, las=2)# cex.axis=2.5, cex.lab=3, cex.main=3)
boxplot(outputData[outputData$sp == 'spHW',]$Schoeners_D_simi ~ outputData[outputData$sp == 'spHW',]$sdmType, ylim=c(0,1), ylab="Schoener's D", main='HW species')
boxplot(outputData[outputData$sp == 'spHW',]$Hellinger_I_simi ~ outputData[outputData$sp == 'spHW',]$sdmType, ylim=c(0,1), ylab="Hellinger", main='HW species')
boxplot(outputData[outputData$sp == 'spCD',]$Schoeners_D_simi ~ outputData[outputData$sp == 'spCD',]$sdmType, ylim=c(0,1), ylab="Schoener's D",  main='CD species')
boxplot(outputData[outputData$sp == 'spCD',]$Hellinger_I_simi ~ outputData[outputData$sp == 'spCD',]$sdmType, ylim=c(0,1), ylab="Hellinger", main='CD species')
dev.off()

##tests (all statistically non-significative)
kruskal.test(outputData[outputData$sp == 'spHW',]$Schoeners_D_simi ~ outputData[outputData$sp == 'spHW',]$sdmType)
kruskal.test(outputData[outputData$sp == 'spHW',]$Hellinger_I_simi ~ outputData[outputData$sp == 'spHW',]$sdmType)
kruskal.test(outputData[outputData$sp == 'spCD',]$Schoeners_D_simi ~ outputData[outputData$sp == 'spCD',]$sdmType)
kruskal.test(outputData[outputData$sp == 'spCD',]$Hellinger_I_simi ~ outputData[outputData$sp == 'spCD',]$sdmType)

##Density for full dataset
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/densidadeDadosTotais.jpeg', width=600, height = 400)
par(mfrow=c(1,2), lwd=2, cex=1)
plot(density(outputData[outputData$sdmType == 'multitemporal',]$Schoeners_D_simi),ylim=c(0,5), lwd=2, col='red', main='', xlab="Schoener's D", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal',]$Schoeners_D_simi), lwd=2)
#
plot(density(outputData[outputData$sdmType == 'multitemporal',]$Hellinger_I_simi),ylim=c(0,5), lwd=2, col='red', main='', xlab='Hellinger', ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal',]$Hellinger_I_simi), lwd=2)
##
legend(x='topright', legend=c('Multitemporal calibration', 'Monotemporal calibration'), lty=1, col=c('red','black'), bty='n')
dev.off()

##Density for full dataset
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/densidade_sps.jpeg')
par(mfrow=c(2,2), mar=c(5,4,3,1), lwd=2, cex=1)
#
plot(density(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi),ylim=c(0,7), lwd=2, col='red', main='HW species', xlab="Schoener's D", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi), lwd=2)
#
plot(density(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi),ylim=c(0,5), lwd=2, col='red', main='HW species', xlab="Hellinger", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi), lwd=2)
legend(x='topright', legend=c('Multitemporal calibration', 'Monotemporal calibration'), lty=1, col=c('red','black'), bty='n', cex=0.8)
#
plot(density(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi),ylim=c(0,5), lwd=2, col='red', main='CD species', xlab="Schoener's D", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi), lwd=2)
#
plot(density(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi),ylim=c(0,5), col='red', lwd=2, main='CD species', xlab="Hellinger", ylab='Density')
lines(density(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi), lwd=2)
dev.off()

##Schoener and Hellinger for time
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Shoener&HellingerXtempo.jpeg',width=600, height=600)
par(mfrow=c(2,2), mar=c(4,4,4,1), cex=1.2)
plot(outputData[outputData$sdmType == 'multitemporal',]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal',]$kyrBP),type='p',ylab="Schoeners' D", xlab="Time (kyr BP)", ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), main='Multitemporal')
#
plot(outputData[outputData$sdmType == 'multitemporal',]$Hellinger_I_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal',]$kyrBP),type='p',ylab="Hellinger",xlab="Time (kyr BP)",ylim=c(0,1),col=rgb(0,0,0,alpha=0.5), main='Multitemporal')
#
plot(outputData[outputData$sdmType == 'monotemporal',]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal',]$kyrBP),type='p',ylab="Schoeners' D",xlab="Time (kyr BP)",ylim=c(0,1),col=rgb(0,0,0,alpha=0.5), main='Monotemporal')
#
plot(outputData[outputData$sdmType == 'monotemporal',]$Hellinger_I_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal',]$kyrBP),type='p',ylab="Hellinger",xlab="Time (kyr BP)",ylim=c(0,1),col=rgb(0,0,0,alpha=0.5), main='Monotemporal')
dev.off()

## Shoener e Hellinger for tempo - species
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Shoener&HellingerXtempo_sps.jpeg', width=1200, height=1200)
par(mfrow=c(2,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='HW species', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$kyrBP),type='p',ylab="Hellinger",xlab="Time (kyr BP)", main='HW species',ylim=c(0,1), col=rgb(0,0,0,alpha=0.5))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$kyrBP),type='p',ylab="Schoeners' D",xlab="Time (kyr BP)", main='CD species', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$kyrBP),type='p',ylab="Hellinger",xlab="Time (kyr BP)", main='CD species',ylim=c(0,1), col=rgb(0,0,0,alpha=0.5))
dev.off()

## Schoener X sample size X time - spHW
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/SchoenerXtempoXsample_spHW.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
dev.off()

## Schoener X sample size X time - spCD
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/SchoenerXtempoXsample_spCD.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$Schoeners_D_simi ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Schoener's D",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
dev.off()


## Hellinger X sample size X time - spHW
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/HellingerXtempoXsample_spHW.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
dev.off()

## Hellinger X sample size X time - spCD
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/HellingerXtempoXsample_spCD.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(7,7,3,3), cex=1.5, cex.lab=2, cex.axis=2, cex.main=2)
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 10,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='10 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 50,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='50 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
#
plot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$Hellinger_I_equiv ~ as.factor(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD' & outputData$sampleSize == 100,]$kyrBP),type='p',ylab="Hellinger's I",xlab="Time (kyr BP)", main='100 pts', ylim=c(0,1), col=rgb(0,0,0,alpha=0.5), yaxt='n')
axis(side=2, seq(0,1,by=0.1), labels=c(0,NA,NA,NA,NA,0.5,NA,NA,NA,NA,1))
dev.off()



## Sample size - full dataset
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplot_sampleSize_dadosTotais.jpeg', width=800, height=900)
par(mfrow=c(2,2), cex=1.5)
boxplot(outputData[outputData$sdmType == 'multitemporal',]$Schoeners_D_simi~ outputData[outputData$sdmType == 'multitemporal',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoeners' D", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'monotemporal',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoeners' D", main='Monotemporal')
boxplot(outputData[outputData$sdmType == 'multitemporal',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'multitemporal',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'monotemporal',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Monotemporal')
dev.off()

## Sample size - spHW
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplot_sampleSize_spHW.jpeg', width=800, height=900)
par(mfrow=c(2,2), cex=1.5)
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoener's D", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoener's D", main='Monotemporal')
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spHW',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Monotemporal')
dev.off()

## Sample size - spCD
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplot_sampleSize_spCD.jpeg', width=800, height=900)
par(mfrow=c(2,2), cex=1.5)
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoeners' D", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Schoener's D", main='Monotemporal')
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Multitemporal')
boxplot(outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'monotemporal' & outputData$sp == 'spCD',]$sampleSize, ylim=c(0,1), xlab='Sample Size', ylab="Hellinger", main='Monotemporal')
dev.off()

## Shoener's D e Hellinger X number of time layers - SDMmulti
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/boxplot_NumberOfTimeLayers.jpeg', height=1000, width=600)
par(mfrow=c(3,2), cex=1.3)
boxplot(outputData[outputData$sdmType == 'multitemporal',]$Schoeners_D_simi ~ outputData[outputData$sdmType == 'multitemporal',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Schoener's D", main='Full dataset')
##
boxplot(outputData[outputData$sdmType == 'multitemporal',]$Hellinger_I_simi ~ outputData[outputData$sdmType == 'multitemporal',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Hellinger", main='Full dataset')
##
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Schoeners_D_simi~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Schoener's D", main='HW species')
##
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$Hellinger_I_simi~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spHW',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Hellinger", main='HW species')
##
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Schoeners_D_simi~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Schoener's D", main='CD species')
##
boxplot(outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$Hellinger_I_simi~ outputData[outputData$sdmType == 'multitemporal' & outputData$sp == 'spCD',]$numbOfTimeLayers, xlab='Number of time layers', ylab="Hellinger", main='CD species')
dev.off()


##graphics for clamping

projectFolder = "/home/anderson/Projetos/Sps artificiais/"
sdmTypes = c("multitemporal", "monotemporal")
spsTypes = c("spHW", "spCD")
sampleSizes = c(10, 50, 100)
numRep = 5
clampList = list()
territory = list()


for(h in 1:length(sdmTypes)){
    for(i in 1:length(spsTypes)){
        for(m in 1:length(sampleSizes)){
            for(n in 1:numRep){
                for(l in 1:24){
                    ##clamping maps
                    sdmClampPath = paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[m],'.replica',n,'/proj_',l-1,'kyr/','proj_',l-1,'kyr_ClampingMask.grd',sep='') #path to suitability mapsfrom SDM
                    clampLayer_i = raster(sdmClampPath)
                    scenName = paste(sdmTypes[h],'_proj_',l-1,'kyr_',spsTypes[i],'.sample',sampleSizes[m],'.replica',n,sep='')
                    clampList[[scenName]] = clampLayer_i
                    clamping = (sum(getValues(clampLayer_i)>0, na.rm=TRUE)/ncell(getValues(clampLayer_i))) * 100
                    ##map of species distribution
                    sdmDistPath = paste(projectFolder,'maxent/',sdmTypes[h],'/',spsTypes[i],'/',spsTypes[i],'.sample',sampleSizes[m],'.replica',n,'/proj_',l-1,'kyr/','proj_',l-1,'kyr_',spsTypes[i],'.sample',sampleSizes[m],'.replica',n,'_TSSbin.grd',sep='') #path to suitability mapsfrom SDM
                    distLayer_i = raster(sdmDistPath)
                    ##computing extent of clamping being performed in modelled distribution
                    distUnderClamp = (clampLayer_i + distLayer_i)==2
                    distUnderClamp = ( freq(distUnderClamp, value=1)/sum(freq(distUnderClamp)[1:2,2]) ) * 100
                    ##output table
                    outputDF = outputData[ which(outputData$sdmType==sdmTypes[h] & outputData$sp==spsTypes[i] & outputData$sampleSize==sampleSizes[m] & outputData$replicate==n & outputData$kyrBP==l-1 ), ]
                    if(nrow(outputDF)>0){
                        territory[[scenName]] = data.frame( outputDF,
                                                           clamping = clamping,
                                                           distUnderClamp = distUnderClamp )
                    }
                }
            }
        }
    }
}




##transforming the list into stack of gridfiles
clampStack = stack(clampList)


##multitemporal

## all clamping maps - multitemporal, spHW, sample size 10 pts (all temporal layers)
scenNames =  grep(pattern='^multitemporal.*spHW.*sample10.*replica1', x=names(clampStack), value=TRUE) #separating the names
scenNames =  grep(pattern='^multitemporal.*spHW.*sample100.*replica1', x=scenNames, value=TRUE, invert=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Multitemporal_SpHW_sample10_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 10',
                     names.attr=c(paste('spHW ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - multitemporal, spCD, sample 10 pts (all temporal layers)
scenNames =  grep(pattern='^multitemporal.*spCD.*sample10.*replica1', x=names(clampStack), value=TRUE) #separating the names
scenNames =  grep(pattern='^multitemporal.*spCD.*sample100.*replica1', x=scenNames, value=TRUE, invert=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Multitemporal_SpCD_sample10_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 10',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - multitemporal, spHW, sample 50 pts (all time layers)
scenNames =  grep(pattern='^multitemporal.*spHW.*sample50.*replica1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Multitemporal_SpHW_sample50_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 50',
                     names.attr=c(paste('HW sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - multitemporal, spCD, sample 50 pts (all temporal layers)
scenNames =  grep(pattern='^multitemporal.*spCD.*sample50.*replica1', x=names(clampStack), value=TRUE) #separando os nomes
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Multitemporal_SpCD_sample50_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 50',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - multitemporal, spHW, sample 100 pts (all temporal layers)
scenNames =  grep(pattern='^multitemporal.*spHW.*sample100.*replica1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Multitemporal_SpHW_sample100_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 100',
                     names.attr=c(paste('HW sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - multitemporal, spCD, sample 100 pts (all temporal layers)
scenNames =  grep(pattern='^multitemporal.*spCD.*sample100.*replica1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Multitemporal_SpCD_sample100_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 100',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()


###monotemporal

## all clamping maps - monotemporal, spHW, sample 10 pts (all time layers)
scenNames =  grep(pattern='^monotemporal.*spHW.*sample10.*replica1', x=names(clampStack), value=TRUE) #separating the names
scenNames =  grep(pattern='^monotemporal.*spHW.*sample100.*replica1', x=scenNames, value=TRUE, invert=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Monotemporal_SpHW_sample10_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 10',
                     names.attr=c(paste('HW sp. ',0:23,'kyr BP',sep='')))
dev.off()

## all clamping maps - monotemporal, spCD, sample 10 pts (all time layers)
scenNames =  grep(pattern='^monotemporal.*spCD.*sample10.*replica1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Monotemporal_SpCD_sample10_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 10',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - monotemporal, spHW, sample 50 pts (all time layers)
scenNames =  grep(pattern='^monotemporal.*spHW.*sample50.*replica1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Monotemporal_SpHW_sample50_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 50',
                     names.attr=c(paste('HW sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - monotemporal, spCD, sample 50 pts (all time layers)
scenNames =  grep(pattern='^monotemporal.*spCD.*sample50.*replica1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Monotemporal_SpCD_sample50_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 50',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - monotemporal, spHW, sample 100 pts (all time layers)
scenNames =  grep(pattern='^monotemporal.*spHW.*sample100.*replica1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Monotemporal_SpHW_sample100_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spHW, sample size = 100',
                     names.attr=c(paste('HW sp. ',0:22,'kyr BP',sep='')))
dev.off()

## all clamping maps - monotemporal, spCD, sample 100 pts (all temporal layers)
scenNames =  grep(pattern='^monotemporal.*spCD.*sample100.*replica1', x=names(clampStack), value=TRUE) #separating the names
clampStack[[scenNames]] = mask(clampStack[[scenNames]], AmSulShape)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/Clamping/clamp_Monotemporal_SpCD_sample100_replica1.jpg', width=1000, height=1100)
rasterVis::levelplot(clampStack[[scenNames[1:23]]],
                     col.regions=colorRampPalette(c("lightgrey","red")),
                     main='spCD, sample size = 100',
                     names.attr=c(paste('CD sp. ',0:22,'kyr BP',sep='')))
dev.off()




### trend graphs ###



##dataClamp = data.frame(clamping=as.numeric(territory), scenario=names(territory))
dataClamp = do.call('rbind', territory)
dataClamp$scenario = names(territory)
rownames(dataClamp) = NULL

######
## dataClamp$sdm = NA
## dataClamp[ grep('monotemporal', dataClamp$scenario), ]$sdm = 'monotemporal'
## dataClamp[ grep('multitemporal', dataClamp$scenario), ]$sdm = 'multitemporal'
## dataClamp$kyr = c(0:23)
## dataClamp$sp = NA
## dataClamp[grep('spHW', dataClamp$scenario),]$sp = 'spHW'
## dataClamp[grep('spCD', dataClamp$scenario),]$sp = 'spCD'
## dataClamp$sample = NA
## dataClamp[grep('10.replica', dataClamp$scenario),]$sample = 10
## dataClamp[grep('50.replica', dataClamp$scenario),]$sample = 50
## dataClamp[grep('100.replica', dataClamp$scenario),]$sample = 100
######


## clamping for America do Sul X Schoener's D (D = overlap between real X modelled niche-based distribution) ##

jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/clampAmSulXschoenerXsampleXsps.jpeg', width=1200, height=1200)
par(mfrow=c(2,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.5, cex.axis=2, cex.main=2)
##
##spCD monotemporal
##
plot(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("CD species - SDM"["mono"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping), ylim=c(0,1), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping), ylim=c(0,1), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping)), lty=2, col='black')
##
##spCD multitemporal
##
plot(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(South America area (in %))", main=expression("CD species - SDM"["multi"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping), ylim=c(0,1), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping), ylim=c(0,1), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping)), lty=2, col='black')
##
##spHW monotemporal
##
plot(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(South America area (in %))", main=expression("HW species - SDM"["mono"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping)), lty=2, col='black')
##
##spHW multitemporal
##
plot(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(South America area (in %))", main=expression("HW species - SDM"["multi"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$clamping>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping)), lty=2, col='black')
##
dev.off()



## clamping across modelled distribution X Schoener's D (D = overlap between real X modelled niche-based distribution) ##

jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/clampSpsDistXschoenerXsampleXsps.jpeg', width=1200, height=1200)
par(mfrow=c(2,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.5, cex.axis=2, cex.main=2)
##
##spCD monotemporal
##
plot(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-3,1), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("CD species - SDM"["mono"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spCD multitemporal
##
plot(dataClamp[which(dataClamp$distUnderClamp>=0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>=0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-6.5,0), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("CD species - SDM"["multi"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp), pch=20, ylab="Schoener's D", col='black')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spHW monotemporal
##
plot(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-10,5), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("HW species - SDM"["mono"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spHW multitemporal
##
plot(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-10,0), ylim=c(0,1), pch=20, ylab="Schoener's D", xlab="log(Clamping in modeled distribution (in %))", main=expression("HW species - SDM"["multi"]), col='red')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp), pch=20, col='blue')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp), pch=20, col='black')
##
abline(lm(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$Schoeners_D_equiv ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
dev.off()



## clamping for America do Sul X clamping across modelled distributions (D = overlap between real X modelled niche-based distribution) ##

jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/clampSpsDistXclampAmSul.jpeg', width=1200, height=1200)
par(mfrow=c(2,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.5, cex.axis=2, cex.main=2)
##
##spCD monotemporal
##
plot(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-3.5,1), ylim=c(-2,3), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="log(South America area (in %))", main=expression("CD species - SDM"["mono"]), col='red')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["mono"]), col='blue')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["mono"]), col='black')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spCD multitemporal
##
plot(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-6.5,-0.5), ylim=c(-3,0.5), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="log(South America area (in %))", main=expression("CD species - SDM"["multi"]), col='red')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["multi"]), col='blue')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["multi"]), col='black')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spHW monotemporal
##
plot(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-3.5,1), ylim=c(-2,3), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="log(South America area (in %))", main=expression("HW species - SDM"["mono"]), col='red')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["mono"]), col='blue')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["mono"]), col='black')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
##spHW multitemporal
##
plot(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp), xlim=c(-6.5,-0.5), ylim=c(-3,0.5), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="log(South America area (in %))", main=expression("HW species - SDM"["multi"]), col='red')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 10),]$distUnderClamp)), lty=2, col='red')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["multi"]), col='blue')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 50),]$distUnderClamp)), lty=2, col='blue')
##
points(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp), ylim=c(0,1), pch=20, ylab="log(Clamping in modeled distribution (in %))", xlab="South America area (in %)", main=expression("SDM"["multi"]), col='black')
##
abline(lm(log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$clamping) ~ log(dataClamp[which(dataClamp$distUnderClamp>0 & dataClamp$sdmType == 'multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize == 100),]$distUnderClamp)), lty=2, col='black')
##
dev.off()



## clamping trends across geological time


jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/clampXtempo.jpeg', width=900)
par(mfrow=c(1,2))
plot(dataClamp[which(dataClamp$sdm=='multitemporal'),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal'),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='South América area (in %)', main='Multitemporal')
plot(dataClamp[which(dataClamp$sdm=='monotemporal'),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal'),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='South América area (in %)', main='Monotemporal')
dev.off()




## clamping trends across geological time, accounting for species, sample size e model calibration (mono and multitemporal)

## graphics for spCD
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/clampXtempoXsample_spCD.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.7, cex.axis=2, cex.main=2)
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$kyr), ylim=c(0,15), xlab='',  ylab='', main=expression('SDM'['mono']), col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$kyr), ylim=c(0,15), xlab='',  ylab='', main=expression('SDM'['multi']), col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$kyr), ylim=c(0,15), xlab='',  ylab='South América area (in %)', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$kyr), ylim=c(0,15), xlab='',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
dev.off()


## graphics for spHW
jpeg('/home/anderson/Projetos/Sps artificiais/graficos - resultados oficiais/clampXtempoXsample_spHW.jpeg', width=1200, height=1200)
par(mfrow=c(3,2), pch=1, mar=c(5,5,3,2), cex=1.5, cex.lab=1.7, cex.axis=2, cex.main=2)
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$kyr), ylim=c(0,15), xlab='',  ylab='', main=expression('SDM'['mono']), col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$kyr), ylim=c(0,15), xlab='',  ylab='', main=expression('SDM'['multi']), col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$kyr), ylim=c(0,15), xlab='',  ylab='South América area (in %)', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$kyr), ylim=c(0,15), xlab='',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
plot(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$clamping ~ as.factor(dataClamp[which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$kyr), ylim=c(0,15), xlab='Time (in kyr BP)',  ylab='', main='', col=rgb(0,0,0,alpha=0.5))
dev.off()




## boxplot clamping (full dataset and sps)
jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/BoxplotClamp.jpeg', width=1200)
par(mfrow=c(1,3), cex=1.3)
boxplot(log(dataClamp$clamping) ~ dataClamp$sdm, ylab='log(South America area (in %))', main='Full dataset')
##
boxplot(log(dataClamp[which(dataClamp$sp == 'spHW'),]$clamping) ~ dataClamp[which(dataClamp$sp == 'spHW'),]$sdm, ylab='log(South America area (in %))', main='HW species')
##
boxplot(log(dataClamp[which(dataClamp$sp == 'spCD'),]$clamping) ~ dataClamp[which(dataClamp$sp == 'spCD'),]$sdm, ylab='log(South America area (in %))', main='CD species')
dev.off()



##correlation between Schoener and clamping across sps' modelled distributions

##monotemporal - spCD
corSpCDmono10 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$distUnderClamp, method='pearson')

corSpCDmono50 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$distUnderClamp, method='pearson')

corSpCDmono100 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$distUnderClamp, method='pearson')

##monotemporal - spHW
corSpHWmono10 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$distUnderClamp, method='pearson')

corSpHWmono50 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$distUnderClamp, method='pearson')

corSpHWmono100 = cor.test(dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='monotemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$distUnderClamp, method='pearson')

##multitemporal - spCD
corSpCDmulti10 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==10),]$distUnderClamp, method='pearson')

corSpCDmulti50 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==50),]$distUnderClamp, method='pearson')

corSpCDmulti100 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spCD' & dataClamp$sampleSize==100),]$distUnderClamp, method='pearson')

##multitemporal - spHW
corSpHWmulti10 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==10),]$distUnderClamp, method='pearson')

corSpHWmulti50 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==50),]$distUnderClamp, method='pearson')

corSpHWmulti100 = cor.test(dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$Schoeners_D_simi,
                         dataClamp[ which(dataClamp$sdm=='multitemporal' & dataClamp$sp=='spHW' & dataClamp$sampleSize==100),]$distUnderClamp, method='pearson')



##table
corTable = data.frame( scenario = c('corFullDataset','corMulti','corMono','corSpCDmono10','corSpCDmono50','corSpCDmono100','corSpHWmono10','corSpHWmono50','corSpHWmono100','corSpCDmulti10','corSpCDmulti50','corSpCDmulti100','corSpHWmulti10','corSpHWmulti50','corSpHWmulti100'),
                      correlation = c(as.numeric(corFullDataset$estimate),as.numeric(corMulti$estimate),as.numeric(corMono$estimate),as.numeric(corSpCDmono10$estimate),as.numeric(corSpCDmono50$estimate),as.numeric(corSpCDmono100$estimate),as.numeric(corSpHWmono10$estimate),as.numeric(corSpHWmono50$estimate),as.numeric(corSpHWmono100$estimate),as.numeric(corSpCDmulti10$estimate),as.numeric(corSpCDmulti50$estimate),as.numeric(corSpCDmulti100$estimate),as.numeric(corSpHWmulti10$estimate),as.numeric(corSpHWmulti50$estimate),as.numeric(corSpHWmulti100$estimate)),
                      p.value = c(as.numeric(corFullDataset$p.value),as.numeric(corMulti$p.value),as.numeric(corMono$p.value),as.numeric(corSpCDmono10$p.value),as.numeric(corSpCDmono50$p.value),as.numeric(corSpCDmono100$p.value),as.numeric(corSpHWmono10$p.value),as.numeric(corSpHWmono50$p.value),as.numeric(corSpHWmono100$p.value),as.numeric(corSpCDmulti10$p.value),as.numeric(corSpCDmulti50$p.value),as.numeric(corSpCDmulti100$p.value),as.numeric(corSpHWmulti10$p.value),as.numeric(corSpHWmulti50$p.value),as.numeric(corSpHWmulti100$p.value)) )

corTable[,'p.value'] = round(corTable[,'p.value'], 3)

write.csv(corTable, paste(projectFolder,'correlationTable.csv'), row.names=FALSE)



### maps to compare climatic conditions between 0 and 22 kyrBP

temp0kyr = raster('/home/anderson/PosDoc/dados_ambientais/dados_projeto/000/bioclim_01.asc')
temp0kyr = mask(temp0kyr, AmSulShape)
preci0kyr = raster('/home/anderson/PosDoc/dados_ambientais/dados_projeto/000/bioclim_12.asc')
preci0kyr = mask(preci0kyr, AmSulShape)

temp22kyr = raster('/home/anderson/PosDoc/dados_ambientais/dados_projeto/022/bioclim_01.asc')
temp22kyr = mask(temp22kyr, AmSulShape)
preci22kyr = raster('/home/anderson/PosDoc/dados_ambientais/dados_projeto/022/bioclim_12.asc')
preci22kyr = mask(preci22kyr, AmSulShape)

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/temp0&22kyr.jpeg', width=800)
par(mfrow=c(1,2), mar=c(5,5,5,6))
plot(temp0kyr, main='0 kyr BP'); grid()
plot(temp22kyr, main='22 kyr BP'); grid()
dev.off()

jpeg('/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/preci0&22kyr.jpeg', width=800)
par(mfrow=c(1,2), mar=c(5,5,5,6))
plot(preci0kyr, main='0 kyr BP'); grid()
plot(preci22kyr, main='22 kyr BP'); grid()
dev.off()

clamp0 = mask(clampStack$multitemporal_proj_0kyr_spHW.sample50.replica1, AmSulShape)
clamp22 = mask(clampStack$multitemporal_proj_22kyr_spHW.sample50.replica1, AmSulShape)

par(mfrow=c(1,2))
plot(clamp0, main=c(paste('spHW ',0,'kyr BP',sep='')), col=c('lightgrey','red'), legend=FALSE)
plot(clamp22, main=c(paste('spHW ',22,'kyr BP',sep='')), col=c('lightgrey','red'), legend=FALSE)

dev.off()



##mapped distributions for present, interglacial and maximum glacial

library(raster)
library(maptools)
library(RColorBrewer)
AmSulShape = readShapePoly("/home/anderson/shapefiles/Am_Sul/borders.shp")

### MULTITEMPORAL ###

## spHW ##

## real distributions
HWcurrentReal = raster(paste(projectFolder,'NichoReal/spHW/000.asc',sep='')) > 0.2
HW22Real = raster(paste(projectFolder,'NichoReal/spHW/022.asc',sep='')) > 0.2

## multitemporal
HWModel_0kyrSample10 = raster(paste(projectFolder,'/maxent/multitemporal/spHW/spHW.sample10.replica1/proj_0kyr/proj_0kyr_spHW.sample10.replica1_TSSbin.grd',sep=''))
HWModel_0kyrSample50 = raster(paste(projectFolder,'/maxent/multitemporal/spHW/spHW.sample50.replica1/proj_0kyr/proj_0kyr_spHW.sample50.replica1_TSSbin.grd',sep=''))
HWModel_0kyrSample100 = raster(paste(projectFolder,'/maxent/multitemporal/spHW/spHW.sample100.replica1/proj_0kyr/proj_0kyr_spHW.sample100.replica1_TSSbin.grd',sep=''))
##
HWModel_22kyrSample10 = raster(paste(projectFolder,'/maxent/multitemporal/spHW/spHW.sample10.replica1/proj_22kyr/proj_22kyr_spHW.sample10.replica1_TSSbin.grd',sep=''))
HWModel_22kyrSample50 = raster(paste(projectFolder,'/maxent/multitemporal/spHW/spHW.sample50.replica1/proj_22kyr/proj_22kyr_spHW.sample50.replica1_TSSbin.grd',sep=''))
HWModel_22kyrSample100 = raster(paste(projectFolder,'/maxent/multitemporal/spHW/spHW.sample100.replica1/proj_22kyr/proj_22kyr_spHW.sample100.replica1_TSSbin.grd',sep=''))

## spCD ##

## real distributions
CDcurrentReal = raster(paste(projectFolder,'NichoReal/spCD/000.asc',sep='')) > 0.2
CD22Real = raster(paste(projectFolder,'NichoReal/spCD/022.asc',sep='')) > 0.2

## SDM multitemporal
CDModel_0kyrSample10 = raster(paste(projectFolder,'maxent/multitemporal/spCD/spCD.sample10.replica1/proj_0kyr/proj_0kyr_spCD.sample10.replica1_TSSbin.grd', sep=''))
CDModel_0kyrSample50 = raster(paste(projectFolder,'maxent/multitemporal/spCD/spCD.sample50.replica1/proj_0kyr/proj_0kyr_spCD.sample50.replica1_TSSbin.grd', sep=''))
CDModel_0kyrSample100 = raster(paste(projectFolder,'/maxent/multitemporal/spCD/spCD.sample100.replica1/proj_0kyr/proj_0kyr_spCD.sample100.replica1_TSSbin.grd', sep=''))
##
CDModel_22kyrSample10 = raster(paste(projectFolder,'/maxent/multitemporal/spCD/spCD.sample10.replica1/proj_22kyr/proj_22kyr_spCD.sample10.replica1_TSSbin.grd', sep='')) 
CDModel_22kyrSample50 = raster(paste(projectFolder,'/maxent/multitemporal/spCD/spCD.sample50.replica1/proj_22kyr/proj_22kyr_spCD.sample50.replica1_TSSbin.grd', sep=''))
CDModel_22kyrSample100 = raster(paste(projectFolder,'/maxent/multitemporal/spCD/spCD.sample100.replica1/proj_22kyr/proj_22kyr_spCD.sample100.replica1_TSSbin.grd', sep=''))

##overlaps spHW

jpeg(filename='/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/sobreposicoesHWmulti.jpg', width = 1400 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,20), mar=c(3,3,5,6))
plot(HWcurrentReal*1+HWModel_0kyrSample10*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample50*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample100*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
##legend
legend("topright",legend=c('Virtual species','Maxent projection','Overlap'),inset=c(-0.7,0),xpd=NA,pch=20,col=c('green','blue','dark green'),cex=2.5)
##
plot(HW22Real*1+HWModel_22kyrSample10*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample50*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample100*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
mtext('HW species',outer=TRUE,cex=4)
dev.off()

##overlap spCD

jpeg(filename='/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/sobreposicoesCDmulti.jpg', width = 1400 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,20), mar=c(3,3,5,6))
plot(CDcurrentReal*1+CDModel_0kyrSample10*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample50*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample100*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
##legend
legend("topright",legend=c('Virtual species','Maxent projection','Overlap'),inset=c(-0.7,0),xpd=NA,pch=20,col=c('green','blue','dark green'),cex=2.5)
##
plot(CD22Real*1+CDModel_22kyrSample10*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample50*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample100*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
mtext('CD species',outer=TRUE,cex=4)
dev.off()


### MONOTEMPORAL ###


## spHW ##

## real distributions
HWcurrentReal = raster(paste(projectFolder,'NichoReal/spHW/000.asc',sep='')) > 0.2
HW22Real = raster(paste(projectFolder,'NichoReal/spHW/022.asc',sep='')) > 0.2

## monotemporal
HWModel_0kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample10.replica1/proj_0kyr/proj_0kyr_spHW.sample10.replica1_TSSbin.grd')
HWModel_0kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample50.replica1/proj_0kyr/proj_0kyr_spHW.sample50.replica1_TSSbin.grd')
HWModel_0kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample100.replica1/proj_0kyr/proj_0kyr_spHW.sample100.replica1_TSSbin.grd')
##
HWModel_22kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample10.replica1/proj_22kyr/proj_22kyr_spHW.sample10.replica1_TSSbin.grd') 
HWModel_22kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample50.replica1/proj_22kyr/proj_22kyr_spHW.sample50.replica1_TSSbin.grd')
HWModel_22kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spHW/spHW.sample100.replica1/proj_22kyr/proj_22kyr_spHW.sample100.replica1_TSSbin.grd')

## spCD ##

## real distributions
CDcurrentReal = raster(paste(projectFolder,'NichoReal/spCD/000.asc',sep='')) > 0.2
CD22Real = raster(paste(projectFolder,'NichoReal/spCD/022.asc',sep='')) > 0.2

## SDM monotemporal
CDModel_0kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample10.replica1/proj_0kyr/proj_0kyr_spCD.sample10.replica1_TSSbin.grd')
CDModel_0kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample50.replica1/proj_0kyr/proj_0kyr_spCD.sample50.replica1_TSSbin.grd')
CDModel_0kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample100.replica1/proj_0kyr/proj_0kyr_spCD.sample100.replica1_TSSbin.grd')
##
CDModel_22kyrSample10 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample10.replica1/proj_22kyr/proj_22kyr_spCD.sample10.replica1_TSSbin.grd') 
CDModel_22kyrSample50 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample50.replica1/proj_22kyr/proj_22kyr_spCD.sample50.replica1_TSSbin.grd')
CDModel_22kyrSample100 = raster('/home/anderson/Documentos/Projetos/Sps artificiais/maxent/monotemporal/spCD/spCD.sample100.replica1/proj_22kyr/proj_22kyr_spCD.sample100.replica1_TSSbin.grd')

##overlaps spHW
jpeg(filename='/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/sobreposicoesHWmono.jpg', width = 1400 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,20), mar=c(3,3,5,6))
plot(HWcurrentReal*1+HWModel_0kyrSample10*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample50*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(HWcurrentReal*1+HWModel_0kyrSample100*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
##legend
legend("topright",legend=c('Virtual species','Maxent projection','Overlap'),inset=c(-0.7,0),xpd=NA,pch=20,col=c('green','blue','dark green'),cex=2.5)
##
plot(HW22Real*1+HWModel_22kyrSample10*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample50*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(HW22Real*1+HWModel_22kyrSample100*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
mtext('HW species.',outer=TRUE,cex=4)
dev.off()

##overlap spCD
jpeg(filename='/home/anderson/Documentos/Projetos/Sps artificiais/graficos - resultados oficiais/sobreposicoesCDmono.jpg', width = 1400 , height = 1100) 
par(mfrow=c(2,3),oma=c(0,0,5,20), mar=c(3,3,5,6))
plot(CDcurrentReal*1+CDModel_0kyrSample10*2,main='(A)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample50*2,main='(B)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(CDcurrentReal*1+CDModel_0kyrSample100*2,main='(C)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 0 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
##legend
legend("topright",legend=c('Virtual species','Maxent projection','Overlap'),inset=c(-0.7,0),xpd=NA,pch=20,col=c('green','blue','dark green'),cex=2.5)
##
plot(CD22Real*1+CDModel_22kyrSample10*2,main='(D)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=2,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 10 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample50*2,main='(E)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 50 records',cex=2)
grid()
plot(CD22Real*1+CDModel_22kyrSample100*2,main='(F)',col=c('white','green','blue','darkgreen'),legend=FALSE,cex.axis=1.7,cex.main=4)
plot(AmSulShape,add=TRUE)
text(-50,-45,'Time: 22 kyr BP',cex=2)
text(-50,-50,'Sample: 100 records',cex=2)
grid()
mtext('CD species',outer=TRUE,cex=4)
dev.off()

