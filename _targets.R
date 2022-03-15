library(targets)
library(tarchetypes)
library(tidyverse)
source("R/functionsDAG.R")
source("R/functionsSignal.R")
source("R/functionsReview.R")
source("R/functionsSimulation.R")
source("R/functionsReplications.R")
options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("brms", "cowplot", "conleyreg", "countrycode", 
                            "dagitty", "geosphere", "ggdag", "ggrepel", "ggtext", 
                            "haven", "lmtest", "MASS", "papaja", "psych", "readxl", 
                            "rstan", "sjlabelled", "tidybayes", "tidyverse"))

# targets for simulation (see below)
simulationTargets <-
  tar_map(
    # return nested list from tar_map()
    unlist = FALSE,
    # map simulation over different covariance matrices and different values of lambda and rho
    values = expand_grid(covMat = rlang::syms(c("simGeoCov", "simLinCov")),
                         lambda = c(0.2, 0.5, 0.8), rho = c(0.2, 0.5, 0.8)),
    # simulate data
    tar_target(simData, simulateData(covMat = covMat, continent, iso, langFam, 
                                     seed = seed, lambda = lambda, rho = rho), 
               pattern = map(seed), iteration = "list"),
    # ols analyses
    tar_target(olsModel1, fitOLSModel(y ~ x,              data = simData), pattern = map(simData)),
    tar_target(olsModel2, fitOLSModel(y ~ x + latitude,   data = simData), pattern = map(simData)),
    tar_target(olsModel3, fitOLSModel(y ~ x + longitude,  data = simData), pattern = map(simData)),
    tar_target(olsModel4, fitOLSModel(y ~ x + continent,  data = simData), pattern = map(simData)),
    tar_target(olsModel5, fitOLSModel(y ~ x + langFamily, data = simData), pattern = map(simData)),
    # conley se analyses
    tar_target(conleyModel1, fitConleyModel(simData, dist_cutoff = 100), pattern = map(simData)),
    tar_target(conleyModel2, fitConleyModel(simData, dist_cutoff = 1000), pattern = map(simData)),
    tar_target(conleyModel3, fitConleyModel(simData, dist_cutoff = 10000), pattern = map(simData)),
    # brms analyses
    tar_target(brmsModel1, fitBrmsModel(brmsInitial1, simData), pattern = map(simData)),
    tar_target(brmsModel2, fitBrmsModel(brmsInitial2, simData), pattern = map(simData)),
    tar_target(brmsModel3, fitBrmsModel(brmsInitial3, simData), pattern = map(simData))
    )

# pipeline
list(
  
  #### Causal model ####
  
  tar_target(dag, plotDAG()),
  
  #### Geographic and linguistic signal ####
  
  # files
  tar_target(fileGeo, "data/networks/1F Population Distance.xlsx", format = "file"),
  tar_target(fileLin, "data/networks/2F Country Distance 1pml adj.xlsx", format = "file"),
  tar_target(fileHDI, "data/hdi/2020_statistical_annex_all.xlsx", format = "file"),
  tar_target(fileISOHDI, "data/hdi/iso.csv", format = "file"),
  tar_target(fileWVS, "data/wvs/Integrated_values_surveys_1981-2021.sav", format = "file"),
  # hdi data
  tar_target(hdi, loadHDIData(fileHDI, fileISOHDI)),
  # wvs data
  tar_target(wvs, read_sav(fileWVS)),
  # covariance matrices
  tar_target(geoCov, loadCovMat(fileGeo, log = TRUE)),
  tar_target(linCov, loadCovMat(fileLin, log = FALSE)),
  # geographic and linguistic signal
  tar_target(hdiSignal, fitHDISignal(hdi, geoCov, linCov)),
  tar_target(tradSignal, fitWVSSignal(wvs, outcome = "trad", geoCov, linCov)),
  tar_target(survSignal, fitWVSSignal(wvs, outcome = "surv", geoCov, linCov)),
  # posterior samples
  tar_target(postHDI,  as_draws_array(hdiSignal , variable = "^sd_", regex = TRUE)),
  tar_target(postTrad, as_draws_array(tradSignal, variable = "^sd_", regex = TRUE)),
  tar_target(postSurv, as_draws_array(survSignal, variable = "^sd_", regex = TRUE)),
  # calculate signal
  tar_target(hdiGeo,  hypothesis(hdiSignal,  "(sd_isoGeo__Intercept^2/(sd_isoGeo__Intercept^2+sd_isoLin__Intercept^2+sd_iso__Intercept^2))=0",  class = NULL)),
  tar_target(hdiLin,  hypothesis(hdiSignal,  "(sd_isoLin__Intercept^2/(sd_isoGeo__Intercept^2+sd_isoLin__Intercept^2+sd_iso__Intercept^2))=0",  class = NULL)),
  tar_target(tradGeo, hypothesis(tradSignal, "(sd_isoGeo__Intercept^2/(sd_isoGeo__Intercept^2+sd_isoLin__Intercept^2+sd_S009__Intercept^2))=0", class = NULL)),
  tar_target(tradLin, hypothesis(tradSignal, "(sd_isoLin__Intercept^2/(sd_isoGeo__Intercept^2+sd_isoLin__Intercept^2+sd_S009__Intercept^2))=0", class = NULL)),
  tar_target(survGeo, hypothesis(survSignal, "(sd_isoGeo__Intercept^2/(sd_isoGeo__Intercept^2+sd_isoLin__Intercept^2+sd_S009__Intercept^2))=0", class = NULL)),
  tar_target(survLin, hypothesis(survSignal, "(sd_isoLin__Intercept^2/(sd_isoGeo__Intercept^2+sd_isoLin__Intercept^2+sd_S009__Intercept^2))=0", class = NULL)),
  # plot signal
  tar_target(plotSignal, plotGeoLinSignal(postHDI, postTrad, postSurv)),
  
  #### Review ####
  
  # files
  tar_target(fileReview, "review/mainReview.xlsx", format = "file"),
  tar_target(fileMethods, "review/methods.csv", format = "file"),
  # load review
  tar_target(review, loadReview(fileReview, fileMethods)),
  # fit random effects models to retrieve corrected analysis-level proportions
  tar_target(rM1, fitReviewMLM(review, outcome = "ControlNI")),
  tar_target(rM2, fitReviewMLM(review, outcome = "method_RegionalFixedEffects")),
  tar_target(rM3, fitReviewMLM(review, outcome = "method_Distance")),
  tar_target(rM4, fitReviewMLM(review, outcome = "method_SharedCulturalHistory")),
  tar_target(rM5, fitReviewMLM(review, outcome = "method_Other")),
  # posterior samples
  tar_target(postRM1, as_draws_array(rM1, variable = "^b_", regex = TRUE)),
  tar_target(postRM2, as_draws_array(rM2, variable = "^b_", regex = TRUE)),
  tar_target(postRM3, as_draws_array(rM3, variable = "^b_", regex = TRUE)),
  tar_target(postRM4, as_draws_array(rM4, variable = "^b_", regex = TRUE)),
  tar_target(postRM5, as_draws_array(rM5, variable = "^b_", regex = TRUE)),
  # fit spline model
  tar_target(spline, fitSpline(review)),
  # plot review summary
  tar_target(plotReview, plotReviewSummary(review, spline, postRM1, postRM2, 
                                           postRM3, postRM4, postRM5)),
  
  #### Simulation ####
  
  # files
  tar_target(fileContinent, "data/countryData/continent.csv", format = "file"),
  tar_target(fileLangFam, "data/countryData/langFamily.xlsx", format = "file"),
  tar_target(fileISO, "data/countryData/countries_codes_and_coordinates.csv", format = "file"),
  # continent data (https://datahub.io/JohnSnowLabs/country-and-continent-codes-list#data)
  tar_target(continent, read.csv(fileContinent, na.strings = "")),
  # language family data (from glottolog)
  tar_target(langFam, read_xlsx(fileLangFam)),
  # isocodes (https://gist.github.com/tadast/8827699)
  tar_target(iso, loadISO(fileISO)),
  # geographic and linguistic covariance matrices
  tar_target(simGeoCov, loadSimCovMat(fileGeo, log = TRUE , continent, iso, langFam)),
  tar_target(simLinCov, loadSimCovMat(fileLin, log = FALSE, continent, iso, langFam)),
  # initialise brms models
  tar_target(brmsInitial1, setupBrms(simData_simLinCov_0.8_0.8[[1]], simLinCov, type = "spatial")),
  tar_target(brmsInitial2, setupBrms(simData_simLinCov_0.8_0.8[[1]], simLinCov, type = "linguistic")),
  tar_target(brmsInitial3, setupBrms(simData_simLinCov_0.8_0.8[[1]], simLinCov, type = "both")),
  # seed for simulation
  tar_target(seed, 1:2),
  # simulation
  simulationTargets,
  # combine
  tar_combine(olsModel1_simGeoCov, simulationTargets$olsModel1[1:9]  , command = dplyr::bind_rows(!!!.x)),
  tar_combine(olsModel1_simLinCov, simulationTargets$olsModel1[10:18], command = dplyr::bind_rows(!!!.x)),
  tar_combine(olsModel2_simGeoCov, simulationTargets$olsModel2[1:9]  , command = dplyr::bind_rows(!!!.x)),
  tar_combine(olsModel2_simLinCov, simulationTargets$olsModel2[10:18], command = dplyr::bind_rows(!!!.x)),
  tar_combine(olsModel3_simGeoCov, simulationTargets$olsModel3[1:9]  , command = dplyr::bind_rows(!!!.x)),
  tar_combine(olsModel3_simLinCov, simulationTargets$olsModel3[10:18], command = dplyr::bind_rows(!!!.x)),
  tar_combine(olsModel4_simGeoCov, simulationTargets$olsModel4[1:9]  , command = dplyr::bind_rows(!!!.x)),
  tar_combine(olsModel4_simLinCov, simulationTargets$olsModel4[10:18], command = dplyr::bind_rows(!!!.x)),
  tar_combine(olsModel5_simGeoCov, simulationTargets$olsModel5[1:9]  , command = dplyr::bind_rows(!!!.x)),
  tar_combine(olsModel5_simLinCov, simulationTargets$olsModel5[10:18], command = dplyr::bind_rows(!!!.x)),
  tar_combine(conleyModel1_simGeoCov, simulationTargets$conleyModel1[1:9]  , command = dplyr::bind_rows(!!!.x)),
  tar_combine(conleyModel1_simLinCov, simulationTargets$conleyModel1[10:18], command = dplyr::bind_rows(!!!.x)),
  tar_combine(conleyModel2_simGeoCov, simulationTargets$conleyModel2[1:9]  , command = dplyr::bind_rows(!!!.x)),
  tar_combine(conleyModel2_simLinCov, simulationTargets$conleyModel2[10:18], command = dplyr::bind_rows(!!!.x)),
  tar_combine(conleyModel3_simGeoCov, simulationTargets$conleyModel3[1:9]  , command = dplyr::bind_rows(!!!.x)),
  tar_combine(conleyModel3_simLinCov, simulationTargets$conleyModel3[10:18], command = dplyr::bind_rows(!!!.x)),
  tar_combine(brmsModel1_simGeoCov, simulationTargets$brmsModel1[1:9]  , command = dplyr::bind_rows(!!!.x)),
  tar_combine(brmsModel1_simLinCov, simulationTargets$brmsModel1[10:18], command = dplyr::bind_rows(!!!.x)),
  tar_combine(brmsModel2_simGeoCov, simulationTargets$brmsModel2[1:9]  , command = dplyr::bind_rows(!!!.x)),
  tar_combine(brmsModel2_simLinCov, simulationTargets$brmsModel2[10:18], command = dplyr::bind_rows(!!!.x)),
  tar_combine(brmsModel3_simGeoCov, simulationTargets$brmsModel3[1:9]  , command = dplyr::bind_rows(!!!.x)),
  tar_combine(brmsModel3_simLinCov, simulationTargets$brmsModel3[10:18], command = dplyr::bind_rows(!!!.x)),
  # plot simulation results
  tar_target(plotSim1, plotSimInd(olsModel1_simGeoCov, olsModel2_simGeoCov, olsModel3_simGeoCov, olsModel4_simGeoCov, 
                                  olsModel5_simGeoCov, conleyModel1_simGeoCov, conleyModel2_simGeoCov, conleyModel3_simGeoCov,
                                  brmsModel1_simGeoCov, brmsModel2_simGeoCov, brmsModel3_simGeoCov, 
                                  type = "spatial", file = "figures/simIndGeo.pdf")),
  tar_target(plotSim2, plotSimInd(olsModel1_simLinCov, olsModel2_simLinCov, olsModel3_simLinCov, olsModel4_simLinCov, 
                                  olsModel5_simLinCov, conleyModel1_simLinCov, conleyModel2_simLinCov, conleyModel3_simLinCov,
                                  brmsModel1_simLinCov, brmsModel2_simLinCov, brmsModel3_simLinCov, 
                                  type = "cultural phylogenetic", file = "figures/simIndLin.pdf")),
  tar_target(plotSim3, plotSimAll(olsModel1_simGeoCov, olsModel2_simGeoCov, olsModel3_simGeoCov, olsModel4_simGeoCov, 
                                  olsModel5_simGeoCov, conleyModel1_simGeoCov, conleyModel2_simGeoCov, conleyModel3_simGeoCov,
                                  brmsModel1_simGeoCov, brmsModel2_simGeoCov, brmsModel3_simGeoCov, 
                                  type = "spatial", file = "figures/simAllGeo.pdf")),
  tar_target(plotSim4, plotSimAll(olsModel1_simLinCov, olsModel2_simLinCov, olsModel3_simLinCov, olsModel4_simLinCov, 
                                  olsModel5_simLinCov, conleyModel1_simLinCov, conleyModel2_simLinCov, conleyModel3_simLinCov,
                                  brmsModel1_simLinCov, brmsModel2_simLinCov, brmsModel3_simLinCov, 
                                  type = "cultural phylogenetic", file = "figures/simAllLin.pdf")),
  
  #### Replications ####
  
  # files
  tar_target(fileAlesina,     "data/replications/Alesina2013/tradPloughUse.xlsx", format = "file"),
  tar_target(fileBeck1,       "data/replications/Beck2003/AppendixTable10.xlsx", format = "file"),
  tar_target(fileBeck2,       "data/replications/Beck2005/Table1.xlsx", format = "file"),
  tar_target(fileBockstette1, "data/replications/Bockstette2002/AppendixTable1.xlsx", format = "file"),
  tar_target(fileBockstette2, "data/replications/Bockstette2002/Data_Extract_From_World_Development_Indicators.xlsx", format = "file"),
  tar_target(fileEasterly1,   "data/replications/Easterly2003/Acemoglu et al. 2001 Appendix Table A2.xlsx", format = "file"),
  tar_target(fileEasterly2,   "data/replications/Easterly2003/wgidataset.xlsx", format = "file"),
  tar_target(fileEasterly3,   "data/replications/Easterly2007/AppendixA.xlsx", format = "file"),
  tar_target(fileEasterly4,   "data/replications/Easterly2007/WIID2C.xls", format = "file"),
  tar_target(fileFincher,     "data/replications/Fincher2008/suppData.xlsx", format = "file"),
  tar_target(fileGelfand1,    "data/replications/Gelfand2011/pathogens.xlsx", format = "file"),
  tar_target(fileGelfand2,    "data/replications/Gelfand2011/tightnessTable1.xlsx", format = "file"),
  tar_target(fileGelfand3,    "data/replications/Gelfand2011/2005-esi-all-countries.xls", format = "file"),
  tar_target(fileGelfand4,    "data/replications/Gelfand2011/WorldBank_API_NY.GNP.PCAP.CD_DS2_en_csv_v2_3358861.csv", format = "file"),
  tar_target(fileGelfand5,    "data/replications/Gelfand2011/tuberculosis.csv", format = "file"),
  tar_target(fileInglehart,   "data/replications/Inglehart2000/API_SL.IND.EMPL.ZS_DS2_en_csv_v2_3478718.csv", format = "file"),
  tar_target(fileKnack,       "data/replications/Knack1997/DataAppendix.xlsx", format = "file"),
  tar_target(fileSkidmore1,   "data/replications/Skidmore2002/emdat_public_2021_11_30_query_uid-kUNKpe.xlsx", format = "file"),
  tar_target(fileSkidmore2,   "data/replications/Skidmore2002/pwt56_forweb.xls", format = "file"),
  tar_target(fileSkidmore3,   "data/replications/Skidmore2002/link.csv", format = "file"),
  tar_target(fileSkidmore4,   "data/replications/Skidmore2002/API_AG.LND.TOTL.K2_DS2_en_csv_v2_3359379.csv", format = "file"),
  # Adamczyk 2009
  tar_target(adamczykData, loadDataAdamczyk2009(wvs, iso)),
  tar_target(adamczykM1a, fitModelAdamczyk2009(adamczykData, linCov, control = "none")),
  tar_target(adamczykM1b, fitModelAdamczyk2009(adamczykData, linCov, control = "spatial")),
  tar_target(adamczykM1c, fitModelAdamczyk2009(adamczykData, linCov, control = "cultural")),
  tar_target(adamczykM1d, fitModelAdamczyk2009(adamczykData, linCov, control = "both")),
  tar_target(adamczykSlope1a, posterior_samples(adamczykM1a, pars = "b_countrySurv")[,1]),
  tar_target(adamczykSlope1b, posterior_samples(adamczykM1b, pars = "b_countrySurv")[,1]),
  tar_target(adamczykSlope1c, posterior_samples(adamczykM1c, pars = "b_countrySurv")[,1]),
  tar_target(adamczykSlope1d, posterior_samples(adamczykM1d, pars = "b_countrySurv")[,1]),
  tar_target(adamczykCond1a, conditional_effects(adamczykM1a)[["countrySurv"]]),
  tar_target(adamczykCond1b, conditional_effects(adamczykM1b)[["countrySurv"]]),
  tar_target(adamczykCond1c, conditional_effects(adamczykM1c)[["countrySurv"]]),
  tar_target(adamczykCond1d, conditional_effects(adamczykM1d)[["countrySurv"]]),
  # Alesina 2013
  tar_target(alesinaData, loadDataAlesina2013(fileAlesina, iso)),
  tar_target(alesinaM1a, fitModelAlesina2013(alesinaData, linCov, control = "none")),
  tar_target(alesinaM1b, fitModelAlesina2013(alesinaData, linCov, control = "spatial")),
  tar_target(alesinaM1c, fitModelAlesina2013(alesinaData, linCov, control = "cultural")),
  tar_target(alesinaM1d, fitModelAlesina2013(alesinaData, linCov, control = "both")),
  tar_target(alesinaSlope1a, posterior_samples(alesinaM1a, pars = "b_tradPloughMean")[,1]),
  tar_target(alesinaSlope1b, posterior_samples(alesinaM1b, pars = "b_tradPloughMean")[,1]),
  tar_target(alesinaSlope1c, posterior_samples(alesinaM1c, pars = "b_tradPloughMean")[,1]),
  tar_target(alesinaSlope1d, posterior_samples(alesinaM1d, pars = "b_tradPloughMean")[,1]),
  tar_target(alesinaCond1a, conditional_effects(alesinaM1a)[["tradPloughMean"]]),
  tar_target(alesinaCond1b, conditional_effects(alesinaM1b)[["tradPloughMean"]]),
  tar_target(alesinaCond1c, conditional_effects(alesinaM1c)[["tradPloughMean"]]),
  tar_target(alesinaCond1d, conditional_effects(alesinaM1d)[["tradPloughMean"]]),
  # Beck 2003
  tar_target(beckData1, loadDataBeck2003(fileBeck1, iso)),
  tar_target(beckM1a, fitModelBeck2003(beckData1, linCov, control = "none")),
  tar_target(beckM1b, fitModelBeck2003(beckData1, linCov, control = "spatial")),
  tar_target(beckM1c, fitModelBeck2003(beckData1, linCov, control = "cultural")),
  tar_target(beckM1d, fitModelBeck2003(beckData1, linCov, control = "both")),
  tar_target(beckSlope1a, posterior_samples(beckM1a, pars = "b_logSettlerMortality")[,1]),
  tar_target(beckSlope1b, posterior_samples(beckM1b, pars = "b_logSettlerMortality")[,1]),
  tar_target(beckSlope1c, posterior_samples(beckM1c, pars = "b_logSettlerMortality")[,1]),
  tar_target(beckSlope1d, posterior_samples(beckM1d, pars = "b_logSettlerMortality")[,1]),
  tar_target(beckCond1a, conditional_effects(beckM1a)[["logSettlerMortality"]]),
  tar_target(beckCond1b, conditional_effects(beckM1b)[["logSettlerMortality"]]),
  tar_target(beckCond1c, conditional_effects(beckM1c)[["logSettlerMortality"]]),
  tar_target(beckCond1d, conditional_effects(beckM1d)[["logSettlerMortality"]]),
  # Beck 2005
  tar_target(beckData2, loadDataBeck2005(fileBeck2, iso)),
  tar_target(beckM2a, fitModelBeck2005(beckData2, linCov, control = "none")),
  tar_target(beckM2b, fitModelBeck2005(beckData2, linCov, control = "spatial")),
  tar_target(beckM2c, fitModelBeck2005(beckData2, linCov, control = "cultural")),
  tar_target(beckM2d, fitModelBeck2005(beckData2, linCov, control = "both")),
  tar_target(beckSlope2a, posterior_samples(beckM2a, pars = "b_SME250")[,1]),
  tar_target(beckSlope2b, posterior_samples(beckM2b, pars = "b_SME250")[,1]),
  tar_target(beckSlope2c, posterior_samples(beckM2c, pars = "b_SME250")[,1]),
  tar_target(beckSlope2d, posterior_samples(beckM2d, pars = "b_SME250")[,1]),
  tar_target(beckCond2a, conditional_effects(beckM2a)[["SME250"]]),
  tar_target(beckCond2b, conditional_effects(beckM2b)[["SME250"]]),
  tar_target(beckCond2c, conditional_effects(beckM2c)[["SME250"]]),
  tar_target(beckCond2d, conditional_effects(beckM2d)[["SME250"]]),
  # Bockstette 2002
  tar_target(bockstetteData, loadDataBockstette2002(fileBockstette1, fileBockstette2, iso)),
  tar_target(bockstetteM1a, fitModelBockstette2002(bockstetteData, linCov, control = "none")),
  tar_target(bockstetteM1b, fitModelBockstette2002(bockstetteData, linCov, control = "spatial")),
  tar_target(bockstetteM1c, fitModelBockstette2002(bockstetteData, linCov, control = "cultural")),
  tar_target(bockstetteM1d, fitModelBockstette2002(bockstetteData, linCov, control = "both")),
  tar_target(bockstetteSlope1a, posterior_samples(bockstetteM1a, pars = "b_Statehist5")[,1]),
  tar_target(bockstetteSlope1b, posterior_samples(bockstetteM1b, pars = "b_Statehist5")[,1]),
  tar_target(bockstetteSlope1c, posterior_samples(bockstetteM1c, pars = "b_Statehist5")[,1]),
  tar_target(bockstetteSlope1d, posterior_samples(bockstetteM1d, pars = "b_Statehist5")[,1]),
  tar_target(bockstetteCond1a, conditional_effects(bockstetteM1a)[["Statehist5"]]),
  tar_target(bockstetteCond1b, conditional_effects(bockstetteM1b)[["Statehist5"]]),
  tar_target(bockstetteCond1c, conditional_effects(bockstetteM1c)[["Statehist5"]]),
  tar_target(bockstetteCond1d, conditional_effects(bockstetteM1d)[["Statehist5"]]),
  # Easterly 2003
  tar_target(easterlyData1, loadDataEasterly2003(fileEasterly1, fileEasterly2, iso)),
  tar_target(easterlyM1a, fitModelEasterly2003(easterlyData1, linCov, control = "none")),
  tar_target(easterlyM1b, fitModelEasterly2003(easterlyData1, linCov, control = "spatial")),
  tar_target(easterlyM1c, fitModelEasterly2003(easterlyData1, linCov, control = "cultural")),
  tar_target(easterlyM1d, fitModelEasterly2003(easterlyData1, linCov, control = "both")),
  tar_target(easterlySlope1a, posterior_samples(easterlyM1a, pars = "b_institIndex")[,1]),
  tar_target(easterlySlope1b, posterior_samples(easterlyM1b, pars = "b_institIndex")[,1]),
  tar_target(easterlySlope1c, posterior_samples(easterlyM1c, pars = "b_institIndex")[,1]),
  tar_target(easterlySlope1d, posterior_samples(easterlyM1d, pars = "b_institIndex")[,1]),
  tar_target(easterlyCond1a, conditional_effects(easterlyM1a)[["institIndex"]]),
  tar_target(easterlyCond1b, conditional_effects(easterlyM1b)[["institIndex"]]),
  tar_target(easterlyCond1c, conditional_effects(easterlyM1c)[["institIndex"]]),
  tar_target(easterlyCond1d, conditional_effects(easterlyM1d)[["institIndex"]]),
  # Easterly 2007
  tar_target(easterlyData2, loadDataEasterly2007(fileEasterly3, fileEasterly4, iso)),
  tar_target(easterlyM2a, fitModelEasterly2007(easterlyData2, linCov, control = "none")),
  tar_target(easterlyM2b, fitModelEasterly2007(easterlyData2, linCov, control = "spatial")),
  tar_target(easterlyM2c, fitModelEasterly2007(easterlyData2, linCov, control = "cultural")),
  tar_target(easterlyM2d, fitModelEasterly2007(easterlyData2, linCov, control = "both")),
  tar_target(easterlySlope2a, posterior_samples(easterlyM2a, pars = "b_LWHEATSUGAR")[,1]),
  tar_target(easterlySlope2b, posterior_samples(easterlyM2b, pars = "b_LWHEATSUGAR")[,1]),
  tar_target(easterlySlope2c, posterior_samples(easterlyM2c, pars = "b_LWHEATSUGAR")[,1]),
  tar_target(easterlySlope2d, posterior_samples(easterlyM2d, pars = "b_LWHEATSUGAR")[,1]),
  tar_target(easterlyCond2a, conditional_effects(easterlyM2a)[["LWHEATSUGAR"]]),
  tar_target(easterlyCond2b, conditional_effects(easterlyM2b)[["LWHEATSUGAR"]]),
  tar_target(easterlyCond2c, conditional_effects(easterlyM2c)[["LWHEATSUGAR"]]),
  tar_target(easterlyCond2d, conditional_effects(easterlyM2d)[["LWHEATSUGAR"]]),
  # Fincher 2008
  tar_target(fincherData, loadDataFincher2008(fileFincher, iso)),
  tar_target(fincherM1a, fitModelFincher2008(fincherData, linCov, control = "none")),
  tar_target(fincherM1b, fitModelFincher2008(fincherData, linCov, control = "spatial")),
  tar_target(fincherM1c, fitModelFincher2008(fincherData, linCov, control = "cultural")),
  tar_target(fincherM1d, fitModelFincher2008(fincherData, linCov, control = "both")),
  tar_target(fincherSlope1a, posterior_samples(fincherM1a, pars = "b_pathPrevHistorical")[,1]),
  tar_target(fincherSlope1b, posterior_samples(fincherM1b, pars = "b_pathPrevHistorical")[,1]),
  tar_target(fincherSlope1c, posterior_samples(fincherM1c, pars = "b_pathPrevHistorical")[,1]),
  tar_target(fincherSlope1d, posterior_samples(fincherM1d, pars = "b_pathPrevHistorical")[,1]),
  tar_target(fincherCond1a, conditional_effects(fincherM1a)[["pathPrevHistorical"]]),
  tar_target(fincherCond1b, conditional_effects(fincherM1b)[["pathPrevHistorical"]]),
  tar_target(fincherCond1c, conditional_effects(fincherM1c)[["pathPrevHistorical"]]),
  tar_target(fincherCond1d, conditional_effects(fincherM1d)[["pathPrevHistorical"]]),
  # Gelfand 2011
  tar_target(gelfandData, loadDataGelfand2011(fileGelfand1, fileGelfand2, fileGelfand3, fileGelfand4, fileGelfand5, iso)),
  tar_target(gelfandM1a, fitModelGelfand2011(gelfandData, linCov, control = "none")),
  tar_target(gelfandM1b, fitModelGelfand2011(gelfandData, linCov, control = "spatial")),
  tar_target(gelfandM1c, fitModelGelfand2011(gelfandData, linCov, control = "cultural")),
  tar_target(gelfandM1d, fitModelGelfand2011(gelfandData, linCov, control = "both")),
  tar_target(gelfandSlope1a, posterior_samples(gelfandM1a, pars = "b_DISCAS")[,1]),
  tar_target(gelfandSlope1b, posterior_samples(gelfandM1b, pars = "b_DISCAS")[,1]),
  tar_target(gelfandSlope1c, posterior_samples(gelfandM1c, pars = "b_DISCAS")[,1]),
  tar_target(gelfandSlope1d, posterior_samples(gelfandM1d, pars = "b_DISCAS")[,1]),
  tar_target(gelfandCond1a, conditional_effects(gelfandM1a)[["DISCAS"]]),
  tar_target(gelfandCond1b, conditional_effects(gelfandM1b)[["DISCAS"]]),
  tar_target(gelfandCond1c, conditional_effects(gelfandM1c)[["DISCAS"]]),
  tar_target(gelfandCond1d, conditional_effects(gelfandM1d)[["DISCAS"]]),
  # Inglehart 2000
  tar_target(inglehartData, loadDataInglehart2000(fileInglehart, wvs, iso)),
  tar_target(inglehartM1a, fitModelInglehart2000(inglehartData, linCov, control = "none")),
  tar_target(inglehartM1b, fitModelInglehart2000(inglehartData, linCov, control = "spatial")),
  tar_target(inglehartM1c, fitModelInglehart2000(inglehartData, linCov, control = "cultural")),
  tar_target(inglehartM1d, fitModelInglehart2000(inglehartData, linCov, control = "both")),
  tar_target(inglehartSlope1a, posterior_samples(inglehartM1a, pars = "b_indust")[,1]),
  tar_target(inglehartSlope1b, posterior_samples(inglehartM1b, pars = "b_indust")[,1]),
  tar_target(inglehartSlope1c, posterior_samples(inglehartM1c, pars = "b_indust")[,1]),
  tar_target(inglehartSlope1d, posterior_samples(inglehartM1d, pars = "b_indust")[,1]),
  tar_target(inglehartCond1a, conditional_effects(inglehartM1a)[["indust"]]),
  tar_target(inglehartCond1b, conditional_effects(inglehartM1b)[["indust"]]),
  tar_target(inglehartCond1c, conditional_effects(inglehartM1c)[["indust"]]),
  tar_target(inglehartCond1d, conditional_effects(inglehartM1d)[["indust"]]),
  # Knack 1997
  tar_target(knackData, loadDataKnack1997(fileKnack, iso)),
  tar_target(knackM1a, fitModelKnack1997(knackData, linCov, control = "none")),
  tar_target(knackM1b, fitModelKnack1997(knackData, linCov, control = "spatial")),
  tar_target(knackM1c, fitModelKnack1997(knackData, linCov, control = "cultural")),
  tar_target(knackM1d, fitModelKnack1997(knackData, linCov, control = "both")),
  tar_target(knackSlope1a, posterior_samples(knackM1a, pars = "b_Trust")[,1]),
  tar_target(knackSlope1b, posterior_samples(knackM1b, pars = "b_Trust")[,1]),
  tar_target(knackSlope1c, posterior_samples(knackM1c, pars = "b_Trust")[,1]),
  tar_target(knackSlope1d, posterior_samples(knackM1d, pars = "b_Trust")[,1]),
  tar_target(knackCond1a, conditional_effects(knackM1a)[["Trust"]]),
  tar_target(knackCond1b, conditional_effects(knackM1b)[["Trust"]]),
  tar_target(knackCond1c, conditional_effects(knackM1c)[["Trust"]]),
  tar_target(knackCond1d, conditional_effects(knackM1d)[["Trust"]]),
  # Skidmore 2002
  tar_target(skidmoreData, loadDataSkidmore2002(fileSkidmore1, fileSkidmore2, fileSkidmore3, fileSkidmore4, iso)),
  tar_target(skidmoreM1a, fitModelSkidmore2002(skidmoreData, linCov, control = "none")),
  tar_target(skidmoreM1b, fitModelSkidmore2002(skidmoreData, linCov, control = "spatial")),
  tar_target(skidmoreM1c, fitModelSkidmore2002(skidmoreData, linCov, control = "cultural")),
  tar_target(skidmoreM1d, fitModelSkidmore2002(skidmoreData, linCov, control = "both")),
  tar_target(skidmoreSlope1a, posterior_samples(skidmoreM1a, pars = "b_logDisastersPerLandArea")[,1]),
  tar_target(skidmoreSlope1b, posterior_samples(skidmoreM1b, pars = "b_logDisastersPerLandArea")[,1]),
  tar_target(skidmoreSlope1c, posterior_samples(skidmoreM1c, pars = "b_logDisastersPerLandArea")[,1]),
  tar_target(skidmoreSlope1d, posterior_samples(skidmoreM1d, pars = "b_logDisastersPerLandArea")[,1]),
  tar_target(skidmoreCond1a, conditional_effects(skidmoreM1a)[["logDisastersPerLandArea"]]),
  tar_target(skidmoreCond1b, conditional_effects(skidmoreM1b)[["logDisastersPerLandArea"]]),
  tar_target(skidmoreCond1c, conditional_effects(skidmoreM1c)[["logDisastersPerLandArea"]]),
  tar_target(skidmoreCond1d, conditional_effects(skidmoreM1d)[["logDisastersPerLandArea"]]),
  # get spatial autocorrelation at 1000km
  tar_target(inglehartSA, getSpatialAutocorrelation1000km(inglehartData, inglehartM1d)),
  tar_target(skidmoreSA, getSpatialAutocorrelation1000km(skidmoreData, skidmoreM1d)),
  # get cultural phylogenetic signal
  tar_target(alesinaSignal, getCulturalPhylogeneticSignal(alesinaM1d)),
  tar_target(gelfandSignal, getCulturalPhylogeneticSignal(gelfandM1d)),
  # plot results of replications
  tar_target(plotReplications1, plotReplicationResults1(list(adamczykSlope1a, adamczykSlope1b, adamczykSlope1c, adamczykSlope1d,
                                                             alesinaSlope1a, alesinaSlope1b, alesinaSlope1c, alesinaSlope1d,
                                                             beckSlope1a, beckSlope1b, beckSlope1c, beckSlope1d,
                                                             beckSlope2a, beckSlope2b, beckSlope2c, beckSlope2d,
                                                             bockstetteSlope1a, bockstetteSlope1b, bockstetteSlope1c, bockstetteSlope1d,
                                                             easterlySlope1a, easterlySlope1b, easterlySlope1c, easterlySlope1d,
                                                             easterlySlope2a, easterlySlope2b, easterlySlope2c, easterlySlope2d,
                                                             fincherSlope1a, fincherSlope1b, fincherSlope1c, fincherSlope1d,
                                                             gelfandSlope1a, gelfandSlope1b, gelfandSlope1c, gelfandSlope1d,
                                                             inglehartSlope1a, inglehartSlope1b, inglehartSlope1c, inglehartSlope1d,
                                                             knackSlope1a, knackSlope1b, knackSlope1c, knackSlope1d,
                                                             skidmoreSlope1a, skidmoreSlope1b, skidmoreSlope1c, skidmoreSlope1d))),
  tar_target(plotReplications2, plotReplicationResults2(list(adamczykM1c$data, alesinaM1c$data, beckM1c$data,
                                                             beckM2c$data, bockstetteM1c$data, easterlyM1c$data,
                                                             easterlyM2c$data, fincherM1c$data, gelfandM1c$data,
                                                             inglehartM1c$data, knackM1c$data, skidmoreM1c$data),
                                                        list(adamczykCond1a, adamczykCond1b, adamczykCond1c, adamczykCond1d,
                                                             alesinaCond1a, alesinaCond1b, alesinaCond1c, alesinaCond1d,
                                                             beckCond1a, beckCond1b, beckCond1c, beckCond1d,
                                                             beckCond2a, beckCond2b, beckCond2c, beckCond2d,
                                                             bockstetteCond1a, bockstetteCond1b, bockstetteCond1c, bockstetteCond1d,
                                                             easterlyCond1a, easterlyCond1b, easterlyCond1c, easterlyCond1d,
                                                             easterlyCond2a, easterlyCond2b, easterlyCond2c, easterlyCond2d,
                                                             fincherCond1a, fincherCond1b, fincherCond1c, fincherCond1d,
                                                             gelfandCond1a, gelfandCond1b, gelfandCond1c, gelfandCond1d,
                                                             inglehartCond1a, inglehartCond1b, inglehartCond1c, inglehartCond1d,
                                                             knackCond1a, knackCond1b, knackCond1c, knackCond1d,
                                                             skidmoreCond1a, skidmoreCond1b, skidmoreCond1c, skidmoreCond1d))),
  tar_target(plotReplications3, plotReplicationResults3(list(adamczykM1d, alesinaM1d, beckM1d, beckM2d, bockstetteM1d,
                                                             easterlyM1d, easterlyM2d, fincherM1d, gelfandM1d, inglehartM1d,
                                                             knackM1d, skidmoreM1d),
                                                        list(adamczykData, alesinaData, beckData1, beckData2, bockstetteData,
                                                             easterlyData1, easterlyData2, fincherData, gelfandData, inglehartData,
                                                             knackData, skidmoreData))),
  tar_target(plotReplications4, plotReplicationResults4(list(adamczykM1d, alesinaM1d, beckM1d, beckM2d, bockstetteM1d,
                                                             easterlyM1d, easterlyM2d, fincherM1d, gelfandM1d, inglehartM1d,
                                                             knackM1d, skidmoreM1d))),
  tar_target(plotReplications5, plotReplicationResults5(list(adamczykM1d, alesinaM1d, beckM1d, beckM2d, bockstetteM1d,
                                                             easterlyM1d, easterlyM2d, fincherM1d, gelfandM1d, inglehartM1d,
                                                             knackM1d, skidmoreM1d),
                                                        list(adamczykM1a, alesinaM1a, beckM1a, beckM2a, bockstetteM1a,
                                                             easterlyM1a, easterlyM2a, fincherM1a, gelfandM1a, inglehartM1a,
                                                             knackM1a, skidmoreM1a),
                                                        list(adamczykData, alesinaData, beckData1, beckData2, bockstetteData,
                                                             easterlyData1, easterlyData2, fincherData, gelfandData, inglehartData,
                                                             knackData, skidmoreData))),
  
  #### Manuscript ####
  
  tar_render(manuscript, "manuscript.Rmd")
)
