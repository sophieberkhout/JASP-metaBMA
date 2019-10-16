#
# Copyright (C) 2018 University of Amsterdam
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# Main function ----

BayesianMetaAnalysis <- function(jaspResults, dataset, options) {
  
  # Ready: variables needed for the analysis (confidence interval missing)
  ready <- options[["effectSize"]] != "" && (options[["standardError"]] != "" || (all(unlist(options$confidenceInterval) != "")  && !is.null(unlist(options[["confidenceInterval"]])))) 
  
  # Dependencies: basically everything
  dependencies <- c("effectSize", "standardError", "confidenceInterval","studyLabels", "modelSpecification",
                    "allPos", "allNeg",
                    "priorH0FE", "priorH1FE", "priorH0RE", "priorH1RE",
                    "priorES", "cauchy", "normal", "t",
                    "informativeCauchyLocation", "informativeCauchyScale",
                    "checkLowerTruncCauchy", "checkUpperTruncCauchy",
                    "lowerTruncCauchy", "upperTruncCauchy",
                    "informativeNormalMean", "informativeNormalStd",
                    "checkLowerTruncNormal", "checkUpperTruncNormal",
                    "lowerTruncNormal", "upperTruncNormal",
                    "informativeTLocation", "informativeTScale", "informativeTDf",
                    "checkLowerTruncT", "checkUpperTruncT",
                    "lowerTruncT", "upperTruncT",
                    "priorSE", "inverseGamma", "inverseGammaShape", "inverseGammaScale",
                    "halfT", "informativehalfTScale", "informativehalfTDf",
                    "BFComputation", "integration", "bridgeSampling", "iterBridge",
                    "iterMCMC", "chainsMCMC")
  
  # Dataset with effectSize, standardError, and studyLabels
  # If data is null stuff is missing
  dataset <- .readData(jaspResults, options)
    
  # Table: Posterior Model Estimates
  if(options$mainTable){
    .bmaMainTable(jaspResults, dataset, options, ready, dependencies)
  }
  
  # Table: Model Probabilities
  if(options$postTable){
    .bmaPostModelTable(jaspResults, dataset, options, ready, dependencies)
  }
  
  # Table: Effect Sizes per Study
  if(options$esTable){
    .bmaEffectSizeTable(jaspResults, dataset, options, ready, dependencies)
  }
  
  # Plot: Prior(s); only when checked  
  if(options$plotPrior){
    .priorPlot(jaspResults, dataset, options, ready)
  }
  
  # Plot: Prior(s) and Posterior(s); only when checked 
  if(options$plotPosterior){  
    .priorAndPosteriorPlot(jaspResults, dataset, options, ready, dependencies)
  }
  
  # Plot: Forest plot; only when checked
  if(options$checkForestPlot || options$plotCumForest){  
    .forestPlot(jaspResults, dataset, options, ready, dependencies)
  }
  
  # Plot: Cumulative forest plot and sequential; only when checked
  if(options$plotSequential || options$plotSeqPM){
    .sequentialPlot(jaspResults, dataset, options, ready, dependencies)
  }
}
  
  # Get dataset
  .readData <- function(jaspResults, options){
    varES <- options[["effectSize"]]
    varSE <- options[["standardError"]]
    CI <- unlist(options$confidenceInterval)
    lower <- CI[[1]]
    upper <- CI[[2]]
    study <- options[["studyLabels"]]
    if(varES == "") varES <- NULL
    if(varSE == "") varSE <- NULL
    if(CI[[1]] == ""  || CI[[2]] == "" || is.null(CI)) {
      lower <- NULL
      upper <- NULL
    }
    if(study == "") study <- NULL
    variables.to.read <- c(varES, varSE, lower, upper, study)
    dataset <- .readDataSetToEnd(columns.as.numeric = variables.to.read, 
                                 exclude.na.listwise = variables.to.read)
    return(dataset)
  }
  
  # Save priors for later use (without data)
  .bmaPriors <- function(jaspResults, options) {
    if (!is.null(jaspResults[["bmaPriors"]])) return(jaspResults[["bmaPriors"]]$object)

    # Effect size prior parameters
    # Lower and upper limits without truncation
    lowerES <- -Inf
    upperES <- Inf
    
    # Cauchy prior
    if(options$priorES == "cauchy"){
      familyES <- "t"
      paramES <- c(options$informativeCauchyLocation, 
                   options$informativeCauchyScale, 1)
      # If truncated is checked
      if(options$checkLowerTruncCauchy){
        lowerES <- options$lowerTruncCauchy
      }
      if(options$checkUpperTruncCauchy){
        upperES <- options$upperTruncCauchy
      }
    }
    # Normal prior
    if(options$priorES == "normal"){
      familyES <- "norm"
      paramES <- c(options$informativeNormalMean, 
                   options$informativeNormalStd)
      # If truncated is checked
      if(options$checkLowerTruncNormal){
        lowerES <- options$lowerTruncNormal
      }
      if(options$checkUpperTruncNormal){
        upperES <- options$upperTruncNormal
      }
    }
    # T prior
    if(options$priorES == "t"){
      familyES <- "t"
      paramES <- c(options$informativeTLocation, 
                   options$informativeTScale, 
                   options$informativeTDf)
      # If truncated is checked
      if(options$checkLowerTruncT){
        lowerES <- options$lowerTruncT
      }
      if(options$checkUpperTruncT){
        upperES <- options$upperTruncT
      }
    }
  
    if (lowerES >= upperES)
      .quitAnalysis("The prior lower bound is not smaller than the upper bound.") 
    
  # Heterogeneity prior parameters
    # Inverse gamma prior
    if(options$priorSE == "inverseGamma"){
      familySE <- "invgamma"
      paramSE <- c(options$inverseGammaShape, 
                   options$inverseGammaScale)
    }
    # Half t prior
    if(options$priorSE == "halfT"){
      familySE <- "t"
      paramSE <- c(0, # location is always zero
                   options$informativehalfTScale, 
                   options$informativehalfTDf)
    }

    # Make priors (probability density functions)
    d <- metaBMA::prior(familyES, paramES, lowerES, upperES)
    tau <- metaBMA::prior(familySE, paramSE, 0)
    
    # Save priors
    jaspResults[["bmaPriors"]] <- createJaspState(
                    object=list(d = d, tau = tau),
                    dependencies=c("priorES", 
                                   "cauchy", "informativeCauchyLocation", "informativeCauchyScale",
                                   "checkLowerTruncCauchy", "lowerTruncCauchy", "checkUpperTruncCauchy", "upperTruncCauchy",
                                   "normal", "informativeNormalMean", "informativeNormalStd",
                                   "checkLowerTruncNormal", "lowerTruncNormal", "checkUpperTruncNormal", "upperTruncNormal",
                                   "t", "informativeTLocation", "informativeTScale","informativeTDf",
                                   "checkLowerTruncT", "lowerTruncT", "checkUpperTruncT", "upperTruncT",
                                   "priorSE", "inverseGamma", "inverseGammaShape", "inverseGammaScale",
                                   "halfT", "informativehalfTScale", "informativehalfTDf"))
    
    return(jaspResults[["bmaPriors"]]$object)
  }
 
  #For state
  .bmaResultsState <- function(jaspResults, dataset, options, dependencies) {
    if(!is.null(jaspResults[["bmaResults"]])) return(jaspResults[["bmaResults"]]$object);
    
    bmaResults                  <- createJaspState(object=.bmaResults(jaspResults, dataset, options), dependencies=dependencies)
    jaspResults[["bmaResults"]] <- bmaResults
    
    return(bmaResults$object)
}
   
  # Save the Bayesian meta-analysis
  .bmaResults <- function(jaspResults, dataset, options) {
    
    varES <- options[["effectSize"]]

    # Get necessary variables
    y <- dataset[, .v(options[["effectSize"]])]
    
    if(all(unlist(options[["confidenceInterval"]]) != "") && !is.null(unlist(options[["confidenceInterval"]]))){
      lower <- dataset[, .v(options$confidenceInterval[[1]][[1]])]
      upper <- dataset[, .v(options$confidenceInterval[[1]][[2]])]
 
      .hasErrors(dataset = dataset,
                 exitAnalysisIfErrors= TRUE,
                 custom = function() {
                   if (!all(lower < upper))
                   return("The 95% CI Lower Bound must be smaller than the Upper Bound.")
                   })
      
      SE <- (upper - lower)/3.92
    }
    if(options$standardError != ""){
      SE <- dataset[, .v(options[["standardError"]])]
      .hasErrors(dataset = dataset,
                 exitAnalysisIfErrors= TRUE,
                 type = "negativeValues",
                 negativeValues.target = options$standardError)
    }
    

    # Advanced: estimation settings
    iter <- options$iterMCMC
    chains <- options$chainsMCMC
    
    # Advanced: bayes factor computation
    if(options$BFComputation == "integration"){
      logml <- "integrate"
      logml_iter <- 5000
    } else if(options$BFComputation == "bridgeSampling"){
      logml <- "stan"
      logml_iter <- options$iterBridge
    }
    
    # Prior model probabilities
    prior <- c(options[["priorH0FE"]], options[["priorH1FE"]], 
                    options[["priorH0RE"]], options[["priorH1RE"]])
    
    # Get priors from jasp state
    .bmaPriors(jaspResults, options)
    
    d   <- jaspResults[["bmaPriors"]]$object[["d"]]
    tau <- jaspResults[["bmaPriors"]]$object[["tau"]]

    # Bayesian meta analysis
    if(options$modelSpecification != "CRE"){
    # Bayesian model averaging (includes fixed and random effects)
      results <- metaBMA::meta_bma(y     = y, 
                                   SE    = SE, 
                                   prior = prior, 
                                   d     = d, 
                                   tau   = tau,
                                   logml   = logml,
                                   logml_iter = logml_iter,
                                   iter     = iter,
                                   chains = chains)
    } else {
    # Ordered effects
      results <- metaBMA::meta_ordered(y = y, 
                                       SE = SE, 
                                       d = d, 
                                       tau = tau,
                                       # logml = logml,
                                       # logml_iter = logml_iter,
                                       iter = 10000 # because of an issue with stored variables, it is not yet possible to make it reactive.
                                       # chains = chains
      )
    }
    
    return(results)
  }
  
  .bmaGetModelName <- function(options) {
    if(options$modelSpecification == "CRE") return("ordered")
    if(options$modelSpecification == "BMA") return("averaged")
    if(options$modelSpecification == "RE")  return("random")
                                            return("fixed")
  }
  
  .bmaSequentialResults <- function(jaspResults, dataset, options, dependencies) {
    
    if(!is.null(jaspResults[["bmaSeqResults"]])) return(jaspResults[["bmaSeqResults"]]$object)
    
    startProgressbar(nrow(dataset)-1)
    
    seqResults <- list(mean=numeric(), lowerMain=numeric(), upperMain=numeric(), 
                       BFs=numeric(1), posterior_models=list()) 
    
    d                       <- .bmaPriors(jaspResults, options)[["d"]]
    seqResults$mean[1]      <- attr(d, "param")[1]
    s                       <- attr(d, "param")[2]
    seqResults$lowerMain[1] <- seqResults$mean[1] - 1.96 * s
    seqResults$upperMain[1] <- seqResults$mean[1] + 1.96 * s
    
    modelName <- .bmaGetModelName(options)
    
    for(i in 2:nrow(dataset)){
      m <- .bmaResults(jaspResults, dataset[1:i, ], options)
      
      seqResults$mean[i]      <- m$estimates[modelName, "mean"]
      seqResults$lowerMain[i] <- m$estimates[modelName, "2.5%"]
      seqResults$upperMain[i] <- m$estimates[modelName, "97.5%"]
    
      if(options$modelSpecification == "BMA") seqResults$BFs[i] <- m$inclusion$incl.BF
      if(options$modelSpecification == "FE")  seqResults$BFs[i] <- m$BF["fixed_H1", "fixed_H0"]
      if(options$modelSpecification == "RE")  seqResults$BFs[i] <- m$BF["random_H1", "random_H0"]
      if(options$modelSpecification == "CRE") seqResults$BFs[i] <- m$BF["ordered", "null"]
      
      seqResults$posterior_models[[i]] <- m$posterior_models
      
      progressbarTick()
    }
    
    jaspResults[["bmaSeqResults"]] <- createJaspState(object=seqResults, dependencies=dependencies)
    
    return(jaspResults[["bmaSeqResults"]]$object)
  }
  
  # Table: Posterior Model Estimates
  .bmaMainTable <- function(jaspResults, dataset, options, ready, dependencies) {
    if (!is.null(jaspResults[["bmaTable"]])) return()
    bmaTable <- createJaspTable(title = "Posterior Model Estimates")
    bmaTable$position <- 1
    
    # Add standard depencies
    bmaTable$dependOn(c(dependencies, "mainTable", "BF"))
    
    if (options$BF == "BF10")
      bfTitle <- "BF\u2081\u2080"
    else if (options$BF == "BF01")
      bfTitle <- "BF\u2080\u2081"
    else
      bfTitle <- "Log(BF\u2081\u2080)"
    
    # Add columns
    bmaTable$addColumnInfo(name = "model", title = "", type = "string", combine = TRUE)
    bmaTable$addColumnInfo(name = "parameter", title = "", type = "string")
    bmaTable$addColumnInfo(name = "ES", title = "Mean", type = "number")
    bmaTable$addColumnInfo(name = "SD", title = "SD", type = "number")
    bmaTable$addColumnInfo(name = "lb", title = "Lower", type = "number",
                           overtitle = "95% CI")
    bmaTable$addColumnInfo(name = "ub", title = "Upper", type = "number",
                           overtitle = "95% CI")    
    bmaTable$addColumnInfo(name = "BF", title = bfTitle, type = "number")

    jaspResults[["bmaTable"]] <- bmaTable
    
    # Check if ready
    if(!ready){
      return()
    }
    
    # Get analysis results
    m <- .bmaResultsState(jaspResults, dataset, options, dependencies)
    
    # Row names (tried to get modelRE idented, but failed)
    modelBMA <- "Averaged"
    modelFE  <- "Fixed effects"   
    modelRE  <- "Random effects"
    modelCRE <- "Ordered effects"
    
    eta <- "\u03B7" 
    tau <- "\u03C4"
    mu <- "\u03BC"

    # Get results per column (different per model)
    if(options$modelSpecification == "BMA"){
      model <- c(modelFE, modelRE, modelRE, modelBMA)
      parameter <- c(eta, eta, tau, eta)
      group <- c(T, T, F, T)
      meanES <- c(m$estimates["fixed", "mean"], 
                  m$estimates["random", "mean"],
                  m$meta$random$estimates["tau", "mean"],
                  m$estimates["averaged", "mean"])
      meanSD <- c(m$estimates["fixed", "sd"], 
                  m$estimates["random", "sd"], 
                  m$meta$random$estimates["tau", "sd"],
                  m$estimates["averaged", "sd"])
      lower <- c(m$estimates["fixed", "2.5%"], 
                 m$estimates["random", "2.5%"], 
                 m$meta$random$estimates["tau", "2.5%"],
                 m$estimates["averaged", "2.5%"])
      upper <- c(m$estimates["fixed", "97.5%"], 
                 m$estimates["random", "97.5%"], 
                 m$meta$random$estimates["tau", "97.5%"],
                 m$estimates["averaged", "97.5%"])
      BF <- c(m$BF["fixed_H1", "fixed_H0"], 
              m$BF["random_H1", "random_H0"], 
              m$BF["random_H1", "fixed_H1"],
              m$inclusion$incl.BF)
    }
    else if(options$modelSpecification == "RE"){
      model <- c(modelRE, modelRE)
      parameter <- c(tau, eta)
      group <- c(T, F)
      meanES <- m$meta$random$estimates[, "mean"]
      meanSD <- m$meta$random$estimates[, "sd"]
      lower <- m$meta$random$estimates[, "2.5%"]
      upper <- m$meta$random$estimates[, "97.5%"]
      BF <- c(m$BF["random_H1", "random_H0"], 
              m$BF["random_H1", "fixed_H1"])
    }
    else if(options$modelSpecification == "FE"){
      model <- modelFE
      parameter <- tau
      group <- T
      meanES <- m$meta$fixed$estimates[, "mean"]
      meanSD <- m$meta$fixed$estimates[, "sd"]
      lower <- m$meta$fixed$estimates[, "2.5%"]
      upper <- m$meta$fixed$estimates[, "97.5%"]
      BF <- m$BF["fixed_H1", "fixed_H0"]
    }
    else if(options$modelSpecification == "CRE"){
      model <- c(modelFE, modelCRE, modelCRE, modelRE, modelRE)
      parameter <- c(eta, eta, tau, eta, tau)
      group <- c(T, T, F, T, F)
      meanES <- c(m$estimates["fixed", "mean"],
                  m$meta$ordered$estimates[c("average_effect", "tau"), "mean"],
                  m$meta$random$estimates[, "mean"])
      meanSD <- c(m$estimates["fixed", "sd"],
                  m$meta$ordered$estimates[c("average_effect", "tau"), "sd"],
                  m$meta$random$estimates[, "sd"])
      lower <- c(m$estimates["fixed", "2.5%"],
                 m$meta$ordered$estimates[c("average_effect", "tau"), "2.5%"],
                 m$meta$random$estimates[, "2.5%"])
      upper <- c(m$estimates["fixed", "97.5%"],
                 m$meta$ordered$estimates[c("average_effect", "tau"), "97.5%"],
                 m$meta$random$estimates[, "97.5%"])
      BF <- c(m$BF["fixed", "null"],
              m$BF["ordered", "null"],
              m$BF["ordered", "fixed"],
              m$BF["random", "null"],
              m$BF["random", "fixed"])
    }

    if(options$modelSpecification == "CRE") creBF <- m$BF["ordered", "random"]
    if(options$BF == "BF01"){
      BF <- 1/BF
      if(options$modelSpecification == "CRE"){
        creBF <- 1/m$BF["ordered", "random"]
      }
    }
    if(options$BF == "logBF10"){
      BF <- log(BF)
      if(options$modelSpecification == "CRE"){
        creBF <- log(m$BF["ordered", "random"])
      }
    }
    # Add results to table
    rows <- data.frame(model = model, 
                       parameter = parameter,
                       ES = meanES, 
                       SD = meanSD, 
                       lb = lower, 
                       ub = upper, 
                       BF = BF,
                       .isNewGroup = group)
    row.names(rows) <- paste0("row", 1:length(model))
    
    bmaTable$setData(rows)
    
    if(options$modelSpecification == "BMA") {
      bmaTable$addFootnote("Averaged over the fixed effects model and the random effects model.",
                           colNames = "model", rowNames="row4") 
      # bmaTable$addFootnote("The estimates for the model averaged \u03C4 are not yet available. It is currently in development.",
      #                    colNames = "parameter", rowNames="row5") 
    }
    if(options$modelSpecification == "CRE"){
      bmaTable$addFootnote(paste0("Bayes Factor of the constrained random effects H\u2081 versus the fixed effects H\u2080 model. The Bayes Factor for the ordered H\u2081 versus the unconstrained (random) effects H\u2081 model is ",
                                 round(creBF, 3),
                                 "."),
                           colNames = "BF", rowNames="row2") 
      
    }
  }
  
  # Table: Model Probabilities
  .bmaPostModelTable <- function(jaspResults, dataset, options, ready, dependencies) {
    if (!is.null(jaspResults[["postTable"]])) return()
    postTable <- createJaspTable(title = "Model Probabilities")
    postTable$dependOn(c(dependencies, "postTable"))
    
    # Add columns
    postTable$addColumnInfo(name = "model", title = "", type = "string")
    postTable$addColumnInfo(name = "priorProb",   title = "Prior",   type = "number")
    postTable$addColumnInfo(name = "postProb",   title = "Posterior",   type = "number")
    
    # Add table to output
    jaspResults[["postTable"]] <- postTable
    jaspResults[["postTable"]]$position <- 2
    
    # Check if ready
    if(!ready){
      return()
    }
    
    # Get results from jasp state
    m <- .bmaResultsState(jaspResults, dataset, options, dependencies)
    
    
    # Get results per column (different per model)
    if(options$modelSpecification == "BMA"){
      model <- c("Fixed H\u2080", "Fixed H\u2081", "Random H\u2080", "Random H\u2081")
      postProb <- m$posterior_models
      priorProb <- m$prior_models
    }
    if(options$modelSpecification == "FE"){
      model <- c("Fixed H\u2080", "Fixed H\u2081")
      postProb <- m$posterior_models[c("fixed_H0", "fixed_H1")]
      priorProb <- m$prior_models[1:2]
    }
    if(options$modelSpecification == "RE"){
      model <- c("Random H\u2080", "Random H\u2081")
      postProb <- m$posterior_models[c("random_H0", "random_H1")]
      priorProb <- m$prior_models[3:4]
    }
    if(options$modelSpecification == "CRE"){
      model <- c("Fixed H\u2080", "Fixed H\u2081", "Ordered H\u2081", "Random H\u2081")
      postProb <- m$posterior_models
      priorProb <- m$prior_models
    }
    
    # Fill table
    row <- data.frame(model = model, priorProb =  priorProb, postProb = postProb)
    postTable$addRows(row)
    
  }
  
  # Table: Effect Sizes per Study
  .bmaEffectSizeTable <- function(jaspResults, dataset, options, ready, dependencies) {
    if (!is.null(jaspResults[["esTable"]])) return()
    esTable <- createJaspTable(title = "Effect Sizes per Study")
    esTable$dependOn(c(dependencies, "esTable"))
    
    # Add standard columns
    esTable$addColumnInfo(name = "study", title = "", type = "string")
    esTable$addColumnInfo(name = "observedES", title = "Observed", type = "number")

    # Add conditional columns
    if(options$modelSpecification != "FE"){
      esTable$addColumnInfo(name = "estimatedES", title = "Mean", type = "number",
                            overtitle = "Estimated")
      esTable$addColumnInfo(name = "estimatedLower", title = "Lower", type = "number",
                            overtitle = "Estimated")
      esTable$addColumnInfo(name = "estimatedUpper", title = "Upper", type = "number",
                            overtitle = "Estimated")
    }
    
    # Add table to output
    jaspResults[["esTable"]] <- esTable
    jaspResults[["esTable"]]$position <- 3
    
    # Check if ready
    if(!ready){
      return()
    }    

    # Get results from jasp state
    m <- .bmaResultsState(jaspResults, dataset, options, dependencies)
    
    # Get effect size variable
    varES <- dataset[, .v(options[["effectSize"]])]
    
    # Create empty vectors
    estimatedES <- rep(NA, length(varES))
    estimatedLower <- rep(NA, length(varES))
    estimatedUpper <- rep(NA, length(varES))
    
    # Fill vectors with estimation variables if not FE
    if(options$modelSpecification != "FE"){
      estimatedES    <- rstan::summary(m$meta$random$stanfit_dstudy)$summary[3:(length(varES) + 2), "mean"]
      estimatedLower <- rstan::summary(m$meta$random$stanfit_dstudy)$summary[3:(length(varES) + 2), "2.5%"]
      estimatedUpper <- rstan::summary(m$meta$random$stanfit_dstudy)$summary[3:(length(varES) + 2), "97.5%"]
    }
    
    if(options$modelSpecification == "CRE"){
      estimatedES <- rstan::summary(m$meta$ordered$stanfit_dstudy)$summary[3:(length(varES) + 2), "mean"]
      estimatedLower <- rstan::summary(m$meta$ordered$stanfit_dstudy)$summary[3:(length(varES) + 2), "2.5%"]
      estimatedUpper <- rstan::summary(m$meta$ordered$stanfit_dstudy)$summary[3:(length(varES) + 2), "97.5%"]

    }

    # Add studylabels when given, otherwise use "Study n"
    if(options[["studyLabels"]] != ""){
      studyLabels <- dataset[, .v(options[["studyLabels"]])]
    } else {
      studyLabels <- paste("Study", 1:length(varES))
    }
    
    # Only show conditional columns for right analysis
    esTable$showSpecifiedColumnsOnly <- TRUE
    
    # Add results to table
    row <- data.frame(study = studyLabels, 
                      observedES = varES, 
                      estimatedES = estimatedES, 
                      estimatedLower = estimatedLower, 
                      estimatedUpper = estimatedUpper)
    esTable$addRows(row)
    
    if(options$modelSpecification != "FE"){
    esTable$addFootnote("Posterior mean and 95% CI estimates from the random effects model.",
                        colNames = c("estimatedES", "estimatedLower", "estimatedUpper"))
    }
  }
  
  # Plot: prior(s)  
  .priorPlot <- function(jaspResults, dataset, options, ready) {
    priorContainer <- createJaspContainer(title = "Prior")
    priorContainer$dependOn("plotPrior")
    jaspResults[["priorContainer"]] <- priorContainer
    jaspResults[["priorContainer"]]$position <- 4
    
    # Create empty plot
    priorPlot <- createJaspPlot(plot = NULL, title = "Effect Size", 
                                width = 450, height = 350)
    
    # Custom dependencies (only dependent on prior settings)
    priorPlot$dependOn(c("priorES", 
                         "cauchy", "informativeCauchyLocation", "informativeCauchyScale",
                         "checkLowerTruncCauchy", "lowerTruncCauchy", "checkUpperTruncCauchy", "upperTruncCauchy",
                         "normal", "informativeNormalMean", "informativeNormalStd",
                         "checkLowerTruncNormal", "lowerTruncNormal", "checkUpperTruncNormal", "upperTruncNormal",
                         "t", "informativeTLocation", "informativeTScale","informativeTDf",
                         "checkLowerTruncT", "lowerTruncT", "checkUpperTruncT", "upperTruncT"))
    
    # Fill plot with effect size prior
    .fillPriorPlot(priorPlot, jaspResults, dataset, options, type = "ES")
    priorContainer[["ES"]] <- priorPlot

    # Make plot hetergeneity prior
    if(options$modelSpecification != "FE"){
      priorPlotSE <- createJaspPlot(plot = NULL, title = "Heterogeneity", 
                                    width = 350, height = 350)
      priorPlotSE$dependOn(c("priorSE", "inverseGamma", "inverseGammaShape", "inverseGammaScale",
                             "halfT", "informativehalfTScale", "informativehalfTDf"))
      .fillPriorPlot(priorPlotSE, jaspResults, dataset, options, type = "SE")
      priorContainer[["SE"]] <- priorPlotSE
    }
  }
  
  # Fill prior plot
  .fillPriorPlot <- function(priorPlot, jaspResults, dataset, options, type){
    # Get priors from jasp state
    if (is.null(jaspResults[["bmaPriors"]]))
      .bmaPriors(jaspResults, options)
    priors <- jaspResults[["bmaPriors"]]$object
    
    # Get parameters and x limits
    if(type == "ES"){
      prior <- priors$d
      mean <- attr(prior, "param")[1]
      s <- attr(prior, "param")[2]
      xlimLeft <- mean - (s * 5)
      xlab <- expression("Effect size "*eta)
    } else if(type == "SE"){
      prior <- priors$tau
      mean <- attr(prior, "param")[1]
      s <- attr(prior, "param")[2]
      xlimLeft <- 0
      xlab <- expression("Heterogeneity "*tau)
    }
    if(options$modelSpecification == "CRE" && options$direction == "allPos"){
      xlimLeft <- 0
    } else if(options$modelSpecification == "CRE" && options$direction == "allNeg"){
      xlimRight <- 0
    }
    xlimRight <- mean + (s * 5)

    xlimLeft <- xlimLeft - 0.05
    xlimRight <- xlimRight + 0.05
    
    # Create dataframe for ggplot
    x <- c(xlimLeft, xlimRight)
    df <- data.frame(x = x)
    
    xBreaks <- getPrettyAxisBreaks(seq(xlimLeft, xlimRight, 0.5))
    
    # Plot density function
    plot <- ggplot2::ggplot(df, ggplot2::aes(x)) +
            ggplot2::stat_function(fun = prior, n = 1000, size = 1) +
            ggplot2::labs(x = xlab, y = "Density") +
            ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
            ggplot2::xlim(xlimLeft, xlimRight) +
            ggplot2::scale_x_continuous(breaks = xBreaks)
    # df2 <- data.frame(x = x, y = x)
    # plot <- ggplot2::ggplot(df2, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
    plot <- themeJasp(plot)
    priorPlot$plotObject <- plot
    return()
  }

  # Plot: Prior and Posterior
  .priorAndPosteriorPlot <- function(jaspResults, dataset, options, ready, dependencies) {
    postContainer <- createJaspContainer(title = "Prior and Posteriors")
    postContainer$dependOn(c(dependencies, "plotPosterior", "shade"))
    jaspResults[["postContainer"]] <- postContainer
    jaspResults[["postContainer"]]$position <- 5
    
    # Create empty plot
    postPlotES <- createJaspPlot(plot = NULL, title = "Effect size", width = 450, height = 350)
    
    # Custom dependencies
    postPlotES$dependOn(c("priorES", "cauchy", "normal", "t",
                           "informativeCauchyLocation", "informativeCauchyScale",
                           "lowerTruncCauchy", "upperTruncCauchy",
                           "informativeNormalMean", "informativeNormalStd",
                           "lowerTruncNormal", "upperTruncNormal",
                           "informativeTLocation", "informativeTScale", "informativeTDf",
                           "lowerTruncT", "upperTruncT"))
    
    # Check if ready
    if(!ready){
      return()
    }
    
    # Fill posterior plot effect size
    .fillPostPlot(postPlotES, jaspResults, dataset, options, type = "ES")
    postContainer[["ES"]] <- postPlotES

    # Make posterior plot heterogeneity
    if(options$modelSpecification != "FE"){
      postPlotSE <- createJaspPlot(plot = NULL, title = "Heterogeneity", width = 450, height = 350)
      postPlotSE$dependOn(c("priorSE", "inverseGamma", "inverseGammaShape", "inverseGammaScale",
                            "halfT", "informativehalfTScale", "informativehalfTDf"))
      postContainer[["SE"]] <- postPlotSE
      .fillPostPlot(postPlotSE, jaspResults, dataset, options, type = "SE")
   }
  }
  
  # Fill prior and posterior plot
  .fillPostPlot <- function(postPlot, jaspResults, dataset, options, type){
    # Get results from jasp state
    m <- jaspResults[["bmaResults"]]$object

    # Get prior and posterior functions, and 95% CI intervals
    mPostFixed <- m$meta$fixed$posterior_d
    mPostRandom <- m$meta$random$posterior_d
    alpha <- 0.2
    est <- m$estimates[1, 1]
    x <- seq(est - 10, est + 10, .001)
    postName <- "Posterior"
    # Effect size priors
    if(type == "ES"){
      mPrior <- m$prior_d$fixed
      xlab <- expression("Effect size "*eta)
      xlim <- c(-4, 4)
      valuesCol <- c("#009E73", "#D55E00", "black", "black")
      valuesLine <- c("solid", "solid", "solid", "dotted")
      if(options$modelSpecification == "BMA"){
        mPost <- m$posterior_d
        int <- c(m$estimates["averaged", "2.5%"], m$estimates["averaged", "97.5%"])
        postName <- "Averaged"
        labelsModel <- c(expression("Fixed H"[1]), expression("Random H"[1]), expression("Averaged H"[1]), expression("Prior H"[1]))
      } else if(options$modelSpecification == "RE"){
        mPost <- m$meta$random$posterior_d
        int <- c(m$estimates["random", "2.5%"], m$estimates["random", "97.5%"])
        valuesCol <- c("black", "black")
        valuesLine <- c("solid", "dotted")
        postName <- "Random"
        labelsModel <- c(expression("Random H"[1]), expression("Prior H"[1]))
      } else if(options$modelSpecification == "FE"){
        mPost <- m$meta$fixed$posterior_d
        int <- c(m$estimates["fixed", "2.5%"], m$estimates["fixed", "97.5%"])
        valuesCol <- c("black", "black")
        valuesLine <- c("solid", "dotted")
        postName <- "Fixed"
        labelsModel <- c(expression("Fixed H"[1]), expression("Prior H"[1]))
      } else if(options$modelSpecification == "CRE"){
        mPost <- m$posterior_d
        int <- c(m$estimates["ordered", "2.5%"], m$estimates["ordered", "97.5%"])
        xlim <- c(0, 4)
        postName <- "Ordered"
        valuesCol <- c("#009E73", "black", "#D55E00", "black")
        labelsModel <- c(expression("Fixed H"[1]), expression("Ordered H"[1]), expression("Random H"[1]), expression("Prior H"[1]))
        
      }
    # Heterogeneity priors
    } else if(type == "SE"){
      mPrior <- m$meta$random$prior_tau
      mPost <- m$meta$random$posterior_tau
      int <- c(m$meta$random$estimates["tau", "2.5%"], m$meta$random$estimates["tau", "97.5%"])
      xlab <- expression("Heterogeneity "*tau)
      xlim <- c(0, 3)
      valuesCol <- c("black", "black")
      valuesLine <- c("solid", "dotted")
      alpha <- 0.3
      postName <- "Random"
      if(options$modelSpecification == "BMA") labelsModel <- c(expression("Random H"[1]), expression("Prior H"[1]))
      if(options$modelSpecification == "CRE") labelsModel <- c(expression("Ordered H"[1]), expression("Random H"[1]), expression("Prior H"[1]))
      if(options$modelSpecification == "FE") labelsModel <- c(expression("Fixed H"[1]), expression("Prior H"[1]))
      if(options$modelSpecification == "RE") labelsModel <- c(expression("Random H"[1]), expression("Prior H"[1]))
    }

    limitsX <- mPost(x)
    x <- x[limitsX > 0.00001]
    x2 <- x

    x[which.min(x)] <- min(x) - 0.05
    x[which.max(x)] <- max(x) + 0.05
    
    if(type == "ES"){
      x2 <- c(x, x)
      yPostSE <- c(mPostFixed(x), mPostRandom(x))
      g2 <- rep(c("Fixed", "Random"), each = length(x))
      
    }    
    yPost <- mPost(x)
    yPrior <- mPrior(x)

      
    df <- data.frame(x = c(x, x), y = c(mPrior(x), mPost(x)), g = rep(c("Prior", postName), each = length(x)))
    df$g <- factor(df$g, levels = c(postName, "Prior"))
    if(type == "ES" && (options$modelSpecification == "BMA" || options$modelSpecification == "CRE")){
      mPostFixed <- m$meta$fixed$posterior_d
      mPostRandom <- m$meta$random$posterior_d
      df2 <- data.frame(x = x2, y = yPostSE, g = g2)
      df <- rbind(df2, df)
    }
    
    if(type == "SE" && options$modelSpecification == "BMA") valuesCol <- c("#D55E00", "black")
    
    if(options$modelSpecification == "CRE" && type == "ES"){
      df$g <- factor(df$g, levels = c("Fixed", "Ordered", "Random", "Prior"))
    }
    if(options$modelSpecification == "CRE" && type == "SE"){
      mPost <- m$meta$ordered$posterior_tau
      mPostRandom <- m$meta$random$posterior_tau
      int <- c(m$meta$ordered$estimates["tau", "2.5%"], m$meta$ordered$estimates["tau", "97.5%"])
      df <- data.frame(x = c(x, x, x), 
                       y = c(mPrior(x), mPost(x), mPostRandom(x)), 
                       g = rep(c("Prior", "Ordered", "Random"), each = length(x)))
      df$g <- factor(df$g, levels = c("Ordered", "Random", "Prior"))
      valuesCol <- c("black", "#D55E00", "black")
      valuesLine <- c("solid", "solid", "dotted")
    }
  
    plot <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = x, y = y, color = g, linetype = g)) +
          JASPgraphs::geom_line() +
          ggplot2::scale_x_continuous(xlab, breaks = getPrettyAxisBreaks(df$x)) +
          ggplot2::scale_y_continuous("Density", breaks = JASPgraphs::getPrettyAxisBreaks(df$y))
    
    plot <- plot + 
      ggplot2::scale_linetype_manual("", values = valuesLine, labels = labelsModel) +
      ggplot2::scale_color_manual("", values = valuesCol, labels = labelsModel) +
      ggplot2::theme(legend.text.align = 0)
    
    if(options$shade){
      plot <- plot + 
        ggplot2::stat_function(fun = mPost,
                      xlim = int,
                      geom = "area", alpha = alpha, show.legend = F, size = 0, fill = "grey") +
        ggplot2::geom_vline(xintercept = 0, linetype = "dotted")
    }
    
    if(any(x == 0)){
      plot <- plot + 
        ggplot2::geom_vline(xintercept = 0, linetype = "dotted")
    }
    
    xr   <- range(df$x)
    idx  <- which.max(df$y)
    xmax <- df$x[idx]
    if (xmax > mean(xr)) {
      legend.position = c(0.2, 0.875)
    } else {
      legend.position = c(0.80, 0.875)
    }
    plot <- JASPgraphs::themeJasp(plot, legend.position = legend.position)
    postPlot$plotObject <- plot
    return()
  }
  
  # Plot: Forest plot
  .forestPlot <- function(jaspResults, dataset, options, ready, dependencies) {
    forestContainer <- createJaspContainer(title = "Forest Plot")
    forestContainer$dependOn(dependencies)
    jaspResults[["forestContainer"]] <- forestContainer
    jaspResults[["forestContainer"]]$position <- 6
    
    # Get studylabels
    if(options[["studyLabels"]] != ""){
      studyLabels <- as.character(dataset[, .v(options[["studyLabels"]])])
    } else {
      studyLabels <- paste("Study", 1:nrow(dataset))
    }
    
    # Scale the height and width of the plot
    height <- nrow(dataset) * 50
    width <- 500 + (nchar(max(studyLabels)) * 5)
    
    # title of plot based on observed/estimated 
    if(options$forestPlot == "plotForestBoth"){
      title <- "Observed and estimated study effects"
      height <- nrow(dataset) * 50 * 1.5
    } else  if(options$forestPlot == "plotForestEstimated"){
      title <- "Estimated study effects"
    } else  if(options$forestPlot == "plotForestObserved" || options$modelSpecification == "FE"){
      title <- "Observed study effects"
    }
    
    # Check if ready    
    if(!ready){
      return()
    } 
    
    # Create empty plot
    if(options$checkForestPlot){
      forestPlot <- createJaspPlot(plot = NULL, title = title, height = height, width = width)
      # Fill plot
      forestPlot$dependOn(c("plotForestObserved", "plotForestEstimated", "plotForestBoth", 
                            "checkForestPlot", "ascendingForest", "labelForest",
                            "orderForest"))
      .fillForestPlot(forestPlot, jaspResults, dataset, options, studyLabels)
      # Add plot to container
      forestContainer[["forestPlot"]] <- forestPlot
    }
    
    
    
    
    if(options$plotCumForest){
      cumForestPlot <- createJaspPlot(plot = NULL, title = "Cumulative forest plot", height = height, width = width)
      cumForestPlot$dependOn("plotCumForest")
      .fillCumForest(cumForestPlot, jaspResults, dataset, options, studyLabels, dependencies)
      forestContainer[["cumForestPlot"]] <- cumForestPlot
    }

    
  }
  
  .fillForestPlot <- function(forestPlot, jaspResults, dataset, options, studyLabels){
    # Get analysis results from jasp state
    m <- jaspResults[["bmaResults"]]$object
    
    # Create effect size and standard error variable and make dataframe
    varES <- dataset[, .v(options[["effectSize"]])]
    
    if(all(unlist(options[["confidenceInterval"]]) != "") && !is.null(unlist(options[["confidenceInterval"]]))){
      lower <- dataset[, .v(options[["confidenceInterval"]][[1]][[1]])]
      upper <- dataset[, .v(options[["confidenceInterval"]][[1]][[2]])]
      varSE <- (upper - lower)/3.92
    }
    if(options[["standardError"]] != ""){
      varSE <- dataset[, .v(options[["standardError"]])]
    }

    # Assign weights for the observed point sizes
    weight <- 1/varSE^2
    weight_scaled <- ((4 - 1)*(weight - min(weight))) / (max(weight) - min(weight)) + 2
    
    # Assign weights for the estimated point sizes
    # Should be different for ordered analysis
    se_estimated <- rstan::summary(m$meta$random$stanfit_dstudy)$summary[3:(length(varES) + 2), "se_mean"]
    
    if(options[["modelSpecification"]] == "CRE"){
      se_estimated <- rstan::summary(m$meta$ordered$stanfit_dstudy)$summary[3:(length(varES) + 2), "se_mean"]
    }
    
    weight_estimated <- 1 / se_estimated^2
    weight_estimated_scaled <- ((4 - 1) * (weight_estimated - min(weight_estimated))) / (
                                max(weight_estimated) - min(weight_estimated)) + 2

    # Create text object for next to the observed points
    ci <- .95
    lower <- varES - qnorm((ci+1)/2) * varSE
    upper <- varES + qnorm((ci+1)/2) * varSE

    text_observed <- paste(sprintf('%.2f', varES),
                           " [",
                           sprintf('%.2f', lower),
                           ", ",
                           sprintf('%.2f', upper),
                           "]",
                           sep = "")

    # Get estimated points and CI's
    if(options$modelSpecification == "BMA" || options$modelSpecification == "RE" || options$modelSpecification == "CRE"){
      mean_estimates <- rstan::summary(m$meta$random$stanfit_dstudy)$summary[3:(length(varES) + 2), "mean"]
      lower_estimates <- rstan::summary(m$meta$random$stanfit_dstudy)$summary[3:(length(varES) + 2), "2.5%"]
      upper_estimates <- rstan::summary(m$meta$random$stanfit_dstudy)$summary[3:(length(varES) + 2), "97.5%"]
    } 
    # The estimates for the ordered analysis are not always saved
    if(options$modelSpecification == "CRE"){
      mean_estimates <- rstan::summary(m$meta$ordered$stanfit_dstudy)$summary[1:length(varES) + 2, "mean"]
      lower_estimates <- rstan::summary(m$meta$ordered$stanfit_dstudy)$summary[1:length(varES) + 2, "2.5%"]
      upper_estimates <- rstan::summary(m$meta$ordered$stanfit_dstudy)$summary[1:length(varES) + 2, "97.5%"]
    }

    # Create text object for estimated points
    if(options$modelSpecification != "FE"){
      text_estimated <- paste(sprintf('%.2f', mean_estimates),
                              " [",
                              sprintf('%.2f', lower_estimates),
                              ", ",
                              sprintf('%.2f', upper_estimates),
                              "]",
                              sep = "")
    }

    
    

    
    # Make index for model diamond
    modelIndex <- .bmaGetModelName(options)

    
    yDiamond <- -0.5
    
    if (options$modelSpecification == "BMA" || options$modelSpecification == "CRE") {
      if (options$forestPlot == "plotForestBoth")
        yDiamond <- c(-0.5, -1.1, -1.7)
      else
        yDiamond <- c(-0.5, -1.5, -2.5)
    }
    
    # Create diamond for averaged or ordered model
    meanMain <- m$estimates[modelIndex, "mean"]
    lowerMain <- m$estimates[modelIndex, "2.5%"]
    upperMain <- m$estimates[modelIndex, "97.5%"]
    if(modelIndex == "ordered"){
      yMain <- yDiamond[2]
    } else if(options$modelSpecification == "BMA"){
      yMain <- yDiamond[3]
    } else yMain <- -0.5

    d <- data.frame(x = c(lowerMain, meanMain,
                          upperMain, meanMain),
                    y = c(yMain, yMain + 0.25, 
                          yMain, yMain - 0.25))

    # Text object for next to model diamond
    textMain <- paste0(sprintf('%.2f', meanMain), " [",
                       sprintf('%.2f', lowerMain), ", ",
                       sprintf('%.2f', upperMain), "]")

    
    # Create diamond for fixed model
    meanFixed <- m$estimates["fixed", "mean"]
    lowerFixed <- m$estimates["fixed", "2.5%"]
    upperFixed <- m$estimates["fixed", "97.5%"]
    yFixed <- yDiamond[1]
    
    d.fixed <- data.frame(x = c(lowerFixed, meanFixed,
                                upperFixed, meanFixed),
                          y = c(yFixed, yFixed + 0.25, 
                                yFixed, yFixed - 0.25))

    text_fixed <- paste0(sprintf('%.2f', meanFixed), " [",
                         sprintf('%.2f', lowerFixed), ", ",
                         sprintf('%.2f', upperFixed), "]")

    # Create diamond for random model
    meanRandom <- m$estimates["random", "mean"]
    lowerRandom <- m$estimates["random", "2.5%"]
    upperRandom <- m$estimates["random", "97.5%"]
    if(options$modelSpecification == "RE"){
      yRandom <- -0.5
    } else if(options$modelSpecification == "BMA"){
      yRandom <- yDiamond[2]
    } else if(options$modelSpecification == "CRE"){
      yRandom <- yDiamond[3]
    } else yRandom <- 0

    d.random <- data.frame(x = c(lowerRandom, meanRandom,
                                 upperRandom, meanRandom),
                           y = c(yRandom, yRandom + 0.25, 
                                 yRandom, yRandom - 0.25))

    text_random <- paste0(sprintf('%.2f', meanRandom), " [",
                          sprintf('%.2f', lowerRandom), ", ",
                          sprintf('%.2f', upperRandom), "]")

    # Get y coordinates, labels, and text for diamonds
    if(options$modelSpecification == "BMA"){
      model <- c("Fixed effects", "Random effects", "Averaged")
      textDiamond <- c(text_fixed, text_random, textMain)
    } else if(options$modelSpecification == "RE"){
      model <- "Random effects"
      textDiamond <- text_random
    } else if(options$modelSpecification == "FE"){
      model <- "Fixed effects"
      textDiamond <- text_fixed
    } else if(options$modelSpecification == "CRE"){
      model <- c("Fixed effects", "Ordered effects", "Random effects")
      textDiamond <- c(text_fixed, textMain, text_random)
    }

    # Shape if only observed points
    shape <- 15

    df <- data.frame(effectSize = varES, y = length(varES):1,
                     studyLabels = studyLabels,
                     weight_scaled = weight_scaled, 
                     lower = lower, upper = upper,
                     text = text_observed)
    
    # Change objects if only estimated points
    if(options$forestPlot == "plotForestEstimated"){
      df <- data.frame(effectSize = mean_estimates, y = length(varES):1, 
                       studyLabels = studyLabels,
                       weight_scaled = weight_estimated_scaled,
                       lower = lower_estimates, upper = upper_estimates,
                       text = text_estimated)
      shape <- 19
    }
    
    # Get y values for the estimated points 
    yEst <- rev(seq(.6, length(varES) - .4, 1))
    
    ranked <- rank(df$effectSize, ties.method="first")
    if(options$orderForest == "ascendingForest"){
      ord <- (length(varES) + 1) - ranked
      df$y <- ord
      yEst <- yEst[ranked]
    } 
    
    if(options$orderForest == "descendingForest"){
      ord <- ranked
      df$y <- ord
      yEst <- yEst[(length(varES) + 1) - ranked]
    } 

    # Create plot
    plot <-  ggplot2::ggplot(df,
                             ggplot2::aes(x = effectSize,
                                          y = y)) +
      # Add dotted vertical line at x = 0
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted")+

      # Add observed or estimated points with CI, and text
      ggplot2::geom_point(shape = shape, size = df$weight_scaled) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = df$lower, xmax = df$upper), height = .2) +
      
      # Add text on the right as secondary axis
      ggplot2::scale_y_continuous(breaks = c(df$y, 
                                             yDiamond),
                                  labels = c(as.character(df$studyLabels), model),
                                  sec.axis = ggplot2::sec_axis(~ .,
                                                               breaks = c(df$y, yDiamond),
                                                               labels = c(as.character(df$text), textDiamond))) +
      
      # Name of x axis
      ggplot2::xlab(expression("Effect size "*eta))
    
    if(options$forestPlot == "plotForestBoth"){
      dfBoth <- data.frame(effectSize = c(varES, mean_estimates),
                           y = c(df$y, yEst),
                           studyLabels = c(studyLabels, studyLabels),
                           weight_scaled = c(weight_scaled, weight_estimated_scaled), 
                           lower = c(lower, lower_estimates), upper = c(upper, upper_estimates),
                           text = c(text_observed, text_estimated),
                           g = rep(c("Observed", "Estimated"), each = length(varES)))
      
      plot <-  ggplot2::ggplot(dfBoth,
                               ggplot2::aes(x = effectSize,
                                            y = y)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
        ggplot2::geom_point(ggplot2::aes(shape = as.factor(dfBoth$g), colour = as.factor(dfBoth$g)), size = dfBoth$weight_scaled) +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = dfBoth$lower, xmax = dfBoth$upper, colour = as.factor(dfBoth$g),), height = .1, show.legend = FALSE) +
        # Add estimated axis labels
        ggplot2::scale_y_continuous(breaks = c(df$y,
                                               yDiamond),
                                    labels = c(as.character(df$studyLabels), model),
                                    sec.axis = ggplot2::sec_axis(~ .,
                                                                 breaks = c(df$y,
                                                                            yEst,
                                                                            yDiamond),
                                                                 labels = c(text_observed,
                                                                            text_estimated,
                                                                            textDiamond))) +
        ggplot2::scale_color_manual("", values = c("slategrey", "black"), labels = c("Estimated", "Observed")) +
        ggplot2::scale_shape_manual("", values = c(16, 15)) +
        ggplot2::guides(shape = ggplot2::guide_legend(reverse=TRUE, override.aes = list(size=3)),
                        colour = ggplot2::guide_legend(reverse=TRUE)) +
        # Adjust colour of secondary axis
        ggplot2::theme(axis.text.y.right = ggplot2::element_text(colour = c(rep(c("black", "slategrey"), each = nrow(df)), rep("black", 3)))) +
        ggplot2::xlab(expression("Effect size "*eta))
    }


    # Add jasp theme to plot
    plot <- themeJasp(plot,
                      yAxis = FALSE)
    
    # Add other theme elements (no y axis and aligning y axis labels)
    plot <- plot +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.line.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(hjust = 0),
                     axis.text.y.right = ggplot2::element_text(hjust = 1),
      )

    if(options$forestPlot == "plotForestBoth"){
      plot <- plot + ggplot2::theme(
        legend.position = c(1, 1),
        legend.justification=c(0, 0),
        plot.margin = ggplot2::unit(c(5, 1, 0.5, 0.5), "lines"),
        legend.title = ggplot2::element_blank()
      )
    }
    # Add the model diamond
    plot <- plot +
      ggplot2::geom_polygon(data = d, ggplot2::aes(x = x, y = y))

    # Add the diamonds of the other models for BMA or ordered analysis
    if(options$modelSpecification == "BMA" || options$modelSpecification == "CRE"){
      plot <- plot +
      ggplot2::geom_polygon(data = d.fixed, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_polygon(data = d.random, ggplot2::aes(x = x, y = y))
    }

    forestPlot$plotObject <- plot
    return()
  }
  
  .fillCumForest <- function(cumForestPlot, jaspResults, dataset, options, studyLabels, dependencies){

    rowResults <- .bmaSequentialResults(jaspResults, dataset, options, dependencies)
    
    meanMain   <- rowResults$mean
    lowerMain  <- rowResults$lowerMain
    upperMain  <- rowResults$upperMain
    
    text <- paste(sprintf('%.2f', meanMain),
                           " [",
                           sprintf('%.2f', lowerMain),
                           ", ",
                           sprintf('%.2f', upperMain),
                           "]",
                           sep = "")
    
    studyLabels[2] <- paste(studyLabels[1], "\n &", studyLabels[2])
    studyLabels    <- paste("+", studyLabels)
    studyLabels[1] <- "Prior"
    
    df <- data.frame(effectSize = meanMain, studyLabels = studyLabels, y = length(meanMain):1)
    
    plot <-  ggplot2::ggplot(df,
                              ggplot2::aes(x = effectSize,
                                           y = y))+
      # Add dotted vertical line at x = 0
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted")+
      
      # Add observed or estimated points with CI, and text
      ggplot2::geom_point(shape = 16, size = 4) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = lowerMain, xmax = upperMain), height = .2) +
      
      # Name of x axis
      ggplot2::xlab(expression("Overall effect size "*eta)) +
      
      ggplot2::scale_y_continuous(breaks = df$y,
                                  labels = studyLabels,
                                  sec.axis = ggplot2::sec_axis(~ .,
                                                               breaks = df$y,
                                                               labels = text))
      
    # Add jasp theme to plot
    plot <- themeJasp(plot,
                      yAxis = FALSE)
    
    # Add other theme elements (no y axis and aligning y axis labels)
    plot <- plot +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.line.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(hjust = 0),
                     axis.text.y.right = ggplot2::element_text(hjust = 1)
      )
    
      
    cumForestPlot$plotObject <- plot
    return()
      
  }
  
  .sequentialPlot <- function(jaspResults, dataset, options, ready, dependencies) {
    # Create empty plot
    seqContainer <- createJaspContainer(title = "Sequential Analysis")
    seqContainer$dependOn(dependencies)
    jaspResults[["seqContainer"]] <- seqContainer
    jaspResults[["seqContainer"]]$position <- 6
    

    # Check if ready
    if(!ready){
      return()
    }
    
    # Fill posterior plot effect size
    if(options$plotSequential){
      seqPlot <- createJaspPlot(plot = NULL, title = "Bayes factors", height = 400, width = 530)
      seqPlot$dependOn(c("plotSequential", "BF")) 
      seqContainer[["seqPlot"]] <- seqPlot
      .fillSeqPlot(seqPlot, jaspResults, dataset, options, dependencies)
    }
    
    if(options$plotSeqPM){
      seqPMPlot <- createJaspPlot(plot = NULL, title = "Posterior model probabilities", height = 400, width = 530)
      seqPMPlot$dependOn("plotSeqPM")
      .fillSeqPM(seqPMPlot, jaspResults, dataset, options, dependencies)
      seqContainer[["seqPMPlot"]] <- seqPMPlot
    }
    
  }
    
  .fillSeqPlot <- function(seqPlot, jaspResults, dataset, options, dependencies){
    
    rowResults <- .bmaSequentialResults(jaspResults, dataset, options, dependencies)
    
    BFs <- rowResults$BFs
    BFs[1] <- 1
    
    if(options$BF == "BF01") {
      BFs    <- 1/BFs
      bfType <- "BF01"
    } else if(options$BF == "logBF10") {
      bfType <- "LogBF10"
    } else {
      bfType <- "BF10" 
    }
    
    if(nrow(dataset) < 40) plotLineOrPoint <- "point" else plotLineOrPoint <- "line"
    
    df <- data.frame(x = 1:nrow(dataset), y = log(BFs))
    # }
    
    plot <- JASPgraphs::PlotRobustnessSequential(dfLines = df,
                                     plotLineOrPoint = plotLineOrPoint,
                                     xName = "Number of studies",
                                     BF = BFs[nrow(dataset)],
                                     bfType = bfType,
                                     hasRightAxis = TRUE
                                     )

    
    # plot + ggplot2::coord_flip()
    seqPlot$plotObject <- plot
    return()
  }
  
  .fillSeqPM <- function(seqPMPlot, jaspResults, dataset, options, dependencies){
    n     <- nrow(dataset)
    x     <- 0:n
    x     <- x[-2]
    dfPMP <- data.frame(prob = 0, g = rep(c("FE0", "FE1", "RE0", "RE1"), each = n))
    m     <- .bmaResultsState(jaspResults, dataset, options, dependencies)
    pM    <- m$prior_models
    
    dfPMP[c(1, 1 + n, 1 + 2*n, 1 + 3*n), 1] <- pM
    
    rowResults <- .bmaSequentialResults(jaspResults, dataset, options, dependencies)
    
    for(i in 2:nrow(dataset)){
      posterior_models <- rowResults$posterior_models[i]
    
      if(options$modelSpecification == "BMA" || options$modelSpecification == "CRE")
        dfPMP[c(i, i + n, i + 2*n, i + 3*n), 1] <- posterior_models
      
      if(options$modelSpecification == "FE")
        dfPMP[c(i, i + n), 1] <- posterior_models[1:2]
      
      if(options$modelSpecification == "RE")
        dfPMP[c(i, i + n), 1] <- posterior_models[3:4]
    }
    
    gridLines <- makeGridLines(x = 0,
                               y = c(0.25, 0.5, 0.75),
                               xend = nrow(dataset),
                               colour = rep("gray", 3),
                               linetype = rep("dashed", 3))
    
    df <- data.frame(x = x, y = dfPMP$prob, g = dfPMP$g)
    plot <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = g)) +
      gridLines + 
      ggplot2::geom_line(size = 1.5) +
      ggplot2::scale_y_continuous(limits = c(0,1.25), breaks = c(0, .25, .5, .75, 1)) +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2)) +
      ggplot2::theme(legend.spacing.x = ggplot2::unit(0.35, 'cm')) +
      ggplot2::labs(x = "Number of studies", y = "Posterior Model \n Probability") + 
      ggplot2::scale_colour_manual(labels = c(expression("Fixed H"[0]),
                                              expression("Fixed H"[1]),
                                              expression("Random H"[0]),
                                              expression("Random H"[1])),
                                   values = c("#56B4E9", "#009E73", "#CC79A7", "#D55E00")) 

    if(nrow(dataset) < 40) plot <- plot + ggplot2::geom_point(size = 3) 
    
    plot <- themeJasp(plot, legend.position = c(0.7, 1), legend.title = "none")

    seqPMPlot$plotObject <- plot
    return()                              
  }
  
  