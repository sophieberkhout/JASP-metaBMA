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
    
    dataset <- .readData(jaspResults, options)
    
    .bmaMainTable(jaspResults, dataset, options)
    
    .bmaEffectSizeTable(jaspResults, dataset, options)
    
    .bmaPostModelTable(jaspResults, dataset, options)
    
    .priorPlotES(jaspResults, dataset, options)
  }
  
  .readData <- function(jaspResults, options){
    varES <- options[["effectSize"]]
    varSE <- options[["standardError"]]
    study <- options[["studyLabels"]]
    if(varES == "") varES <- NULL
    if(varSE == "") varSE <- NULL
    if(study == "") study <- NULL
    variables.to.read <- c(varES, varSE, study)
    dataset <- .readDataSetToEnd(columns.as.numeric = variables.to.read)
    return(dataset)
  }
  
  .bmaResults <- function(jaspResults, dataset, options) {
    
    bmaResults <- createJaspState()
    jaspResults[["bmaResults"]] <- bmaResults
    bmaResults$dependOn(c("effectSize", "standardError", "studyLabels", "modelSpecification",
                          "priorH0FE", "priorH1FE", "priorH0RE", "priorH1RE",
                          "priorES", "cauchy", "normal", "t",
                          "informativeCauchyLocation", "informativeCauchyScale",
                          "lowerTruncCauchy", "upperTruncCauchy",
                          "informativeNormalMean", "informativeNormalStd",
                          "lowerTruncNormal", "upperTruncNormal",
                          "informativeTLocation", "informativeTScale", "informativeTDf",
                          "lowerTruncT", "upperTruncT"))
    
    varES <- dataset[, .v(options[["effectSize"]])]
    varSE <- dataset[, .v(options[["standardError"]])]
    
  
    # effect size prior
    lower <- -Inf
    upper <- Inf
    
    if(options$priorES == "cauchy"){
      family <- "t"
      param <- c(options$informativeCauchyLocation, options$informativeCauchyScale, 1)
      if(options$truncCauchy){
        lower <- options$lowerTruncCauchy
        upper <- options$upperTruncCauchy
      }
    }
    if(options$priorES == "normal"){
      family <- "norm"
      param <- c(options$informativeNormalMean, options$informativeNormalStd)
      if(options$truncCauchy){
        lower <- options$lowerTruncNormal
        upper <- options$upperTruncNormal
      }
    }
    if(options$priorES == "t"){
      family <- "t"
      param <- c(options$informativeTLocation, options$informativeTScale, options$informativeTDf)
      if(options$truncCauchy){
        lower <- options$lowerTruncT
        upper <- options$upperTruncT
      }
    }
    
    # heterogeneity prior
    if(options$priorSE == "inverseGamma"){
      family <- "invgamma"
      param <- c(options$inverseGammaShape, options$inverseGammaScale)
      lower <- -Inf
    }
    if(options$priorSE == "halfT"){
      family <- "t"
      param <- c(0, options$informativehalfTScale, options$informativehalfTDf)
      lower <- 0
    }
    
   
    d <- metaBMA::prior(family, param, lower, upper)
    tau <- metaBMA::prior(family, param, lower)
    
    modelPrior <- c(options[["priorH0FE"]], options[["priorH1FE"]], 
                options[["priorH0RE"]], options[["priorH1RE"]])
      
    results <- metaBMA::meta_bma(varES, varSE, prior = modelPrior, d = d, tau = tau)

    bmaResults$object <- results
    
    return()
  }
  
  
  .bmaMainTable <- function(jaspResults, dataset, options) {
    if (!is.null(jaspResults[["bmaTable"]])) return()
    
    bmaTable <- createJaspTable(title = "Bayesian Meta Analysis")
    jaspResults[["bmaTable"]] <- bmaTable
    bmaTable$dependOn(c("effectSize", "standardError", "confidenceInterval",
                        "modelSpecification", "priorH0FE", "priorH1FE", "priorH0RE", "priorH1RE",
                        "priorES", "cauchy", "normal", "t",
                        "informativeCauchyLocation", "informativeCauchyScale",
                        "lowerTruncCauchy", "upperTruncCauchy",
                        "informativeNormalMean", "informativeNormalStd",
                        "lowerTruncNormal", "upperTruncNormal",
                        "informativeTLocation", "informativeTScale", "informativeTDf",
                        "lowerTruncT", "upperTruncT"))
    
    bmaTable$addColumnInfo(name = "model", title = "", type = "string")
    bmaTable$addColumnInfo(name = "ES",   title = "Estimate",   type = "number", format = "dp:3")
    bmaTable$addColumnInfo(name = "SE",      title = "Standard Error",      type = "number", format = "dp:3")
    bmaTable$addColumnInfo(name = "BF", title = "BF\u2081\u2080", type = "number", format = "dp:3")
    bmaTable$addColumnInfo(name = "lb",      title = "Lower",      type = "number", format = "dp:3",
                           overtitle = "95% HDI")
    bmaTable$addColumnInfo(name = "ub",      title = "Upper",      type = "number", format = "dp:3",
                           overtitle = "95% HDI")
    
    ready <- options[["effectSize"]] != "" && options[["standardError"]] != ""
    
    if(!ready){
      return()
    }
    
    
    if (is.null(jaspResults[["bmaResults"]]))
      .bmaResults(jaspResults, dataset, options)
    
    m <- jaspResults[["bmaResults"]]$object
    
    if(options$modelSpecification == "BMA"){
      model <- c("Averaged \u03B7", "Fixed effects \u03B7", "Random effects \u03B7", "\u03C4")
      meanES <- c(m$estimates[, 1], m$meta$random$estimates[2, 1])
      meanSE <- c(m$estimates[, 2], m$meta$random$estimates[2, 2])
      lower <- c(m$estimates[, 6], m$meta$random$estimates[2, 6])
      upper <- c(m$estimates[, 7], m$meta$random$estimates[2, 7])
      BF <- c(m$inclusion$incl.BF, m$BF["fixed_H1", "fixed_H0"], 
              m$BF["random_H1", "random_H0"], m$BF["random_H1", "fixed_H1"])
    }
    else if(options$modelSpecification == "RE"){
      model <- c("Random effects \u03B7", "\u03C4")
      meanES <- m$meta$random$estimates[, 1]
      meanSE <- m$meta$random$estimates[, 2]
      lower <- m$meta$random$estimates[, 6]
      upper <- m$meta$random$estimates[, 7]
      BF <- c(m$BF["random_H1", "random_H0"], m$BF["random_H1", "fixed_H1"])
    }
    else if(options$modelSpecification == "FE"){
      model <- "Fixed effects \u03B7"
      meanES <- m$meta$fixed$estimates[, 1]
      meanSE <- m$meta$fixed$estimates[, 2]
      lower <- m$meta$fixed$estimates[, 6]
      upper <- m$meta$fixed$estimates[, 7]
      BF <- m$BF["fixed_H1", "fixed_H0"]
    }

    
    row <- data.frame(model = model, ES = meanES, SE = meanSE, BF = BF, lb = lower, ub = upper)
    bmaTable$addRows(row)
    
  }
  
  .bmaEffectSizeTable <- function(jaspResults, dataset, options) {
    if (!is.null(jaspResults[["esTable"]])) return()
    
    esTable <- createJaspTable(title = "Effect Sizes per Study")
    jaspResults[["esTable"]] <- esTable
    esTable$dependOn(c("effectSize", "standardError", "confidenceInterval",
                        "modelSpecification", "studyLabels",
                       "priorES", "cauchy", "normal", "t",
                       "informativeCauchyLocation", "informativeCauchyScale",
                       "lowerTruncCauchy", "upperTruncCauchy",
                       "informativeNormalMean", "informativeNormalStd",
                       "lowerTruncNormal", "upperTruncNormal",
                       "informativeTLocation", "informativeTScale", "informativeTDf",
                       "lowerTruncT", "upperTruncT"))
    
    esTable$addColumnInfo(name = "study", title = "", type = "string")
    esTable$addColumnInfo(name = "observedES",   title = "Observed",   type = "number", format = "dp:3",
                           overtitle = "Effect Sizes")
    
    ready <- options[["effectSize"]] != ""

    if(!ready){
      return()
    }

    m <- jaspResults[["bmaResults"]]$object
    varES <- dataset[, .v(options[["effectSize"]])]
    estimatedES <- rep(NA, length(varES))
   
    if(options$modelSpecification == "BMA" || options$modelSpecification == "RE"){
      esTable$addColumnInfo(name = "estimatedES",      title = "Estimated",      type = "number", format = "dp:3",
                            overtitle = "Effect Sizes")
      estimatedES <- rstan::summary(m$meta$random$stanfit_dstudy)$summary[3:(length(varES) + 2), "mean"]
    }
    
    
    
    

    if(options[["studyLabels"]] != ""){
      studyLabels <- dataset[, .v(options[["studyLabels"]])]
    } else {
      studyLabels <- 1:length(varES)
    }
    
    esTable$showSpecifiedColumnsOnly <- TRUE
    
    row <- data.frame(study = studyLabels, observedES = varES, estimatedES = estimatedES)
    esTable$addRows(row)
    
  }
  
  .bmaPostModelTable <- function(jaspResults, dataset, options) {
    if (!is.null(jaspResults[["postTable"]])) return()

    postTable <- createJaspTable(title = "Model Probabilities")
    jaspResults[["postTable"]] <- postTable
    postTable$dependOn(c("effectSize", "standardError", "confidenceInterval",
                       "modelSpecification",
                       "priorH0FE", "priorH1FE", "priorH0RE", "priorH1RE",
                       "priorES", "cauchy", "normal", "t",
                       "informativeCauchyLocation", "informativeCauchyScale",
                       "lowerTruncCauchy", "upperTruncCauchy",
                       "informativeNormalMean", "informativeNormalStd",
                       "lowerTruncNormal", "upperTruncNormal",
                       "informativeTLocation", "informativeTScale", "informativeTDf",
                       "lowerTruncT", "upperTruncT"))

    postTable$addColumnInfo(name = "model", title = "", type = "string")
    postTable$addColumnInfo(name = "priorProb",   title = "Prior",   type = "number", format = "dp:3")
    postTable$addColumnInfo(name = "postProb",   title = "Posterior",   type = "number", format = "dp:3")
    
    m <- jaspResults[["bmaResults"]]$object

    ready <- options[["effectSize"]] != "" && options[["standardError"]] != ""

    if(!ready){
      return()
    }

    if (is.null(jaspResults[["bmaResults"]]))
      .bmaResults(jaspResults, dataset, options)


    postTable$showSpecifiedColumnsOnly <- TRUE
    
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



    row <- data.frame(model = model, priorProb =  priorProb, postProb = postProb)
    postTable$addRows(row)

  }
  
  .priorPlotES <- function(jaspResults, dataset, options) {
    priorPlot <- createJaspPlot(plot = NULL, title = "Prior")
    
    priorPlot$dependOn(c("priorES", "cauchy", "normal", "t",
                          "informativeCauchyLocation", "informativeCauchyScale",
                          "lowerTruncCauchy", "upperTruncCauchy",
                          "informativeNormalMean", "informativeNormalStd",
                          "lowerTruncNormal", "upperTruncNormal",
                          "informativeTLocation", "informativeTScale", "informativeTDf",
                          "lowerTruncT", "upperTruncT"))
    
    jaspResults[["priorPlot"]] <- priorPlot
   # priorPlot$plotObject <- ggplot2::ggplot(<code>)
    
    .fillPriorPlotES(priorPlot, jaspResults, dataset, options)
   # return()
  }
  
  .fillPriorPlotES <- function(priorPlot, jaspResults, dataset, options){
    m <- jaspResults[["bmaResults"]]$object
    x <- seq(0, 1, .001)
    df <- data.frame(x = x)
    plot <- ggplot2::ggplot(df, ggplot2::aes(x)) +
            ggplot2::stat_function(fun = m$prior_d$fixed, n = 1000, size = 1) +
            ggplot2::labs(x = expression(eta), y = "Density")
    # df2 <- data.frame(x = x, y = x)
    # plot <- ggplot2::ggplot(df2, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
    plot <- themeJasp(plot)
    priorPlot$plotObject <- plot
    
    return()
  }
  
  
  
  