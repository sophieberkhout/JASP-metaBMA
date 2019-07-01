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
    
    if(options$plotPrior){
      .priorPlot(jaspResults, dataset, options)
    }
    
    .priorAndPosteriorPlot(jaspResults, dataset, options)
    
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
                          "lowerTruncT", "upperTruncT",
                          "priorSE", "inverseGamma", "inverseGammaShape", "inverseGammaScale",
                          "halfT", "informativehalfTScale", "informativehalfTDf"))
    
    varES <- dataset[, .v(options[["effectSize"]])]
    varSE <- dataset[, .v(options[["standardError"]])]
    
  
    # effect size prior
    lowerES <- -Inf
    upperES <- Inf
    
    if(options$priorES == "cauchy"){
      familyES <- "t"
      paramES <- c(options$informativeCauchyLocation, options$informativeCauchyScale, 1)
      if(options$truncCauchy){
        lowerES <- options$lowerTruncCauchy
        upperES <- options$upperTruncCauchy
      }
    }
    if(options$priorES == "normal"){
      familyES <- "norm"
      paramES <- c(options$informativeNormalMean, options$informativeNormalStd)
      if(options$truncCauchy){
        lowerES <- options$lowerTruncNormal
        upperES <- options$upperTruncNormal
      }
    }
    if(options$priorES == "t"){
      familyES <- "t"
      paramES <- c(options$informativeTLocation, options$informativeTScale, options$informativeTDf)
      if(options$truncCauchy){
        lowerES <- options$lowerTruncT
        upperES <- options$upperTruncT
      }
    }
    
    # heterogeneity prior
    if(options$priorSE == "inverseGamma"){
      familySE <- "invgamma"
      paramSE <- c(options$inverseGammaShape, options$inverseGammaScale)
      lowerSE <- -Inf
    }
    if(options$priorSE == "halfT"){
      familySE <- "t"
      paramSE <- c(0, options$informativehalfTScale, options$informativehalfTDf)
      lowerSE <- 0
    }
    
    d <- metaBMA::prior(familyES, paramES, lowerES, upperES)
    tau <- metaBMA::prior(familySE, paramSE, lowerSE)
    
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
                        "lowerTruncT", "upperTruncT",
                        "priorSE", "inverseGamma", "inverseGammaShape", "inverseGammaScale",
                        "halfT", "informativehalfTScale", "informativehalfTDf"))
    
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
                       "lowerTruncT", "upperTruncT",
                       "priorSE", "inverseGamma", "inverseGammaShape", "inverseGammaScale",
                       "halfT", "informativehalfTScale", "informativehalfTDf"))
    
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
                       "lowerTruncT", "upperTruncT",
                       "priorSE", "inverseGamma", "inverseGammaShape", "inverseGammaScale",
                       "halfT", "informativehalfTScale", "informativehalfTDf"))

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
  
  .priorPlot <- function(jaspResults, dataset, options) {
    priorContainer <- createJaspContainer(title = "Priors")
    priorContainer$dependOn(c("plotPrior"))
    jaspResults[["priorContainer"]] <- priorContainer
    
    priorPlotES <- createJaspPlot(plot = NULL, title = "Effect size")
    priorPlotSE <- createJaspPlot(plot = NULL, title = "Heterogeneity")
    
    
    priorPlotES$dependOn(c("priorES", "cauchy", "normal", "t",
                          "informativeCauchyLocation", "informativeCauchyScale",
                          "lowerTruncCauchy", "upperTruncCauchy",
                          "informativeNormalMean", "informativeNormalStd",
                          "lowerTruncNormal", "upperTruncNormal",
                          "informativeTLocation", "informativeTScale", "informativeTDf",
                          "lowerTruncT", "upperTruncT"))
    
    priorPlotSE$dependOn(c("priorSE", "inverseGamma", "inverseGammaShape", "inverseGammaScale",
                           "halfT", "informativehalfTScale", "informativehalfTDf"))
    
   # jaspResults[["priorPlot"]] <- priorPlot
   # priorPlot$plotObject <- ggplot2::ggplot(<code>)
    
    .fillPriorPlotES(priorPlotES, jaspResults, dataset, options)
    .fillPriorPlotSE(priorPlotSE, jaspResults, dataset, options)
    
    
    priorContainer[["ES"]] <- priorPlotES
    priorContainer[["SE"]] <- priorPlotSE
    
   # return()
  }
  
  .fillPriorPlotES <- function(priorPlotES, jaspResults, dataset, options){
    m <- jaspResults[["bmaResults"]]$object
    if(options$priorES == "cauchy"){
      mean <- options$informativeCauchyLocation
    } else if(options$priorES == "normal"){
      mean <- options$informativeNormalMean
    } else {
      mean <- options$informativeTLocation
    }
    
    xlimLeft <- mean - 3
    xlimRight <- mean + 3
    x <- seq(xlimLeft, xlimRight, .001)
    df <- data.frame(x = x)
    plot <- ggplot2::ggplot(df, ggplot2::aes(x)) +
            ggplot2::stat_function(fun = m$prior_d$fixed, n = 1000, size = 1) +
            ggplot2::labs(x = expression(eta), y = "Density")
    # df2 <- data.frame(x = x, y = x)
    # plot <- ggplot2::ggplot(df2, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
    plot <- themeJasp(plot)
    priorPlotES$plotObject <- plot
    
    return()
  }
  
  .fillPriorPlotSE <- function(priorPlotSE, jaspResults, dataset, options){
    m <- jaspResults[["bmaResults"]]$object
    x <- seq(0, 3, .001)
    df <- data.frame(x = x)
    plot <- ggplot2::ggplot(df, ggplot2::aes(x)) +
      ggplot2::stat_function(fun = m$meta$random$prior_tau, n = 1000, size = 1) +
      ggplot2::labs(x = expression(tau), y = "Density")
    # df2 <- data.frame(x = x, y = x)
    # plot <- ggplot2::ggplot(df2, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
    plot <- themeJasp(plot)
    priorPlotSE$plotObject <- plot
    
    return()
  }
  
  .priorAndPosteriorPlot <- function(jaspResults, dataset, options) {
    postContainer <- createJaspContainer(title = "Posteriors")
    postContainer$dependOn(c("plotPosterior"))
    jaspResults[["postContainer"]] <- postContainer
    
    postPlotES <- createJaspPlot(plot = NULL, title = "Effect size")
    postPlotSE <- createJaspPlot(plot = NULL, title = "Heterogeneity")
    
    
    postPlotES$dependOn(c("priorES", "cauchy", "normal", "t",
                           "informativeCauchyLocation", "informativeCauchyScale",
                           "lowerTruncCauchy", "upperTruncCauchy",
                           "informativeNormalMean", "informativeNormalStd",
                           "lowerTruncNormal", "upperTruncNormal",
                           "informativeTLocation", "informativeTScale", "informativeTDf",
                           "lowerTruncT", "upperTruncT"))
    
    postPlotSE$dependOn(c("priorSE", "inverseGamma", "inverseGammaShape", "inverseGammaScale",
                           "halfT", "informativehalfTScale", "informativehalfTDf"))
    
    .fillPostPlot(postPlotES, jaspResults, dataset, options, type = "ES")
    .fillPostPlot(postPlotSE, jaspResults, dataset, options, type = "SE")
    
    
    postContainer[["ES"]] <- postPlotES
    postContainer[["SE"]] <- postPlotSE
    
    # return()
  }
  
  .fillPostPlot <- function(postPlot, jaspResults, dataset, options, type){
    m <- jaspResults[["bmaResults"]]$object

    df <- data.frame(x = c(-1, 1), l = c("Prior", "Posterior"))
    
    if(type == "ES"){
      mPrior <- m$prior_d$fixed
      mPost <- m$posterior_d
      int <- c(m$estimates["averaged", "2.5%"], m$estimates["averaged", "97.5%"])
      e <- m$estimates["averaged", "mean"]
      xlab <- expression(eta)
    } else if(type == "SE"){
      mPrior <- m$meta$random$prior_tau
      mPost <- m$meta$random$posterior_tau
      int <- c(m$meta$random$estimates["tau", "2.5%"], m$meta$random$estimates["tau", "97.5%"])
      e <- m$meta$random$estimates["tau", "mean"]
      xlab <- expression(tau)
    }
    
    plot <- ggplot2::ggplot(df, ggplot2::aes(x)) +
      ggplot2::stat_function(fun = mPost, n = 1000, size = 1, ggplot2::aes(linetype = "Posterior")) +
      ggplot2::stat_function(fun = mPrior, n = 1000, size = 1, ggplot2::aes(linetype = "Prior")) +
      ggplot2::scale_linetype_manual(values = c("solid", "dotted")) +
      ggplot2::stat_function(fun = mPost, 
                    xlim = int,
                    geom = "area", alpha = 0.2, show.legend = F) +
      ggplot2::geom_segment(ggplot2::aes(x = e, xend = e, y = 0, yend = mPost(e)-.01), size = 0.5) +
      # ggplot2::geom_vline(ggplot2::aes(xintercept = mode(x))) + 
      # theme_classic() +
      # ggplot2::theme(legend.position = "top",
      #       legend.spacing.x = ggplot2::unit(0.3, "cm"),
      #       legend.title = ggplot2::element_blank()) +
      ggplot2::labs(x = xlab, y = "Density")
    
    mPostFixed <- m$meta$fixed$posterior_d
    mPostRandom <- m$meta$random$posterior_d
    
    plot <- plot + ggplot2::stat_function(fun = mPostFixed, n = 1000, linetype = "dashed", size = 1, ggplot2::aes(colour = "Fixed")) +
      ggplot2::stat_function(fun = mPostRandom, n = 1000, linetype = "dashed", size = 1, ggplot2::aes(colour = "Random"))
    
    plot <- themeJasp(plot)
    
    postPlot$plotObject <- plot
    
    return()
  }
  

  
  
  