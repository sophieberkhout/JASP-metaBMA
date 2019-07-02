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
  
  ready <- options[["effectSize"]] != "" && options[["standardError"]] != "" 
  
  dependencies <- c("effectSize", "standardError", "studyLabels", "modelSpecification",
                    "priorH0FE", "priorH1FE", "priorH0RE", "priorH1RE",
                    "priorES", "cauchy", "normal", "t",
                    "informativeCauchyLocation", "informativeCauchyScale",
                    "lowerTruncCauchy", "upperTruncCauchy",
                    "informativeNormalMean", "informativeNormalStd",
                    "lowerTruncNormal", "upperTruncNormal",
                    "informativeTLocation", "informativeTScale", "informativeTDf",
                    "lowerTruncT", "upperTruncT",
                    "priorSE", "inverseGamma", "inverseGammaShape", "inverseGammaScale",
                    "halfT", "informativehalfTScale", "informativehalfTDf",
                    "BFComputation", "integration", "bridgeSampling", "iterBridge",
                    "iterMCMC", "chainsMCMC")
  
  dataset <- .readData(jaspResults, options)
    
  .bmaMainTable(jaspResults, dataset, options, ready, dependencies)
    
  .bmaEffectSizeTable(jaspResults, dataset, options, ready, dependencies)
    
  .bmaPostModelTable(jaspResults, dataset, options, ready, dependencies)
    
  if(options$plotPrior){
    .priorPlot(jaspResults, dataset, options, ready)
  }
  
  if(options$plotPosterior){  
    .priorAndPosteriorPlot(jaspResults, dataset, options, ready)
  }
  
  .forestPlot(jaspResults, dataset, options, ready, dependencies)
  
}
  
  .readData <- function(jaspResults, options){
    varES <- options[["effectSize"]]
    varSE <- options[["standardError"]]
    #varCI <- options[["confidenceInterval"]]
    study <- options[["studyLabels"]]
    if(varES == "") varES <- NULL
    if(varSE == "") varSE <- NULL
    #if(varCI == "") varCI <- NULL
    if(study == "") study <- NULL
    variables.to.read <- c(varES, varSE, study)
    dataset <- .readDataSetToEnd(columns.as.numeric = variables.to.read)
    return(dataset)
  }
  
  .bmaResults <- function(jaspResults, dataset, options, dependencies) {
    
    bmaResults <- createJaspState()
    jaspResults[["bmaResults"]] <- bmaResults
    bmaResults$dependOn(dependencies)
    
    varES <- dataset[, .v(options[["effectSize"]])]
    varSE <- dataset[, .v(options[["standardError"]])]
    
    # if(options[["standardError"]] != ""){
    #   varSE <- dataset[, .v(options[["standardError"]])]
    # } else if(options[["confidenceInterval"]] != ""){
    #   confidenceInterval <- dataset[, .v(options[["confidenceInterval"]])] 
    #   varSE <- (confidenceInterval)/3.92
    # }
    
  
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
    
    # estimation settings
    iter <- options$iterMCMC
    chains <- options$chainsMCMC
    
    # bayes factor computation
    if(options$BFComputation == "integration"){
      logml <- "integrate"
      logml_iter <- 5000
    } else if(options$BFComputation == "bridgeSampling"){
      logml <- "stan"
      logml_iter <- options$iterBridge
    }
    
    d <- metaBMA::prior(familyES, paramES, lowerES, upperES)
    tau <- metaBMA::prior(familySE, paramSE, lowerSE)
    
    modelPrior <- c(options[["priorH0FE"]], options[["priorH1FE"]], 
                options[["priorH0RE"]], options[["priorH1RE"]])
      
    results <- metaBMA::meta_bma(varES, varSE, 
                                 prior = modelPrior, 
                                 d = d, 
                                 tau = tau,
                                 logml = logml,
                                 logml_iter = logml_iter,
                                 iter = iter,
                                 chains = chains)

    bmaResults$object <- results
    
    return()
  }
  
  
  .bmaMainTable <- function(jaspResults, dataset, options, ready, dependencies) {
    if (!is.null(jaspResults[["bmaTable"]])) return()
    
    bmaTable <- createJaspTable(title = "Bayesian Meta Analysis")
    jaspResults[["bmaTable"]] <- bmaTable
    bmaTable$dependOn(dependencies)
    
    bmaTable$addColumnInfo(name = "model", title = "", type = "string")
    bmaTable$addColumnInfo(name = "ES",   title = "Estimate",   type = "number", format = "dp:3")
    bmaTable$addColumnInfo(name = "SE",      title = "Standard Error",      type = "number", format = "dp:3")
    bmaTable$addColumnInfo(name = "BF", title = "BF\u2081\u2080", type = "number", format = "dp:3")
    bmaTable$addColumnInfo(name = "lb",      title = "Lower",      type = "number", format = "dp:3",
                           overtitle = "95% HDI")
    bmaTable$addColumnInfo(name = "ub",      title = "Upper",      type = "number", format = "dp:3",
                           overtitle = "95% HDI")
    
    if(!ready){
      return()
    }
    
    
    if (is.null(jaspResults[["bmaResults"]]))
      .bmaResults(jaspResults, dataset, options, dependencies)
    
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
  
  .bmaEffectSizeTable <- function(jaspResults, dataset, options, ready, dependencies) {
    if (!is.null(jaspResults[["esTable"]])) return()
    
    esTable <- createJaspTable(title = "Effect Sizes per Study")
    jaspResults[["esTable"]] <- esTable
    esTable$dependOn(dependencies)
    
    esTable$addColumnInfo(name = "study", title = "", type = "string")
    esTable$addColumnInfo(name = "observedES",   title = "Observed",   type = "number", format = "dp:3",
                           overtitle = "Effect Sizes")
    
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
      studyLabels <- paste("Study", 1:length(varES))
    }
    
    esTable$showSpecifiedColumnsOnly <- TRUE
    
    row <- data.frame(study = studyLabels, observedES = varES, estimatedES = estimatedES)
    esTable$addRows(row)
    
  }
  
  .bmaPostModelTable <- function(jaspResults, dataset, options, ready, dependencies) {
    if (!is.null(jaspResults[["postTable"]])) return()

    postTable <- createJaspTable(title = "Model Probabilities")
    jaspResults[["postTable"]] <- postTable
    postTable$dependOn(dependencies)

    postTable$addColumnInfo(name = "model", title = "", type = "string")
    postTable$addColumnInfo(name = "priorProb",   title = "Prior",   type = "number", format = "dp:3")
    postTable$addColumnInfo(name = "postProb",   title = "Posterior",   type = "number", format = "dp:3")
    
    m <- jaspResults[["bmaResults"]]$object

    if(!ready){
      return()
    }

    if (is.null(jaspResults[["bmaResults"]]))
      .bmaResults(jaspResults, dataset, options, dependencies)


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
  
  .priorPlot <- function(jaspResults, dataset, options, ready) {
    priorContainer <- createJaspContainer(title = "Prior")
    priorContainer$dependOn(c("plotPrior"))
    jaspResults[["priorContainer"]] <- priorContainer
    
    priorPlotES <- createJaspPlot(plot = NULL, title = "Effect size")
    
    priorPlotES$dependOn(c("priorES", "cauchy", "normal", "t",
                          "informativeCauchyLocation", "informativeCauchyScale",
                          "lowerTruncCauchy", "upperTruncCauchy",
                          "informativeNormalMean", "informativeNormalStd",
                          "lowerTruncNormal", "upperTruncNormal",
                          "informativeTLocation", "informativeTScale", "informativeTDf",
                          "lowerTruncT", "upperTruncT"))
    
    
   # jaspResults[["priorPlot"]] <- priorPlot
   # priorPlot$plotObject <- ggplot2::ggplot(<code>)
    
    if(!ready){
      return()
    }
    
    .fillPriorPlotES(priorPlotES, jaspResults, dataset, options)
    
    priorContainer[["ES"]] <- priorPlotES

    if(options$modelSpecification == "BMA" || options$modelSpecification == "RE"){
      priorPlotSE <- createJaspPlot(plot = NULL, title = "Heterogeneity")
      priorPlotSE$dependOn(c("priorSE", "inverseGamma", "inverseGammaShape", "inverseGammaScale",
                             "halfT", "informativehalfTScale", "informativehalfTDf"))
      .fillPriorPlotSE(priorPlotSE, jaspResults, dataset, options)
      priorContainer[["SE"]] <- priorPlotSE
    }
    
    
    
    
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
            ggplot2::labs(x = expression(eta), y = "Density") +
            ggplot2::geom_vline(xintercept = 0, linetype = "dotted")
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
      ggplot2::labs(x = expression(tau), y = "Density") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted")
    # df2 <- data.frame(x = x, y = x)
    # plot <- ggplot2::ggplot(df2, ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
    plot <- themeJasp(plot)
    priorPlotSE$plotObject <- plot
    
    return()
  }
  
  .priorAndPosteriorPlot <- function(jaspResults, dataset, options, ready) {
    postContainer <- createJaspContainer(title = "Prior and Posteriors")
    postContainer$dependOn(c("plotPosterior"))
    jaspResults[["postContainer"]] <- postContainer
    
    postPlotES <- createJaspPlot(plot = NULL, title = "Effect size")
    
    
    postPlotES$dependOn(c("priorES", "cauchy", "normal", "t",
                           "informativeCauchyLocation", "informativeCauchyScale",
                           "lowerTruncCauchy", "upperTruncCauchy",
                           "informativeNormalMean", "informativeNormalStd",
                           "lowerTruncNormal", "upperTruncNormal",
                           "informativeTLocation", "informativeTScale", "informativeTDf",
                           "lowerTruncT", "upperTruncT"))
    
    
    if(!ready){
      return()
    }
    
    .fillPostPlot(postPlotES, jaspResults, dataset, options, type = "ES")
    
    postContainer[["ES"]] <- postPlotES

    if(options$modelSpecification == "BMA" || options$modelSpecification == "RE"){
      postPlotSE <- createJaspPlot(plot = NULL, title = "Heterogeneity")
      postPlotSE$dependOn(c("priorSE", "inverseGamma", "inverseGammaShape", "inverseGammaScale",
                            "halfT", "informativehalfTScale", "informativehalfTDf"))
      postContainer[["SE"]] <- postPlotSE
      .fillPostPlot(postPlotSE, jaspResults, dataset, options, type = "SE")
   }
    
    
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
      ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
      # ggplot2::geom_vline(ggplot2::aes(xintercept = mode(x))) + 
      # theme_classic() +
      # ggplot2::theme(legend.position = "top",
      #       legend.spacing.x = ggplot2::unit(0.3, "cm"),
      #       legend.title = ggplot2::element_blank()) +
      ggplot2::labs(x = xlab, y = "Density") +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol=1,nrow=2,byrow=TRUE)) +
      ggplot2::guides(linetype = ggplot2::guide_legend(ncol=1,nrow=2,byrow=TRUE))
    
    
    mPostFixed <- m$meta$fixed$posterior_d
    mPostRandom <- m$meta$random$posterior_d
    
    plot <- plot + ggplot2::stat_function(fun = mPostFixed, n = 1000, size = 1, ggplot2::aes(colour = "Fixed")) +
      ggplot2::stat_function(fun = mPostRandom, n = 1000, size = 1, ggplot2::aes(colour = "Random"))
    
    plot <- themeJasp(plot, legend.position = "bottom", legend.title = "none")
    
    postPlot$plotObject <- plot
    
    return()
  }
  
  .forestPlot <- function(jaspResults, dataset, options, ready, dependencies) {
    forestContainer <- createJaspContainer(title = "Forest Plot")
    forestContainer$dependOn(c("plotForestObserved"))
    jaspResults[["forestContainer"]] <- forestContainer
    
    height <- (nrow(dataset) + 2) * 40
    
    forestPlot <- createJaspPlot(plot = NULL, title = "Forest Plot", height = height, width = 400)
    
    
    forestPlot$dependOn(dependencies)
    
    
    if(!ready){
      return()
    }
    
    .fillForestPlot(forestPlot, jaspResults, dataset, options)
    
    forestContainer[["forestPlot"]] <- forestPlot
    
  }
  
  .fillForestPlot <- function(forestPlot, jaspResults, dataset, options){
    m <- jaspResults[["bmaResults"]]$object
    
    varES <- dataset[, .v(options[["effectSize"]])]
    varSE <- dataset[, .v(options[["standardError"]])]
    if(options[["studyLabels"]] != ""){
      studyLabels <- dataset[, .v(options[["studyLabels"]])]
    } else {
      studyLabels <- paste("Study", 1:length(varSE))
    }
    df <- data.frame(effectSize = varES, studyLabels = studyLabels)
    
    # assign weights for the point sizes
    weight <- 1/varSE^2
    weight_scaled <- ((4-1)*(weight - min(weight)))/(max(weight) - min(weight)) + 2
    
    # text for next to the observed points
    ci <- .95
    lower <- varES - qnorm((ci+1)/2) * varSE
    upper <- varES + qnorm((ci+1)/2) * varSE
    
    text_observed <- format(paste(format(round(varES, 2),
                                         nsmall = 2, trim = T), " [",
                                  format(round(lower, 2), nsmall = 2, trim = T), ", ",
                                  format(round(upper, 2), nsmall = 2, trim = T), "]",
                                  sep = ""),
                            width = 20, justify = "right") # super ugly code
    
    xlim <- c(min(lower), max(upper))
    ylim <- c(0, nrow(df))
    
    # number of characters left of plot (required for scaling. right side: nchar = 20)
    nchar_labels <- max(nchar(as.character(studyLabels)))
    shift_right <- max(xlim)+ diff(xlim)/2 * sqrt(nchar_labels / 20) + 3 # DANIEL: just a heuristic for scaling

    # create dataframe for model diamond
    p.left <- m$estimates[1, 3]
    p.right <- m$estimates[1, 5]
    p.top <- -0.5 + m$estimates[1, 1] - m$estimates[1, 3]
    p.bottom <- -0.5 + m$estimates[1, 1] - m$estimates[1, 5]
    
    d <- data.frame(x = c(p.left, m$estimates[1, 1],
                          p.right, m$estimates[1, 1]),
                    y = c(-0.5, p.top, -0.5, p.bottom))
    
    # text for next to model diamond
    text_overall <- format(paste(format(round(m$estimates[1, 1], 2),
                                        nsmall = 2, trim = T), " [",
                                 format(round(m$estimates[1, 3], 2),
                                        nsmall = 2, trim = T), ", ",
                                 format(round(m$estimates[1, 5], 2),
                                        nsmall = 2, trim = T), "]",
                                 sep = ""),
                           width = 20, justify = "right")
    
    # get data & model specific x limits
    # if(class(m) == "meta_bma" || m$model == "random"){
    #   xlim <- c(min(c(lower, lower_estimates)), max(c(upper, upper_estimates)))
    # } else {
      xlim <- c(min(lower), max(upper))
    # }
    
    # get data & model specific y limits
    # if(class(m) != "m_bma"){
      ylim <- c(-0.5, nrow(df))
    # } else {
    #   ylim <- c(-1.5, nrow(data))
    # }
      
    plot <-  ggplot2::ggplot(df,
                             ggplot2::aes(x = effectSize,
                                          y = as.numeric(reorder(studyLabels, -effectSize))))+
      # add dotted vertical line at x = 0
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted")+
      
      # add observed points and text
      ggplot2::geom_point(shape = 15, size = weight_scaled) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper), height = .2) +
      # ggplot2::annotate("text", label = text_observed,
      #          x = shift_right, y = reorder(df$studyLabels, -df$effectSize),
      #          hjust = "right", size = 6) +
      ggplot2::scale_y_continuous(breaks = 1:length(studyLabels),
                                  labels = studyLabels,
                                  sec.axis = ggplot2::sec_axis(~ .,
                                                               breaks = 1:length(text_observed),
                                                               labels = text_observed)) +
      ggplot2::geom_segment(x = -99, xend = 99, y = 0, yend = 0) +
      # ggplot2::theme(axis.title.y = ggplot2::element_blank(),
      #       axis.line.y = ggplot2::element_blank(),
      #       axis.ticks.y = ggplot2::element_blank(),
      #       plot.margin = ggplot2::unit(c(1,10,1,1), "lines"), # margin on the right for text
      #       panel.grid.major = ggplot2::element_blank(),
      #       panel.grid.minor = ggplot2::element_blank(),
      #       axis.ticks.x = ggplot2::element_line(size = .3),
      #       axis.ticks.length = ggplot2::unit(-1.4, "mm"), # ticks on inside of plot
      #       axis.text.x = ggplot2::element_text(margin = ggplot2::unit(c(2.5, 0, 0, 0), "mm")),
      #       axis.text.y = ggplot2::element_text(hjust = 1)
      # ) +
      ggplot2::xlab("Effect Size") +
      
      # focus x and y axis on range of interest and clip = 'off' to add the text on the right
      ggplot2::coord_cartesian(xlim = xlim,
                      ylim = ylim,
                      clip = 'off')
    
    
    
    plot <- themeJasp(plot,
                      yAxis = FALSE)
    
    plot <- plot +       
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
            axis.line.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            # plot.margin = ggplot2::unit(c(1,10,1,1), "lines"), # margin on the right for text
            # panel.grid.major = ggplot2::element_blank(),
            # panel.grid.minor = ggplot2::element_blank(),
            # axis.ticks.x = ggplot2::element_line(size = .3),
            # axis.ticks.length = ggplot2::unit(-1.4, "mm"), # ticks on inside of plot
            # axis.text.x = ggplot2::element_text(margin = ggplot2::unit(c(2.5, 0, 0, 0), "mm")),
            axis.text.y = ggplot2::element_text(hjust = 1)
      )
      
    plot <- plot +
      ggplot2::geom_polygon(data = d, ggplot2::aes(x = x, y = y)) +
      ggplot2::annotate("text", label = text_overall,
               x = Inf, y = -0.5, hjust = 0, size = 6)+
      ggplot2::geom_text(ggplot2::aes(x = -Inf, y = -0.5, label = "FE model"),
                hjust = 1, size = 6)
    
    
    forestPlot$plotObject <- plot
    
    return()
  }
    
  
  
  