#
# Copyright (C) 2013-2015 University of Amsterdam
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

BayesianMetaAnalysis <- function(jaspResults, dataset, options, state=NULL) {
  jaspResults$title <- "Bayesian Meta-Analysis"
  
  dataset <- .readData(jaspResults, options)
  
  .bmaMainTable(jaspResults, dataset, options)
}

.readData <- function(jaspResults, options){
  varES <- options[["effectSize"]]
  varSE <- options[["standardError"]]
  if(varES == "") varES <- NULL
  if(varSE == "") varSE <- NULL
  variables.to.read <- c(varES, varSE)
  dataset <- .readDataSetToEnd(columns.as.numeric = variables.to.read)
  return(dataset)
}

.bmaMainTable <- function(jaspResults, dataset, options) {
  if (!is.null(jaspResults[["bmaTable"]])) return()
  
  bmaTable <- createJaspTable(title = "BMA test")
  jaspResults[["bmaTable"]] <- bmaTable
  bmaTable$dependOnOptions(c("effectSize", "standardError", "confidenceInterval",
                                  "modelSpecification"))
  
  bmaTable$addColumnInfo(name = "ES",   title = "Effect Size",   type = "number", format = "dp:3")
  bmaTable$addColumnInfo(name = "SE",      title = "Standard Error",      type = "number", format = "dp:3")
  
  ready <- options[["effectSize"]] != "" && options[["standardError"]] != ""
  
  if(!ready){
    return()
  }
  
  varES <- dataset[, .v(options[["effectSize"]])]
  varSE <- dataset[, .v(options[["standardError"]])]
  
  meanES <- mean(varES)
  meanSE <- mean(varSE)
  
  row <- data.frame(ES = meanES, SE = meanSE)
  bmaTable$addRows(row)
  
}

