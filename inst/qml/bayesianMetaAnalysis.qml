//
// Copyright (C) 2013-2018 University of Amsterdam
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public
// License along with this program.  If not, see
// <http://www.gnu.org/licenses/>.
//
import QtQuick 2.8
import QtQuick.Layouts 1.3
import JASP.Controls 1.0

Form
{
    id: form

    VariablesForm
    {
      AvailableVariablesList {name: "variablesList"}
      AssignedVariablesList
      {
          name: "effectSize"
          title: qsTr("Effect Size")
          singleVariable: true
      }
      AssignedVariablesList
      {
          enabled: confidenceInterval.count < 1
          id: standardError
          name: "standardError"
          title: qsTr("Effect Size Standard Error")
          singleVariable: true
      }
      AssignedVariablesList
      {
          enabled: standardError.count == 0
          id: confidenceInterval
          name: "confidenceInterval"
          title: qsTr("95% CI Lower and Upper Bound")
          singleVariable: true
          listViewType: "Pairs"
      }
      AssignedVariablesList
      {
          name: "studyLabels"
          title: qsTr("Study Labels")
          singleVariable: true
      }
    }

    GridLayout
    {
      columns: 2
        RadioButtonGroup
        {
          name: "modelSpecification"
          title: qsTr("Model")
            RadioButton
            {
              value: "FE"
              label: qsTr("Fixed effects")
              id: checkFE
            }
            RadioButton
            {
              value: "RE"
              label: qsTr("Random effects")
              id: checkRE
            }
            RadioButton
            {
              value: "BMA"
              label: qsTr("Model averaging")
              id: checkBMA
              checked: true
            }
            RadioButton
            {
              value: "CRE"
              label: qsTr("Constrained random effects")
              id: checkCRE
                RadioButtonGroup
                {
                  name: "direction"
                  columns: 2
                    RadioButton
                    {
                      value: "allPos"
                      label: qsTr("All positive")
                      checked: true
                    }
                    RadioButton
                    {
                      value: "allNeg"
                      label: qsTr("All negative")
                    }
                }
            }
        }

        Group
        {
          enabled: checkFE.checked || checkRE.checked || checkBMA.checked
          title: qsTr("Prior model probability")
            Group
            {
              enabled: checkFE.checked || checkBMA.checked
              title: qsTr("Fixed effects")
                DoubleField
                {
                  name: "priorH0FE"
                  label: "H\u2080"
                  id: "priorH0FE"
                  defaultValue:
                  {
                    if (checkFE.checked) 0.5
                    else if (checkRE.checked) 0
                    else if (checkBMA.checked) 0.25
                    else 0
                  }
                }
                DoubleField
                {
                  name: "priorH1FE"
                  label: "H\u2081"
                  id: "priorH1FE"
                  defaultValue:
                  {
                    if(checkFE.checked){0.5}
                    else if(checkRE.checked){0}
                    else if(checkBMA.checked){0.25}
                    else {0}
                  }
                }
            }
            Group
            {
              title: qsTr("Random effects")
              enabled: checkRE.checked || checkBMA.checked
                DoubleField
                {
                  name: "priorH0RE"
                  label: "H\u2080"
                  id: "priorH0RE"
                  defaultValue:
                  {
                    if(checkFE.checked){0}
                    else if(checkRE.checked){0.5}
                    else if(checkBMA.checked){0.25}
                    else {0}
                  }
                }
                DoubleField
                {
                  name: "priorH1RE"
                  label: "H\u2081"
                  id: "priorH1RE"
                  defaultValue:
                  {
                    if(checkFE.checked){0}
                    else if(checkRE.checked){0.5}
                    else if(checkBMA.checked){0.25}
                    else {0}
                  }
                }
            }
        }

        Section {
            title: qsTr("Prior")
            columns: 1

            RadioButtonGroup
            				{
                      title: qsTr("Effect size")
            					name: "priorES"
            					RadioButton
            					{
            						label: qsTr("Cauchy"); name: "cauchy"; checked: true; childrenOnSameRow: true; id: cauchyInformative
            						DoubleField { label: qsTr("location:"); name: "informativeCauchyLocation"; visible: cauchyInformative.checked;
                                      defaultValue: 0; negativeValues: true }
            						DoubleField { label: qsTr("scale:"); name: "informativeCauchyScale"; visible: cauchyInformative.checked;
                                      defaultValue: 0.707; fieldWidth: 50}
                        CheckBox
                        {
                          name: "truncCauchy"; label: qsTr("truncation:"); visible: cauchyInformative.checked; childrenOnSameRow: true
                          DoubleField { label: qsTr("from"); name: "lowerTruncCauchy"; visible: cauchyInformative.checked; id: lowerTC;
                                        fieldWidth: 50; negativeValues: true; defaultValue: -1; max: upperTC.defaultValue - 1}
                          DoubleField { label: qsTr("to"); name: "upperTruncCauchy"; visible: cauchyInformative.checked; id: upperTC;
                                        fieldWidth: 50; negativeValues: true; defaultValue: 1; min: lowerTC.defaultValue + 1}
                        }
                      }
            					RadioButton
            					{
            						label: qsTr("Normal"); name: "normal"; childrenOnSameRow: true; id: normalInformative
            						DoubleField { label: qsTr("mean:"); name: "informativeNormalMean"; visible: normalInformative.checked; defaultValue: 0 }
            						DoubleField { label: qsTr("std:"); name: "informativeNormalStd"; visible: normalInformative.checked; defaultValue: 0.707; fieldWidth: 50 }
                        CheckBox
                        {
                          name: "truncNormal"; label: qsTr("truncation:"); visible: normalInformative.checked; childrenOnSameRow: true
                          DoubleField { label: qsTr("from"); name: "lowerTruncNormal"; visible: normalInformative.checked; id: lowerTN;
                                        fieldWidth: 50; negativeValues: true; defaultValue: -1; max: upperTN.defaultValue - 1}
                          DoubleField { label: qsTr("to"); name: "upperTruncNormal"; visible: normalInformative.checked; id: upperTN;
                                        fieldWidth: 50; negativeValues: true; defaultValue: 1; min: lowerTN.defaultValue + 1}
                        }
            					}
            					RadioButton
            					{
            						label: qsTr("t"); name: "t"; childrenOnSameRow: true; id: tInformative
            						DoubleField { label: qsTr("location:"); name: "informativeTLocation"; visible: tInformative.checked; defaultValue: 0; negativeValues: true }
            						DoubleField { label: qsTr("scale:"); name: "informativeTScale"; visible: tInformative.checked; defaultValue: 0.707; fieldWidth: 50 }
            						IntegerField { label: qsTr("df:"); name: "informativeTDf"; visible: tInformative.checked; defaultValue: 1 }
                        CheckBox
                        {
                          name: "truncT"; label: qsTr("truncation:"); visible: tInformative.checked; childrenOnSameRow: true
                          DoubleField { label: qsTr("from"); name: "lowerTruncT"; visible: tInformative.checked; id: lowerTT;
                                        fieldWidth: 50; negativeValues: true; defaultValue: -1; max: upperTT.defaultValue - 1}
                          DoubleField { label: qsTr("to"); name: "upperTruncT"; visible: tInformative.checked; id: upperTT;
                                        fieldWidth: 50; negativeValues: true; defaultValue: 1; min: lowerTT.defaultValue + 1}
                        }
            					}
            				}
              RadioButtonGroup
              {
              enabled: checkRE.checked || checkCRE.checked || checkBMA.checked
              title: qsTr("Heterogeneity (Between study SD)")
              name: "priorSE"
                      RadioButton {
                      name: "inverseGamma"; label: qsTr("Inverse gamma"); childrenOnSameRow: true; id: "igInformative"; checked: true
                      DoubleField { label: qsTr("shape:"); name: "inverseGammaShape"; visible: igInformative.checked; defaultValue: 1; fieldWidth: 50 }
                      DoubleField { label: qsTr("scale:"); name: "inverseGammaScale"; visible: igInformative.checked; defaultValue: 0.15; fieldWidth: 50 }
                      }
                      RadioButton {
                      name: "halfT"; label: qsTr("Half t"); childrenOnSameRow: true; id: halfTInformative
                      DoubleField { label: qsTr("scale:"); name: "informativehalfTScale"; visible: halfTInformative.checked; defaultValue: 0.707; fieldWidth: 50 }
                      IntegerField { label: qsTr("df:"); name: "informativehalfTDf"; visible: halfTInformative.checked; defaultValue: 1 }
                      }
                }
                Group {
                title: " "
                CheckBox {
                name: "plotPrior"
                label: qsTr("Plot prior(s)")
                }
                CheckBox {
                name: "plotRobustCheck"
                label: qsTr("Bayes factor robustness check")
                }
              }
          }


    Section {
      columns: 1
      title: qsTr("Plots")
        CheckBox {
        name: "forestPlot"
        label: qsTr("Forest plot")
        checked: true
          CheckBox {
          name: "plotForestObserved"
          label: qsTr("Observed")
          checked: true
          }
          CheckBox {
          name: "plotForestEstimated"
          label: qsTr("Estimated")
          }
          }
          CheckBox {
          name: "plotPosterior"
          label: qsTr("Prior and posterior")
          }

  }

    Section {
      columns: 1
        title: qsTr("Advanced")
        Group
        {
          title: qsTr("Estimation settings (MCMC)")
          columns: 2
          IntegerField { label: qsTr("iterations:"); name: "iterMCMC";  defaultValue: 2000; max: 1000000; fieldWidth: 50 }
          IntegerField { label: qsTr("chains:"); name: "chainsMCMC";  defaultValue: 4; max: 10; fieldWidth: 50 }
          }
          Group {
          title: qsTr("Bayes factor computation")
            RadioButtonGroup {
            name: "bfcomputation"
              RadioButton {
              name: "integration"; label: qsTr("Integration")
              checked: true
              }
              RadioButton {
              name: "bridgeSampling"; label: qsTr("Bridge sampling"); childrenOnSameRow: true; id: "bridge"; visible: bfComp.checked
              IntegerField { label: qsTr("iterations:"); name: "iterBridge"; visible: bridge.checked; defaultValue: 5000; max: 1000000; fieldWidth: 50 }
              }
            }

          }
        }
    }
}
