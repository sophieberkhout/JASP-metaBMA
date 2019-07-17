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

//// Variable inputs ////
  VariablesForm
  {
    AvailableVariablesList {name: "variablesList"}
// Effect size
    AssignedVariablesList
    {
      name: "effectSize"
      title: qsTr("Effect Size")
      singleVariable: true
      allowedColumns: ["scale"]
    }
// Standard error
    AssignedVariablesList
    {
      enabled: confidenceInterval.count < 1 // Only if no confidence interval input
      id: standardError
      name: "standardError"
      title: qsTr("Effect Size Standard Error")
      singleVariable: true
      allowedColumns: ["scale"]
    }
// Confidence interval, only if no standard error input (only one of the two is necessary)
    AssignedVariablesList
    {
      enabled: standardError.count == 0
      id: confidenceInterval
      name: "confidenceInterval"
      title: qsTr("95% CI Lower and Upper Bound")
      singleVariable: true
      listViewType: "Pairs"
      allowedColumns: ["scale"]
    }
// Study labels are optional. Nominal is not working somehow.
    AssignedVariablesList
    {
      name: "studyLabels"
      title: qsTr("Study Labels")
      singleVariable: true
      // allowedColumns: ["nominal"]
    }
  }
//// End variable inputs ////

  GridLayout
  {
    columns: 2
//// Analysis choices ////
    RadioButtonGroup
    {
      name: "modelSpecification"
      title: qsTr("Model")
// Fixed effect
      RadioButton
      {
        value: "FE"
        label: qsTr("Fixed effects")
        id: checkFE
      }
// Random effects
      RadioButton
      {
        value: "RE"
        label: qsTr("Random effects")
        id: checkRE
      }
// Model averaging
      RadioButton
      {
        value: "BMA"
        label: qsTr("Model averaging")
        id: checkBMA
        checked: true
      }
// Ordered effects
      RadioButton
      {
        value: "CRE"
        label: qsTr("Constrained random effects")
        id: checkCRE
    // Constrain effect sizes to be all positive or all negative
        RadioButtonGroup
        {
          name: "direction"
          columns: 2
          RadioButton
          {
            value: "allPos"
            label: qsTr("All positive")
            checked: true
            id: checkPos
          }
          RadioButton
          {
            value: "allNeg"
            label: qsTr("All negative")
            id: checkNeg
          }
        }
      }
    }
//// End analysis choices ////

//// Prior model probabilities (defaults differ per model) ////
    Group
    {
      enabled: checkFE.checked || checkRE.checked || checkBMA.checked
      title: qsTr("Prior model probability")
      Group
      {
        enabled: checkFE.checked || checkBMA.checked
        title: qsTr("Fixed effects")
    // Fixed effect null hypothesis
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
    // Fixed effects alternative
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
    // Random effects null
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
    // Random effects alternative
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
//// End prior model probabilities ////

//// Priors ////
    Section
    {
      title: qsTr("Prior")
      columns: 1

// Prior for effect size //
      RadioButtonGroup
			{
        title: qsTr("Effect size")
				name: "priorES"

// Cauchy prior
				RadioButton
				{
					label: qsTr("Cauchy");
          name: "cauchy";
          checked: true;
          childrenOnSameRow: true;
          id: cauchyInformative
    // Cauchy location
					DoubleField
          {
            label: qsTr("location:");
            name: "informativeCauchyLocation";
            visible: cauchyInformative.checked;
            defaultValue: 0;
            negativeValues: true
          }
    // Cauchy scale
					DoubleField
          {
            label: qsTr("scale:");
            name: "informativeCauchyScale";
            visible: cauchyInformative.checked;
            defaultValue: 0.707;
            fieldWidth: 50
          }
    // Cauchy truncation NOT WORKING PROPERLY YET
          CheckBox
          {
            name: "truncCauchy";
            label: qsTr("truncation:");
            visible: cauchyInformative.checked;
            childrenOnSameRow: false;
            checked: if(checkCRE.checked){true} else {false}
            DoubleField
            {
              label: qsTr("from");
              name: "lowerTruncCauchy";
              visible: cauchyInformative.checked;
              id: lowerTC;
              fieldWidth: 50;
              negativeValues:
              {
                if(checkCRE.checked && checkPos.checked){false}
                else {true}
              }
              defaultValue:
              {
                if(checkCRE.checked && checkPos.checked){0}
                else {-1}
              }
              max:
              {
                if(checkCRE.checked && checkNeg.checked){-0.001}
              //  else {upperTC.defaultValue - 1}
              }
            // This is the lower limit of the prior.
            // It must be smaller than the upper limit.
            // Also, for the ordered effects analysis,
            // the truncation needs to be either all positive or all negative.
            }
            DoubleField
            {
              label: qsTr("to");
              name: "upperTruncCauchy";
              visible: cauchyInformative.checked;
              id: upperTC;
              fieldWidth: 50;
              negativeValues:
              {
                if(checkCRE.checked && checkPos.checked){false}
                else {true}
              }
              defaultValue:
              {
                if(checkCRE.checked && checkNeg.checked){0}
                else {1}
              }
              // min: lowerTC.defaultValue + 1;
              max:
              {
              //  if(checkCRE.checked && checkNeg.checked){0}
              }
            }
          }
        }

// Normal prior
				RadioButton
				{
					label: qsTr("Normal");
          name: "normal";
          childrenOnSameRow: true;
          id: normalInformative
    // Normal mean
					DoubleField
          {
            label: qsTr("mean:");
            name: "informativeNormalMean";
            visible: normalInformative.checked;
            defaultValue: 0
          }
    // Normal SD
					DoubleField
          {
            label: qsTr("std:");
            name: "informativeNormalStd";
            visible: normalInformative.checked;
            defaultValue: 0.707;
            fieldWidth: 50
          }
    // Normal truncation NOT WORKING PROPERLY YET
          CheckBox
          {
            name: "truncNormal";
            label: qsTr("truncation:");
            visible: normalInformative.checked;
            childrenOnSameRow: false
            DoubleField
            {
              label: qsTr("from");
              name: "lowerTruncNormal";
              visible: normalInformative.checked;
              id: lowerTN;
              fieldWidth: 50;
              negativeValues: true;
              defaultValue: -1;
              max: upperTN.defaultValue - 1
            }
            DoubleField
            {
              label: qsTr("to");
              name: "upperTruncNormal";
              visible: normalInformative.checked;
              id: upperTN;
              fieldWidth: 50;
              negativeValues: true;
              defaultValue: 1;
              min: lowerTN.defaultValue + 1
            }
          }
			  }

// T prior
				RadioButton
				{
					label: qsTr("t");
          name: "t";
          childrenOnSameRow: true;
          id: tInformative;
    // T location
					DoubleField
          {
            label: qsTr("location:");
            name: "informativeTLocation";
            visible: tInformative.checked;
            defaultValue: 0;
            negativeValues: true
          }
    // T scale
					DoubleField
          {
            label:
            qsTr("scale:");
            name: "informativeTScale";
            visible: tInformative.checked;
            defaultValue: 0.707;
            fieldWidth: 50
          }
    // T DF
					IntegerField
          {
            label: qsTr("df:");
            name: "informativeTDf";
            visible: tInformative.checked;
            defaultValue: 1
          }
    // T truncation NOT WORKING PROPERLY YET
          CheckBox
          {
            name: "truncT";
            label: qsTr("truncation:");
            visible: tInformative.checked;
            childrenOnSameRow: false
            DoubleField
            {
              label: qsTr("from");
              name: "lowerTruncT";
              visible: tInformative.checked;
              id: lowerTT;
              fieldWidth: 50;
              negativeValues: true;
              defaultValue: -1;
              max: upperTT.defaultValue - 1
            }
            DoubleField
            {
              label: qsTr("to");
              name: "upperTruncT";
              visible: tInformative.checked;
              id: upperTT;
              fieldWidth: 50;
              negativeValues: true;
              defaultValue: 1;
              min: lowerTT.defaultValue + 1
            }
          }
				}
			}
// End prior effect size //

// Prior heterogeneity //
      RadioButtonGroup
      {
        enabled: checkRE.checked || checkCRE.checked || checkBMA.checked
        title: qsTr("Heterogeneity (Between study SD)")
        name: "priorSE"
// Inverse gamma prior
        RadioButton
        {
          name: "inverseGamma";
          label: qsTr("Inverse gamma");
          childrenOnSameRow: true;
          id: "igInformative";
          checked: true
    // Inverse gamma shape
          DoubleField
          {
            label: qsTr("shape:");
            name: "inverseGammaShape";
            visible: igInformative.checked;
            defaultValue: 1;
            fieldWidth: 50
          }
    // Inverse gamma scale
          DoubleField
          {
            label: qsTr("scale:");
            name: "inverseGammaScale";
            visible: igInformative.checked;
            defaultValue: 0.15;
            fieldWidth: 50
          }
        }
// Half t prior
        RadioButton
        {
          name: "halfT";
          label: qsTr("Half t");
          childrenOnSameRow: true;
          id: halfTInformative
    // Half t scale
          DoubleField
          {
            label: qsTr("scale:");
            name: "informativehalfTScale";
            visible: halfTInformative.checked;
            defaultValue: 0.707;
            fieldWidth: 50
          }
    // Half t DF
          IntegerField
          {
            label: qsTr("df:");
            name: "informativehalfTDf";
            visible: halfTInformative.checked;
            defaultValue: 1
          }
        }
      }
// Option to plot the priors. An extra empty line would be nice.
      CheckBox
      {
        name: "plotPrior"
        label: qsTr("Plot prior(s)")
      }
    }
//// End priors ////

//// Plots section ////
    Section
    {
      columns: 1
      title: qsTr("Plots")
// Forest plot with obtion to show observed, estimated or both effect sizes.
      CheckBox
      {
        name: "forestPlot"
        label: qsTr("Forest plot")
        checked: true
        CheckBox
        {
        name: "plotForestObserved"
        label: qsTr("Observed")
        checked: true
        }
        CheckBox
        {
        enabled: !checkFE.checked
        name: "plotForestEstimated"
        label: qsTr("Estimated")
        }
      }
// Prior and posterior plot
      CheckBox
      {
      name: "plotPosterior"
      label: qsTr("Prior and posterior")
      }
    }
//// End plots section ////

//// Advanced section for sampling settings ////
    Section
    {
      columns: 1
      title: qsTr("Advanced")
// MCMC estimation settings
      Group
      {
        title: qsTr("Estimation settings (MCMC)")
        columns: 2
    // MCMC number of iterations
        IntegerField
        {
          label: qsTr("iterations:");
          name: "iterMCMC";
          defaultValue: 2000;
          max: 1000000;
          fieldWidth: 50
        }
    // MCMC number of chains
        IntegerField
        {
          label: qsTr("chains:");
          name: "chainsMCMC";
          defaultValue: 4;
          max: 10;
          fieldWidth: 50
        }
      }
// Choice between integation or bridge sampling BF computation
      Group
      {
        title: qsTr("Bayes factor computation")
        RadioButtonGroup
        {
          name: "BFComputation"
    // Integation
          RadioButton
          {
            name: "integration";
            label: qsTr("Integration")
            checked: true
          }
    // Bridge sampling
          RadioButton
          {
            name: "bridgeSampling";
            label: qsTr("Bridge sampling");
            childrenOnSameRow: true;
            id: "bridge";
            visible: bfComp.checked
        // Bridge sampling iterations
            IntegerField
            {
              label: qsTr("iterations:");
              name: "iterBridge";
              visible: bridge.checked;
              defaultValue: 5000;
              max: 1000000;
              fieldWidth: 50
            }
          }
        }
      }
    }
//// End advanced section ////
  }
}
