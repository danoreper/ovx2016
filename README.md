# ovx2016

## Installation:
Assumes R>=3.2. Assumes Jags has been installed. As of this writing all other dependencies are R packages, installed from biomaRt and CRAN.

## Running the code:
All main scripts assume the working directory is ovxgold/src

To run permutation testing:
ovxfold/src$ R CMD BATCH --no-save --no-restore ./ovx/perm.testing.main.R

To run permutation testing with a yaml property override file:
ovxgold/src$ R CMD BATCH --no-save --no-restore '--args OVERRIDEFILENAME.yaml'./ovx/perm.testing.main.R 
e.g.,
ovxgold/src$ R CMD BATCH --no-save --no-restore '--args ../config/defaultBayes.yaml' ./ovx/perm.testing.main.R 

To run on the cluster, use bsub; e.g.,
ovxgold/src$ bsub -M 32 R CMD BATCH --no-save --no-restore ./ovx/perm.testing.main.R '--args ../config/defaultBayes.yaml' 


To run estimation, and generate bayesian p-values, caterpillar plots, etc:
ovxgold/src$ R CMD BATCH --no-save --no-restore ./ovx/estimation.main.R
Properties can be similarly overriden as for permutation testing, and script can be run on the cluster


2) YAML propery files:
The default properties have been set assuming a single small machine, and also assuming the user expects the code to run in less than an hour. As such, the default sampling iterations, number of permutations, number of random imputations, etc, are all low, and so the results may be somewhat innacurate. If in possession of a more powerful machine, or running on a compute cluster administered by LSF, property override files may be appropriate.

Existing property override files include:
ovxgold/config/defaultBayes.yaml:   Property settings assuming a single 48 processor machine.
ovxgold/config/defaultCluster.yaml: Property settings assuming the code is run on an LSF cluster (e.g. killdevil at UNC)

The default property file, which includes description of all properties, is:
ovxgold/config/defaultParams.yaml


3) Output files:
ovxgold/output/ovx/totalstrain: figures describing the overall (including intercept) strain effects on the sham population.

ovxgold/output/ovx/totalstrain/cat_*_ranef.strain.pdf: caterpillar plots of the overall strain effect  on various phenotypes. The point in each caterpillar corresponds to the mean effect, the thick line segment represents a 68% credible interval (.16, .84), and the thin line segment represents a 95% credible interval (.025, .975). Values within the caterpillar plot include the intercept; i.e.,
The confidence interval includes the average effect of a typical strain plus the particular effect of a particular strain.

ovxgold/output/ovx/strain: figures describing just the offset effect from the strain mean of each strain, on the sham population.
ovxgold/output/ovx/strain/cat_*_ranef.strain.pdf: Similar to total strain caterpillar plots, but not including intercept
ovxgold/output/ovx/strain/compare.big.pdf: A figure presenting the matrix of bayesian p-values for all strain1-strain2 effect contrasts.
Each panel shows the p-values for the contrasts for a particular phenotype. 
Within a given cell in a given panel, the color intensity represents the p-value of the difference 
between the strain2 effect and strain1 effect, where strain1  
is specified on the x-axis, and strain2 is specified on the y-axis. 
The color intensity is capped at an intensity corresponding to a p-value of .001 
(purely for vizualization purposes; otherwise, the presence of a handful of nearly 0 p-values would drown out all other significant p-values) 
The hue represents the direction of the difference between strain1 and strain2.
Contrasts whose p-value<.05 are denoted by a star in the coresponding cell. 
ovxgold/output/ovx/strain/compare.small.pdf: the same as compare.big.pdf, but zoomed in on just a handful of interesting phenotypes.
ovxgold/output/ovx/strain/heritability_herit.empirical.pdf: caterpillar plots of the heritability of the strain effect on all the phenotypes, where
heritability is computed as Var(Zu)/Var(Zu + epsilon).

ovxgold/output/ovx/strain/heritability_herit.standard.pdf: same as empirical, but now heritability is computed using the typical method:
heritability is computed as Tau^2/(Tau^2 + sigma^2). 

ovxgold/output/ovx/strain/heritability_herit.varp.pdf: same as empirical, but now heritability is computed using the varp method from BayesDiallel.
(Will has explained this in methods. 


ovxgold/output/ovx/strain.by.treatment: figures describing strain-by-treatment effects, not including the mean treatment effect,in the whole population.
ovxgold/output/ovx/strain.by.treatment/*: figures of an identical form as ./strain/*, but now describing strain-by-treatment effects per strain, in the whole population.


ovxgold/output/ovx/totaltreatment: figures describing strain-by-treatment effects, including the mean treatment effect, in the whole population.
ovxgold/output/ovx/strain.by.treatment/*: figures of an identical form as ./totalstrain/*, but now describing strain-by-treatment effects per strain, in the whole population.


ovxgold/output/ovx/compare.p.values.csv: the data used to generate the contrast plots; the contrast information in csv form. 

ovxgold/output/ovx/effects.summary.csv: file describing the mean and 95 percent confidence interval of all the effect types of interest.
effect.name is the name of a strain for which an effect is computed.
effect.type is the type of effect we are considering.
 Possibilities for effect.type include:
1) strain: the effect of strain on some phenotype less the overall effect of any strain. 
2) totalstrain: the total effect of strain, including its intercept.
3) strain.by.treatment: the effect of strain on the treatment, less the overall treatment effect.
4) totaltreatment: the effect of the treatment for a particular strain, taking into accoun the overall treatment effect.
5) treatment: the average treatment effect.  

ovxgold/output/ovx/effectsPlotByPhen.pdf: A scatterplot of scaled strain effect vs scaled strain-by-treatment effect, broken down by phenotype. Reveals that strain and strain-by-treatment effects arent particularlly correlated, suggesting that our addititve model-- rather than a multiplicative model -- is appropriate.

ovxgold/output/ovx/herit.summary.csv: a summary akin to effects.summary.csv, 
but now the statistics are per phenotype, rather than per phentype and effect.name.
Includes all the various type of heritiability.

ovxgold/output/ovx/perm.p.values: median pvalues from permutation testing describing the significance of strain and strain by diet effects.

ovxgold/output/ovx/treatment.pvalues.csv: median pvalues from multiple imputations (no need for perms) describing the significance of treatment effects 




4) Source files:

Parallelizing code, on single and multiple machines:
ovxgold/src/bsub.R
ovxgold/src/bsubScript.R

Loading yaml files, properties, command line args:
ovxgold/src/loadParams.R
ovxgold/src/loadParamsFunc.R

General permutation testing functionality:
ovxgold/src/permutationTesting.R

Miscellaneous generally useful methods:
ovxgold/utils.R

Running Jags and processing gibbs samples
ovxgold/src/bayes
ovxgold/src/bayes/fitjags.R
ovxgold/src/bayes/processSamples.R

Running linear models, building linear formulas, pulling results out of linear models
ovxgold/src/lm
ovxgold/srs/fitBoxCoxModels.R
ovxgold/src/formulaWrapper.R
ovxgold/src/lm.parsing.R


Methods for matching pairs of samples, creating types of difference phenotypes, evaluating the matched phenotypes
ovxgold/src/matching.R
ovxgold/src/generate.R
ovxgold/src/eval.R

OVX project specific functionality
ovxgold/src/ovx
ovxgold/src/ovx/estimation.main.R   -- script for estimation analysis
ovxgold/src/ovx/perm.testing.main.R -- script for permutation testing
ovxgold/src/ovx/ovxGibbs.R          -- functions used by estimation analysis
ovxgold/src/ovx/modelFunc.R         -- functions used throughout ovx; e.g. parsing inputs
ovxgold/src/ovx/model.delta.bug     -- jags model for estimation of strain-by-treatment effects
ovxgold/src/ovx/model.strain.bug    -- jags model for estimation of strain effects


5) data files
Forced swim data:
ovxgold/data/ovx/FINAL\_FST\_forR\_v2.csv

Open field data:
ovxgold/data/ovx/FINAL\_OF\_forR\_v3.csv
