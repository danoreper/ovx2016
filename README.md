# ovx2016

## Installation:
Assumes R>=3.2. Assumes Jags has been installed. As of this writing all other dependencies are R packages, installed from biomaRt and CRAN.

## Running the code:
All main scripts assume the working directory is ovx2016/src

To run permutation testing:
ovxfold/src$ R CMD BATCH --no-save --no-restore ./ovx/perm.testing.main.R

To run permutation testing with a yaml property override file:
```ovx2016/src$ R CMD BATCH --no-save --no-restore '--args OVERRIDEFILENAME.yaml'./ovx/perm.testing.main.R ```
e.g.,
ovx2016/src$ R CMD BATCH --no-save --no-restore '--args ../config/defaultBayes.yaml' ./ovx/perm.testing.main.R 

To run on the cluster, use bsub; e.g.,
ovx2016/src$ bsub -M 32 R CMD BATCH --no-save --no-restore ./ovx/perm.testing.main.R '--args ../config/defaultBayes.yaml' 


To run estimation, and generate bayesian p-values, caterpillar plots, etc:
ovx2016/src$ R CMD BATCH --no-save --no-restore ./ovx/estimation.main.R
Properties can be similarly overriden as for permutation testing, and script can be run on the cluster


## YAML propery files:
The default properties have been set assuming a single small machine, and also assuming the user expects the code to run in less than an hour. As such, the default sampling iterations, number of permutations, number of random imputations, etc, are all low, and so the results may be somewhat innacurate. If in possession of a more powerful machine, or running on a compute cluster administered by LSF, property override files may be appropriate.

Existing property override files include:
ovx2016/config/defaultBayes.yaml:   Property settings assuming a single 48 processor machine.
ovx2016/config/defaultCluster.yaml: Property settings assuming the code is run on an LSF cluster (e.g. killdevil at UNC)

The default property file, which includes description of all properties, is:
ovx2016/config/defaultParams.yaml


## Output files:
ovx2016/output/ovx/totalstrain: figures describing the overall (including intercept) strain effects on the sham population.

ovx2016/output/ovx/totalstrain/cat_*_ranef.strain.pdf: caterpillar plots of the overall strain effect  on various phenotypes. The point in each caterpillar corresponds to the mean effect, the thick line segment represents a 68% credible interval (.16, .84), and the thin line segment represents a 95% credible interval (.025, .975). Values within the caterpillar plot include the intercept; i.e.,
The confidence interval includes the average effect of a typical strain plus the particular effect of a particular strain.

ovx2016/output/ovx/strain: figures describing just the offset effect from the strain mean of each strain, on the sham population.
ovx2016/output/ovx/strain/cat_*_ranef.strain.pdf: Similar to total strain caterpillar plots, but not including intercept
ovx2016/output/ovx/strain/compare.big.pdf: A figure presenting the matrix of bayesian p-values for all strain1-strain2 effect contrasts.
Each panel shows the p-values for the contrasts for a particular phenotype. 
Within a given cell in a given panel, the color intensity represents the p-value of the difference 
between the strain2 effect and strain1 effect, where strain1  
is specified on the x-axis, and strain2 is specified on the y-axis. 
The color intensity is capped at an intensity corresponding to a p-value of .001 
(purely for vizualization purposes; otherwise, the presence of a handful of nearly 0 p-values would drown out all other significant p-values) 
The hue represents the direction of the difference between strain1 and strain2.
Contrasts whose p-value<.05 are denoted by a star in the coresponding cell. 
ovx2016/output/ovx/strain/compare.small.pdf: the same as compare.big.pdf, but zoomed in on just a handful of interesting phenotypes.
ovx2016/output/ovx/strain/heritability_herit.empirical.pdf: caterpillar plots of the heritability of the strain effect on all the phenotypes, where
heritability is computed as Var(Zu)/Var(Zu + epsilon).

ovx2016/output/ovx/strain/heritability_herit.standard.pdf: same as empirical, but now heritability is computed using the typical method:
heritability is computed as Tau^2/(Tau^2 + sigma^2). 

ovx2016/output/ovx/strain/heritability_herit.varp.pdf: same as empirical, but now heritability is computed using the varp method from BayesDiallel.
(Will has explained this in methods. 


ovx2016/output/ovx/strain.by.treatment: figures describing strain-by-treatment effects, not including the mean treatment effect,in the whole population.
ovx2016/output/ovx/strain.by.treatment/*: figures of an identical form as ./strain/*, but now describing strain-by-treatment effects per strain, in the whole population.


ovx2016/output/ovx/totaltreatment: figures describing strain-by-treatment effects, including the mean treatment effect, in the whole population.
ovx2016/output/ovx/strain.by.treatment/*: figures of an identical form as ./totalstrain/*, but now describing strain-by-treatment effects per strain, in the whole population.


ovx2016/output/ovx/compare.p.values.csv: the data used to generate the contrast plots; the contrast information in csv form. 

ovx2016/output/ovx/effects.summary.csv: file describing the mean and 95 percent confidence interval of all the effect types of interest.
effect.name is the name of a strain for which an effect is computed.
effect.type is the type of effect we are considering.
 Possibilities for effect.type include:
1) strain: the effect of strain on some phenotype less the overall effect of any strain. 
2) totalstrain: the total effect of strain, including its intercept.
3) strain.by.treatment: the effect of strain on the treatment, less the overall treatment effect.
4) totaltreatment: the effect of the treatment for a particular strain, taking into accoun the overall treatment effect.
5) treatment: the average treatment effect.  

ovx2016/output/ovx/effectsPlotByPhen.pdf: A scatterplot of scaled strain effect vs scaled strain-by-treatment effect, broken down by phenotype. Reveals that strain and strain-by-treatment effects arent particularlly correlated, suggesting that our addititve model-- rather than a multiplicative model -- is appropriate.

ovx2016/output/ovx/herit.summary.csv: a summary akin to effects.summary.csv, 
but now the statistics are per phenotype, rather than per phentype and effect.name.
Includes all the various type of heritiability.

ovx2016/output/ovx/perm.p.values: median pvalues from permutation testing describing the significance of strain and strain by diet effects.

ovx2016/output/ovx/treatment.pvalues.csv: median pvalues from multiple imputations (no need for perms) describing the significance of treatment effects 




## Source files:

Parallelizing code, on single and multiple machines:
ovx2016/src/bsub.R
ovx2016/src/bsubScript.R

Loading yaml files, properties, command line args:
ovx2016/src/loadParams.R
ovx2016/src/loadParamsFunc.R

General permutation testing functionality:
ovx2016/src/permutationTesting.R

Miscellaneous generally useful methods:
ovx2016/utils.R

Running Jags and processing gibbs samples
ovx2016/src/bayes
ovx2016/src/bayes/fitjags.R
ovx2016/src/bayes/processSamples.R

Running linear models, building linear formulas, pulling results out of linear models
ovx2016/src/lm
ovx2016/srs/fitBoxCoxModels.R
ovx2016/src/formulaWrapper.R
ovx2016/src/lm.parsing.R


Methods for matching pairs of samples, creating types of difference phenotypes, evaluating the matched phenotypes
ovx2016/src/matching.R
ovx2016/src/generate.R
ovx2016/src/eval.R

OVX project specific functionality
ovx2016/src/ovx
ovx2016/src/ovx/estimation.main.R   -- script for estimation analysis
ovx2016/src/ovx/perm.testing.main.R -- script for permutation testing
ovx2016/src/ovx/ovxGibbs.R          -- functions used by estimation analysis
ovx2016/src/ovx/modelFunc.R         -- functions used throughout ovx; e.g. parsing inputs
ovx2016/src/ovx/model.delta.bug     -- jags model for estimation of strain-by-treatment effects
ovx2016/src/ovx/model.strain.bug    -- jags model for estimation of strain effects


## data files
Forced swim data:
ovx2016/data/ovx/FINAL\_FST\_forR\_v2.csv

Open field data:
ovx2016/data/ovx/FINAL\_OF\_forR\_v3.csv
