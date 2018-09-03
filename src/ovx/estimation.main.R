# TODO: Add comment
# 
# Author: doreper
###############################################################################

library(ggplot2)
library(data.table)

library(lme4)
library(utils)
library(RLRsim)

source("./ovx/modelFunc.R")
source("./utils.R")
source("./lm/fitBoxCoxModels.R")
source("./bayes/processSamples.R")
source("./loadParams.R")

rebuildData = prop$ovx$rebuildData
pracma::tic()
if(rebuildData)
{
    delta.inputs  = ovx$get.delta.inputs()
    strain.inputs = ovx$get.strain.inputs()


    rangePerPhen   = list()
    lambdasPerPhen = list()
    for(i in 1:length(strain.inputs$originals))
    {
        trainingData = strain.inputs$originals[[i]]
        phenz        = strain.inputs$phens[[i]]
        m.string     = strain.inputs$mA.string
        checkAnova   = F
        lambdasToTry = prop$ovx$lambdasToTry
        normalizeBeforeTransform = F
        normalizeAfterTransform = F
        uselme = F
        gurka = F
        for(phen in phenz)
        {
            model = fit.model(exp.mat = as.matrix(trainingData[[phen]], ncol=1),
                              cov.data = trainingData, 
                             covariateModelString = m.string,
                              checkAnova = checkAnova,
                              lambdasToTry = lambdasToTry,
                              normalizeBeforeTransform = normalizeBeforeTransform,
                              normalizeAfterTransform  = normalizeAfterTransform,
                              uselme = uselme,
                              gurka=gurka)[[1]]
            lambd = model$selectedLambda
            
            lambdasPerPhen[[phen]] = lambd
            if(lambd!=0)
            {
                rangePerPhen[[phen]] = range(trainingData[[phen]]^lambd)
            } else {
                rangePerPhen[[phen]] = range(log(trainingData[[phen]]))
            }
        }
    }

    print(lambdasPerPhen)
    print(rangePerPhen)
    
    print("running gibbs delta")
    contr.delta = ovx$runDelta(originals = delta.inputs$originals,
                               phens       = delta.inputs$phens,
                               lambdasPerPhen = lambdasPerPhen,
                               forDelta = T,
                               mA.string   = delta.inputs$mA.string,
                               matchings   = delta.inputs$matchings,
                               shuffles    = delta.inputs$shuffles,
                               observeInfo = delta.inputs$gibbsObserveInfo,
                               jagsModel   = delta.inputs$jagsModel,
                               effectToInspect = delta.inputs$effectToInspect,
                               colsToIndex     = delta.inputs$colsToIndex,
                               ##gibbsBatchSize  = delta.inputs$gibbsBatchSize,
                               froot           = delta.inputs$froot)

    save(file =fp(prop$output, "ovx/ovx4.RData"), list = ls())

    print("running gibbs strain")
    contr.strain = ovx$runDelta(originals   = strain.inputs$originals,
                                phens       = strain.inputs$phens,
                                lambdasPerPhen = lambdasPerPhen,
                                forDelta = F,
                                mA.string   = strain.inputs$mA.string,
                                matchings   = strain.inputs$matchings,
                                shuffles    = strain.inputs$shuffles,
                                observeInfo = strain.inputs$gibbsObserveInfo,
                                jagsModel   = strain.inputs$jagsModel,
                                effectToInspect = strain.inputs$effectToInspect,
                                colsToIndex     = strain.inputs$colsToIndex,
                                ##gibbsBatchSize  = strain.inputs$gibbsBatchSize,
                                froot           = strain.inputs$froot)
    save(file =fp(prop$output, "ovx/ovx3.RData"), list = ls())
pracma::toc()
    
} else {
    print("loading")
    load(file =fp(prop$output, "ovx/ovx3.RData"))
    print("done loading")
    source("./ovx/modelFunc.R")
    source("./utils.R")
    source("./lm/fitBoxCoxModels.R")
    source("./multipleTesting.R")
    source("./bayes/processSamples.R")
    source("./loadParams.R")
}


pracma::tic()
contr.delta = ovx$adjustDeltas(contr.delta)



#### Heritability estimates
print("processing heritability")
herit.summary = rbind(ovx$processHerit(contr.strain, "strain"),
                      ovx$processHerit(contr.delta, "strain.by.treatment"))
print(herit.summary)
write.table(file = fp(prop$output, "ovx", "herit.summary.csv"),
            herit.summary, row.names = F, sep = ",")

ovx$plotHerits(contr.strain, "strain")
ovx$plotHerits(contr.delta, "strain.by.treatment")


#### Effect size estimates

##strain
print("processing effect sizes")
effects.strain.summary = ovx$processEffectSizes(contr.strain, "ranef.strain", exact= F, lambdasPerPhen = lambdasPerPhen)
effects.strain.summary$effect.type = "strain"

contr.totalstrain           = ovx$getTreatmentPlusStrain(contr.strain)
effects.totalstrain.summary = ovx$processEffectSizes(contr.totalstrain, "ranef.strain", exact = F, lambdasPerPhen=lambdasPerPhen)
effects.totalstrain.summary$effect.type = "totalstrain"


##strain by treatment
effects.strain.by.treatment.summary = ovx$processEffectSizes(contr.delta,  "ranef.strain", exact =F, lambdasPerPhen=NULL)
effects.strain.by.treatment.summary$effect.type = "strain.by.treatment"

##treat + strain.by.treat
contr.totaltreatment   = ovx$getTreatmentPlusStrain(contr.delta)
effects.totaltreatment.summary = ovx$processEffectSizes(contr.totaltreatment, "ranef.strain", exact = F, lambdasPerPhen=NULL)
effects.totaltreatment.summary$effect.type = "totaltreatment"

##Treatment
treatment.summary = ovx$processEffectSizes(contr.delta, "intercept", exact = T, lambdasPerPhen = NULL)
treatment.summary$effect.type = "treatment"

##combined
effect.summary = rbind(effects.strain.summary,
                       effects.totalstrain.summary,
                       effects.strain.by.treatment.summary,
                       effects.totaltreatment.summary,
                       treatment.summary)


##print(effects.summary)
write.table(file = fp(prop$output, "ovx", "effect.summary.csv"), effect.summary,
            row.names = F, sep = ",")

ovx$plotEffects(contr.strain,         "ranef.strain", F, "strain", lambdasPerPhen)
ovx$plotEffects(contr.totalstrain,    "ranef.strain", F, "totalstrain", lambdasPerPhen)
ovx$plotEffects(contr.delta,          "ranef.strain", F, "strain.by.treatment", lambdasPerPhen = NULL)
ovx$plotEffects(contr.totaltreatment, "ranef.strain", F, "totaltreatment", lambdasPerPhen=NULL)


ovx$plotEffectsSmall(contr.totalstrain,     "totalstrain", lambdasPerPhen)
ovx$plotEffectsSmall(contr.totaltreatment,  "totaltreatment", lambdasPerPhen = NULL)

## Contrast p-value
dfplot.small.strain = contr.strain$dfplots.compare[contr.strain$dfplots.compare$phen %in% c("PctCenter", "CtrDist"),]


dfplot.small.delta = contr.delta$dfplots.compare[contr.delta$dfplots.compare$phen %in% c("PctCenter", "CtrDist"),]


##panel plots of various contrasts
## STRAIN-STRAIN contrasts
ovx$getPlot.compare.big(contr.strain$dfplots.compare,
                        effect = "ranef.strain",
                        fname = fp(prop$output, "ovx", "strain", "compare.big.pdf"),
                        textsize = 4,
                        themargin = .5,
                        reallybig = T)

ovx$getPlot.compare.big(dfplot.small.strain,
                        effect = "ranef.strain",
                        fname = fp(prop$output, "ovx", "strain", "compare.small.pdf"),
                        textsize = NULL,
                        themargin = 2,
                        reallybig = F)

## STRAIN.BY.TREATMENT-STRAIN.BY.TREATMENT contrasts
ovx$getPlot.compare.big(contr.delta$dfplots.compare,
                        effect = "ranef.strain",
                        fname = fp(prop$output, "ovx", "strain.by.treatment", "compare.big.pdf"),
                        textsize = 4,
                        themargin = .5,
                        reallybig = T)

ovx$getPlot.compare.big(dfplot.small.delta,
                        effect = "ranef.strain",
                        fname = fp(prop$output, "ovx", "strain.by.treatment", "compare.small.pdf"),
                        textsize = NULL,
                        themargin = 2,
                        reallybig = F)


contr.strain$dfplots.compare$effect.type = "strain"
contr.delta$dfplots.compare$effect.type = "strain.by.treatment"
write.table(file = fp(prop$output, "ovx", "compare.p.values.csv"),
            rbind(contr.strain$dfplots.compare, contr.delta$dfplots.compare),
            row.names = F,
            sep= ",")

pracma::toc()


plot.df = (dcast(effect.summary, effect.name+phen~effect.type, value.var = "mean"))
plot.df = plot.df[!is.na(plot.df$strain),]

plot.df = data.table(plot.df)
plot.df[,strain.scaled := scale(strain), "phen"]
plot.df[,strain.by.treatment.scaled := scale(strain.by.treatment), "phen"]
plot.df$phen = factor(x = plot.df$phen, ordered = T, levels=c("TotDist", "PctCenter", "VMovNo", "MovTime", "MrgDist", "CtrDist","PctImb"))

cor.p.values = plot.df[,j=list(cor.p = cor.test(strain.scaled, strain.by.treatment.scaled)$p.value), by = "phen"]


pdf(fp(prop$output, "ovx", "effectsPlotByPhen.pdf"))
aplot = ggplot(plot.df, aes(x = strain.scaled, y=strain.by.treatment.scaled))
aplot = aplot + facet_wrap(~phen)
aplot = aplot + geom_point()
print(aplot)
dev.off()


## number of discarded per imputation, for both experiments.
discarded.openfield = nrow(delta.inputs$originals[[1]]) - 2*nrow(delta.inputs$matchings[[1]][[1]]("PctCenter")$df)
print("open field")
print(paste0("discarded ", discarded.openfield, " out of ", nrow(delta.inputs$originals[[1]]))) 



discarded.swim = nrow(delta.inputs$originals[[2]]) - 2*nrow(delta.inputs$matchings[[2]][[1]]("PctImb")$df)
print("forced swim")
print(paste0("discarded ", discarded.swim, " out of ", nrow(delta.inputs$originals[[2]]))) 

