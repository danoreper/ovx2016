## TODO: Add comment
## 
## Author: doreper
###############################################################################

library(ggplot2)
library(data.table)
library(lme4)
library(utils)
library(RLRsim)

##library(runjags)

source("./ovx/modelFunc.R")
source("./utils.R")
source("./lm/fitBoxCoxModels.R")
source("./bayes/processSamples.R")
source("./loadParams.R")


rebuildData = prop$ovx$rebuildData

processPermOutcomes <- function(permOutcomes.strain, effectType)
{

    permOutcomes.strain              = do.call(rbind, permOutcomes.strain)
    permOutcomes.strain$p.value.fdr  = p.adjust(permOutcomes.strain$p.value, method = "fdr")
    permOutcomes.strain$p.value.bonf = p.adjust(permOutcomes.strain$p.value, method = "bonferroni")
    permOutcomes.strain$effect.type  = effectType
    return(permOutcomes.strain)
}


if(rebuildData)
{
    delta.inputs  = ovx$get.delta.inputs()
    strain.inputs = ovx$get.strain.inputs()
    
    print("running perm delta")
    permOutcomes.delta = ovx$getPermOutcomes(delta.inputs$originals,
                                             delta.inputs$phens,
                                             delta.inputs$shuffles,
                                             delta.inputs$matchings,
                                             delta.inputs$shuffleCol,
                                             delta.inputs$mA.string,
                                             delta.inputs$mN.string)

    save(file =fp(prop$output, "ovx/ovx1.RData"), list = ls())


    print("running perm strain")
    permOutcomes.strain = ovx$getPermOutcomes(dataset.all     = strain.inputs$originals,
                                              phens.all       = strain.inputs$phens,
                                              shuffles.all    = strain.inputs$shuffles,
                                              matchings.all   = strain.inputs$matchings,
                                              colnameToShuffle= strain.inputs$shuffleCol,
                                              mA.string       = strain.inputs$mA.string,
                                              mN.string       = strain.inputs$mN.string)
    save(file =fp(prop$output, "ovx/ovx2.RData"), list = ls())
} else {

    load(file =fp(prop$output, "ovx/ovx2.RData"))
}


##### p-values
perm.p.values = (rbind(processPermOutcomes(permOutcomes.strain, "strain"),
                       processPermOutcomes(permOutcomes.delta, "strain.by.treatment")))
print(perm.p.values)
write.table(file = fp(prop$output, "ovx", "perm.p.values.csv"), perm.p.values, row.names = F, sep = ",")


##Treatment p.value evaluated by taking median over imputations
out = ovx$testTreatment(delta.inputs )
write.table(file = fp(prop$output, "ovx","treatment.p.values.csv"), row.names = F, x=out) 
print("treatment p values")
print(out)
