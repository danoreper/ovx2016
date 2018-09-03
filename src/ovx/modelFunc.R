ovx = new.env(hash=T )
fp = file.path

source("./bayes/fitjags.R")
source("./bsub.R")
source("./lm/fitBoxCoxModels.R")
source("./matching/matching.R")
source("./ovx/ovxGibbs.R")
source("./permutationTesting.R")

ovx$get.delta.inputs <- function()
{
    ##browser()
    originals = ovx$getDatas()
    phens = ovx$getPhens()
    
    matchings = list()
    shuffles  = list()
    print("getting delta inputs")
    for(ori in 1:length(originals))
    {
        if(ori%%100==0)
        {
            print(ori)
        }
        matchings[[ori]] = ovx$getDeltaMatchings(originals[[ori]], prop$ovx$numImp, prop$ovx$idcol)
        
        amatching = matchings[[ori]][[1]]
        aphen    = phens[[ori]][1]
        delta.n = nrow(amatching(aphen)$df)  
        shuffles[[ori]] = util$generateShuffles(1:delta.n, prop$ovx$numPerm) 
    }

    mA.string = "~ 1 + (1|STRAIN.1)"
##    mA.sham.string = " ~ 1 + (1_
    ##mA.string = "~ 1 + STRAIN.1"
    mN.string = "~ 1"
    shuffleCol = "STRAIN.1"
    gibbsObserveInfo = ovx$getObserveInfoDelta()
    jagsModel = "./ovx/model.delta.bug"
    effectToInspect = "ranef.strain"
    froot = fp(prop$output, "ovx", "strain.by.treatment") 
    colsToIndex = "STRAIN.1"
    dir.create(froot, showWarnings = F, recursive = T)

    
    return(list(originals  = originals,
                phens     = phens,
                matchings = matchings,
                shuffles  = shuffles,
                mA.string = mA.string,
                mN.string = mN.string,
                shuffleCol = shuffleCol,
                gibbsObserveInfo = gibbsObserveInfo,
                jagsModel = jagsModel,
                effectToInspect = effectToInspect,
                colsToIndex = colsToIndex,
                gibbsBatchSize = 10*prop$ovx$mc.cores,
                froot = froot))
}

ovx$get.strain.inputs <- function()
{
    originals = ovx$getDatas()
    print("getting strain inputs")
    for(i in 1:length(originals))
    {
        originals[[i]] = originals[[i]][originals[[i]]$TREATMENT == "SHAM"]
    }
    phens = ovx$getPhens()
    
    matchings = list()
    shuffles  = list()
    for(ori in 1:length(originals))
    {
        if(ori%%100==0)
        {
            print(ori)
        }
        matchings[[ori]] = ovx$getStrainMatchings(originals[[ori]], prop$ovx$numImp)
        
        amatching = matchings[[ori]][[1]]
        aphen    = phens[[ori]][1]
        delta.n = nrow(amatching(aphen)$df)  
        shuffles[[ori]] = util$generateShuffles(1:delta.n, prop$ovx$numPerm) 
    }

    mA.string = "~ 1 + (1|BATCH) + (1|STRAIN)"
    mN.string = "~ 1 + (1|BATCH)"
    shuffleCol = "STRAIN"    
    froot = fp(prop$output, "ovx", "strain") 
    dir.create(froot, showWarnings = F, recursive = T)
    gibbsObserveInfo = ovx$getObserveInfo()
    jagsModel = "./ovx/model.strain.bug"
    effectToInspect = "ranef.strain"
    colsToIndex = c("STRAIN", "BATCH")
    
    return(list(originals  = originals,
                phens     = phens,
                matchings = matchings,
                shuffles  = shuffles,
                mA.string = mA.string,
                mN.string = mN.string,
                shuffleCol = shuffleCol,
                gibbsObserveInfo = gibbsObserveInfo,
                jagsModel = jagsModel,
                effectToInspect = effectToInspect,
                colsToIndex = colsToIndex,
                gibbsBatchSize = 100*prop$ovx$mc.cores,
                froot = froot))
}

ovx$getDatas <- function()
{
    originals = list()
    originals[[1]] = (fread(fp(prop$data, "ovx/FINAL_OF_forR_v3.csv")))
    originals[[2]] = (fread(fp(prop$data, "ovx/FINAL_FST_forR_v2.csv"))) 
    return(originals)
}

    
ovx$getPhens <- function()
{
    phens = list()
    phens[[1]] = c("TotDist", "PctCenter", "VMovNo", "MovTime", "MrgDist", "CtrDist")
##    phens[[1]] = c("PctCenter", "CtrDist")
##    phens[[1]] = c("PctCenter", "CtrDist", "TotDist")
##    phens[[1]] = c("PctCenter")
##    phens[[2]] = c("ImbTime", "MobileTime","PctImb")
    phens[[2]] = "PctImb"
    return(phens)
}

ovx$getDeltaMatchings <- function(original, numImp, idcol)
{
    matchon  = c("STRAIN", "BATCH")
    matchoff = c("TREATMENT")
    original$TREATMENT = factor(original$TREATMENT, levels = c("SHAM", "OVX"))
    
    matched = list()
    upperlimit = ovx$getUpperLimit(numImp)
    for(r in 1:upperlimit)
    {
        print(r)

        matched[[r]]     = matching$generatePairSubsetMatching(original,
                                                               idcol = idcol,
                                                               matchon = matchon, 
                                                               matchoff = matchoff)
    }
    return(matched)
}

ovx$getStrainMatchings <- function(original, numImp)
{
    matched = list()
    ## upperlimit = ovx$getUpperLimit(numImp)
    ## for(r in 1:upperlimit)
    ## {
    ##     matched[[r]]     = matching$generateIdentityMatching(original)
    ## }
    matched[[1]] = matching$generateIdentityMatching(original)
    return(matched)
}


ovx$callFit <- function(lambdasToTry, trainingData, m.string, phen)
{
    uselme = F
    checkAnova = F
    mA = fit.model(exp.mat = as.matrix(trainingData[[phen]], ncol=1),
                   cov.data = trainingData,
                   covariateModelString = m.string,
                   checkAnova = checkAnova,
                   lambdasToTry = lambdasToTry,
                   normalizeBeforeTransform = prop$ovx$normalizeBeforeTransform,
                   normalizeAfterTransform = prop$ovx$normalizeAfterTransform,
                   uselme = uselme,
                   gurka=prop$ovx$gurka)[[1]]
 
    return(mA)
}


ovx$getUpperLimit <- function(numImp)
{
    out = (as.integer(ceiling(numImp * prop$ovx$extra.imp.multiple)))
    return(out)
}

ovx$getAccumulator <- function()
{
    if(prop$onCluster)
    {
        batchSize = 1
        bsubCommand = bsub$get.default.killdevil.bsub(numProcessPerNode = 1, memoryLimit.GB = prop$ovx$clusterMemLim,  queue=prop$ovx$clusterQueue)
        accum = bsub$get.bsub.accumulator("./ovx/modelFunc.R", bsubCommand, batchSize=batchSize)
    } else {

        accum = bsub$get.mc.accumulator(mc.cores= prop$ovx$mc.cores)
    }
    return(accum)
}

ovx$getPermOutcomes <- function(dataset.all,
                                phens.all,
                                shuffles.all,
                                matchings.all,
                                colnameToShuffle,
                                mA.string,
                                mN.string)
{
    lambdasToTry = prop$ovx$lambdasToTry
       
    getStatFunc <- function()
    {
        statFunc = permutationTesting$get.lr.statfunc(mA.string = mA.string,
                                                      mN.string = mN.string,
                                                      lambdasToTry = lambdasToTry,
                                                      normalizeBeforeTransform = prop$ovx$normalizeBeforeTransform,
                                                      normalizeAfterTransform = prop$ovx$normalizeAfterTransform,
                                                      uselme =F,
                                                      prop$ovx$gurka)
        return(statFunc)
    }
    out = matching$getPermOutcomesAll(dataset.all = dataset.all,
                                      shuffles.all = shuffles.all,
                                      matchings.all = matchings.all,
                                      phens.all = phens.all,
                                      colnameToShuffle = colnameToShuffle,
                                      computeStatFunc = getStatFunc(),
                                      matchingsAccumulator  = ovx$getAccumulator())
                                                                                            
    return(out)
}

