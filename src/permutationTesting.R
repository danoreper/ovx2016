source("./bsub.R")
source("./utils.R")

permutationTesting = new.env(hash=T)

## datasets: a data frame of covariates with a phenotype
## phen.mat: a matrix of phenotypes with rows equal to dataset, where each column is a phenotype; for example,each phenotype could be a an entirely distinct phenotype, or  the same phenotype under a different random imputation
##shuffles: a matrix of permutations indexes, in which each column is a permutation. Should have rows equal to dataset
##computeStatFunc: a function which takes a dataset, a phenotype name, and an accumulator object as an argument. Will be
## passed through to compute a statistic for every phenotype and permutation.
##accum: an accumulator object, whose key methods are the "addCall" method and the runAll method; batches up a set of functions calls using the addCall method, and then runs them in parallel with the runAll method
permutationTesting$getPermOutcomes <- function(dataset,
                                               phen.mat,
                                               shuffles,
                                               colnameToShuffle,
                                               computeStatFunc, 
                                               accum = bsub$get.mc.accumulator(mc.cores=1))
{
    accum$init(computeStatFunc)
    cnamez = colnames(phen.mat)
    for(r in 1:ncol(phen.mat))
    {
        dataset[[cnamez[r]]] = phen.mat[,r]
        accum$addCall(
            list(trainingData=dataset,
                 phen = cnamez[r],
                 colnameToShuffle = colnameToShuffle,
                 shufflez = shuffles))
    }
    outs = accum$runAll()
    T.stats = permutationTesting$collate.permutations(outs, isFile = accum$outputs.files)
    return(T.stats)
}

permutationTesting$collate.permutations <- function(outs, isFile)
{
    upperLimit = length(outs)
    T.ident = rep(NA, upperLimit)
    
    if(isFile)
    {
        outs = bsub$getOutputs(outs)
    }
    for(r in 1:length(outs))
    {
        out = outs[[r]]
        if(r==1)
        {
            T.perm = matrix(NA, nrow = length(out$T.perm), ncol = upperLimit)
        }
        
        T.perm[,r]  = out$T.perm
        T.ident[r]  = out$T.ident
    }
    return(list(T.perm = T.perm, T.ident = T.ident))
}

    
## returns a function which, takes a data set, a phenotype column label, and a set of shuffles as input.
## The function in turn returns the likelihood ratio between the alternate and null lmm models of the phenotype, for every
## permutation of the dataset, and for the identity permutation. Each element of the output T.perm corresponds to a different
## permutation.
##
## Intent of returning a function is to satisfy the getPermOutcomes interface, and to make getPermOutcomes
## general rather than tied in to a specific statistic; the likelihood ratio is just one potential stat we may want to consider
permutationTesting$get.lr.statfunc <- function(mA.string,
                                               mN.string,
                                               lambdasToTry,
                                               normalizeBeforeTransform,
                                               normalizeAfterTransform,
                                               uselme = F,
                                               gurka)
{
    x = force(mA.string)
    x = force(mN.string)
    x = force(lambdasToTry)
    x = force(normalizeBeforeTransform)
    x = force(normalizeAfterTransform)
    x = force(uselme)
    x = force(gurka)
    callFit <- function(trainingData,phen, m.string, lambdasToTry)
    {

        checkAnova = F
        model = fit.model(exp.mat = as.matrix(trainingData[[phen]], ncol=1),
                          cov.data = trainingData,
                          covariateModelString = m.string,
                          checkAnova = checkAnova,
                          lambdasToTry = lambdasToTry,
                          normalizeBeforeTransform = normalizeBeforeTransform,
                          normalizeAfterTransform  = normalizeAfterTransform,
                          uselme = uselme,
                          gurka=gurka)[[1]]
        
        return(model)
    }

    func <- function(trainingData, phen, colnameToShuffle, shufflez)
    {
        
        identityModelAlt  = callFit(trainingData,  phen = phen, m.string = mA.string, lambdasToTry = lambdasToTry)
        lambdasToTry      = identityModelAlt$selectedLambda
        identityModelNull = callFit(trainingData,  phen = phen, m.string = mN.string, lambdasToTry = lambdasToTry)
        
        T.ident = lm.parsing$get.lik.ratio.stat(identityModelAlt$fit, identityModelNull$fit)
        
        T.perm = rep(NA, ncol(shufflez))
        for(p in 1:ncol(shufflez))
        {
            if(p%%100==0)
            {
                print(p)
            }
            trainingDataCopy = trainingData
            if(!(colnameToShuffle %in% colnames(trainingDataCopy)))
            {
                stop(paste0("shuffle col ", colnameToShuffle, "not in colnames: ", colnames(trainingDataCopy)))
            }
            trainingDataCopy[[colnameToShuffle]] = trainingDataCopy[[colnameToShuffle]][shufflez[,p]]
            permModelAlt  = callFit(trainingDataCopy, phen, mA.string, lambdasToTry = identityModelAlt$selectedLambda)
            permModelNull = callFit(trainingDataCopy, phen, mN.string, lambdasToTry = identityModelAlt$selectedLambda)
            T.perm[p] = lm.parsing$get.lik.ratio.stat(permModelAlt$fit, permModelNull$fit)
        }
        return(list(T.ident = T.ident, T.perm = T.perm))
    }
    return(func)
}

##T.perm is a matrix in which each cell is a statistic computed for a given permutation (row) and phenotype (column).
##T.ident is a vector in which every cell is a statistic computed for the identity permutation.
permutationTesting$compute.p.values <- function(T.perm, T.ident, pfunc = util$empirical.one.tail.p.value)
{
    ##TODO fix this
    direction = 1 
    if(direction == 1)
    {
        globfunc = max
    } else {
        globfunc = min
    }
    T.perm.global     = apply(MARGIN = 1, FUN = globfunc, X = T.perm, na.rm = T)

    dfs = list()
    for(j in 1:ncol(T.perm))
    {
        phen = colnames(T.perm)[j]
        p.value              = pfunc(T.ident[j], T.perm[,j],        direction)
        p.value.global       = pfunc(T.ident[j], T.perm.global, direction) 
        if(p.value>p.value.global)
        {
            browser()
        }
        df = data.frame(phen = phen,
                        p.value = p.value,
                        p.value.global = p.value.global)
        
        dfs = util$appendToList(dfs, df)
    }
    dfs        = do.call(rbind, dfs)
    ##dfs$p.fdr  = p.adjust(dfs$p.value, method = "fdr")
    
    return(dfs)
}
