source("./permutationTesting.R")

##dataset.all      a list of datasets, one for each pipeline/distinct collection of phenotypes
##shuffles.all     a list of shuffle matrixes of length equal to dataset.all
##matchings.all    a list of matching functions, which each, when given a phenotype labee, return a delta Y; corresponds to a set of random imputations that are phenotype agnostic.
##phens.all:       a list of phenotype label vectors, one vector of phenotype labels per dataset
##computeStatFunc: a function which computes a vector of statistics (one per permutation), given a dataset, a phenotype label,
## and a set of permutations. Also computes the identity statistic.
## an accumulator object for running calls in parallel
matching$getPermOutcomesAll <- function(dataset.all,
                                        shuffles.all,
                                        matchings.all,
                                        phens.all,
                                        colnameToShuffle,
                                        computeStatFunc,
                                        matchingsAccumulator  = bsub$get.mc.accumulator(mc.cores=1),
                                        phenotypesAccumulator = bsub$get.mc.accumulator(mc.cores=1))
{

    T.perm.all = list()
##    T.ident.all = list()
    for(ori in 1:length(dataset.all))
    {
        dataset   = dataset.all[[ori]]
        phenz     = phens.all[[ori]]
        shufflez  = shuffles.all[[ori]]
        matchingz = matchings.all[[ori]]
        p.values  = matching$getPermOutcomes(dataset = dataset,
                                           shuffles = shufflez,
                                           matchings = matchingz,
                                           phenz = phenz,
                                           colnameToShuffle = colnameToShuffle,
                                           computeStatFunc = computeStatFunc,
                                           matchingsAccumulator = matchingsAccumulator,
                                           phenotypesAccumulator = phenotypesAccumulator)

        T.perm.all[[ori]]  = p.values
        ##T.ident.all[[ori]] = T.stats$T.ident
    }
    
    out = T.perm.all
    return(out)
}

##TODO incorporate direction, or include it in pfunc?
matching$getPermOutcomes <- function(dataset,
                                     shuffles,
                                     matchings,
                                     phenz,
                                     colnameToShuffle,
                                     computeStatFunc,
                                     imputationMergeFunc = median,
                                     pfunc = util$get.empirical.one.tail.p.value,
                                     matchingsAccumulator  = bsub$get.mc.accumulator(mc.cores=1),
                                     phenotypesAccumulator = bsub$get.mc.accumulator(mc.cores=1))
{
    print("getting perm outcomes")
    matchingsAccumulator$init(matching$getPermOutcomes.single.matching,
                              otherGlobals = list(accum=phenotypesAccumulator,
                                                  shufflez = shuffles,
                                                  computeStatFunc = computeStatFunc,
                                                  colnameToShuffle = colnameToShuffle,
                                                  pfunc = pfunc,
                                                  phenz = phenz))
    for(amatching in matchings)
    {
        matchingsAccumulator$addCall(list(amatching = amatching))
    }


    T.stats = matchingsAccumulator$runAll()
    
    isFile = matchingsAccumulator$outputs.files
    if(isFile)
    {
        T.stats = bsub$getOutputs(T.stats)
    }    

    p.values = list()
    for(imp in 1:length(T.stats))
    {
        p.value = T.stats[[imp]]$p.values
        p.values = util$appendToList(p.values, p.value)
##        p.value$imp = imp 
    }

    p.values = data.table(do.call(rbind, p.values))

    p.values = p.values[,list(p.value = imputationMergeFunc(p.value),
                              p.value.global = imputationMergeFunc(p.value.global)),
                        by = "phen"]
    
    return(p.values)
}

##    matchings$collateMatchings(T.stats, imputationMergeFunc)

## Computes identity and permutation statistics for a given matching
##
##matchingz: a matching function that takes a phenotype label
##phen:      a phenotype label
##shufflez:  a matrix of permutations, each column is a permutation.
##computeStatFunc: takes a dataset, a phenotype label, and permutations, and computes statistics.
##accum: an accumulator for running calls in parallel
##
matching$getPermOutcomes.single.matching <- function(amatching,
                                                     phenz,
                                                     shufflez,
                                                     colnameToShuffle,
                                                     computeStatFunc,
                                                     pfunc,
                                                     accum)

{
    x = matching$generate.y.mat(amatching, phenz)
    Y.mat = x$Y.mat
    df = x$df
    Ts = permutationTesting$getPermOutcomes(dataset=df,
                                            phen.mat = Y.mat,
                                            shuffles = shufflez,
                                            colnameToShuffle = colnameToShuffle,
                                            computeStatFunc = computeStatFunc,
                                            accum = accum)

    colnames(Ts$T.perm) = phenz

    p.values = permutationTesting$compute.p.values(Ts$T.perm, Ts$T.ident, pfunc)

    return(list(Ts = Ts, p.values = p.values))
}

