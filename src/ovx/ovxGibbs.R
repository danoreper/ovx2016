source("./bayes/processSamples.R")

ovx$runDelta <- function(originals,
                         phens,
                         lambdasPerPhen,
                         forDelta,
                         mA.string,
                         matchings,
                         shuffles,
                         observeInfo,
                         jagsModel,
                         effectToInspect,
                         colsToIndex,
                         gibbsBatchSize, 
                         froot)
{
    pracma::tic()
    effect = effectToInspect
    dfplots.compare = list()
    herits          = list()
    gibbs.samplez   = list()
    
    
    for(ori in 1:length(originals))
    {
        original = originals[[ori]]
        phenz = phens[[ori]]
        for(phen in phenz)
        {
            lambda = lambdasPerPhen[[phen]]##median(lmds)
            gibbs.samples = ovx$get.gibbs.samples.delta(original = original,
                                                        phen = phen,
                                                        lambda = lambda,
                                                        forDelta = forDelta,
                                                        matchings = matchings[[ori]],
                                                        modelname = jagsModel,
                                                        mA.string = mA.string,
                                                        colsToIndex = colsToIndex,
                                                        observeInfo = observeInfo)
            tic()
            gibbs.samples = mcmc(do.call(rbind, gibbs.samples))
            gibbs.samplez [[phen]] = gibbs.samples


            n.treevarp = 10
            herit = ovx$getHerit(gibbs.samples, n.treevarp)
            herits[[phen]] = herit
            
##            plot(herit)
            
            toc()

##            ovx$caterplot(gibbs.samples, effect, phen, froot)
            
            dfplot.compare = processSamples$getComparisonFrame( gibbs.samples, effect, shortname = T)
            dfplot.compare = ovx$getPlot.compare(dfplot.compare, phen, prop$ovx$discreteBuckets,effect)
            dfplots.compare = util$appendToList(dfplots.compare, dfplot.compare)
        }
    }
    dfplots.compare = do.call(rbind, dfplots.compare)
    pracma::toc()
##    ovx$getPlot.compare.big(dfplots.compare, effect, fname = fp(froot, "contrasts.pdf"),T )
    return(list(gibbs.samples =gibbs.samplez, dfplots.compare = dfplots.compare, herits = herits))
}

ovx$getHerit <- function(gibbs.samples, n.treevarp)
{
    herit.standard = (gibbs.samples[,"precision.strain"])^(-1)/((gibbs.samples[,"precision.strain"])^(-1)+(gibbs.samples[,"precision.epsilon"])^-1)
    herit.standard = as.vector(herit.standard)
    
    strain.cols             = processSamples$getColsForEffect(gibbs.samples, "ranef.strain")
    strain.plus.eps.cols    = processSamples$getColsForEffect(gibbs.samples, "ranef.strain.plus.eps")
    precision.eps.cols      = "precision.epsilon"
    
    num = apply(gibbs.samples[,strain.cols],1, var)  
    denom = apply(gibbs.samples[,strain.plus.eps.cols],1,var)
    herit.empirical = as.vector(num/denom)


    func = function(ro, sigma.inv.2)
    {
        streff = rep(ro, each = n.treevarp)
        eps = rnorm(length(streff), mean = 0, sd = (sigma.inv.2)^-.5)
        tot = streff + eps
        out = var(streff)/var(tot)
        return(out)
    }

    herit.varp = rep(NA, nrow(gibbs.samples))
    for(i in 1:nrow(gibbs.samples))
    {
        ro          = gibbs.samples[i, strain.cols]
        sigma.inv.2 = gibbs.samples[i,precision.eps.cols]
        herit.varp[i] = func(ro, sigma.inv.2)
    }
    
    df = data.frame(herit.standard  = herit.standard,
                    herit.empirical = herit.empirical,
                    herit.varp      = herit.varp)
    herit.mcmc = mcmc(df)
    return(herit.mcmc)
}

ovx$getObserveInfoDelta <- function() 
{

    observeInfo = data.frame(variable = c("precision.strain", "precision.epsilon",
                                          "ranef.strain", "ranef.strain.plus.eps", 
                                          "intercept"),
                             factorForIndexConversion = c(NA, NA,
                                                          "STRAIN.1", NA,
                                                          NA),
                             keep = T)
    observeInfo$variable = as.character(observeInfo$variable)
    observeInfo$factorForIndexConversion = as.character(observeInfo$factorForIndexConversion)
    return(observeInfo)
}

ovx$getObserveInfo <- function() 
{
    observeInfo = data.frame(variable = c("precision.strain", "precision.batch", "precision.epsilon",
                                          "ranef.strain", "ranef.batch", "ranef.strain.plus.eps",
                                          "intercept"),
                             factorForIndexConversion = c(NA, NA, NA,
                                                          "STRAIN", "BATCH", NA,
                                                          NA),
                             keep = T)
    observeInfo$variable = as.character(observeInfo$variable)
    observeInfo$factorForIndexConversion = as.character(observeInfo$factorForIndexConversion)
    return(observeInfo)
}

ovx$get.gibbs.samples.delta <- function(original,
                                        phen,
                                        lambda,
                                        forDelta,
                                        matchings,
                                        modelname,
                                        mA.string,
                                        colsToIndex,
                                        observeInfo)
                   
{
    gibbs.samples = list()
    y.samples     = list()
    print(paste0(phen,":", lambda)) 
    commandList = list()


    accum = ovx$getAccumulator()
    accum$init(fitjags$runJags)    

    for (r in 1:length(matchings))
    {

        matched = matchings[[r]]
        trainingData = matched(phen)$df

##        trainingData$y = trainingData[[phen]]
     #   mA = ovx$callFit(prop$ovx$lambdasToTry, trainingData, mA.string, phen)
        

        ## trainingData[[phen]] =transform.by.lambda(trainingData$y2,lambdaPerGene,F,F) -
        ##     transform.by.lambda(trainingData$y1,lambdaPerGene,F,F)
        if(forDelta)
        {
            if(!any(trainingData$y1==0))
            {
                smallest = 0
            } else {
                smallest = min(trainingData$y1[trainingData$y1!=0])
            }
            ## print(smallest)
            ## print(phen)
            ## print(smallest)

            ##            trainingData[[phen]] = trainingData$y2^lambda - trainingData$y1^lambdas
            ##trainingData[[phen]] = log(   (smallest+trainingData$y2)/(smallest + trainingData$y1))

            trainingData[[phen]] = trainingData$y2 - trainingData$y1
            

        } else{
            trainingData[[phen]] = trainingData[[phen]]^lambda
        }
        
        ## trainingData[[phen]] = transform.by.lambda(trainingData[[phen]],
        ##                                            lambdaPerGene,
        ##                                            prop$ovx$normalizeBeforeTransform,
        ##                                            prop$ovx$normalizeAfterTransform)
        
        ##mA = ovx$callFit(1, trainingData, mA.string, phen)
        
        
        iter=prop$ovx$iter
        burninFrac = prop$ovx$burninfrac
        nchains = 1
        thin = prop$ovx$thin

        
        ## mG = do.call(fitjags$runJags, list(trainingData= trainingData,
        ##                               modelname   = modelname,
        ##                               observeInfo = observeInfo,
        ##                               iter        = iter,
        ##                               burninFrac  = burninFrac,
        ##                               nchains     = nchains,
        ##                               thin        = thin,
        ##                               phen        = phen,
        ##                               colsToIndex = colsToIndex))

        
        accum$addCall(
                      list(trainingData= trainingData,
                           modelname   = modelname,
                           observeInfo = observeInfo,
                           iter        = iter,
                           burninFrac  = burninFrac,
                           nchains     = nchains,
                           thin        = thin,
                           phen        = phen,
                           colsToIndex = colsToIndex))
    }
    outfiles = accum$runAll()
    if(accum$outputs.files)
    {
        collated = ovx$collateGibbs(outfiles)
    } else {
        collated = outfiles
    }
    return(collated)
 }       

ovx$collateGibbs <- function(outfiles)
{
    upperLimit = length(outfiles)


    gibbs.samples.all = list()
    counter = 1
    for(r in 1:upperLimit)
    {
        outfile = outfiles[r]
        x = try(load(file = outfile))
        if(class(x)=="try-error")
        {
            print(r)
        } else { 
            gibbs.samples.all[[r]] = clusterOut
        }
    }
    return((gibbs.samples = gibbs.samples.all))
}

ovx$getPlot.compare.big <- function(dfplot.compare, effect, fname, textsize, themargin, reallybig)
{
    
    dfplot.compare$phen = factor(x = dfplot.compare$phen, ordered = T, levels=c("TotDist", "PctCenter", "VMovNo", "MovTime", "MrgDist", "CtrDist","PctImb"))

    aplot = ovx$getPlot.compare.helper(dfplot.compare, effect)

    if(reallybig)
    {
        aplot = aplot + facet_wrap(~phen, nrow = 2)
    } else {
        aplot = aplot + facet_wrap(~phen)
    }
    aplot = aplot + theme(aspect.ratio=1,
                          panel.margin = unit(themargin, "lines"))

    aplot = aplot + xlab("Strain 1")
    aplot = aplot + ylab("Strain 2")
    ##browser()
    if(!is.null(textsize))
    {
        aplot = aplot + theme(axis.text=element_text(size = textsize))
    }
    print(fname)
    pdf(fname, width=15, height=9)
    print(aplot)
    dev.off()
    
}

ovx$getPlot.compare <- function(dfplot.compare, phen, discreteBuckets, effect)
{
    ## dfplot.compare$phen = factor(x = dfplot.compare$phen, ordered = T, levels=c("TotDist", "PctCenter", "VMovNo", "MovTime", "MrgDist", "CtrDist","PctImb"))

    dfplot.compare$discreteValue = cut(dfplot.compare$pval, discreteBuckets, include.lowest=T)
    dfplot.compare$phen = phen
    aplot = ovx$getPlot.compare.helper(dfplot.compare, effect)
    aplot = aplot + ggtitle(label=phen)
##    print(aplot)
    return(dfplot.compare)
}

ovx$getPlot.compare.helper <- function(dfplot, effectName)
{

    limz = c(-3,3)
    dfplot$pval[-log10(dfplot$pval)>(limz[2]-.01)] = 10^(-(max(limz)-.01))##.00001#.05##min(dfplot$pval)/2
    ##()
    dfplot$dir.string = ""
    dfplot$dir.discrete = factor(paste0(dfplot$direction, dfplot$discreteValue))
##    levels(dfplot$dir.discrete) = bounds
    ## dfplot$col = NA
    ## dfplot$col = blues[as.integer(dfplot$discreteValue)]
    ## dfplot$col[dfplot$direction == -1] = reds[as.integer(dfplot$discreteValue)]
    dfplot$p.with.direction = -log10(dfplot$pval) * dfplot$direction
    dfplot$stringg = 32
    dfplot$stringg[dfplot$pval<.05] = 8

    effect1 = paste0(effectName, ".1")
    effect2 = paste0(effectName, ".2")

    aplot = ggplot(dfplot, aes(x=get(effect1), y=get(effect2), fill = p.with.direction, label=stringg, shape = stringg))##dir.discrete))

    aplot = aplot + geom_tile()#geom_tile(aes(color = stringg))
    aplot = aplot + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
    
    ##aplot = aplot + scale_fill_manual(values = allcolors, name = "P-value (with direction)", drop = FALSE)
    aplot = aplot + scale_shape_identity()
    ##aplot = aplot + scale_fill_identity()
    aplot = aplot + xlab("Strain1" ) + ylab("Strain2")
    aplot = aplot + scale_fill_gradient2(low = "red", high = "blue",
                                         name = expression('-Log'[10]*'(p)'%.%'EffectDirection'),
                                         limits =limz)
    aplot = aplot + geom_point(size = 2.5, color="white")
    return(aplot)
}

ovx$caterplot <- function(gibbs.samples, effect, exact = F, thetitle=phen, limz = NULL)
{
    if(!exact)
    {
        colz = processSamples$getColsForEffect(gibbs.samples, effect)
    } else {
        colz = effect
    }
    toPlot = gibbs.samples[,colz]
    if(!exact)
    {
        colnames(toPlot) = processSamples$getShortName(effect, colnames(toPlot))
    }
    
    if(!is.null(limz))
    {
        caterplot(toPlot, col="black", val.lim = limz)
    } else
    {
        caterplot(toPlot, col="black")
    }
##    caterplot(toPlot,style = "plain", col="black")
    title(thetitle)
}

ovx$processHerit <- function(contr.strain, effect.type)
{
    dfs = list()
    for(phen in names(contr.strain$herits))
    {
        samples.phen = contr.strain$herits[[phen]]
        int.for.phen = processSamples$computeCriticalIntervalInfo(samples.phen)
        int.for.phen$phen = phen
        int.for.phen$effect.type = effect.type
        int.for.phen$herit.type = int.for.phen$effect.type
        int.for.phen$effect.type = NULL
        dfs = util$appendToList(dfs, int.for.phen)
    }
    dfs = do.call(rbind, dfs)
    rownames(dfs) = NULL
    return(dfs)
}

ovx$adjustDeltas <- function(contr.delta)
{
    for(phen in names(contr.delta$gibbs.samples))
    {
        samplez = ovx$getSamples.for.effect.type(contr.delta, phen, "ranef.strain", exact = F)
        offset = rowSums(samplez)/ncol(samplez)
        
        samplez = samplez - offset

        contr.delta$gibbs.samples[[phen]][,processSamples$getColsForEffect(contr.delta$gibbs.samples, "ranef.strain")] = samplez
        
        samplez.treat = ovx$getSamples.for.effect.type(contr.delta, phen, "intercept", exact = T)
        samplez.treat = samplez.treat +  offset
        
        contr.delta$gibbs.samples[[phen]][,"intercept"] = samplez.treat
    }
    return(contr.delta)
}

ovx$getTreatmentPlusStrain <- function(contr.delta)
{
    treat.plus.strain = list()
    treat.plus.strain$gibbs.samples = list()
    for(phen in names(contr.delta$gibbs.samples))
    {
        
        samplez = ovx$getSamples.for.effect.type(contr.delta, phen, "ranef.strain", exact = F)
        samplez.treat = ovx$getSamples.for.effect.type(contr.delta, phen, "intercept", exact = T)
        df = samplez + samplez.treat
        colnames(df) = processSamples$getColsForEffect(contr.delta$gibbs.samples[[phen]], "ranef.strain")
        treat.plus.strain$gibbs.samples[[phen]] = df
    }
    return(treat.plus.strain)
}

ovx$plotHerits <- function(contr.strain,  effect.type)
{
    froot = fp(prop$output, "ovx", effect.type)
    dir.create(froot, recursive = T, showWarnings=F)
    herits = contr.strain$herits

    herit.types = colnames(herits[[1]])
    for(herit.type in herit.types)
    {
        colz = list()
        for(phen in names(herits))
        {
            herit.phen = herits[[phen]]
            colz = util$appendToList(colz, herit.phen[,herit.type])
        }
        anout = mcmc(do.call(cbind, colz))
        colnames(anout) = names(herits)
        afile = fp(froot, paste0("heritability_", herit.type,".pdf"))
        print(afile)
        pdf(afile)
        ovx$caterplot(anout, colnames(anout), paste0("Heritability_",effect.type), exact = T,
                      limz = c(0,1))
        dev.off()
    }
}


ovx$plotEffects <- function(contr.strain,  effect, exact, atitle, lambdasPerPhen)
{
    
    froot = fp(prop$output, "ovx", atitle)
    dir.create(froot, recursive = T, showWarnings=F)
    gibbs.samples = contr.strain$gibbs.samples
    
    for(phen in names(gibbs.samples))
    {

            
        samplez = ovx$getSamples.for.effect.type(contr.strain, phen, effect, exact = exact, lambdasPerPhen)
        afile = fp(froot, paste0("cat_", phen, "_",effect,".pdf"))
        print(afile)
        pdf(afile)
        ovx$caterplot(samplez, colnames(samplez), paste0(atitle,":", phen), exact = T)
        dev.off()
    }
}

ovx$plotEffectsSmall <- function(contr.strain,  atitle, lambdasPerPhen)
{
    
    froot = fp(prop$output, "ovx", atitle)
    dir.create(froot, recursive = T, showWarnings=F)
    gibbs.samples = contr.strain$gibbs.samples

    afile = fp(froot, paste0("cat_combined_small.pdf"))
    print(afile)
    pdf(afile)
    par(mfrow = c(1,2))
    for(phen in c("PctCenter", "CtrDist"))
    {
        samplez = ovx$getSamples.for.effect.type(contr.strain, phen, "ranef.strain", exact = F, lambdasPerPhen)
        ovx$caterplot(samplez, colnames(samplez), paste0(atitle,":", phen), exact = T)
    }
    dev.off()
    
}

##TODO move into bayes
ovx$getSamples.for.effect.type <- function(contr.strain, phen, effect.type, exact, lambdasPerPhen = NULL)
{
    contr.strain.phen = contr.strain$gibbs.samples[[phen]]
    if(!exact)
    {
        colz  = processSamples$getColsForEffect(contr.strain.phen, effect.type)
        samplez = contr.strain.phen[,colz]

        if(!is.null(lambdasPerPhen))
        {
            lambd = lambdasPerPhen[[phen]]
            correct = 1/lambd
            if(any(samplez<1))
            {
                correct = as.integer(correct)
            }
                
            samplez = samplez^(correct)
        }

        colnames(samplez) = processSamples$getShortName(effect.type, colnames(samplez))
    } else {
        samplez = contr.strain.phen[,effect.type]
    }
    return(samplez)
}

    
ovx$processEffectSizes <- function(contr.strain, effect.type, exact=F, lambdasPerPhen)
{
##    ()
    dfs = list()
    for(phen in names(contr.strain$gibbs.samples))
    {
        samplez = ovx$getSamples.for.effect.type(contr.strain, phen, effect.type, exact = exact, lambdasPerPhen)
        ## lambd = lambdasPerPhen[[phen]]
        ## samplez = samples^(1/lambd)

        if(is.null(dim(samplez)))
        {
            asum = summary(samplez)
            aqt   = (asum$quantiles)
            ast   = (asum$statistics)

            df = data.frame(effect.name = effect.type,
                            mean = ast["Mean"],
                            l2.5 = aqt["2.5%"],
                            u97.5= aqt["97.5%"])
            
        } else {
            df = processSamples$computeCriticalIntervalInfo(samplez)
        }
        df$phen = phen
        dfs = util$appendToList(dfs, df)

    }
    dfs = do.call(rbind, dfs)
    rownames(dfs) = NULL
    
    return(dfs)
}

ovx$testTreatment <- function(delta.inputs)
{
    
    empirical.p.values = list()
    for(i in 1:2) 
    {
        phenz = delta.inputs$phens[[i]]
        for(phen in phenz)
        {
            print(phen)
            for(j in 1:prop$ovx$numImp)
            {
                if(j%%100==0)
                {
                    print(j)
                }
                df = delta.inputs$matchings[[i]][[j]](phen)$df
                ##            browser()
                
                mA = fit.model(exp.mat = as.matrix(df[[phen]], ncol=1),
                               cov.data = df,
                               covariateModelString = delta.inputs$mA.string,
                               checkAnova = T,
                               lambdasToTry = 1,
                               normalizeBeforeTransform = F,
                               normalizeAfterTransform = F,
                               uselme = T,
                               gurka= F)[[1]]
                an = anova(mA[[1]])
                p.value = an["(Intercept)","p-value"]
                
                out.df = data.frame(phen = phen, p.value = p.value)
                empirical.p.values = util$appendToList(empirical.p.values, out.df)
            }
        }
    }
    
    empirical.p.values = do.call(rbind, empirical.p.values)
    empirical.p.values = data.table(empirical.p.values)
    out = empirical.p.values[, list(p.value = median(p.value)), by=phen]
    
    return(out)
}
