source("./lm/formulaWrapper.R")
source("./lm/lm.parsing.R")
source("./utils.R")

##used in batching mclapply jobs- memory blows up otherwise, since garbage collection apparently never happens untill a given call to mclapply finishes
##fits a model, walking over possible box cox constants, using the best one. Allows the user to pass in a surrogate variable matrix optionally.
fit.model <- function(js = NULL,
                      exp.mat,
                      cov.data,
                      lambdaPerGene=NULL,
                      covariateModelString,
                      checkAnova = T,
                      batchSize = ncol(exp.mat), ##TODO fix this bug right here.
                      mc.cores = 1,
                      method = "REML",
                      lambdasToTry = c(1),
                      uselme = T,
                      normalizeBeforeTransform = T,
                      normalizeAfterTransform = T,
                      gurka = T)
{
    runsinglefit <- function( jj)
    {
        if (jj %% 100 == 0) { cat(sep="","[",jj,"]") } 
        
        y <- exp.mat[ , jj]
        pheno = "y"
        cov.data[[pheno]]=y
        
        if(!is.null(lambdaPerGene))
        {
            lambdasToTry = lambdaPerGene[jj]
        }

        ##The arguments to getBestLambda are defined in the overall enclosing function.
        fit.bestLambda = try( getBestLambda(lambdasToTry,
                                            pheno = pheno,
                                            covariateModelString,
                                            cov.data,
                                            checkAnova = checkAnova,
                                            normalizeBeforeTransform = normalizeBeforeTransform,
                                            normalizeAfterTransform = normalizeAfterTransform,
                                            method = method,
                                            uselme = uselme,
                                            gurka  = gurka))

        if (class(fit.bestLambda)=="try-error")
        {
            fit.bestLambda = list(fit=NULL, selectedLambda=NULL, anovaWrapper=NULL)
        } 
        
        fit.bestLambda$Probe.Set.ID = colnames(exp.mat)[jj]
        return(fit.bestLambda)
    }

    if(is.null(js)) { js = 1:ncol(exp.mat) }
    ##    pracma::tic()

    ## if(length(js)==1 & js==T)
    ## {
    ##     browser()
    ## }
    fits = bsub$lapply.wrapper(FUN  = runsinglefit,
                               inds = js,
                               batchSize = batchSize,
                               mc.cores = mc.cores)
    ##pracma::toc()
    ##TODO check that this is ok.
    ##fitsAll = do.call(c, fits)
    return(fits)
}

getBestLambda <- function(lambdas,
                          pheno,
                          covariateModelString,
                          dataSet,
                          checkAnova=T,
                          method="REML",
                          uselme = F,
                          normalizeBeforeTransform = T,
                          normalizeAfterTransform = T,
                          extremelb = -3,
                          extremeub = 3,
                          gurka=T) 
{

    ##private function to store identical calling arguments over multiple invocations in code,
    ##prevening us from mistakenly calling with different arguments in another place in this method.
    callFit = function(lambda)
    {
        fit = fit.model.at.lambda(pheno = pheno,
                                  lambda = lambda,
                                  dataSet = dataSet,
                                  covariateModelString = covariateModelString,
                                  method = method,
                                  normalizeBeforeTransform = normalizeBeforeTransform,
                                  normalizeAfterTransform = normalizeAfterTransform,
                                  uselme = uselme,
                                  checkAnova = checkAnova)
        return(fit)
    }

    
    pvals = rep(NA, length(lambdas))
    for(l in 1:length(lambdas))
    {
        lambda = lambdas [l]
        fit = callFit(lambda)
        if(class(fit)=="try-error")
        {
            next
        }

        noiseVector = lm.parsing$getNoiseVector(fit, gurka=gurka)
        
        anovaWrapper = NULL
        if(checkAnova)
        {
            anovaWrapper = try(lm.parsing$getAnovaWrapper(fit))
            if(class(anovaWrapper)=="try-error"|
               !anovaWrapper$pvalueCol %in% colnames(anovaWrapper$an))
            {
                print("anova failure")
                next
            }
        }	

        
        pval = try(shapiro.test(noiseVector)$p.value)
        
        
        ## print(lambda)
        ## print(pval)
        ## hist(noiseVector)
        ## browser()
        
        if(class(pval)=="try-error")
        {
            browser()
        }
        pvals[l] = pval
    }
    
    selectedLambda  = lambdas[which.max(pvals)]
    if(length(selectedLambda)==0|is.na(selectedLambda))
    {
        return(list(fit=NULL, selectedLambda = NULL, anovaWrapper = NULL))
    } 
    ## In principle we could cache the model for the best lambda, and avoid some computation,
    ## but really, most of the time it turns out that we require
    ## recomputing with the inverse normal transform anyway.
    ## Also for simplicity, just recompute
    if((selectedLambda>=extremeub) || (selectedLambda<=extremelb))
    {
        selectedLambda = "inverse_normal"
    }
    ##print(paste0("selected:", selectedLambda))
    fit = callFit(lambda = selectedLambda)
    if(class(fit)=="try-error")
    {
        browser()
    }
    if(checkAnova) { anovaWrapper = try(lm.parsing$getAnovaWrapper(fit)) }

    return(list(fit=fit, selectedLambda = selectedLambda, anovaWrapper = anovaWrapper))
}


fit.model.at.lambda <- function(pheno, lambda, dataSet, covariateModelString,
                                method = "REML",
                                normalizeBeforeTransform,
                                normalizeAfterTransform,
                                uselme,
                                checkAnova) 
{
    ##print(lambda)
    if(is.na(lambda))
    {
        browser()
        stop("lambda is na??")
    }


    y = dataSet[[pheno]]
    y = transform.by.lambda(y, lambda, normalizeBeforeTransform, normalizeAfterTransform)
 
    ##TODO consider passing this as a formula, rather than modifying the dataframe itself,
    ## where the string is "transform.by.lambda(pheno, lambda, ..., etc)"
    dataSet[[pheno]] = y
   
    fit = try(fit.model.string(pheno, covariateModelString, dataSet, uselme, checkAnova, method = method))
    if(class(fit)=="try-error")
    {
        print("erroring")
    }
    return(fit)
}

##TODO move to separate file
transform.by.lambda <- function(y, lambda, normalizeBeforeTransform, normalizeAfterTransform)
{
    if(normalizeBeforeTransform)
    {
        y = scale(y,center = T)
    }

    if(lambda == "inverse_normal")
    {
        y = inverseNormal(y)
    } else {
        y = boxcox.for.lambda(y, lambda)
    }

    if(normalizeAfterTransform)
    {
        y = scale(y, center = T)
    }

    y = y[1:length(y)]
    return(y)
}

inverseNormal <- function(y)
{
    allRanks = rank(y, na.last="keep")
    maxRank = max(allRanks, na.rm=TRUE);
    return(qnorm(allRanks/(maxRank+1)));
}


boxcox.for.lambda <- function(y, lambda1)
{
    offset = .1  #we cant have lambda2==-y_i, so we will make sure lambda2>=-y_i+offset
    lambda2 = max(0, max(-y, na.rm=T)+offset) #dont bother having a lambda2 unless there is at least one zero value


    if(lambda1 == 0)
    {
        out = log(y + lambda2)
    }
    else
    {
        out = (((y + lambda2)^lambda1)-1)/lambda1
    }

    return(out)
}

fit.model.string <- function(pheno, covariateModelString, dataSet, uselme, checkAnova, method="REML")
{
    ## if(uselme & checkAnova)
    ## {
        
    ##     stop("lme doesn't provide anova functionality right now (if ever)")
    ## }
    
    covInfo = formulaWrapper$parseCovariateString(covariateModelString)

    if(length(covInfo$ranef)==0)
    {
        fit = (lm(as.formula(paste0(pheno, paste0("~ ",paste(covInfo$fixef, collapse=" + ")))), data=dataSet))
    } else if(!uselme) {
        formla = as.formula(paste0(pheno, covariateModelString))
        ##TODO: may need to use substitute and eval if environment screws up internal fields
        if(checkAnova)
        {
            fit = (lmerTest::lmer(formla, data=dataSet))
        } else {
            fit = (lme4::lmer(formla, data=dataSet))
        }
    } else {
        ## we are doing this rigamarole with substitute so that the "fixed" and "random"
        ## call fields of the resulting lme object make sense in the environment of dataSet,
        ## and not just in this functions scope.
        ## NB anova comparing 2 lme models needs the fixed and random fields.
        fit = (eval(substitute(
            lme (data   = dataSet,
                 fixed  = as.formula(paste0(pheno, paste0("~ ",paste(covInfo$fixef, collapse=" + ")))),
                 random = as.formula(paste0("~",gsub("\\(|\\)", "", covInfo$ranef))),
                 method = method))))

        
        if(checkAnova)
        {
            an = anova(fit)
        }
        ##TODO: why on earth is this necessary?? data is supposed to be stored
        ##fit$data = dataSet
    }
    return(fit)
}


compareModels <- function(modelfull, modelrestricted)
{
    if(length(modelfull)!=length(modelrestricted))
    {
        stop("mismatched lengths")
    }
    comparisons = rep(NA, length(modelfull))
    
    for(i in 1:length(modelfull))
    {
        pval = try(lm.parsing$compareModels(modelfull[[i]]$fit, modelrestricted[[i]]$fit))
        
        if(class(pval)=="try-error")
        {
            pvals[i] = NA
        } else {
            pvals[i] = pval
        }
    }
    return(pvals)
}

