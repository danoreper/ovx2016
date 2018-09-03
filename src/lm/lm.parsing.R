library("nlme")
library("lmerTest")
source("./utils.R")


lm.parsing = new.env(hash=T)


lm.parsing$get.Z.U <- function(afit, frame = lm.parsing$getCovFrame(afit), sum=F)
{

    if(class(afit)=="lme")
    {
        randomformulas = formula(afit$modelStruct$reStruct)
        Z = list()
        U = list()
        for(i in 1:length(randomformulas))
        {
            formName = names(randomformulas)[i]
            
            randomformula = randomformulas[[i]]
            randomformula = as.character(randomformula)
            if(length(randomformula)>2)
            {
                stop("unimplemented")
            }
            randomformula = randomformula[2]
            termz = strsplit(randomformula, split = "\\+")[[1]]
            termz = str_trim(termz)
            for(j in 1:length(termz))
            {
                term = termz[j]

                if(term =="1")
                {
                    randomformula = formName
                } else {
                    randomformula = paste0(formName,":",term)
                }
                
                randomformula = paste("~", randomformula)
                contrs = list()
                contrs[[formName]] = contrasts(eval(parse(text=formName), envir = frame), contrasts=F)
                if(term!="1")
                {
                    contrs[[term]]     = contrasts(eval(parse(text=term),     envir = frame), contrasts=F)
                }
                ##modelmat = model.matrix(formula(randomformula),  frame, contrasts.arg = contrs)
                modelmat = model.matrix(formula(randomformula),  frame, contrasts.arg = contrs)
                modelmat = modelmat[,-which(colnames(modelmat) %in% "(Intercept)")]

                ##remove unnecessary levels from model matrix
                if(term!="1")
                {
                    modelmatorig = modelmat
                    modelmat = list()

                    for(k in 2:length(colnames(ranef(afit))))
                    {
                        colm = paste0(formName, rownames(ranef(afit)), ":", colnames(ranef(afit))[k]) 
                        modelmat = util$appendToList(modelmat, modelmatorig[,colm])
                    }
                    modelmat = do.call(cbind, modelmat)
                }
                
                fullname = paste0(formName, ".", names(ranef(afit))[j])
                U[[fullname]] = ranef(afit)[,j]
                Z[[fullname]] = modelmat
            }
        }

        if(length(names(randomformulas))>1)
        {
            stop("unimplemented")
        }
    }
    else {
##        browser()
        Z     = lme4::getME(afit, "Ztlist")
        Z     = lapply(FUN=t, Z)
        ## Znew  = list()
        ## for(zname in names(Z))
        ## {
        ##     nesting = strsplit(zname, "\\.")
        ##     nest.1 = nesting[[1]]
        ##     nest.2 = nesting[[2]]
        ##     if(!(nest.1 %in% names(Znew)))
        ##     {
        ##         Znew[[nest.1]] = list()
        ##     }
        ##     Znew[[nest.1]][[nest.2]] = Z[[zname]]
        ## }
        ## Z = Znew

        U     = ranef(afit)
        Unew = list()
        for(u.name in names(U))
        {
            usub = U[[u.name]]
            
            for(u.name.sub in colnames(usub))
            {
                Unew[[paste0(u.name, ".", u.name.sub)]] = usub[[u.name.sub]]
            }
        }
        U = Unew
    }


    
    out = (list(Z=Z, U=U))
    if(sum)
    {
        firstname = names(U)[1]
        total = 0*(Z[[firstname]] %*% U[[firstname]]) 
        for(u.name in names(U))
        {
            total = total + Z[[u.name]] %*% U[[u.name]]
        }
        out$total = as.vector(total)
    }


    return(out)
}

lm.parsing$getLHS <- function(afit, fixedterms = lm.parsing$getFixedTerms(afit))
{
    out = stringr::str_trim(as.character(fixedterms)[2])    
    return(out)
}

lm.parsing$getFixedTerms <- function(afit)
{
    theclass = class(afit)
    
    if(theclass=="lme")
    {
        termz = afit$terms
    } else if(theclass %in% c("lmerMod", "merModLmerTest")){
        termz = terms(afit, fixed.only=T)
    } else {
        print(theclass)
        stop("unimplemented")
    }
    return(termz)
}


lm.parsing$getCovFrame <- function(afit)
{
    theclass = class(afit)
    if(theclass=="lme")
    {
        frame = afit$data
    } else if(theclass %in% c("lmerMod", "merModLmerTest")){
        frame = model.frame(afit)
    } else {
        print(theclass)
        stop("unimplemented")
    }
    return(frame)

}

lm.parsing$getX <- function(afit, frame = lm.parsing$getCovFrame(afit))
{
    theclass = class(afit)
    
    if(theclass=="lme")
    {
        termz = afit$terms
        X = model.matrix(as.formula(termz), frame)
    } else if(theclass %in% c("lmerMod", "merModLmerTest")){
        X     = lme4::getME(afit, "X")
    } else {
        print(theclass)
        stop("unimplemented")
    }
    return(X)
}


##doreper version: works on lme or lmer object
lm.parsing$varexp <- function(afit)
{
    frame = lm.parsing$getCovFrame(afit)
    Z.U   = lm.parsing$get.Z.U(afit)
    Z = Z.U[["Z"]]
    U = Z.U[["U"]]

    termz = lm.parsing$getFixedTerms(afit)
    X     = lm.parsing$getX(afit)
    lhs   = lm.parsing$getLHS(afit, termz)
    Beta  = fixef(afit)
  
    resids = resid(afit)
    fixedTerms = unlist(attr(termz, "term.labels"))

    effectsPerSample = X*matrix(rep(as.vector(Beta),nrow(X)), byrow=T,nrow=nrow(X))

    newdf = data.frame(ID = rownames(effectsPerSample), intercept=Beta[["(Intercept)"]])
    for (fixt in fixedTerms)
    {
        newTokens = list()
        interactTokens = unlist(strsplit(fixt, ":"))
        for(interactToken in interactTokens)
        {
            acol = try(frame[[interactToken]])
            if(class(acol)=="try-error"|is.null(acol))
            {
                acol = eval(parse(text=interactToken), envir=frame)
            }
            if(is.factor(acol))
            {
                #levs = levels(frame[[interactToken]])
                levs = levels(acol)
                newTokens = util$appendToList(newTokens, paste0(interactToken, levs))
            } else
            {
                newTokens = util$appendToList(newTokens, interactToken)
            }
        
        }

        interactions = levels(do.call(interaction, c(newTokens, sep=":"))) 
        interactions = intersect(interactions, colnames(effectsPerSample))
        groupedCol = effectsPerSample[,interactions]
        if(length(interactions)>1)
        {
            groupedCol = rowSums(groupedCol)
        }
        groupedCol = matrix(groupedCol,ncol=1)
        colnames(groupedCol) = fixt

        newdf = cbind(newdf, groupedCol)
    }
    rownames(newdf) = newdf$ID
    newdf$ID = NULL

    newraneff = data.frame(ID = rownames(newdf))
    for(aranef in names(U))
    { 
        ranefmat = U[[aranef]]
        zmat     = Z[[aranef]]

        persample = try(matrix(zmat%*%ranefmat,ncol=1))
        if(class(persample)=="try-error")
        {
            browser()
        }
            
        colnames(persample) = aranef
        newraneff = cbind(newraneff, persample)
    }
    rownames(newraneff)=newraneff$ID
    newraneff$ID = NULL

    newdf = cbind(newdf, newraneff)
    newdf$resid = resids

    ##    totalVar = var(frame[[lhs]])
    response = eval(parse(text=lhs), envir=frame)

    vars  = cov(newdf)    
    vars  = vars/var(as.vector(response))
    if(any(diag(vars)<0))
    {
        browser()
    }
    return(list(response = response, components = newdf))
}


##WV version.. doesn't seem to work right
lm.parsing$lmer.vartable <- function(fit.lmer, pctvar=TRUE)
{

     y      <- fitted(fit.lmer)+resid(fit.lmer)
     n      <- length(y)
     totvar <- var(y)
     tss    <- totvar*(n-1)
                                         # variance components
     vc    <- VarCorr(fit.lmer)
     out.df <- NULL
     for (comp.name in names(vc))
     {
         comp      <- vc[[comp.name]]
         comp.sd   <- attr(comp, "stddev") #comp@factors$correlation@sd
         comp.var  <- comp.sd^2
         comp.type <- colnames(comp) #@factors$correlation)
        
         df <- data.frame(
             Name     = rep(comp.name,length(comp.var)),
             Type     = comp.type,
             Variance = comp.var)
         out.df <- rbind(out.df, df)
     }
     scale <- attr(vc, "sc")
     df <- data.frame(
         Name     = "Residual",
         Type     = "(Intercept)",
         Variance = scale^2)
     out.df <- rbind(out.df, df)
     out.df$PctVar <- 100 * out.df$Variance / totvar
    
                                         # fixed components
     an <- anova(fit.lmer, type=1)
     if (0!=nrow(an))
 	{
             fix.df <- data.frame(
                 Name     = rownames(an),
                 Type     = rep("Fixed", nrow(an)),
                 Variance = an$"Sum Sq"/n,
                 PctVar   = 100 * an$"Sum Sq"/tss)
             out.df <- rbind(out.df, fix.df)
 	}
     rownames(out.df) <- 1:nrow(out.df)
    
     return (out.df)
} 


##TODO, implment for case where var2 is also a factor, and var1 is not. Make things more general
##TODO, make this work for lme
lm.parsing$form.interaction.contrast.mat <- function( fit.with.interaction, var1, var2)
{
    var1.levelz = levels(attr(fit.with.interaction, "frame")[[var1]])
    var2.levelz = ""
    
    alleffectnames = names(fixef(fit.with.interaction))
    pairz = paste0(var1, paste(var1.levelz, var2, sep =":"))
    refpair = setdiff(pairz, alleffectnames)
    pairz = intersect(pairz, alleffectnames)
    pairz.inds = match(pairz, alleffectnames)
    ## allDeltas = t(combn(contrast.inds,2))
    allDeltaStrings = t(combn(pairz,2))
    allDeltas = cbind(match(allDeltaStrings[,1], alleffectnames),
                      match(allDeltaStrings[,2], alleffectnames))
    nvar = length(names(fixef(fit.with.interaction)))
    ncontrast = length(pairz) + nrow(allDeltas)
    contrast.mat = matrix(0, nrow = ncontrast, ncol = nvar)
    rownames(contrast.mat) = rep("", nrow(contrast.mat))
    
    counter = 1
    for(c1 in pairz.inds)
    {
        contrast.mat[counter,c1]=1
        rownames(contrast.mat)[counter] = paste(pairz[counter], "-", refpair)
        counter = counter+1
    }
    for(i in 1:nrow(allDeltas))
    {
        contrast.mat[counter, allDeltas[i,2]] = 1
        contrast.mat[counter, allDeltas[i,1]] = -1
        rownames(contrast.mat)[counter] = paste(allDeltaStrings[i,2], "-",
                                                allDeltaStrings[i,1])

        counter = counter+1
    }
    return(contrast.mat)
}


lm.parsing$getAnovaWrapper <- function(fit, type=1)
{

    fit.lmer           = (class(fit) =="merModLmerTest")
    if(!fit.lmer)
    {
        ##Default anova does NOT allow specification of type... not sure what to do here
        return(list(an = stats::anova(fit),
                    pvalueCol = "p-value"))
    }
    else
    {
        return(list(an = lmerTest::anova(fit, type=type),
                    pvalueCol = "Pr(>F)"))
    }
}


lm.parsing$getTsWrapper <- function(fit)
{
    fit.lmer = (class(fit) =="merModLmerTest")

    
    if(!fit.lmer)
    {
        return(list(
            ts        = summary(fit)$tTable,
            estCol    = "Value",
            tvalueCol = "t-value",
            pvalueCol = "p-value",
            seCol    = "Std.Error",
            dfCol    = "DF"
        ))
    }
    else
    {
        return(list(
            ts        = coefficients(summary(fit)),
            estCol    = "Estimate",
            tvalueCol = "t value",
            pvalueCol = "Pr(>|t|)",
            seCol     = "Std. Error",
            dfCol     = "df"))
    }
}


lm.parsing$compareModels <- function(modelfull, modelrestricted)
{
    an.compare = try(anova(modelfull, modelrestricted))
    out = try(list(p.value = an.compare["modelrestricted", "p-value"],
                   L.Ratio = an.compare["modelrestricted", "L.Ratio"],
                   num.df  = an.compare["modelfull",       "df"],
                   denom.df= an.compare["modelrestricted", "df"]))
                   
    return(out)
}


lm.parsing$getNoiseVector <- function(fit.with.interaction, gurka)
{
    noiseVector           = try(resid(fit.with.interaction))
    if(class(noiseVector)=="try-error")
    {
        browser()
    }
    if(gurka & class(fit.with.interaction) %in% c("lmerMod", "merModLmerTest", "lme"))
    {
        noiseVector = noiseVector + lm.parsing$get.Z.U(fit.with.interaction, sum = T)$total
    }
    return(noiseVector)
}

lm.parsing$get.lik.ratio.stat <- function(modelAlt, modelNull)
{
    out = as.numeric(2*(logLik(modelAlt)-logLik(modelNull)))
    return(out)
}


lm.parsing$getY <- function(fitWrapper)
{
    ##TODO implement
    return(0)
}
