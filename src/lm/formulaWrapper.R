library(stringr)

formulaWrapper = new.env(hash=T)
##expects a covariate string (without the outcome part) of lmer form, e.g., "1 + X1 +X2 + X1:X2 + (Z|X1)"
formulaWrapper$parseCovariateString <- function(covariateString)
{
    covariateInfo = list()
    covariateInfo$original.string = covariateString

    tokens = gsub(pattern = "~", replacement = "", covariateInfo$original.string)
    ##pull out the random effects terms
    ranefs = gregexpr(tokens, pattern="\\(.*?\\|.*?\\)", perl=T)
    starts = as.vector(ranefs[[1]])
    ends = starts + attr(ranefs[[1]],"match.length")-1

    
    if(any(starts>0))
    {
        covariateInfo$ranef = rep("", length(starts))
        
        for(i in seq(1, length.out = length(starts)))
        {
            ## x = try(
            ## {
            ##     print("tokens:")
            ##     print(tokens)
            ##     print("subset:")
                
            ##     print(substr(tokens, starts[i], ends[i]))
            ##     print(traceback())
            ## })
            ## if(class(x)=="try-error")
            ## {
            ##     x = 5
            ##     browser()
            ## }
            covariateInfo$ranef[[i]] = substr(tokens, starts[i], ends[i])
        }
    } else {
        covariateInfo$ranef = c()
    }
        
    ##remove the random effects before getting the fixed effects
    tokens = gsub(tokens, pattern="\\(.*?\\|.*?\\)", replacement = "", perl=T)
    tokens = stringr::str_trim(unlist(strsplit(tokens, "\\+")))
    covariateInfo$fixef = tokens[tokens!=""]

    covariateInfo = formulaWrapper$.parseCovariateStringHelper(covariateInfo)

    return(covariateInfo)
}

formulaWrapper$.parseCovariateStringHelper <- function(covariateInfo)
{
    covariateInfo$modified.string = paste0("~ ",paste(covariateInfo$fixef, collapse = " + "))
    if(length(covariateInfo$ranef)>0)
    {
        if(length(covariateInfo$fixef>0))
        {
            covariateInfo$modified.string = paste0(covariateInfo$modified.string,
                                                   " + ",
                                                   paste(covariateInfo$ranef, collapse = " + "))
        }else
        {
            covariateInfo$modified.string = paste0(covariateInfo$modified.string,
                                                   paste(covariateInfo$ranef, collapse = " + "))
        }
    }
    return(covariateInfo)
}

formulaWrapper$onlyKeepEffectAndInteractions <- function(effect.string, covariateString)
{
    effect.string = c("1", effect.string)

    return(formulaWrapper$.setOperationHelper(effect.string, covariateString, keepMatchingTerms = T))
}

##TODO rewrite to allow removal of just interactions?
#effect.string is one (or several) of the main effect terms (fixed or random) in covariateInfo
formulaWrapper$removeEffectAndInteractions <- function(effect.string, covariateString)
{
    return(formulaWrapper$.setOperationHelper(effect.string, covariateString, keepMatchingTerms = F))
}

formulaWrapper$.setOperationHelper <- function(effect.string, covariateString, keepMatchingTerms)
{
    checkOutcome <- function(tokens, effect.string)
    {
        return(any(tokens %in% effect.string))
    }

    tinytokensplit = "\\*|\\:|\\+"

    covariateInfo = formulaWrapper$parseCovariateString(covariateString)
    ranef = c()
    for(i in seq(1, length.out = length(covariateInfo$ranef)))
    {
        arand = covariateInfo$ranef
        arand = gsub(arand, pattern = "\\(|\\)", replacement = "")
        left.right = strsplit(arand, "\\|")[[1]]
        left  = left.right[1]
        right = left.right[2]


        ##Remove any terms with even a partially irrelevant right
        right = strsplit(right, "\\+")[[1]]
        tinytokens = strsplit(right, split  = tinytokensplit)
        
        if(any(unlist(lapply(tinytokens, checkOutcome, effect.string))))
        {
            if(!keepMatchingTerms)
            {
                next
            } 
        }

        ##Remove irrelevant terms from left
        left = stringr::str_trim(strsplit(left,"\\+")[[1]])
        tinytokens = strsplit(left, split  = tinytokensplit)
        tokeep = (lapply(tinytokens, checkOutcome, effect.string) == keepMatchingTerms)
        left = left[tokeep]
        if(length(left)==0)
        {
            next
        }

        arand = paste("(", paste(left, collapse = "+"), "|", paste(right, collapse="+", ")", sep=""))
        ranef = c(ranef, arand)
    }

    fixef = c()
    for(i in seq(1,length.out =length(covariateInfo$fixef)))
    {
        afix = covariateInfo$fixef[i]
        left = stringr::str_trim(strsplit(afix,"\\+")[[1]])
    
        tinytokens = strsplit(left, split  = "\\*|\\:|\\+")
        tokeep = (lapply(tinytokens, checkOutcome, effect.string) == keepMatchingTerms)
        left = left[tokeep]
        if(length(left)==0)
        {
            next
        }
        fixef = c(fixef, left)
    }


    covariateInfo$fixef = fixef
    covariateInfo$ranef = ranef
    covariateInfo = formulaWrapper$.parseCovariateStringHelper(covariateInfo)
    return(covariateInfo)
}

