## matching$generate.y.mat <- function(match.gen.funcs, phen)
## {
##     for(j in 1:length( match.gen.funcs))
##     {
##         match.gen.func = match.gen.funcs[[j]]
##         out = match.gen.func(phen)$y
##         if(j==1)
##         {
##             Y.mat = matrix(nrow = length(out), ncol = length(match.gen.funcs))
##         }
##         Y.mat[,j] = out
##     }
##     return(list(Y.mat = Y.mat, df = match.gen.funcs[[1]](phen)$df))
## }


## matching$generate.y.mat <- function(match.gen.funcs, phen)
## {
##     for(j in 1:length( match.gen.funcs))
##     {
##         match.gen.func = match.gen.funcs[[j]]
##         out = match.gen.func(phen)$y
##         if(j==1)
##         {
##             Y.mat = matrix(nrow = length(out), ncol = length(match.gen.funcs))
##         }
##         Y.mat[,j] = out
##     }
##     return(list(Y.mat = Y.mat, df = match.gen.funcs[[1]](phen)$df))
## }

matching$generate.y.mat <- function(match.gen.func, phens)
{
    for(j in 1:length(phens))
    {
        phen = phens[j]
        out = match.gen.func(phen)$y
        if(j==1)
        {
            Y.mat = matrix(nrow = length(out), ncol = length(phens))
        }
        Y.mat[,j] = out
    }
    colnames(Y.mat) = phens
    return(list(Y.mat = Y.mat, df = match.gen.func(phens[1])$df))
}


matching$generateIdentityMatching <- function(df)
{
    getDeltaForPhen <- function(phen)
    {
        return(list(y = df[[phen]], df = df))
    }
    return(getDeltaForPhen)
}


## TODO: pass in ordered offlevels
## matchon: fields which have to be identical in order to match a pair
## matchoff: a field which has to be different in order to match a pair
matching$generatePairSubsetMatching <- function(df, matchon, matchoff, idcol)
{
##    browser()
    dfbad = list()
    df = data.table(df)
    ##    dfnew = df[,j=list(ID = paste(get(idcol), collapse=",")), by = c(matchon, matchoff)]
    dfnew = df[,j=list(ID = paste(get(idcol), collapse=",")), by = c(matchon)]

    l1.ID = c()
    l2.ID = c()
    offlevels = levels(df[[matchoff]])
##    offlevels = orderedOffLevels
##    offlevels = unique(df[[matchoff]])
    for(i in 1:nrow(dfnew))
    {
        row = dfnew[i,]
        IDs = unlist(strsplit(row$ID, ","))
        matchingOffValues = df[[matchoff]][match(IDs, df[[idcol]])]
        
        level1.IDs  = IDs[which(matchingOffValues == offlevels[1])]
        level2.IDs  = IDs[which(matchingOffValues == offlevels[2])]

        
        numpairs = min(length(level1.IDs), length(level2.IDs))
        if(numpairs == 0)
        {
            dfbad = util$appendToList(dfbad, dfnew[i,])
            next
        }

        
        level1.IDs = sample(x = level1.IDs, replace = F,size = numpairs)
        level2.IDs = sample(x = level2.IDs, replace = F,size = numpairs)

        l1.ID = c(l1.ID, level1.IDs)
        l2.ID = c(l2.ID, level2.IDs)
    }
    dfbad = do.call(rbind, dfbad)
    ## print(nrow(dfbad))
    ## browser()
    
    matched = data.frame(ID.1 = l1.ID, ID.2 = l2.ID)
    for(cname in matchon)
    {
        for(j in 1:2)
        {
            matched[[paste0(cname,".",j)]] = df[[cname]][match(l1.ID, df[[idcol]])]
        }
    }
    
    getDeltaForPhen <- function(phen)
    {
        ##browser()
        matchedCopy = matched
        y1 = df[[phen]][match(matchedCopy$ID.1, df[[idcol]])]
        y2 = df[[phen]][match(matchedCopy$ID.2, df[[idcol]])]
        out = y2 - y1
        matchedCopy$y = out
        matchedCopy$y1 = y1
        matchedCopy$y2 = y2
        matchedCopy[[phen]] = out
        return(list(y = out, df = matchedCopy))
    }
    
    return(getDeltaForPhen)
}

