# FUNCTIONS to dissect names
# -----------------------------------------------

f1 <- function(x) {
    
    y <- strsplit(as.character(x), "_", fixed=T)
    z <- y[[1]][1]
    return(z)

}


f1.2 <- function(x) {
    
    y <- strsplit(as.character(x), "_", fixed=T)
    z <- paste(y[[1]][1],"_",y[[1]][2],sep="")
    return(z)

}


f2 <- function(x) {
    
    y <- strsplit(as.character(x), "_", fixed=T)
    z <- y[[1]][2]
    return(z)

}


f3 <- function(x) {
    
    y <- strsplit(as.character(x), "_", fixed=T)
    z <- y[[1]][3]
    return(z)

}


f4 <- function(x) {
    
    y <- strsplit(as.character(x), "|", fixed=T)
    z <- y[[1]][2]
    return(z)
    
}

f5 <- function(x) {
    
    y <- strsplit(as.character(x), ".", fixed=T)
    z <- y[[1]][1]
    return(z)

}

# -----------------------------------------------