# discordanceMatrix.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


discordanceMatrix <- function(
  GT                          = geno(loadAndGenotypePvVcf())[["GT"]],
  returnValue                 = c("count", "ratio", "proportion")
) {
#  GTList <- split(GT, col(GT))
#  names(GTList) <- dimnames(GT)[[2]]
  if(returnValue[1] == "ratio") {
    discordance <- function(x, y) length(which(x!=y)) / length(which(x==y))
    GTDiscordanceMatrix <- matrix(numeric(), dim(GT)[2], dim(GT)[2], dimnames=list(dimnames(GT)[[2]], dimnames(GT)[[2]]))
  } else if(returnValue[1] == "count"){
    discordance <- function(x, y) length(which(x!=y))
    GTDiscordanceMatrix <- matrix(integer(), dim(GT)[2], dim(GT)[2], dimnames=list(dimnames(GT)[[2]], dimnames(GT)[[2]]))
  } else if(returnValue[1] == "proportion"){
    discordance <- function(x, y) length(which(x!=y)) / length(which(x==y | x!=y))
    GTDiscordanceMatrix <- matrix(numeric(), dim(GT)[2], dim(GT)[2], dimnames=list(dimnames(GT)[[2]], dimnames(GT)[[2]]))
  } else {
    stop("invalid returnValue:", returnValue)
  }
#  GTDiscordanceMatrix <- mapply(discordance, split(GT, col(GT)), split(GT, col(GT)))
  for(i in seq(dim(GT)[2])) {
    cat(sprintf("\n%3.0f", i))
    for(j in seq(dim(GT)[2])) {
      cat(".")
      GTDiscordanceMatrix[i, j] <- discordance(GT[, i], GT[, j])
    }
    gc()
  }
#  browser()
#  vecDiscordance <- Vectorize(discordance)
#  GTDiscordanceMatrix <- outer(GTList, GTList, vecDiscordance)
#  diag(GTDiscordanceMatrix) <- NA
  return(GTDiscordanceMatrix)
}
