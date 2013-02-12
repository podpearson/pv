# discordanceMatrix.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


discordanceMatrix <- function(
  GT                          = geno(loadAndGenotypePvVcf())[["GT"]],
  returnRatio                 = FALSE
) {
  GTList <- split(GT, col(GT))
  names(GTList) <- dimnames(GT)[[2]]
  if(returnRatio) {
    discordance <- function(x, y) length(which(x!=y)) / length(which(x==y))
  } else {
    discordance <- function(x, y) length(which(x!=y))
  }
  vecDiscordance <- Vectorize(discordance)
  GTDiscordanceMatrix <- outer(GTList, GTList, vecDiscordance)
#  diag(GTDiscordanceMatrix) <- NA
  return(GTDiscordanceMatrix)
}
