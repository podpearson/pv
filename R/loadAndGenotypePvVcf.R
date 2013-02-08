# loadAndGenotypePvVcf.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadAndGenotypePvVcf <- function(
  originalVcfFilename         = "data/genotypes/pv_02.vcf.gz",
  typableVcfFilename          = sub("\\.vcf\\.gz", "\\.typable\\.vcf", originalVcfFilename),
  typableRdaFilename          = paste(typableVcf, "rda", sep="."),
  reload                      = FALSE
) {
  if(!reload && file.exists(typableRdaFilename)) {
    load(typableRdaFilename)
  } else {
    if(reload || !file.exists(typableVcfFilename)) {
      cleanVcfCommand <- paste("zcat", originalVcfFilename, "| scripts/perl/cleanAndGenotypeMagnusVcf.pl >", typableVcfFilename)
      system(cleanVcfCommand)
    }
    typableVcf <- readVcf(typableVcfFilename, genome="P. vivax reference, PlasmoDB V6.0")
    save(typableVcf, file=typableRdaFilename)
  }
  return(typableVcf)
}
