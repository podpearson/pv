# loadAndGenotypePvVcf.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


loadAndGenotypePvVcf <- function(
  vcfFilename                 = "data/genotypes/pv_02.vcf.gz",
  rdaFilename                 = sub("\\.gz", "\\.rda", vcfFilename),
  reload                      = TRUE
) {
  if(!reload && file.exists(rdaFilename)) {
    load(rdaFilename)
  } else {
    vcf <- readVcf(vcfFilename, genome="P. vivax reference, PlasmoDB V6.0")
    save(vcf, file=rdaFilename)
  }
  return(vcf)
}

