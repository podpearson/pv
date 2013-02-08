# countryCodes.R
# 
# Package: pv
# Author: Richard D Pearson
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


countryCodes <- function(
  sampleIDs
) {
  require(RMySQL)
#  Requires the following to be run on the same machine first:
#  ssh -L 11122:malsrv1b:22 -L 33309:mcs10:3309 -N -t -x rp7@ssh.sanger.ac.uk
  mydb = dbConnect(MySQL(), user='solaris_ro', dbname='solaris', host='127.0.0.1', port=33309)
  rs <- dbSendQuery(mydb, "select temp_ox_code.sample_id, ox_code, ox_src_code, country from (select sample_id, value as ox_code from vw_sample_ids where tag=\"ox_code\") as temp_ox_code, (select sample_id, value as ox_src_code from vw_sample_ids where tag=\"ox_src_code\") as temp_ox_src_code, vw_sample_geodata where temp_ox_code.sample_id=temp_ox_src_code.sample_id and temp_ox_code.sample_id=vw_sample_geodata.sample_id")
  solarisSamples <- fetch(rs, n=-1)
  
  solarisPvSamples <- subset(solarisSamples, ox_code %in% sampleIDs)
  solarisPvDistinctSamples <- solarisPvSamples[!duplicated(solarisPvSamples[["sample_id"]]),]
  row.names(solarisPvDistinctSamples) <- solarisPvDistinctSamples[["ox_code"]]
  
  sampleCountries <- solarisPvDistinctSamples[sampleIDs, "country"]
  names(sampleCountries) <- sampleIDs
  return(sampleCountries)
}
