#NOTE! this is overwritten by some notebooks, so update everywhere, if adding lines.
sp_id = 'acytrx0'
summary_name = 'canyon'
gbif_req_id = 'GBIFr14'
gbif_filter_id = 'GBIFf4'
workDir = 'T:/Occurrence_Records/'
codeDir = 'T:/Scripts/occurrence-record-wrangler/'
inDir = workDir + 'Inputs/'
outDir = workDir + 'Outputs/'
default_coordUncertainty = 500
SRID_dict = {'WGS84': 4326, 'AlbersNAD83': 102008} # Used in file names for output.
spdb = outDir + sp_id + gbif_req_id + gbif_filter_id + '.sqlite'
