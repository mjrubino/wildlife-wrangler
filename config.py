#NOTE! this is overwritten by some notebooks, so update everywhere, if adding lines.
sp_id = 'bbtnwx1'
summary_name = 'green_warbler1'
gbif_req_id = 'GBIFr12'
gbif_filter_id = 'GBIFf4'
ebird_req_id = None
ebird_filter_id = None
evaluation = 'eval_gbif1'
workDir = '/Users/nmtarr/Documents/RANGES/'
codeDir = '/Users/nmtarr/Code/occurrence-record-wrangler/'
inDir = workDir + 'Inputs/'
outDir = workDir + 'Outputs/'
default_coordUncertainty = 1000
SRID_dict = {'WGS84': 4326, 'AlbersNAD83': 102008} # Used in file names for output.
spdb = outDir + sp_id + gbif_req_id + gbif_filter_id + '.sqlite'
