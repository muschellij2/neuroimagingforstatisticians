######################################
# Example of Hippocampus Segmentation
# Using MALF 
######################################
rm(list=ls())
library(fslr)
library(extrantsr)
library(plyr)
hn = system("hostname", 
            intern = TRUE)
rootdir = NULL
if (grepl("compute", hn)){
  rootdir = file.path(
    "/dcl01/smart/data", 
    "structural", 
    "Templates")    
}
if (grepl("taki", hn) |
    grepl("scisub", hn)){
  rootdir = file.path("/project", 
                      "taki2", 
                      "Templates")
}


## Data from 
# https://masi.vuse.vanderbilt.edu/workshop2012/index.php/Workshop_Program
# SS images were skull stripped using
# fslmaths image -mas label image_SS
# need 
# https://masi.vuse.vanderbilt.edu/workshop2012/images/e/e6/MICCAI-Challenge-2012-Label-Information.xlsx
# for mapping of labels
datadir = file.path(rootdir, 
                    paste0("MICCAI-2012-Multi-", 
                           "Atlas-Challenge-Data"))
template_dir = file.path(datadir, 
                         "all-images")
## all images is all files in one folder

#######################################
# Just getting the template data together
#######################################
niis = list.files(
  path = template_dir,
  pattern = ".nii.gz", 
  full.names = TRUE)

bases = nii.stub(niis, bn = TRUE)
templates = niis[grep("_3$", bases)]

df = data.frame(
  template = templates,
  stringsAsFactors = FALSE)
df$ss_template = paste0(
  nii.stub(df$template),
  "_SS.nii.gz")
df$label = paste0(
  nii.stub(df$template),
  "_glm.nii.gz")
stopifnot(all(file.exists(unlist(df))))

#######################################
# Just keeping a few for demonstration
# More templates = more computation, but
# More templates = better segmenation
#######################################
xdf = df
df = df[1:5,]

hippo_inds = c(47:48, 170:171)

#######################################
# Making a list of labeled files 
# 47 - Right hippocampus
# 48 - left hippocampus
# 170 Right PHG   parahippocampal gyrus
# 171 Left PHG   parahippocampal gyrus
#######################################
lab_list = llply(
  df$label,
  function(x) {
    img = readnii(x)
    hippo = niftiarr(img, 
                     # img %in% c(38:41, 71:73)
                     img %in% hippo_inds
    )
    hippo
  }, .progress = "text")


temp_list = llply(
  df$ss_template, 
  readnii,
  .progress = "text")

################################
# Let's try an image to test
################################
ind = 6
infile = xdf$ss_template[ind]
label_file = xdf$label[ind]
outfile = file.path("~",
                    paste0(nii.stub(infile, bn = TRUE),
                           "_Labeled.nii.gz"))

################################
# Run multi-atlas label fusion (malf)
################################
if (file.exists(outfile)){
  malfer = readnii(outfile)
} else {
  malfer = malf(
    infile = infile,
    template.images = temp_list,
    template.structs = lab_list,
    outfile = outfile,
    typeofTransform = "SyN",
    interpolator = "NearestNeighbor",    
    keep_images = FALSE)    
}


################################
# read in the images
################################
img = readnii(infile)
label = readnii(label_file)
hippo = niftiarr(label, 
                 label %in% hippo_inds)

################################
# Compare to gold standard
################################
tab = table(c(hippo), c(malfer))
lab_tab = table(c(label), c(malfer))

################################
# get center (for plotting)
################################
xyz = xyz(hippo)

################################
# Show differences
################################
ortho_diff(img, 
           pred = malfer, 
           roi = hippo, xyz = xyz)



ortho2(img, malfer, xyz = xyz)
ortho2(img, hippo, xyz = xyz)

# 
# jointLabelFusion(targetI, 
#                  targetIMask, atlasList, 
#                  beta = 4, rad = NA,
#                  labelList = NA, doscale = TRUE, 
#                  doNormalize = TRUE,
#                  maxAtlasAtVoxel = c(1, Inf), 
#                  rho = 0.01, usecor = FALSE,
#                  boundary.condition = "image", 
#                  rSearch = 2, segvals = NA,
#                  probimgs = NA, jifImage = NA)
