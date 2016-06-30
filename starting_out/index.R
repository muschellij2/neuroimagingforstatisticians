## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE)

## ---- cache = FALSE------------------------------------------------------
rm(list = ls())
library(dcm2niir)
library(oro.dicom)
library(fslr)
library(plyr)
library(dplyr)
library(data.table) # just for rbindlist
library(ggplot2)
library(extrantsr)
library(scales)

## ----unzip---------------------------------------------------------------
# Unzipping and creating a df
file = "BRAINIX.zip"
out = unzip(file, exdir = tempdir())
df = data.frame(file = out,
                stringsAsFactors = FALSE)
df$fol = dirname(df$file)

## ----dcm2nii, dependson="unzip"------------------------------------------
#################################
# Convert to NIfTI for each unique path
#################################
paths = unique(df$fol)
p = lapply(paths, function(x) {
  dcm2nii(x)
})
res = sapply(p, `[[`, "result")
stopifnot(all(res == 0))

## ----check_dcm2nii, dependson="dcm2nii"----------------------------------
##########################
# Make sure only one file was produced and return it
##########################
files = sapply(p, check_dcm2nii)
stopifnot(length(files) == length(p))

## ----readnii, results="hide"---------------------------------------------
##########################
# read files as nifti
##########################
niis = llply(files, readnii, .progress = "text")
n = nii.stub(files, bn = TRUE)
n = gsub("-", "_", n) # don't want "-" in names
names(niis) = n

## ---- results="hide", warning=FALSE--------------------------------------
################################
# Snapshots of each DICOM
################################
first_image = ddply(df, .(fol),
                    function(x) x[1,])
dicom_hdrs = llply(first_image$file,
                    function(x) {
                      readDICOMFile(x)$hdr
                      }, .progress = "text")
dcmtab = dicomTable(dicom_hdrs)
rownames(dcmtab) = n

## ----sub_imgs------------------------------------------------------------
inames = c("T1_3D_FFE_C_CLEAR_20061201141645_801" = "T1_withneck",
           "T1_SE_extrp_CLEAR_20061201141645_701" = "T1",
"sT2W_FLAIR_SENSE_20061201141645_401" =   "FLAIR" ,
  "sT2_TSE_T_SENSE_20061201141645_301"= "T2" )
imgs = niis[names(inames)]
names(imgs) = inames
imgs

## ------------------------------------------------------------------------
mat = lapply(inames,
             function(x){
               vals = c(imgs[[x]]) # vector
               vals = vals[ vals > median(vals) ] # rm "background"
               data.frame(value = vals,
                          type = x,
                          stringsAsFactors = FALSE)
             })
mat = rbindlist(mat) # data.table
mat = as.data.frame(mat) # don't want a data.table
# make plot
g = ggplot(aes(x = value), 
           data = mat) + 
  geom_line(stat = "density") +
  facet_wrap(~type, scales = "free_x")

## ----hist, echo = FALSE--------------------------------------------------
print(g)

## ----plot_t1, cache=TRUE-------------------------------------------------
ortho2(imgs$T1)

## ----wplot_t1, cache=TRUE------------------------------------------------
ortho2(robust_window(imgs$T1))

## ----windowing-----------------------------------------------------------
wimgs = llply(imgs, robust_window, .progress = "text")

## ----plot_t2, cache=TRUE, dependson="windowing"--------------------------
ortho2(wimgs$T2)

## ----plot_flair, cache=TRUE, dependson="windowing"-----------------------
ortho2(wimgs$FLAIR)

## ----n4, dependson="sub_imgs", cache = TRUE------------------------------
n4 = llply(imgs, bias_correct, correction = "N4", .progress = "text")

## ----make_ratio----------------------------------------------------------
rat = finite_img(imgs$T1 / n4$T1)
img_cut = function(img, breaks, ...){
  cuts = cut(img, breaks = breaks, ...)
  levs = levels(cuts)
  cuts = as.numeric(cuts)
  # res.p[ rs > ncut ] = cuts
  img = niftiarr(img, array(cuts, dim = dim(img)))
  return(list(img=img, levs = levs))
}

breaks = unique(
  quantile(rat[rat > 0], probs = seq(0, 1, by = 0.2))
)

col.cut = div_gradient_pal(low="blue", 
                                 mid="red", 
                                 high="yellow")
col.cut = col.cut(seq(0, 1, length = length(breaks)-1) )
col.cut = alpha(col.cut, 0.5)
cut_img = img_cut(rat, 
                  breaks = breaks, include.lowest=FALSE)

## ----plot_ratio, cache=TRUE, dependson="make_ratio"----------------------
ortho2(imgs$T1, cut_img$img, col.y = col.cut, 
       ybreaks = 0:length(cut_img$levs), 
       ycolorbar = TRUE,
       clabels = cut_img$levs)

## ----slice_hist----------------------------------------------------------
slices = c(4, 8, 12, 16, 20)
df = lapply(slices, 
            function(slice) {
              vals = c(imgs$T1[,,slice])
              data.frame(
                value = c(vals, n4$T1[,,slice]),
                type = rep(c("T1", "N4"), each = length(vals)),
                slice = slice,
                stringsAsFactors = FALSE)
            })
df = do.call("rbind", df)
df = df[ df$value > 0, ]
g = ggplot(aes(x = value, colour = factor(slice)), 
           data = df) +
  geom_line(stat = "density") + facet_wrap(~type)

## ----plot_slice_hist, cache=TRUE, dependson = "slice_hist"---------------
print(g)

## ----bet,dependson="n4", cache = TRUE, results= 'hide', warning= FALSE----
bet = fslbet_robust(n4$T1, correct = FALSE)

## ----plot_bet------------------------------------------------------------
ortho2(bet)

