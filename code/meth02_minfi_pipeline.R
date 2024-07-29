.libPaths( c('.Rpackages', .libPaths() ) )

#' `minfi` is one of the most commonly used R packages for working with Infinium
#' Methylation BeadChips. In the following we import and preprocess the same dataset
#' as in the first script but using `minfi`, and point out some differences to
#' `ewastools`.

options(warn=-1)
library(data.table)
library(magrittr)
options(warn=0)

pheno = fread("data/pheno_clean.csv")

#' Reminder: this is the `ewastools` pipeline (without QC)

library(ewastools)

beta =
	paste0("data/", pheno$gsm) %>% # idat file paths
	read_idats(quiet=TRUE) %>%    # import fluorescence intensities
	detectionP %>%                # compute detection p-values
	mask(0.01) %>%                # set undetected probes to missing
	correct_dye_bias %>%
	dont_normalize                # calculate beta-values

dim(beta)

#' There are some name collisions between `minfi` and `ewastools` (function 
#' definitions with the same name). We therefore detach `ewastools` (the
#' inverse of `library()`), before we proceed.
detach("package:ewastools")

#' In order to still be able to call functions from a package that is not loaded
#' we can use the namespace operator `::`.
#' The syntax is `package_name::function_name()`

# ------------------------------------------------------------------------------
#' ## Pre-processing with `minfi`
options(warn=-1)
suppressMessages(library(minfi))
suppressMessages(library(ENmix))
suppressMessages(library(IlluminaHumanMethylation450kmanifest))
suppressMessages(library(IlluminaHumanMethylation450kanno.ilmn12.hg19))
suppressMessages(library(IlluminaHumanMethylationEPICmanifest))
suppressMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
suppressMessages(library(FlowSorted.Blood.450k))
suppressMessages(library(FlowSorted.Blood.EPIC))
options(warn=0)
idat_files = paste0("data/", pheno$gsm)

#' importing idat files, result is a RGChannelSet
rgset = read.metharray(basenames=idat_files)
rgset

#' `preprocessRaw` computes beta-value based on raw fluorescence intensities
methylset = preprocessRaw(rgset)
methylset

#' Many other preprocessing methods are available
#' 
#' * `preprocessNoob`: background subtraction/correction
#' * `preprocessFunnorm`: optionally includes preprocessNoob
#' * `preprocessIllumina`: normalization as in GenomeStudio software
#' * `preprocessQuantile`: quantile normalization stratified by U/M signal and
#'       probe type design
#' * `preprocessSWAN`: quantile normalization stratified by probe type design and
#'       number of CpG sites in probe sequence
#' * `preprocessENmix`: background correction (`ENmix` package)

processed = preprocessIllumina(rgset[,1:2])

#' Due to their design, Type II probes feature more background noise, shifting the
#' peaks for completely (un)methylated Cpg sites inwards
plotBetasByType(getBeta(processed)[,1], probeTypes=getAnnotation(processed))

#' `rcp()` (regression on correlated probes) uses Type I probes to shift the
#' distribution of beta-values of Type II probes. Compared to above plot, the peaks
#' should be aligned now.
plotBetasByType(    rcp(processed)[,1], probeTypes=getAnnotation(processed))


#' ## Cell type prediction
#' `minfi` can estimate the cell proportions of three different tissues: "Blood",
#' "CordBlood", or "DLPFC" (frontal cortex). Estimates for "Blood" are based on
#' the Reinius reference dataset. You have to use an `rgset` to estimate cell
#' proportions with `minfi` as the user-provided dataset is normalized together
#' with the reference dataset of purified cell types.  
#' Note: the first time running estimateCellCounts2 locally, you will be prompted 
#' to create a directory for ExperimentHub. 
#' Enter "yes" in the console to proceed 
 
minfi.LC = estimateCellCounts2(rgset, compositeCellType="Blood")
minfi.LC = minfi.LC$prop

#' **QUESTION:** What normalization is used by default?
#'
#' **QUESTION:** Proportions for which cell types can be estimated
#' when `compositeCellType=="CordBlood"`?
#'
#' Comparison with `ewastools` shows that the estimates are not identical but
#' highly correlated
ewastools.LC = ewastools::estimateLC(beta, ref="Reinius")

cor(ewastools.LC$CD4, minfi.LC[,"CD4T"] )
cor(ewastools.LC$GR , minfi.LC[,"Neu"]  )
cor(ewastools.LC$MO , minfi.LC[,"Mono"] )
cor(ewastools.LC$B  , minfi.LC[,"Bcell"])

par(mfrow = c(2, 2))
plot(ewastools.LC$CD4, minfi.LC[,"CD4T"] )
plot(ewastools.LC$GR , minfi.LC[,"Neu"]  )
plot(ewastools.LC$MO , minfi.LC[,"Mono"] )
plot(ewastools.LC$B  , minfi.LC[,"Bcell"])
dev.off()

#' In case you have samples measured on both the 450K and EPIC chips, `minfi` can
#' virtually convert them in both directions.  
convertArray(rgset, outType="IlluminaHumanMethylation450k")

#' clear environment
rm(list = ls()); gc()
