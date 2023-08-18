suppressMessages(library(Rsubread))

setwd('/scratch/gds44474/MIMULUS/RNAseq_endosperm_til/data/star_pass2')

bam_names_caes <- c('TWN36xUTC1_1_STAR_LVR_v1_MD_Split_Q60.bam','UTC1xTWN36_1_STAR_LVR_v1_MD_Split_Q60.bam','TWN36xUTC1_2_STAR_LVR_v1_MD_Split_Q60.bam','UTC1xTWN36_2_STAR_LVR_v1_MD_Split_Q60.bam','TWN36xUTC1_3_STAR_LVR_v1_MD_Split_Q60.bam','UTC1xTWN36_3_STAR_LVR_v1_MD_Split_Q60.bam') 

# to convert gff3 to gtf I used gffread/0.11.6-GCCcore-8.3.0
# gffread MtilingiivarLVRv1.1.primaryTrs.gff3  -T -o my.gtf
feature_counts_caes <- featureCounts(
  files = bam_names_caes,
  annot.ext = 'my.gtf',
  isGTFAnnotationFile = TRUE,  
  useMetaFeatures = TRUE,
  countMultiMappingReads = FALSE,
  ignoreDup = FALSE,
  strandSpecific = 2,
  isPairedEnd = TRUE,
  requireBothEndsMapped = TRUE,
  countChimericFragments = FALSE,
  autosort = TRUE,
  nthreads = 6,
  verbose = TRUE
)

saveRDS(feature_counts_caes,'endosperm_counts_caes_object')
bam_names_til <- c('LVR1xSOP12_1_STAR_LVR_v1_MD_Split_Q60.bam','SOP12xLVR1_1_STAR_LVR_v1_MD_Split_Q60.bam','LVR1xSOP12_2_STAR_LVR_v1_MD_Split_Q60.bam','SOP12xLVR1_2_STAR_LVR_v1_MD_Split_Q60.bam','LVR1xSOP12_3_STAR_LVR_v1_MD_Split_Q60.bam','SOP12xLVR1_3_STAR_LVR_v1_MD_Split_Q60.bam')

feature_counts_til <- featureCounts(
  files = bam_names_til,
  annot.ext = 'my.gtf',  
  isGTFAnnotationFile = TRUE,  
  useMetaFeatures = TRUE,
  countMultiMappingReads = FALSE,
  ignoreDup = FALSE,
  strandSpecific = 2,
  isPairedEnd = TRUE,
  requireBothEndsMapped = TRUE,
  countChimericFragments = FALSE,
  autosort = TRUE,
  nthreads = 6,
  verbose = TRUE
)

saveRDS(feature_counts_til,'endosperm_counts_til_object')
