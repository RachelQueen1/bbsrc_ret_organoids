library(rmarkdown)
library(dplyr)

allData <- "/data/rachel/Linda_Lako/Spatial/spaceranger_outs_2021_102"
tissue <- "AD3"
projDir <- paste0("/data/rachel/Linda_Lako/Spatial/", tissue)
outputDir <- file.path(projDir, "/rObjects/")
## find markers
sampleNames <- paste0("AD3_D", c(10,20,35,60,90,150,210))



for (sampleName in sampleNames){
  thetitle <- paste(sampleName)
  day <- gsub("AD3_D", "", sampleName) %>% as.numeric()
  render('individual_results.rmd', output_file = paste0("AD3_results_", sampleName))
}
