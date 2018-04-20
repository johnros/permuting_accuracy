ksource <- function(x, ...) {
  library(knitr)
  source(purl(x, output = tempfile()), ...)
}


setwd('~/workspace/permuting_accuracy/Round3/')
files <- list.files(pattern = '.Rmd')
base::sort(files)
for(.file in files){
  # .file <- files[[2]]
  ksource(.file)
}
