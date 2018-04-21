ksource <- function(x, ...) {
  library(knitr)
  source(purl(x, output = tempfile()), ...)
}

list.files("~/workspace/permuting_accuracy/Round3/Output/", full.names = TRUE) %>% 
  list %>%
  do.call(file.remove,.)

files <- list.files('~/workspace/permuting_accuracy/Round3', pattern = '.Rmd', full.names = TRUE)
for(.file in files){
  # .file <- files[[2]]
  ksource(.file)
}
