ksource <- function(x, ...) {
  library(knitr)
  source(purl(x, output = tempfile()), ...)
}

library(magrittr)

# Automate #####

# Remove existing files
list.files("~/workspace/permuting_accuracy/Round4/Output/", full.names = TRUE) %>%
  list %>%
  do.call(file.remove,.)

# Run with exceptions
files <- list.files('~/workspace/permuting_accuracy/Round4', pattern = '.Rmd', full.names = TRUE)
exceptions <- c('file1.Rmd', 'file39.Rmd')
files <- files[!grepl(exceptions, files, fixed = TRUE)]

for(.file in files){
  try({ksource(.file))
}
pbPost(type = "note", body = "Done", "All files finished")
