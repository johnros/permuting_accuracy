ksource <- function(x, ...) {
  library(knitr)
  source(purl(x, output = tempfile()), ...)
}

# Automate #####

# Remove existing files
# list.files("~/workspace/permuting_accuracy/Round3/Output/", full.names = TRUE) %>% 
#   list %>%
#   do.call(file.remove,.)

# Run with exceptions
files <- list.files('~/workspace/permuting_accuracy/Round3', pattern = '.Rmd', full.names = TRUE)
exceptions <- c('file1.Rmd')
files <- files[!grepl(exceptions, files, fixed = TRUE)]

for(.file in files){
  # .file <- files[[13]]
  try(ksource(.file))
  }
