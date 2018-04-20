setwd('~/workspace/permuting_accuracy/Round3/')
files <- list.files(pattern = '.Rmd')
base::sort(files)
for(.file in files){
  # .file <- files[[2]]
  source(.file)
}
