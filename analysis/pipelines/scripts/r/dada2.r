#!/usr/bin/env Rscript4.3.1

msg <-"the script gets executed"

#library (tidyverse)
#library(devtools)
#devtools::install_github("benjjneb/dada2")
#library(dada2)
library(rmarkdown)


print (msg)

version$version.string

arguments <- commandArgs(TRUE)

output.folder<-arguments[1]

fastq.folder<-arguments[2]

scripts.folder <-arguments[3]

hashing <-arguments[4]

continuing <- arguments[5:8]


#list.files()
setwd (scripts.folder)
getwd ()
render("r/dada2.Rmd", output_file = paste0(output.folder,"/dada2_report.html"),
 params = list(folder = output.folder,
               fastqs = fastq.folder,
                 hash = hashing,
             # original = source,
                 cont = continuing))

