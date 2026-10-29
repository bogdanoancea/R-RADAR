#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(RADAR))
tryCatch(radar_cli(commandArgs(trailingOnly = TRUE)), error = function(e) {
  message("RADAR: ", conditionMessage(e))
  quit(status = 1L)
})
