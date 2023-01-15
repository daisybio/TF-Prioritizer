#!/usr/bin/env Rscript

library(DiffLogo)
library(seqLogo)
library(argparse)

parser <- ArgumentParser()
parser$add_argument("--input_file", required = TRUE)
parser$add_argument("--plot_file", required = TRUE)
parser$add_argument("--output_file", required = TRUE)
args <- parser$parse_args()

if (file.size(args$input_file) == 0L) {
    print(paste0("Empty file: ", args$input_file))
} else {
    currentPwm <- getPwmFromFastaFile(args$input_file)
    currentPwm <- seqLogo::makePWM(currentPwm)
    png(file=args$plot_file)
    seqLogo::seqLogo(currentPwm,fill=c(A="#4daf4a", C="#377eb8", G="#ffd92f", T="#e41a1c"))
    dev.off()
    currentPwm <- t(pwm(currentPwm))
    write.table(currentPwm, args$output_file, row.names = FALSE, col.names=FALSE, quote = F, sep="\t")
}
