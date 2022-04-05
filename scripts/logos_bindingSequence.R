library(DiffLogo)
library(seqLogo)

generate <- function(inputFile, plotFile, outputFile) {
    currentPwm <- getPwmFromFastaFile(inputFile)
    currentPwm <- seqLogo::makePWM(currentPwm)
    png(file=plotFile)
    seqLogo::seqLogo(currentPwm,fill=c(A="#4daf4a", C="#377eb8", G="#ffd92f", T="#e41a1c"))
    dev.off()
    currentPwm <- t(pwm(currentPwm))
    write.table(currentPwm, outputFile, row.names = FALSE, col.names=FALSE, quote = F, sep="\t")
}

{CALLS}
