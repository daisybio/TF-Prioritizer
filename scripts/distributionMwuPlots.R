library(ggplot2)
library(plyr)
library(dplyr)

background <- read.csv("{BACKGROUNDFILE}", sep = "\t", comment.char = "#")
background$POPULATION <- rep("BACKGROUND", nrow(background))
background_scores <- background %>% select(TF_TG_SCORE)
background_ecdf <- ecdf(background_scores$TF_TG_SCORE)

generate <- function(inputFile, name) {
  distribution <- read.csv(inputFile, sep = "\t", comment.char = "#")

  distribution$POPULATION <- rep(name, nrow(distribution))
  cummulatedDistribution <- rbind(distribution, background)
  mu <- ddply(cummulatedDistribution, "POPULATION", summarise, grp.mean = mean(TF_TG_SCORE))
  p <- ggplot(cummulatedDistribution, aes(x = TF_TG_SCORE, color = POPULATION, fill = POPULATION,)) +
    geom_histogram(position = "dodge", binwidth = 1000, alpha = 0.4) +
    geom_vline(data = mu, aes(xintercept = grp.mean, color = POPULATION),
               linetype = "dashed", size = 1.5) +
    labs(title = paste("Distribution of TF-TG-Scores for background and ", name, "\n", sep=""))
  p +
    theme(panel.background = element_rect(fill = "#F7F7F7", colour = "#D0D0D0",
                                          size = 2, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "#D0D0D0"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "#D0D0D0"),
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    xlim(-1, 200000) +
    ylim(0, 15000) +
    xlab("TF-TG Score") +
    scale_fill_discrete(name = "Population", labels = c("Background", name)) +
    scale_color_discrete(name = "Population", labels = c("Background", name))

  scores <- distribution %>% select(TF_TG_SCORE)

  plot(scores, col = "#F8766D", main = "Distribution of scores from PPI and random network populations", xlab = "score", ylab = "probability")
  plot(background_ecdf, col = "#00BFC4", add = TRUE, lty = 2)
  legend(1, 95, legend = c("PPI", "RANDOM"),
         col = c("#F8766D", "#00BFC4"), lty = 1:2, cex = 0.8)
}

{ CALLS }
