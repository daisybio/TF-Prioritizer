package util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class MANN_WHITNEYU_PLOTS_FILES {

    public Options_intern options_intern;

    public File background;
    public File tf_data;

    public File f_output_script;
    public File f_output_plot;

    public boolean hm = false;

    public void write_and_execute_scripts() throws IOException {
        StringBuilder script = new StringBuilder();
        script.append("library(ggplot2)\n" +
                "library(plyr)\n" +
                "library(dplyr)\n\n");
        script.append("setwd(\""+f_output_plot.getAbsolutePath()+"\")\n");
        script.append("background <- read.csv(\""+background.getAbsolutePath()+"\",sep=\"\\t\",comment.char = \"#\")\n");
        script.append("background$POPULATION <- rep(\"BACKGROUND\",nrow(background))\n");
        script.append("background_scores <- background %>% select(TF_TG_SCORE)\n");
        script.append("background_ecdf <- ecdf(background_scores$TF_TG_SCORE)\n");
        script.append("\n");
        for(File tf : tf_data.listFiles())
        {
            String tf_name = tf.getName().split("\\.")[0].split("_")[0];

            script.append("#######################################################################\n");
            script.append("#####START#############"+tf_name+"#####################################\n");
            script.append("#######################################################################\n");

            script.append(tf_name+"_distribution<-read.csv(\""+tf.getAbsolutePath()+"\", sep=\"\\t\", comment.char = \"#\")\n");
            script.append(tf_name+"_distribution$POPULATION<-rep(\""+tf_name+"\",nrow("+tf_name+"_distribution))\n");
            script.append(tf_name+"_distribution_cumm <- rbind("+tf_name+"_distribution,background)\n");
            script.append(tf_name+"_mu <- ddply("+tf_name+"_distribution_cumm, \"POPULATION\", summarise, grp.mean=mean(TF_TG_SCORE))\n");
            script.append(tf_name+"_p <- ggplot("+tf_name+"_distribution_cumm, aes(x=TF_TG_SCORE, color = POPULATION, fill=POPULATION,)) + \n" +
                    "  geom_histogram( position=\"dodge\",binwidth=1000,alpha=0.4) +\n" +
                    "  geom_vline(data="+tf_name+"_mu, aes(xintercept=grp.mean, color=POPULATION),\n" +
                    "             linetype=\"dashed\", size = 1.5) +\n" +
                    "  labs(title = \"Distribution of TF-TG-Scores for background and "+tf_name+"\\n\")\n");
            script.append(tf_name+"_p + theme(panel.background = element_rect(fill = \"#F7F7F7\", colour = \"#D0D0D0\",\n" +
                    "                                          size = 2, linetype = \"solid\"),\n" +
                    "          panel.grid.major = element_line(size = 0.5, linetype = 'solid',\n" +
                    "                                          colour = \"#D0D0D0\"), \n" +
                    "          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',\n" +
                    "                                          colour = \"#D0D0D0\"),\n" +
                    "          axis.text.x = element_text(angle=30, hjust=1)) + xlim(-1,200000) + ylim(0,15000) + xlab(\"TF-TG Score\") + scale_fill_discrete(name = \"Population\", labels = c(\"Background\", \""+tf_name+"\")) + scale_color_discrete(name = \"Population\", labels = c(\"Background\", \""+tf_name+"\"))");
            script.append("\n");

            script.append("\n" +
                    tf_name+"_scores <- "+tf_name+"_distribution %>% select(TF_TG_SCORE)\n" +
                    tf_name+"_quantiles <- quantile("+tf_name+"_scores$TF_TG_SCORE, probs = c(.25,.75))\n" +
                    tf_name+"_iqr <-("+tf_name+"_scores$TF_TG_SCORE)\n" +
                    tf_name+"_up <- "+tf_name+"_quantiles[2]+1.5*"+tf_name+"_iqr\n" +
                    tf_name+"_down <- "+tf_name+"_quantiles[1]-1.5*"+tf_name+"_iqr\n");

            script.append("\n" +
                    tf_name+"_ecdf <- ecdf("+tf_name+"_scores$TF_TG_SCORE)\n");

            script.append("plot("+tf_name+"_scores, col = \"#F8766D\", main = \"Distribution of scores from PPI and random network populations\", xlab = \"score\", ylab = \"probability\")\n" +
                    "plot(background_ecdf, col = \"#00BFC4\", add = TRUE, lty=2)\n" +
                    "legend(1, 95, legend=c(\"PPI\", \"RANDOM\"),\n" +
                    "       col=c(\"#F8766D\", \"#00BFC4\"), lty=1:2, cex=0.8)\n");

            script.append("#######################################################################\n");
            script.append("#####END##############"+tf_name+"#####################################\n");
            script.append("#######################################################################\n");
        }



        File f_script = new File(f_output_script.getAbsolutePath()+File.separator+options_intern.file_suffix_distribution_analysis_mann_whitneyU_plot_scripts);
        BufferedWriter bw_script = new BufferedWriter(new FileWriter(f_script));
        bw_script.write(script.toString());
        bw_script.close();

    }

}
