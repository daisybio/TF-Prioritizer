package util.Configs.Tools;

import util.Configs.Config;

import java.io.File;

public class General extends AbstractTool
{
    public final Config<Boolean> fileLogging = new Config<>(true);
    public final Config<Boolean> ensgMappingEnabled = new Config<>(true);
    public final Config<Boolean> calculateTpmLengthsEnabled = new Config<>(true);
    public final Config<Boolean> calculateGenePositionsEnabled = new Config<>(true);

    public final Config<String> differentTps = new Config<>("DIFFERENT_TPS");
    public final Config<String> sameTps = new Config<>("SAME_TPS");

    public final Config<File> d_ext =
            new Config<>(new File(sourceDirectory.getAbsolutePath() + File.separator + "ext"));


    public General(File workingDirectory, File sourceDirectory) throws IllegalAccessException
    {
        super(workingDirectory, sourceDirectory);
        registerEntries();
    }
}
