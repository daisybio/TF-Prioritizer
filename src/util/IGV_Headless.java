package util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class IGV_Headless
{
    private final StringBuilder commandBuilder = new StringBuilder();
    private final Options_intern options_intern;
    private final String name;

    public IGV_Headless(Options_intern options_intern, String name)
    {
        this.options_intern = options_intern;
        this.name = name;
        addCommand("new");
    }

    public void addCommand(String command)
    {
        commandBuilder.append(command).append("\n");
    }

    private void save(File file) throws IOException
    {
        BufferedWriter writer = new BufferedWriter(new FileWriter(file));
        writer.write(commandBuilder.toString());
        writer.close();
    }

    public void run(File workingDirectory) throws Exception
    {
        addCommand("exit");

        File batchFile = new File(workingDirectory.getAbsolutePath() + File.separator + name + ".bat");

        save(batchFile);

        String command = "xvfb-run --auto-servernum --server-num=1 " + options_intern.igv_path_to_igv + "/igv.sh -b " +
                batchFile.getAbsolutePath();

        Process child = Runtime.getRuntime().exec(command);
        int code = child.waitFor();
        switch (code)
        {
            case 0:
                break;
            case 1:
                String message = child.getErrorStream().toString();
                throw new Exception(message);
        }
    }
}
