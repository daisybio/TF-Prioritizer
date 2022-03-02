package util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class IGV_Headless implements Runnable
{
    private final StringBuilder commandBuilder = new StringBuilder();
    private final Options_intern options_intern;
    private File batchFile;

    public IGV_Headless(Options_intern options_intern)
    {
        this.options_intern = options_intern;
        addCommand("new");
        addCommand("maxPanelHeight 1080");
    }

    public void addCommand(String command)
    {
        commandBuilder.append(command).append("\n");
    }

    public void save(File file) throws IOException
    {
        addCommand("exit");
        batchFile = file;
        BufferedWriter writer = new BufferedWriter(new FileWriter(batchFile));
        writer.write(commandBuilder.toString());
        writer.close();
    }

    @Override public void run()
    {
        String command = "xvfb-run --auto-servernum --server-num=1 " + options_intern.igv_path_to_igv + "/igv.sh -b " +
                batchFile.getAbsolutePath();

        Process child = null;
        try
        {
            child = Runtime.getRuntime().exec(command);
        } catch (IOException e)
        {
            e.printStackTrace();
        }
        int code = 0;
        try
        {
            code = child.waitFor();
        } catch (InterruptedException e)
        {
            e.printStackTrace();
        }
        switch (code)
        {
            case 0:
                break;
            case 1:
                String message = child.getErrorStream().toString();
                System.out.println(message);
        }
    }
}
