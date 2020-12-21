package util;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;

public class Logger
{
    boolean write_to_file;
    File logfile;

    public Logger(boolean write_to_file, String working_dir) throws IOException {
        this.write_to_file=write_to_file;

        if(write_to_file)
        {
            logfile = new File(working_dir+File.separator+"logfile.txt");
            BufferedWriter bw = new BufferedWriter(new FileWriter(logfile));
            bw.write("");
            bw.close();
        }

    }

    public void logLine(String message) throws IOException {
        SimpleDateFormat formatter = new SimpleDateFormat("dd/MM/yyyy HH:mm:ss");
        Date date = new Date();

        String line = "["+formatter.format(date)+"]\t"+message;
        System.out.println(line);

        if(write_to_file)
        {
            BufferedWriter bw = new BufferedWriter(new FileWriter(logfile,true));
            bw.write(line);
            bw.newLine();
            bw.close();
        }
    }
}
