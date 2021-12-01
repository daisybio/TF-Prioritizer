package conversion_helpers;

import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ENdb_conversion
{
    public static void main(String[] args) throws IOException
    {
        File f_input_download_endb = new File("/home/markus/Downloads/ENdb_enhancer.txt");
        File f_output_endb = new File("/home/markus/IdeaProjects/COM2POSE/ext/Enhancers_DB/ENdb_enhancer.bed");

        //convert into bedfile
        StringBuilder sb_output = new StringBuilder();
        sb_output.append("chrom\tchromStart\tchromEnd\tname\treference_genome\n");

        BufferedReader br_input = new BufferedReader(new FileReader(f_input_download_endb));
        String line_input = br_input.readLine();
        String[] split_header = line_input.split("\t");
        while((line_input=br_input.readLine())!=null)
        {
            String[] split = line_input.split("\t");

            if(split.length<7 || !line_input.startsWith("E") || !split[4].startsWith("chr"))
            {
                continue;
            }

            //check for special characters
            Pattern p = Pattern.compile("[^a-z0-9 ]", Pattern.CASE_INSENSITIVE);
            Matcher m = p.matcher(split[2]);
            boolean b = m.find();

            sb_output.append(split[4]);
            sb_output.append("\t");
            sb_output.append(split[5]);
            sb_output.append("\t");
            sb_output.append(split[6]);
            sb_output.append("\t");
            if(b)
            {
                sb_output.append("ENDb:"+split[0]+";E_Symbol:--");

            }
            else
            {
                sb_output.append("ENDb:"+split[0]+";E_Symbol:"+split[2]);
            }
            sb_output.append("\t");
            sb_output.append(split[3]);
            sb_output.append("\n");
        }
        br_input.close();

        BufferedWriter bw_out = new BufferedWriter(new FileWriter(f_output_endb));
        bw_out.write(sb_output.toString());
        bw_out.close();

    }
}
