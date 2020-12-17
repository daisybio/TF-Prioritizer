package util;

public class Options_intern
{
    /**
     * REQUIRED OPTIONS
     */
    public String input_dir_peaks ="";
    public String input_reference_genome ="";
    public String output_dir = "";
    public String gene_annotation = "";

    /**
     * OPTIONAL OPTIONS
     */
    public int window = 50000;
    public int cores = 1;

    /**
     * FLAGS
     */
    public boolean run_all_subdirectories = false;
}
