package com2pose;

import util.Options_intern;

import java.io.File;

public class COM2POSE_lib
{
    Options_intern options_intern;

    public COM2POSE_lib(Options_intern options_intern)
    {
        this.options_intern = options_intern;
    }

    public void run()
    {


        //process peak file or if directory all files in directory
        File folder = new File(options_intern.input_dir_peaks);
        run_all_files_subdirs(folder);

    }

    /**
     * run COM2POSE for file
     * @param file_in
     */
    private void run_single(File file_in)
    {
        //The main program

    }

    /**
     * recursive method - runs all files in a directory
     * @param folder input peaks file
     */
    private void run_all_files_subdirs(File folder)
    {
        if(folder.isDirectory())
        {
            for(File fileDir : folder.listFiles())
            {
                if(fileDir.isDirectory()&& options_intern.run_all_subdirectories)
                {
                    run_all_files_subdirs(fileDir);
                }
                else if (!fileDir.isDirectory())
                {
                    run_single(fileDir);
                }
            }
        }
        else
        {
            run_single(folder);
        }

    }



}
