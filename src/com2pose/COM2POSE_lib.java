package com2pose;

import org.apache.commons.compress.compressors.lz4.BlockLZ4CompressorOutputStream;
import org.apache.commons.compress.utils.IOUtils;
import util.*;

import javax.net.ssl.HttpsURLConnection;
import javax.net.ssl.SSLContext;
import javax.net.ssl.TrustManager;
import javax.net.ssl.X509TrustManager;
import java.io.*;
import java.lang.annotation.Target;
import java.net.Socket;
import java.net.URL;
import java.net.URLConnection;
import java.nio.Buffer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.*;

public class COM2POSE_lib {
    public Options_intern options_intern;
    public Logger logger;

    //only used in create_overview_website()!!
    ArrayList<File> threshold_folders = new ArrayList<>();
    boolean threshold_folders_filled = false;

    /**
     * constructor for analysis programms, you cannot run the pipeline here
     *
     * @param options_intern
     * @param logger
     */
    public COM2POSE_lib(Options_intern options_intern, Logger logger) {
        this.options_intern = options_intern;
        this.logger = logger;
    }

    /**
     * constructor for the pipeline
     *
     * @param options_intern
     * @throws IOException
     */
    public COM2POSE_lib(Options_intern options_intern) throws IOException {
        this.options_intern = options_intern;

        logger = new Logger(options_intern.write_to_logfile, options_intern.com2pose_working_directory);

        logger.logLine("#########################################");
        logger.logLine("############## COM2POSE #################");
        logger.logLine("#########################################");
        logger.logLine("Working directory set to: " + options_intern.com2pose_working_directory);
        logger.logLine("COM2POSE path set to: " + options_intern.path_to_COM2POSE);
    }

    public void run_igv_on_important_loci() throws Exception {
        logger.logLine("[IGV] Run IGV on important loci.");

        check_tepic_input_with_options();

        HashMap<String, String> gene_to_coordinates = new HashMap<>();

        HashMap<String, String> ensg_gene_symbol_map = new HashMap<>();
        BufferedReader br_ensg_gene_symbol =
                new BufferedReader(new FileReader(new File(options_intern.tepic_ensg_symbol)));
        String line_ensg_gene_symbol = br_ensg_gene_symbol.readLine();
        while ((line_ensg_gene_symbol = br_ensg_gene_symbol.readLine()) != null) {
            String[] split = line_ensg_gene_symbol.split("\t");
            if (split.length > 1) {
                ensg_gene_symbol_map.put(split[0].toUpperCase(), split[1].toUpperCase());
            }
        }
        br_ensg_gene_symbol.close();

        File f_input_gene_to_coordinates = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.folder_name_deseq2_preprocessing_gene_positions + File.separator +
                options_intern.file_suffix_deseq2_preprocessing_gene_positions_data);

        BufferedReader br_gene_coordinates = new BufferedReader(new FileReader(f_input_gene_to_coordinates));
        String line_gene_coordinates = br_gene_coordinates.readLine();
        String[] header_gene_coordinates = line_gene_coordinates.split("\t");
        while ((line_gene_coordinates = br_gene_coordinates.readLine()) != null) {
            String[] split = line_gene_coordinates.split("\t");

            if (split[1].startsWith("CHR")) {
                continue;
            }

            String gene_name = split[0];
            //chr10:95,381,903-95,419,253
            //chr10:95,221,224-95,253,042
            //chr10:95,217,767-95,255,176
            //chr10:95,217,765-95,255,115
            if(gene_name.equals("Socs2"))
            {
                System.out.println("X");
            }

            String chr = "chr" + split[1];
            int begin = Integer.parseInt(split[2]);
            int end = Integer.parseInt(split[3]);

            begin -= 50000;
            end += 50000;

            String position_string = chr + ":" + begin + "-" + end;
            gene_to_coordinates.put(gene_name.toUpperCase(), position_string);
        }

        //create ouput folder
        File f_output_root =
                new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_igv+
                        File.separator+ options_intern.folder_out_igv_important_loci);
        f_output_root.mkdir();

        File root_copy =
                new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_target_genes_dcg);
        for(File f_tp : root_copy.listFiles())
        {
            if(f_tp.isDirectory())
            {
                String key_tp = f_tp.getName();

                for(File f_hm : f_tp.listFiles())
                {
                    if(f_hm.isDirectory())
                    {
                        String key_hm = f_hm.getName();

                        File f_output_tf_hm = new File(f_output_root.getAbsolutePath() + File.separator + key_hm);
                        f_output_tf_hm.mkdir();

                        File f_output_tf_hm_tp =
                                new File(f_output_tf_hm.getAbsolutePath() + File.separator + key_tp);
                        f_output_tf_hm_tp.mkdir();

                        String load_tf_chip_seq = "load ";
                        boolean added_prediction_data=false;

                        //include predictive HMs
                        if(options_intern.igv_include_prediction_data.size()>0)
                        {
                            File f_input =
                                    new File(options_intern.tepic_input_directory+File.separator + key_tp);

                            for(String including_hm : options_intern.igv_include_prediction_data)
                            {
                                File f_hm_input = new File(f_input.getAbsolutePath()+File.separator+including_hm);
                                if(f_hm_input.exists())
                                {
                                    for(File f_file : f_hm_input.listFiles())
                                    {
                                        if(added_prediction_data)
                                        {
                                            load_tf_chip_seq+="," + f_file.getAbsolutePath();
                                        }
                                        else
                                        {
                                            added_prediction_data=true;
                                            load_tf_chip_seq+= f_file.getAbsolutePath();
                                        }
                                    }
                                }
                            }
                        }

                        //add all own TF ChIP-seq if available
                        if(!options_intern.igv_path_to_tf_chip_seq.equals(""))
                        {
                            File f_input_tf_chip_seq = new File(options_intern.igv_path_to_tf_chip_seq);
                            File f_input_tf_chip_seq_for_now =
                                    new File(f_input_tf_chip_seq.getAbsolutePath()+File.separator+ key_tp);
                            if(f_input_tf_chip_seq_for_now.exists())
                            {
                                if(f_input_tf_chip_seq_for_now.isDirectory())
                                {
                                    for(File f_tf_chipseq : f_input_tf_chip_seq_for_now.listFiles())
                                    {
                                        if(f_tf_chipseq.isDirectory())
                                        {

                                            for(File f_tf_chipseq_TF_file:f_tf_chipseq.listFiles() )
                                            {
                                                if(f_tf_chipseq_TF_file.isFile())
                                                {
                                                    if(added_prediction_data)
                                                    {
                                                        load_tf_chip_seq+="," + f_tf_chipseq_TF_file.getAbsolutePath();
                                                    }
                                                    else
                                                    {
                                                        added_prediction_data=true;
                                                        load_tf_chip_seq+= f_tf_chipseq_TF_file.getAbsolutePath();
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        //add all ChIP-ATLAS ChIP-seq if available
                        File f_chip_atlas_tf_chipseq_data_root =
                                new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_chip_atlas+File.separator+options_intern.folder_out_chip_atlas_peak_files);

                        if(f_chip_atlas_tf_chipseq_data_root.exists())
                        {
                            for(File f_tf : f_chip_atlas_tf_chipseq_data_root.listFiles())
                            {
                                if(f_tf.isDirectory())
                                {

                                    for(File f_tf_file : f_tf.listFiles())
                                    {
                                        if(f_tf_file.isFile())
                                        {
                                            String file_name = f_tf_file.getName();
                                            if(file_name.endsWith("idx"))
                                                continue;

                                            if(added_prediction_data)
                                            {
                                                load_tf_chip_seq+="," + f_tf_file.getAbsolutePath();
                                            }
                                            else
                                            {
                                                added_prediction_data=true;
                                                load_tf_chip_seq+= f_tf_file.getAbsolutePath();
                                            }

                                        }
                                    }

                                }
                            }
                        }

                        if(!added_prediction_data)
                            continue;

                        //save session
                        File f_save_session =
                                new File(f_output_tf_hm_tp.getAbsolutePath()+File.separator+options_intern.file_suffix_session);
                        igv_save_sessions(f_save_session,load_tf_chip_seq);

                        Socket socket_session = new Socket("127.0.0.1", options_intern.igv_port_number);
                        PrintWriter out_session = new PrintWriter(socket_session.getOutputStream(), true);
                        BufferedReader in_session = new BufferedReader(new InputStreamReader(socket_session.getInputStream()));

                        out_session.println("genome " + options_intern.igv_species_ref_genome);
                        String response_session = in_session.readLine();

                        //logger.logLine("[IGV] " + load_tf_chip_seq);
                        out_session.println(load_tf_chip_seq);
                        String response_load_sessoin = in_session.readLine();

                        //do pics of important loci
                        for(String locus: options_intern.igv_important_locus_all_prio_tf)
                        {
                            if (!gene_to_coordinates.containsKey(locus))
                            {
                                continue;
                            }

                            String snapshot_name = locus + ".png";
                            String response ="";

                            File f_output_shot = new File(f_output_tf_hm_tp.getAbsolutePath());
                            f_output_shot.mkdir();
                            //logger.logLine("[IGV] "+ "snapshotDirectory "+f_output_shot.getAbsolutePath());
                            out_session.println("snapshotDirectory " + f_output_shot.getAbsolutePath());
                            response = in_session.readLine();
                            //logger.logLine("[IGV] " + response);
                            //logger.logLine("[IGV] "+ "goto " + gene_to_coordinates.get(targets));
                            out_session.println("goto " + gene_to_coordinates.get(locus));
                            response = in_session.readLine();
                            //logger.logLine("[IGV] " + response);
                            //logger.logLine("[IGV] "+ "sort");
                            //out.println("sort");
                            //response=in.readLine();
                            //logger.logLine("[IGV] " + response);
                            //logger.logLine("[IGV] "+ "collapse");
                            //out.println("collapse");
                            //response=in.readLine();
                            //logger.logLine("[IGV] " + response);
                            //logger.logLine("[IGV] "+ "snapshot");
                            out_session.println("snapshot");
                            response = in_session.readLine();
                            //logger.logLine("[IGV] " + response);

                            File f_snapshot_made = new File("NO");

                            //search png name
                            for (File f_snapshot : f_output_shot.listFiles()) {
                                if (f_snapshot.isFile()) {
                                    if (f_snapshot.getName().startsWith("chr")) {
                                        f_snapshot_made = f_snapshot;
                                        break;
                                    }

                                }
                            }

                            if (!f_snapshot_made.exists()) {
                                //its the all overview
                                //logger.logLine("[ERROR]: trying to move file which does not exist!");
                            } else {

                                String command_mv = "mv " + f_snapshot_made.getAbsolutePath() + " " +
                                        f_snapshot_made.getParentFile().getAbsolutePath() + File.separator +
                                        snapshot_name;


                                Process child = Runtime.getRuntime().exec(command_mv);
                                int code = child.waitFor();
                                switch (code) {
                                    case 0:
                                        break;
                                    case 1:
                                        String message = child.getErrorStream().toString();
                                        throw new Exception(message);
                                }


                            }
                        }


                        out_session.println("new");
                        response_session = in_session.readLine();

                        socket_session.close();
                    }
                }
            }
        }




        logger.logLine("[IGV] Finished run IGV on important loci.");


    }

    /**
     * run IGV on own TF ChIP-seq data.
     *
     * @throws Exception
     */
    public void run_igv_own_data() throws Exception {
        logger.logLine(
                "[IGV] start taking screenshots for top TFs and their corresponding top target genes (OWN TF DATA)");

        check_tepic_input_with_options();

        HashMap<String, String> gene_to_coordinates = new HashMap<>();

        HashMap<String, String> ensg_gene_symbol_map = new HashMap<>();
        BufferedReader br_ensg_gene_symbol =
                new BufferedReader(new FileReader(new File(options_intern.tepic_ensg_symbol)));
        String line_ensg_gene_symbol = br_ensg_gene_symbol.readLine();
        while ((line_ensg_gene_symbol = br_ensg_gene_symbol.readLine()) != null) {
            String[] split = line_ensg_gene_symbol.split("\t");
            if (split.length > 1) {
                ensg_gene_symbol_map.put(split[0].toUpperCase(), split[1].toUpperCase());
            }
        }
        br_ensg_gene_symbol.close();

        File f_input_gene_to_coordinates = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.folder_name_deseq2_preprocessing_gene_positions + File.separator +
                options_intern.file_suffix_deseq2_preprocessing_gene_positions_data);
        BufferedReader br_gene_coordinates = new BufferedReader(new FileReader(f_input_gene_to_coordinates));
        String line_gene_coordinates = br_gene_coordinates.readLine();
        String[] header_gene_coordinates = line_gene_coordinates.split("\t");
        while ((line_gene_coordinates = br_gene_coordinates.readLine()) != null) {
            String[] split = line_gene_coordinates.split("\t");

            if (split[1].startsWith("CHR")) {
                continue;
            }

            String gene_name = split[0];
            String chr = "chr" + split[1];
            int begin = Integer.parseInt(split[2]);
            int end = Integer.parseInt(split[3]);

            begin -= 50000;
            end += 50000;

            String position_string = chr + ":" + begin + "-" + end;
            gene_to_coordinates.put(gene_name.toUpperCase(), position_string);
        }

        /*
        //read gtf for genome coordinates
        HashMap<String,String> gene_to_coordinates = new HashMap<>();
        BufferedReader br_gene_coordinates = new BufferedReader(new FileReader(new File(options_intern.tepic_gene_annot)));
        String line_gene_coordinates = br_gene_coordinates.readLine();
        while((line_gene_coordinates=br_gene_coordinates.readLine())!=null)
        {
            if(line_gene_coordinates.startsWith("#"))
            {
                continue;
            }

            String[] split = line_gene_coordinates.split("\t");

            if(!split[2].equals("gene"))
            {
                continue;
            }

            String chr = split[0];
            String position = split[3]+"-"+split[4];
            String gene_name = "";

            String[] gene_name_prep = split[8].split(";");
            for(String key : gene_name_prep)
            {
                key=key.trim();
                if(key.startsWith("gene_name"))
                {
                    String[] split_intern = key.split("\"");
                    gene_name=split_intern[1].toUpperCase();
                }
            }
            gene_to_coordinates.put(gene_name,chr+":"+position);
        }
        */

        //get most significant TFs
        ArrayList<String> tfs_in_order = new ArrayList<>();
        BufferedReader br_tfs = new BufferedReader(new FileReader(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution +
                        File.separator + options_intern.folder_out_distribution_dcg + File.separator +
                        options_intern.file_suffix_distribution_analysis_dcg));
        String line_tfs = br_tfs.readLine();
        while ((line_tfs = br_tfs.readLine()) != null) {
            String[] split = line_tfs.split("\t");
            tfs_in_order.add(split[1]);
        }

        //get target genes for tfs
        HashMap<String, HashMap<String, HashMap<String, ArrayList<String>>>> hm_timepoint_tf_target_genes =
                new HashMap<>();
        File f_root_target_genes = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_out_target_genes_dcg);

        //prefill map
        for (File f_target_genes_tp : f_root_target_genes.listFiles()) {
            if (f_target_genes_tp.isDirectory()) {
                for (File f_target_genes_tp_hm : f_target_genes_tp.listFiles()) {
                    if (f_target_genes_tp_hm.isDirectory()) {
                        HashMap<String, HashMap<String, ArrayList<String>>> tp_tf_target_genes = new HashMap<>();
                        hm_timepoint_tf_target_genes.put(f_target_genes_tp_hm.getName(), tp_tf_target_genes);
                    }
                }
            }
        }
        for (File f_target_genes_tp : f_root_target_genes.listFiles()) {
            String hm = "";
            if (f_target_genes_tp.isDirectory()) {
                for (File f_target_genes_tp_hm : f_target_genes_tp.listFiles()) {
                    if (f_target_genes_tp_hm.isDirectory()) {
                        hm = f_target_genes_tp_hm.getName();
                        HashMap<String, ArrayList<String>> tf_target_genes = new HashMap<>();
                        hm_timepoint_tf_target_genes.get(hm).put(f_target_genes_tp.getName(), tf_target_genes);

                    }
                }
            }
        }

        for (File f_target_genes_tp : f_root_target_genes.listFiles()) {
            if (f_target_genes_tp.isDirectory()) {
                for (File f_target_genes_tp_hm : f_target_genes_tp.listFiles()) {
                    if (f_target_genes_tp_hm.isDirectory()) {
                        HashMap<String, ArrayList<String>> tf_target_genes = new HashMap<>();
                        hm_timepoint_tf_target_genes.get(f_target_genes_tp_hm.getName())
                                .put(f_target_genes_tp.getName(), tf_target_genes);

                        for (File f_target_genes_tp_hm_tf : f_target_genes_tp_hm.listFiles()) {
                            if (f_target_genes_tp_hm_tf.isDirectory()) {
                                String tf_name = f_target_genes_tp_hm_tf.getName();

                                ArrayList<String> target_genes = new ArrayList<>();

                                tf_target_genes.put(tf_name, target_genes);

                                BufferedReader br = new BufferedReader(new FileReader(new File(
                                        f_target_genes_tp_hm_tf.getAbsolutePath() + File.separator + tf_name +
                                                ".csv")));
                                String line = br.readLine();
                                int count = 0;
                                while ((line = br.readLine()) != null) {
                                    String[] split = line.split("\t");
                                    target_genes.add(split[0]);

                                    if (count == options_intern.plot_top_k_genes - 1) {
                                        break;
                                    }

                                    count++;
                                }
                                br.close();
                            }
                        }
                    }
                }
            }
        }


        /*
        File f_root_target_genes = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_target_genes);
        for(File f_hm : f_root_target_genes.listFiles())
        {
            String hm = f_hm.getName();
            HashMap<String,HashMap<String,ArrayList<String>>> hm_hm;
            if(hm_timepoint_tf_target_genes.containsKey(hm))
            {
                hm_hm=hm_timepoint_tf_target_genes.get(hm);
            }
            else
            {
                hm_hm = new HashMap<>();
                hm_timepoint_tf_target_genes.put(hm,hm_hm);
            }

            if(f_hm.isDirectory())
            {
                File f_hm_cutoff = new File(f_hm.getAbsolutePath()+File.separator+"0.1");
                if(f_hm_cutoff.exists() && f_hm_cutoff.isDirectory())
                {
                    File f_hm_cutoff_all_diff_data = new File(f_hm_cutoff.getAbsolutePath()+File.separator+options_intern.folder_out_target_genes_all_different);
                    if(f_hm_cutoff_all_diff_data.exists() && f_hm_cutoff_all_diff_data.isDirectory())
                    {
                        for(File f_hm_cutoff_all_diff_data_tp : f_hm_cutoff_all_diff_data.listFiles())
                        {
                            String tp = f_hm_cutoff_all_diff_data_tp.getName();

                            HashMap<String,ArrayList<String>> hm_tp;
                            if(hm_hm.containsKey(tp))
                            {
                                hm_tp=hm_hm.get(tp);
                            }
                            else
                            {
                                hm_tp=new HashMap<>();
                                hm_hm.put(tp,hm_tp);
                            }

                            if(f_hm_cutoff_all_diff_data_tp.isDirectory())
                            {
                                for(String tfs: tfs_in_order)
                                {
                                    File f_tf_target_genes_input = new File(f_hm_cutoff_all_diff_data_tp.getAbsolutePath()+File.separator+tfs+".csv");
                                    if(f_tf_target_genes_input.exists())
                                    {
                                        ArrayList<String> target_genes_ordered = new ArrayList<>();
                                        hm_tp.put(tfs,target_genes_ordered);
                                        which_tfs.add(tfs);

                                        BufferedReader br_tg = new BufferedReader(new FileReader(f_tf_target_genes_input));
                                        String line_tg = br_tg.readLine();
                                        while((line_tg=br_tg.readLine())!=null)
                                        {
                                            String[] split = line_tg.split("\t");
                                            target_genes_ordered.add(split[1].toUpperCase());
                                        }

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }*/

        //read in all important tfs per hm
        File f_important_tfs_per_hm = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution +
                        File.separator + options_intern.folder_out_distribution_stats + File.separator +
                        options_intern.folder_out_distribution_stats_HM);

        HashMap<String, HashSet<String>> hm_to_important_tfs = new HashMap<>();
        for (File fileDir_hm : f_important_tfs_per_hm.listFiles()) {
            if (fileDir_hm.isDirectory()) {
                HashSet<String> important_tfs = new HashSet<>();
                hm_to_important_tfs.put(fileDir_hm.getName(), important_tfs);

                BufferedReader br = new BufferedReader(new FileReader(new File(
                        fileDir_hm.getAbsolutePath() + File.separator +
                                options_intern.file_suffix_distribution_analysis_plot_stats)));
                String line = br.readLine();
                line = br.readLine();
                while ((line = br.readLine()) != null) {
                    String[] split = line.split("\t");
                    important_tfs.add(split[1]);
                }
            }
        }


        logger.logLine("[IGV] taking pictures ...");
        /*
        http://software.broadinstitute.org/software/igv/automation
         */

        //fill which tfs from all dcg
        HashSet<String> which_tfs = new HashSet<>();
        HashMap<String, Integer> tf_to_dcg_place = new HashMap<>();

        File f_input_all_dcg = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution +
                        File.separator + options_intern.folder_out_distribution_dcg + File.separator +
                        options_intern.file_suffix_distribution_analysis_dcg);
        BufferedReader br_input_all_dcg = new BufferedReader(new FileReader(f_input_all_dcg));
        String line_input_all_dcg = br_input_all_dcg.readLine();
        while ((line_input_all_dcg = br_input_all_dcg.readLine()) != null) {
            String[] split = line_input_all_dcg.split("\t");
            which_tfs.add(split[1]);

            tf_to_dcg_place.put(split[1], Integer.parseInt(split[0]));
        }
        br_input_all_dcg.close();

        File f_output_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_igv +
                        File.separator + options_intern.folder_out_igv_own_data);
        f_output_root.mkdirs();

        HashSet<Double> already_done_percentage = new HashSet<>();

        for (int i = 0; i < tfs_in_order.size(); i++) {
            String tf = tfs_in_order.get(i);
            if (which_tfs.contains(tf)) {
                int rank = i + 1;
                File f_output_tf = new File(f_output_root.getAbsolutePath() + File.separator + rank + "_" + tf);
                f_output_tf.mkdir();
                File f_out_session = new File(
                        f_output_tf.getAbsolutePath() + File.separator + options_intern.folder_out_igv_session);
                f_out_session.mkdir();

                for (String key_hm : hm_timepoint_tf_target_genes.keySet()) {
                    HashMap<String, HashMap<String, ArrayList<String>>> hm = hm_timepoint_tf_target_genes.get(key_hm);

                    for (String key_tp : hm.keySet()) {
                        File f_out_session_tp = new File(f_out_session.getAbsolutePath() + File.separator + key_tp);
                        f_out_session_tp.mkdir();

                        /*
                        String command = "sh ";
                        if(options_intern.igv_path_to_igv.endsWith("sh"))
                        {
                            command+=options_intern.igv_path_to_igv;
                        }
                        else
                        {
                            command+=options_intern.igv_path_to_igv+"igv.sh";
                        }

                        logger.logLine("[IGV] Open IGV to start local server ...");
                        Process child = Runtime.getRuntime().exec("/bin/bash -c "+command);*/
                        /*int code = child.waitFor();
                        switch (code){
                            case 0:
                                break;
                            case 1:
                                String message = child.getErrorStream().toString();
                                throw new Exception(message);
                        }*/
                        HashMap<String, ArrayList<String>> hm_tp = hm.get(key_tp);

                        if (hm_tp.containsKey(tf)) {

                            String load_tf_chip_seq = "load ";
                            ArrayList<String> tf_chip_seq_files = new ArrayList<>();

                            boolean added_prediction_data = false;

                            //include predictive HMs
                            if(options_intern.igv_include_prediction_data.size()>0)
                            {
                                File f_input =
                                        new File(options_intern.tepic_input_directory+File.separator+key_tp);

                                for(String including_hm : options_intern.igv_include_prediction_data)
                                {
                                    File f_hm_input = new File(f_input.getAbsolutePath()+File.separator+including_hm);
                                    if(f_hm_input.exists())
                                    {
                                        for(File f_file : f_hm_input.listFiles())
                                        {
                                            if(added_prediction_data)
                                            {
                                                load_tf_chip_seq+="," + f_file.getAbsolutePath();
                                            }
                                            else
                                            {
                                                added_prediction_data=true;
                                                load_tf_chip_seq+= f_file.getAbsolutePath();
                                            }
                                        }
                                    }
                                }
                            }


                            File f_input_evaluation_lines_root = new File(options_intern.igv_path_to_tf_chip_seq);
                            for (File f_evaluation_tp : f_input_evaluation_lines_root.listFiles()) {
                                if (f_evaluation_tp.getName().equals(key_tp)) {
                                    for (File f_evaluation_tp_tf : f_evaluation_tp.listFiles()) {
                                        if (f_evaluation_tp_tf.isDirectory()) {
                                            for (File f_evaluation_tp_tf_file : f_evaluation_tp_tf.listFiles()) {
                                                if (f_evaluation_tp_tf_file.isFile()) {
                                                    tf_chip_seq_files.add(f_evaluation_tp_tf_file.getAbsolutePath());
                                                }
                                            }

                                        }
                                    }
                                }
                            }



                            for (int j = 0; j < tf_chip_seq_files.size(); j++) {
                                if(tf_chip_seq_files.get(j).endsWith(".idx"))
                                {
                                    continue;
                                }
                                if (!added_prediction_data) {
                                    load_tf_chip_seq += tf_chip_seq_files.get(j);
                                    added_prediction_data=true;
                                } else {
                                    load_tf_chip_seq += "," + tf_chip_seq_files.get(j);
                                }
                            }

                            //check if ChIP-ATLAS data is available
                            File f_chip_atlas_data =
                                    new File(options_intern.com2pose_working_directory + File.separator+options_intern.folder_out_chip_atlas+ File.separator + options_intern.folder_out_chip_atlas_peak_files+ File.separator + rank + "_" + tf);
                            if(f_chip_atlas_data.exists())
                            {
                                for(File f_input : f_chip_atlas_data.listFiles())
                                {
                                    if(f_input.getName().endsWith(".idx"))
                                    {
                                        continue;
                                    }

                                    if(added_prediction_data)
                                    {
                                        load_tf_chip_seq+="," + f_input.getAbsolutePath();
                                    }
                                    else
                                    {
                                        added_prediction_data=true;
                                        load_tf_chip_seq+= f_input.getAbsolutePath();
                                    }
                                }
                            }

                            //save sessions
                            File f_session_tp_name = new File(f_out_session_tp.getAbsolutePath() + File.separator +
                                    options_intern.file_suffix_session);
                            igv_save_sessions(f_session_tp_name, load_tf_chip_seq);

                            /*out.println("genome "+options_intern.igv_species_ref_genome);
                            String response = in.readLine();
                            if(!response.equals("OK"))
                            {
                                logger.logLine("[IGV] igv_species_ref_genome not OK!");
                            }*/

                            File f_output_tf_hm = new File(f_output_tf.getAbsolutePath() + File.separator + key_hm);
                            f_output_tf_hm.mkdir();

                            File f_output_tf_hm_tp =
                                    new File(f_output_tf_hm.getAbsolutePath() + File.separator + key_tp);
                            f_output_tf_hm_tp.mkdir();


                            ArrayList<String> target_genes_for_photo_session = hm_tp.get(tf);

                            //read in regions to target gene table
                            HashMap<String, HashSet<String>> gene_to_regions_all = new HashMap<>();
                            File f_input_regions_to_target_genes = new File(
                                    options_intern.com2pose_working_directory + File.separator +
                                            options_intern.folder_name_tepic_output_raw + File.separator + key_tp +
                                            File.separator + key_hm);
                            File f_input_regions_to_target_genes_search = new File("");
                            if (!f_input_regions_to_target_genes.exists()) {
                                if (!key_hm.equals(options_intern.distribution_analysis_all_name)) {
                                    continue;
                                } else {
                                    //could be unnecessary
                                    for (String key_hm_to_important : hm_to_important_tfs.keySet()) {
                                        hm_to_important_tfs.get(key_hm_to_important).add(tf);
                                    }
                                    continue;
                                }
                            }

                            for (File f_search : f_input_regions_to_target_genes.listFiles()) {
                                f_input_regions_to_target_genes_search = f_search;
                            }
                            f_input_regions_to_target_genes_search = new File(
                                    f_input_regions_to_target_genes_search.getAbsolutePath() + File.separator +
                                            options_intern.file_suffix_tepic_output_regions_to_target_genes);

                            BufferedReader br_regions_to_target_genes =
                                    new BufferedReader(new FileReader(f_input_regions_to_target_genes_search));
                            String line_regions_to_target_genes = br_regions_to_target_genes.readLine();
                            String[] header_regions_to_target_genes = line_regions_to_target_genes.split("\t");
                            while ((line_regions_to_target_genes = br_regions_to_target_genes.readLine()) != null) {
                                String[] split = line_regions_to_target_genes.split("\t");
                                String region = split[0];
                                String[] split_target_genes = split[1].split(";");

                                for (String key_target_gene : split_target_genes) {
                                    key_target_gene = key_target_gene.split("\\.")[0];
                                    key_target_gene = ensg_gene_symbol_map.get(key_target_gene);

                                    if (key_target_gene == null) {
                                        continue;
                                    }

                                    if (gene_to_regions_all.containsKey(key_target_gene)) {
                                        gene_to_regions_all.get(key_target_gene).add(region);
                                    } else {
                                        HashSet<String> x = new HashSet<>();
                                        x.add(region);
                                        gene_to_regions_all.put(key_target_gene, x);
                                    }
                                }
                            }
                            br_regions_to_target_genes.close();

                            //create wides range for each possible target gene
                            HashMap<String, String> gene_to_region = new HashMap<>();
                            for (String key_gene : gene_to_regions_all.keySet()) {
                                HashSet<String> regions_into_one = gene_to_regions_all.get(key_gene);

                                String chr = "";
                                int begin = Integer.MAX_VALUE;
                                int end = Integer.MIN_VALUE;

                                for (String position : regions_into_one) {
                                    String[] split = position.split(":");

                                    if (chr.equals("")) {
                                        chr = split[0];
                                    } else {
                                        if (!chr.equals(split[0])) {
                                            //logger.logLine("[ERROR] DIFFERENT CHROMOSOME FOR LOOKUP?!");
                                            continue;
                                        }
                                    }

                                    String[] split_pos = split[1].split("-");
                                    int begin_n = Integer.parseInt(split_pos[0]);
                                    int end_n = Integer.parseInt(split_pos[1]);

                                    if (begin_n < begin) {
                                        begin = begin_n;
                                    }
                                    if (end < end_n) {
                                        end = end_n;
                                    }

                                }

                                begin += 1000;
                                end += 1000;

                                String position = "chr" + chr + ":" + begin + "-" + end;
                                gene_to_region.put(key_gene, position);
                            }

                            //create new gene ranges
                            HashMap<String, String> gene_to_coordinates_this_run = new HashMap<>(gene_to_coordinates);
                            for (String key_gene_to_region : gene_to_region.keySet()) {
                                if (gene_to_coordinates_this_run.containsKey(key_gene_to_region)) {
                                    String position = gene_to_coordinates_this_run.get(key_gene_to_region);
                                    String region_tf = gene_to_region.get(key_gene_to_region);

                                    String[] split_p = position.split(":");
                                    String[] split_p_position = split_p[1].split("-");

                                    String[] split_r = region_tf.split(":");
                                    String[] split_r_position = split_r[1].split("-");

                                    String chr_p = split_p[0];
                                    String chr_r = split_r[0];


                                    int begin_p = Integer.parseInt(split_p_position[0]);
                                    int end_p = Integer.parseInt(split_p_position[1]);

                                    int begin_r = Integer.parseInt(split_r_position[0]);
                                    int end_r = Integer.parseInt(split_r_position[1]);

                                    if (!chr_p.equals(chr_r)) {
                                        //logger.logLine("[ERROR]: chromosome do not match");
                                        begin_p -= 100000;
                                        end_p += 100000;

                                        String final_string = chr_p + ":" + begin_p + "-" + end_p;

                                        gene_to_coordinates_this_run.put(key_gene_to_region.toUpperCase(),
                                                final_string);
                                        continue;
                                    }


                                    int begin_final = Integer.MAX_VALUE;
                                    int end_final = Integer.MIN_VALUE;

                                    if (begin_r < begin_p) {
                                        begin_final = begin_r;
                                    } else {
                                        begin_final = begin_p;
                                    }
                                    if (end_r < end_p) {
                                        end_final = end_p;
                                    } else {
                                        end_final = end_r;
                                    }

                                    begin_final -= 100000;
                                    end_final += 100000;

                                    String final_string = chr_p + ":" + begin_final + "-" + end_final;

                                    gene_to_coordinates_this_run.put(key_gene_to_region.toUpperCase(), final_string);

                                }
                            }

                            Socket socket = new Socket("127.0.0.1", options_intern.igv_port_number);
                            PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
                            BufferedReader in = new BufferedReader(new InputStreamReader(socket.getInputStream()));


                            //logger.logLine("[IGV] " + load_tf_chip_seq);
                            out.println(load_tf_chip_seq);
                            String response_load = in.readLine();
                            //logger.logLine("[IGV] " + response_load);

                            //logger.logLine("[IGV] "+ "genome "+options_intern.igv_species_ref_genome);
                            out.println("genome " + options_intern.igv_species_ref_genome);
                            String response = in.readLine();
                            //logger.logLine("[IGV] "+ response);

                            int count = 1;
                            for (String targets : target_genes_for_photo_session) {
                                String snapshot_name = count + "_" + tf + "_" + targets + ".png";

                                File f_output_shot = new File(f_output_tf_hm_tp.getAbsolutePath());
                                f_output_shot.mkdir();
                                //logger.logLine("[IGV] "+ "snapshotDirectory "+f_output_shot.getAbsolutePath());
                                out.println("snapshotDirectory " + f_output_shot.getAbsolutePath());
                                response = in.readLine();
                                //logger.logLine("[IGV] " + response);
                                //logger.logLine("[IGV] "+ "goto " + gene_to_coordinates.get(targets));
                                out.println("goto " + gene_to_coordinates_this_run.get(targets));
                                response = in.readLine();
                                //logger.logLine("[IGV] " + response);
                                //logger.logLine("[IGV] "+ "sort");
                                //out.println("sort");
                                //response=in.readLine();
                                //logger.logLine("[IGV] " + response);
                                //logger.logLine("[IGV] "+ "collapse");
                                //out.println("collapse");
                                //response=in.readLine();
                                //logger.logLine("[IGV] " + response);
                                //logger.logLine("[IGV] "+ "snapshot");
                                out.println("snapshot");
                                response = in.readLine();
                                //logger.logLine("[IGV] " + response);

                                File f_snapshot_made = new File("NO");

                                //search png name
                                for (File f_snapshot : f_output_shot.listFiles()) {
                                    if (f_snapshot.isFile()) {
                                        if (f_snapshot.getName().startsWith("chr")) {
                                            f_snapshot_made = f_snapshot;
                                            break;
                                        }

                                    }
                                }

                                if (!f_snapshot_made.exists()) {
                                    //its the all overview
                                    //logger.logLine("[ERROR]: trying to move file which does not exist!");
                                } else {

                                    String command_mv = "mv " + f_snapshot_made.getAbsolutePath() + " " +
                                            f_snapshot_made.getParentFile().getAbsolutePath() + File.separator +
                                            snapshot_name;


                                    Process child = Runtime.getRuntime().exec(command_mv);
                                    int code = child.waitFor();
                                    switch (code) {
                                        case 0:
                                            break;
                                        case 1:
                                            String message = child.getErrorStream().toString();
                                            throw new Exception(message);
                                    }


                                }


                                count++;
                            }

                            out.println("new");
                            response = in.readLine();
                            //logger.logLine("[IGV] " + response);

                            socket.close();

                            //new session -> clear tracks
                            /*logger.logLine("[IGV] "+ "new");
                            out.println("new");
                            String response = in.readLine();
                            logger.logLine("[IGV] " + response);*/
                        }
                    }
                }
            }

            double percentage_done = ((i * 1.0) / tfs_in_order.size());
            DecimalFormat df = new DecimalFormat("0.0");
            String percentage_string = df.format(percentage_done);
            percentage_done = Double.parseDouble(percentage_string) * 100;

            if (percentage_done % 10 == 0 && !already_done_percentage.contains(percentage_done)) {
                logger.logLine("[IGV] percentage of tasks done: " + percentage_done + "%");
                already_done_percentage.add(percentage_done);
            }
        }


        logger.logLine("[IGV] finished taking screenshots for top TFs and their corresponding top target genes");
    }

    /**
     * run IGV and take screenshots chr wide
     */
    public void run_igv_chr_wide_data() throws IOException {
        logger.logLine("[IGV] Take screenshots of peaks chr wide and genome wide");

        if (options_intern.mix_level.equals("SAMPLE_LEVEL")) {
            File root_mix_working_dir = new File(
                    options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option);
            File f_sample_mix_output = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                    options_intern.folder_name_mix_option_sample_mix);
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = f_sample_mix_output.getAbsolutePath();
        }

        if (options_intern.mix_level.equals("HM_LEVEL")) {
            File root_mix_working_dir = new File(
                    options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option);
            File f_output_hm = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                    options_intern.folder_name_mix_option_hm_mix);
            f_output_hm.mkdir();

            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = f_output_hm.getAbsolutePath();
        }
        if (!options_intern.black_list_dir.equals("")) {
            File output_folder = new File(options_intern.com2pose_working_directory + File.separator +
                    options_intern.folder_name_blacklisted_regions);
            File output_folder_new_input = new File(output_folder.getAbsolutePath() + File.separator +
                    options_intern.folder_name_blacklisted_regions_new_input);
            output_folder_new_input.mkdir();

            //set new folder directory for tepic input and save old one
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = output_folder_new_input.getAbsolutePath();
        }
        if (options_intern.mix_mutually_exclusive) {
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory =
                    options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option +
                            File.separator + options_intern.folder_name_mix_option_mutually_exclusive + File.separator +
                            options_intern.folder_name_mix_options_mutually_exclusive_input;
        }

        HashSet<String> chromosomes = new HashSet<>();

        File f_input_userInput_bedfiles = new File(File.separator + options_intern.tepic_input_directory);
        File f_input_chipATLAS_data = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_chip_atlas +
                        File.separator + options_intern.folder_out_chip_atlas_peak_files);

        File f_output_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_igv +
                        File.separator + options_intern.folder_out_igv_chipAtlas_chrWide_genomeWide_views);
        f_output_root.mkdirs();

        Socket socket = new Socket("127.0.0.1", options_intern.igv_port_number);
        PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
        BufferedReader in = new BufferedReader(new InputStreamReader(socket.getInputStream()));

        out.println("new");
        String response = in.readLine();

        String load_chipatlas_overall = "load ";

        //load ChIP-seq CHIPATLAS data
        for (File f_dir_atlas : f_input_chipATLAS_data.listFiles()) {
            if (f_dir_atlas.isDirectory()) {
                for (File f_dir_bed : f_dir_atlas.listFiles()) {
                    if (f_dir_bed.isFile()) {
                        if (f_dir_bed.getName().contains("idx")) {
                            continue;
                        }

                        out.println("load " + f_dir_bed.getAbsolutePath());
                        String response_load = in.readLine();

                        load_chipatlas_overall += " " + f_dir_bed.getAbsolutePath();
                    }
                }
            }
        }

        //load input ChIP-seq data
        for (File f_tp : f_input_userInput_bedfiles.listFiles()) {
            if (f_tp.isDirectory()) {
                for (File f_hm : f_tp.listFiles()) {
                    if (f_hm.isDirectory()) {
                        for (File f_bed : f_hm.listFiles()) {
                            if (f_bed.isFile()) {
                                if (f_bed.getName().contains("idx")) {
                                    continue;
                                }
                                ;

                                out.println("load " + f_bed.getAbsolutePath());
                                String response_load = in.readLine();

                                load_chipatlas_overall += " " + f_bed.getAbsolutePath();
                            }
                        }
                    }
                }
            }

        }

        File f_fasta_for_chromosomes = new File(options_intern.tepic_input_ref_genome);
        BufferedReader br_fasta = new BufferedReader(new FileReader(f_fasta_for_chromosomes));
        String line_fasta = "";
        while ((line_fasta = br_fasta.readLine()) != null) {
            if (line_fasta.startsWith(">")) {
                String[] split = line_fasta.split(">");
                if (split[1].length() <= 3) {
                    chromosomes.add(split[1]);
                }
            }
        }

        //take screenshots of all chromosomes and genome wide

        for (String key_chr : chromosomes) {
            if (!key_chr.contains("chr")) {
                key_chr = "chr" + key_chr;
            }

            out.println("genome " + options_intern.igv_species_ref_genome);
            String response_genome = in.readLine();

            File f_output_shot = new File(f_output_root.getAbsolutePath());
            //logger.logLine("[IGV] "+ "snapshotDirectory "+f_output_shot.getAbsolutePath());
            out.println("snapshotDirectory " + f_output_shot.getAbsolutePath());
            String response_next = in.readLine();
            //logger.logLine("[IGV] " + response);
            //logger.logLine("[IGV] "+ "goto " + gene_to_coordinates.get(targets));
            out.println("goto " + key_chr);
            response_next = in.readLine();
            //logger.logLine("[IGV] " + response);
            //logger.logLine("[IGV] "+ "sort");
            //out.println("sort");
            //response=in.readLine();
            //logger.logLine("[IGV] " + response);
            //logger.logLine("[IGV] "+ "collapse");
            out.println("collapse");
            response_next = in.readLine();
            //logger.logLine("[IGV] " + response);
            //logger.logLine("[IGV] "+ "snapshot");
            out.println("snapshot");
            response_next = in.readLine();
        }

        //take screenshot of genome view

        out.println("genome " + options_intern.igv_species_ref_genome);
        String response_genome = in.readLine();

        File f_output_shot = new File(f_output_root.getAbsolutePath());
        //logger.logLine("[IGV] "+ "snapshotDirectory "+f_output_shot.getAbsolutePath());
        out.println("snapshotDirectory " + f_output_shot.getAbsolutePath());
        String response_next = in.readLine();
        //logger.logLine("[IGV] " + response);
        //logger.logLine("[IGV] "+ "goto " + gene_to_coordinates.get(targets));
        out.println("goto genome");
        response_next = in.readLine();
        //logger.logLine("[IGV] " + response);
        //logger.logLine("[IGV] "+ "sort");
        //out.println("sort");
        //response=in.readLine();
        //logger.logLine("[IGV] " + response);
        //logger.logLine("[IGV] "+ "collapse");
        out.println("collapse");
        response_next = in.readLine();
        //logger.logLine("[IGV] " + response);
        //logger.logLine("[IGV] "+ "snapshot");
        out.println("snapshot");
        response_next = in.readLine();


        out.println("new");
        String response_end = in.readLine();
        //logger.logLine("[IGV] " + response);
        socket.close();

        //save session
        File f_output_session =
                new File(f_output_shot.getAbsolutePath() + File.separator + options_intern.file_suffix_session);
        igv_save_sessions(f_output_session, load_chipatlas_overall);

        logger.logLine("[IGV] Finished taking screenshots of peaks chr wide");

    }

    /**
     * run igv on chip atlas data
     */
    public void run_igv_chip_atlas_data() throws Exception {
        logger.logLine("[ChIP-ATLAS-IGV] Start running ChIP Atlas evaluation data. (ChIP-ATLAS data)");
        //logger.logLine("[ChIP-ATLAS-IGV] This is as far automated as possible. If a .bed file is very big a user entry is prompted at IGV, please click GO there!");

        check_tepic_input_with_options();

        File f_output_igv_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_igv +
                        File.separator + options_intern.folder_out_igv_chip_atlas_data);
        f_output_igv_root.mkdirs();

        File f_bed_data_input_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_chip_atlas +
                        File.separator + options_intern.folder_out_chip_atlas_peak_files);
        File f_target_genes_input_root = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_out_target_genes_dcg);
        File f_important_tfs_per_hm = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution +
                        File.separator + options_intern.folder_out_distribution_stats + File.separator +
                        options_intern.folder_out_distribution_stats_HM);

        //retrieve ordered tfs
        File f_input_dcg_result = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution +
                        File.separator + options_intern.folder_out_distribution_dcg + File.separator +
                        options_intern.file_suffix_distribution_analysis_dcg);
        ArrayList<String> ordered_tfs_dcg = new ArrayList<>();

        HashMap<String, String> gene_to_coordinates = new HashMap<>();

        File f_input_gene_to_coordinates = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.folder_name_deseq2_preprocessing_gene_positions + File.separator +
                options_intern.file_suffix_deseq2_preprocessing_gene_positions_data);
        BufferedReader br_gene_coordinates = new BufferedReader(new FileReader(f_input_gene_to_coordinates));
        String line_gene_coordinates = br_gene_coordinates.readLine();
        String[] header_gene_coordinates = line_gene_coordinates.split("\t");
        while ((line_gene_coordinates = br_gene_coordinates.readLine()) != null) {
            String[] split = line_gene_coordinates.split("\t");

            if (split[1].startsWith("CHR")) {
                continue;
            }

            String gene_name = split[0];
            String chr = "chr" + split[1];
            String begin = split[2];
            String end = split[3];

            String position_string = chr + ":" + begin + "-" + end;
            gene_to_coordinates.put(gene_name.toUpperCase(), position_string);
        }

        /*

        //read gtf for genome coordinates
        //with this kind of thing it got messed up... use the biomart version!
        BufferedReader br_gene_coordinates = new BufferedReader(new FileReader(new File(options_intern.tepic_gene_annot)));
        String line_gene_coordinates = br_gene_coordinates.readLine();
        while((line_gene_coordinates=br_gene_coordinates.readLine())!=null)
        {
            if(line_gene_coordinates.startsWith("#"))
            {
                continue;
            }

            String[] split = line_gene_coordinates.split("\t");

            if(!split[2].equals("gene"))
            {
                continue;
            }

            String chr = split[0];
            String position = split[3]+"-"+split[4];
            String gene_name = "";

            String[] gene_name_prep = split[8].split(";");
            for(String key : gene_name_prep)
            {
                key=key.trim();
                if(key.startsWith("gene_name"))
                {
                    String[] split_intern = key.split("\"");
                    gene_name=split_intern[1].toUpperCase();
                }
            }

            gene_to_coordinates.put(gene_name,chr+":"+position);
        }*/

        HashMap<String, String> ensg_gene_symbol_map = new HashMap<>();
        BufferedReader br_ensg_gene_symbol =
                new BufferedReader(new FileReader(new File(options_intern.tepic_ensg_symbol)));
        String line_ensg_gene_symbol = br_ensg_gene_symbol.readLine();
        while ((line_ensg_gene_symbol = br_ensg_gene_symbol.readLine()) != null) {
            String[] split = line_ensg_gene_symbol.split("\t");
            if (split.length > 1) {
                ensg_gene_symbol_map.put(split[0].toUpperCase(), split[1].toUpperCase());
            }
        }
        br_ensg_gene_symbol.close();


        BufferedReader br_dcg = new BufferedReader(new FileReader(f_input_dcg_result));
        String line_dcg = br_dcg.readLine();
        while ((line_dcg = br_dcg.readLine()) != null) {
            String[] split = line_dcg.split("\t");
            ordered_tfs_dcg.add(split[1]);
        }
        br_dcg.close();

        //read in all important tfs per hm
        HashMap<String, HashSet<String>> hm_to_important_tfs = new HashMap<>();
        for (File fileDir_hm : f_important_tfs_per_hm.listFiles()) {
            if (fileDir_hm.isDirectory()) {
                HashSet<String> important_tfs = new HashSet<>();
                hm_to_important_tfs.put(fileDir_hm.getName(), important_tfs);

                BufferedReader br = new BufferedReader(new FileReader(new File(
                        fileDir_hm.getAbsolutePath() + File.separator +
                                options_intern.file_suffix_distribution_analysis_plot_stats)));
                String line = br.readLine();
                line = br.readLine();
                while ((line = br.readLine()) != null) {
                    String[] split = line.split("\t");
                    important_tfs.add(split[1]);
                }
            }
        }

        //read in TFs from ALL analysis
        File f_input_all_hm_to_important_tfs = new File(
                f_important_tfs_per_hm.getParentFile().getAbsolutePath() + File.separator +
                        options_intern.folder_out_distribution_stats_ALL + File.separator +
                        options_intern.file_suffix_distribution_analysis_plot_stats);
        BufferedReader br_input_all_hm_to_important_tfs =
                new BufferedReader(new FileReader(f_input_all_hm_to_important_tfs));
        String line_input_all_hm_to_important_tfs = br_input_all_hm_to_important_tfs.readLine();
        HashSet<String> hs_input_all_hm_to_important_tfs = new HashSet<>();
        while ((line_input_all_hm_to_important_tfs = br_input_all_hm_to_important_tfs.readLine()) != null) {
            String[] split = line_input_all_hm_to_important_tfs.split("\t");
            hs_input_all_hm_to_important_tfs.add(split[1]);

            //retrieve all TF targets from ALL analysis
            for (String key_hm : hm_to_important_tfs.keySet()) {
                hm_to_important_tfs.get(key_hm).add(split[1]);

            }

        }
        br_input_all_hm_to_important_tfs.close();
        hm_to_important_tfs.put(options_intern.distribution_analysis_all_name, hs_input_all_hm_to_important_tfs);


        //go over all dcg ordered tfs
        int i = 1;
        HashSet<Double> already_done_percentage = new HashSet<>();

        for (String key_tf : ordered_tfs_dcg) {
            File f_input_bed_test =
                    new File(f_bed_data_input_root.getAbsolutePath() + File.separator + i + "_" + key_tf);
            if (!f_input_bed_test.exists()) {
                i++;
                continue;
            }

            File f_output_igv_shots = new File(f_output_igv_root.getAbsolutePath() + File.separator + i + "_" + key_tf);
            f_output_igv_shots.mkdirs();

            for (File file_dir_tp : f_target_genes_input_root.listFiles()) {
                if (file_dir_tp.isDirectory()) {
                    File f_output_igv_shots_tp =
                            new File(f_output_igv_shots.getAbsolutePath() + File.separator + file_dir_tp.getName());
                    f_output_igv_shots_tp.mkdirs();

                    for (String key_hm : hm_to_important_tfs.keySet()) {
                        if (hm_to_important_tfs.get(key_hm).contains(key_tf)) {

                            //read in regions to target gene table
                            HashMap<String, HashSet<String>> gene_to_regions_all = new HashMap<>();
                            File f_input_regions_to_target_genes = new File(
                                    options_intern.com2pose_working_directory + File.separator +
                                            options_intern.folder_name_tepic_output_raw + File.separator +
                                            file_dir_tp.getName() + File.separator + key_hm);
                            File f_input_regions_to_target_genes_search = new File("");
                            if (!f_input_regions_to_target_genes.exists()) {
                                if (!key_hm.equals(options_intern.distribution_analysis_all_name)) {
                                    continue;
                                } else {
                                    //could be unnecessary
                                    for (String key_hm_to_important : hm_to_important_tfs.keySet()) {
                                        hm_to_important_tfs.get(key_hm_to_important).add(key_tf);
                                    }
                                    continue;
                                }
                            }

                            for (File f_search : f_input_regions_to_target_genes.listFiles()) {
                                f_input_regions_to_target_genes_search = f_search;
                            }
                            f_input_regions_to_target_genes_search = new File(
                                    f_input_regions_to_target_genes_search.getAbsolutePath() + File.separator +
                                            options_intern.file_suffix_tepic_output_regions_to_target_genes);

                            BufferedReader br_regions_to_target_genes =
                                    new BufferedReader(new FileReader(f_input_regions_to_target_genes_search));
                            String line_regions_to_target_genes = br_regions_to_target_genes.readLine();
                            String[] header_regions_to_target_genes = line_regions_to_target_genes.split("\t");
                            while ((line_regions_to_target_genes = br_regions_to_target_genes.readLine()) != null) {
                                String[] split = line_regions_to_target_genes.split("\t");
                                String region = split[0];
                                String[] split_target_genes = split[1].split(";");

                                for (String key_target_gene : split_target_genes) {
                                    key_target_gene = key_target_gene.split("\\.")[0];
                                    key_target_gene = ensg_gene_symbol_map.get(key_target_gene);

                                    if (key_target_gene == null) {
                                        continue;
                                    }

                                    if (gene_to_regions_all.containsKey(key_target_gene)) {
                                        gene_to_regions_all.get(key_target_gene).add(region);
                                    } else {
                                        HashSet<String> x = new HashSet<>();
                                        x.add(region);
                                        gene_to_regions_all.put(key_target_gene, x);
                                    }
                                }
                            }
                            br_regions_to_target_genes.close();

                            //create wides range for each possible target gene
                            HashMap<String, String> gene_to_region = new HashMap<>();
                            for (String key_gene : gene_to_regions_all.keySet()) {
                                HashSet<String> regions_into_one = gene_to_regions_all.get(key_gene);

                                String chr = "";
                                int begin = Integer.MAX_VALUE;
                                int end = Integer.MIN_VALUE;

                                for (String position : regions_into_one) {
                                    String[] split = position.split(":");

                                    if (chr.equals("")) {
                                        chr = split[0];
                                    } else {
                                        if (!chr.equals(split[0])) {
                                            //logger.logLine("[ERROR] DIFFERENT CHROMOSOME FOR LOOKUP?!");
                                            continue;
                                        }
                                    }

                                    String[] split_pos = split[1].split("-");
                                    int begin_n = Integer.parseInt(split_pos[0]);
                                    int end_n = Integer.parseInt(split_pos[1]);

                                    if (begin_n < begin) {
                                        begin = begin_n;
                                    }
                                    if (end < end_n) {
                                        end = end_n;
                                    }

                                }

                                begin += 1000;
                                end += 1000;

                                String position = "chr" + chr + ":" + begin + "-" + end;
                                gene_to_region.put(key_gene, position);
                            }

                            //create new gene ranges
                            HashMap<String, String> gene_to_coordinates_this_run = new HashMap<>(gene_to_coordinates);
                            for (String key_gene_to_region : gene_to_region.keySet()) {
                                if (gene_to_coordinates_this_run.containsKey(key_gene_to_region)) {
                                    String position = gene_to_coordinates_this_run.get(key_gene_to_region);
                                    String region_tf = gene_to_region.get(key_gene_to_region);

                                    String[] split_p = position.split(":");
                                    String[] split_p_position = split_p[1].split("-");

                                    String[] split_r = region_tf.split(":");
                                    String[] split_r_position = split_r[1].split("-");

                                    String chr_p = split_p[0];
                                    String chr_r = split_r[0];


                                    int begin_p = Integer.parseInt(split_p_position[0]);
                                    int end_p = Integer.parseInt(split_p_position[1]);

                                    int begin_r = Integer.parseInt(split_r_position[0]);
                                    int end_r = Integer.parseInt(split_r_position[1]);

                                    if (!chr_p.equals(chr_r)) {
                                        //logger.logLine("[ERROR]: chromosome do not match");
                                        begin_p -= 100000;
                                        end_p += 100000;

                                        String final_string = chr_p + ":" + begin_p + "-" + end_p;

                                        gene_to_coordinates_this_run.put(key_gene_to_region.toUpperCase(),
                                                final_string);
                                        continue;
                                    }


                                    int begin_final = Integer.MAX_VALUE;
                                    int end_final = Integer.MIN_VALUE;

                                    if (begin_r < begin_p) {
                                        begin_final = begin_r;
                                    } else {
                                        begin_final = begin_p;
                                    }
                                    if (end_r < end_p) {
                                        end_final = end_p;
                                    } else {
                                        end_final = end_r;
                                    }

                                    begin_final -= 100000;
                                    end_final += 100000;

                                    String final_string = chr_p + ":" + begin_final + "-" + end_final;

                                    gene_to_coordinates_this_run.put(key_gene_to_region.toUpperCase(), final_string);

                                }
                            }

                            File f_output_igv_shots_tp_hm =
                                    new File(f_output_igv_shots_tp.getAbsolutePath() + File.separator + key_hm);
                            f_output_igv_shots_tp_hm.mkdirs();

                            //read in all target genes
                            ArrayList<String> target_genes = new ArrayList<>();

                            File f_input_target_genes_file = new File(
                                    file_dir_tp.getAbsolutePath() + File.separator + key_hm + File.separator + key_tf +
                                            File.separator + key_tf + ".csv");

                            if (!f_input_target_genes_file.exists()) {
                                continue;
                            }

                            BufferedReader br = new BufferedReader(new FileReader(f_input_target_genes_file));
                            String line = br.readLine();
                            int count_target_genes = 0;
                            while ((line = br.readLine()) != null) {
                                String[] split = line.split("\t");
                                target_genes.add(split[0]);

                                if (count_target_genes == options_intern.plot_top_k_genes - 1) {
                                    break;
                                }
                                count_target_genes++;
                            }
                            br.close();

                            //shoot igv now
                            int count = 0;

                            File f_input_bed = new File(
                                    f_bed_data_input_root.getAbsolutePath() + File.separator + i + "_" + key_tf);
                            if (!f_input_bed.exists()) {
                                continue;
                            }

                            String load_for_session_intern = "load ";
                            String load_tf_chip_seq_1 ="load ";
                            boolean added_prediction_data=false;

                            //include predictive HMs
                            if(options_intern.igv_include_prediction_data.size()>0)
                            {
                                File f_input =
                                        new File(options_intern.tepic_input_directory+File.separator+file_dir_tp.getName());

                                for(String including_hm : options_intern.igv_include_prediction_data)
                                {
                                    File f_hm_input = new File(f_input.getAbsolutePath()+File.separator+including_hm);
                                    if(f_hm_input.exists())
                                    {
                                        for(File f_file : f_hm_input.listFiles())
                                        {
                                            if(added_prediction_data)
                                            {
                                                load_tf_chip_seq_1+="," + f_file.getAbsolutePath();
                                                load_for_session_intern+="," + f_file.getAbsolutePath();
                                            }
                                            else
                                            {
                                                added_prediction_data=true;
                                                load_tf_chip_seq_1+= f_file.getAbsolutePath();
                                                load_for_session_intern+=f_file.getAbsolutePath();
                                            }
                                        }
                                    }
                                }
                            }

                            //out.println(load_tf_chip_seq_1);
                            //String response_load_1 = in.readLine();

                            String load_bed_files_overall = "load";

                            for (File f_bed_files : f_input_bed.listFiles()) {
                                if (f_bed_files.getName().endsWith("idx")) {
                                    continue;
                                }
                                String load_tf_chip_seq = "load " + f_bed_files.getAbsolutePath();

                                //out.println(load_tf_chip_seq);
                                //String response_load = in.readLine();

                                load_bed_files_overall += " " + f_bed_files.getAbsolutePath();

                                if(added_prediction_data)
                                {
                                    load_for_session_intern+="," + f_bed_files.getAbsolutePath();
                                }
                                else
                                {
                                    added_prediction_data=true;
                                    load_for_session_intern+=f_bed_files.getAbsolutePath();
                                }

                            }

                            //save session
                            File f_session_intern =
                                    new File(f_output_igv_shots.getAbsolutePath() + File.separator +file_dir_tp.getName()+File.separator+key_hm+File.separator+
                                    options_intern.file_suffix_session);
                            igv_save_sessions(f_session_intern, load_for_session_intern);

                            Socket socket = new Socket("127.0.0.1", options_intern.igv_port_number);
                            PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
                            BufferedReader in = new BufferedReader(new InputStreamReader(socket.getInputStream()));

                            //load again
                            out.println(load_for_session_intern);
                            String response_load_2 = in.readLine();


                            for (String key_target_gene : target_genes) {
                                String snapshot_name = count + "_" + key_tf + "_" + key_target_gene + ".png";

                                File f_output_tf_hm_tp = new File(
                                        f_output_igv_shots_tp_hm.getAbsolutePath() + File.separator + count + "_" +
                                                key_target_gene);

                                //logger.logLine("[IGV] " + load_tf_chip_seq);
                                //logger.logLine("[IGV] " + response_load);
                                //logger.logLine("[IGV] "+ "genome "+options_intern.igv_species_ref_genome);
                                //logger.logLine("[IGV] "+ response);
                                out.println("genome " + options_intern.igv_species_ref_genome);
                                String response_genome = in.readLine();

                                //File f_output_shot= new File(f_output_tf_hm_tp.getAbsolutePath());
                                //f_output_shot.mkdir();

                                File f_output_shot = new File(f_output_igv_shots_tp_hm.getAbsolutePath());
                                f_output_shot.mkdir();


                                //logger.logLine("[IGV] "+ "snapshotDirectory "+f_output_shot.getAbsolutePath());
                                out.println("snapshotDirectory " + f_output_shot.getAbsolutePath());
                                String response = in.readLine();
                                //logger.logLine("[IGV] " + response);
                                //logger.logLine("[IGV] "+ "goto " + gene_to_coordinates.get(targets));
                                //out.println("goto " + gene_to_coordinates.get(key_target_gene));
                                out.println("goto " + gene_to_coordinates_this_run.get(key_target_gene));
                                response = in.readLine();
                                //logger.logLine("[IGV] " + response);
                                //logger.logLine("[IGV] "+ "sort");
                                //out.println("sort");
                                //response=in.readLine();
                                //logger.logLine("[IGV] " + response);
                                //logger.logLine("[IGV] "+ "collapse");
                                //out.println("collapse");
                                //response=in.readLine();
                                //logger.logLine("[IGV] " + response);
                                //logger.logLine("[IGV] "+ "snapshot");
                                out.println("snapshot");
                                response = in.readLine();
                                //logger.logLine("[IGV] " + response);

                                File f_snapshot_made = new File("NO");

                                //search png name
                                for (File f_snapshot : f_output_igv_shots_tp_hm.listFiles()) {
                                    if (f_snapshot.isFile()) {
                                        if (f_snapshot.getName().startsWith("chr")) {
                                            f_snapshot_made = f_snapshot;
                                            break;
                                        }

                                    }
                                }

                                if (!f_snapshot_made.exists()) {
                                    //its the all overview
                                    //logger.logLine("[ERROR]: trying to move file which does not exist!");
                                } else {

                                    String command_mv = "mv " + f_snapshot_made.getAbsolutePath() + " " +
                                            f_snapshot_made.getParentFile().getAbsolutePath() + File.separator +
                                            snapshot_name;


                                    Process child = Runtime.getRuntime().exec(command_mv);
                                    int code = child.waitFor();
                                    switch (code) {
                                        case 0:
                                            break;
                                        case 1:
                                            String message = child.getErrorStream().toString();
                                            throw new Exception(message);
                                    }


                                }


                                count++;

                            }
                            out.println("new");
                            String response = in.readLine();
                            //logger.logLine("[IGV] " + response);
                            socket.close();

                            //save session
                            File f_session = new File(f_output_igv_shots.getAbsolutePath() + File.separator +
                                    options_intern.file_suffix_session);
                            igv_save_sessions(f_session, load_bed_files_overall);

                        }
                    }
                }
            }

            double percentage_done = ((i * 1.0) / ordered_tfs_dcg.size());
            DecimalFormat df = new DecimalFormat("0.0");
            String percentage_string = df.format(percentage_done);
            percentage_done = Double.parseDouble(percentage_string) * 100;
            already_done_percentage.add(-1.0);

            if (percentage_done % 10 == 0 && !already_done_percentage.contains(percentage_done)) {
                logger.logLine("[ChIP-ATLAS-IGV] percentage of tasks done: " + percentage_done + "%");
                already_done_percentage.add(percentage_done);
            }

            i++;
        }


        logger.logLine("[ChIP-ATLAS-IGV] Finished evaluation with ChIP Atlas data.");
    }

    public void get_chip_atlas_data() throws Exception {

        File f_output_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_chip_atlas);
        f_output_root.mkdir();

        File f_output_list =
                new File(f_output_root.getAbsolutePath() + File.separator + options_intern.folder_out_chip_atlas_list);
        f_output_list.mkdir();

        File f_output_peak_files = new File(
                f_output_root.getAbsolutePath() + File.separator + options_intern.folder_out_chip_atlas_peak_files);
        f_output_peak_files.mkdir();

        File f_output_list_csv = new File(
                f_output_list.getAbsolutePath() + File.separator + options_intern.file_suffix_chip_atlas_list_csv);


        logger.logLine("[ChIP-ATLAS] Reading available TFs in ChIP-Atlas ...");

        HashMap<String, String> tf_to_url = new HashMap<>();

        BufferedReader br = new BufferedReader(new FileReader(f_output_list_csv.getAbsolutePath()));
        String line = br.readLine();
        String[] split_header = line.split(",");

        int column_gene_version = -1;
        int column_antigen_class = -1;
        int column_antigen = -1;
        int column_cell_type_class = -1;
        int column_url = -1;

        String regex_tissue_type = "";


        String[] split_types = options_intern.chip_atlas_tissue_type.split(";");

        ArrayList<String> all_regex = new ArrayList<>();
        for (String key_type : split_types) {
            String regex_for_this_type = "";

            String[] split_tissue_type = key_type.split(" ");
            for (String s : split_tissue_type) {
                regex_for_this_type += ".*" + s;
            }
            regex_for_this_type += ".*";
            regex_for_this_type = regex_for_this_type.toUpperCase();

            all_regex.add(regex_for_this_type);
        }

        if (all_regex.size() > 1) {
            regex_tissue_type += "(";
        }

        for (int i = 0; i < all_regex.size(); i++) {
            if (i > 0) {
                regex_tissue_type += "|" + all_regex.get(i);
            } else {
                regex_tissue_type += all_regex.get(i);
            }
        }

        if (all_regex.size() > 1) {
            regex_tissue_type += ")";
        }


        for (int i = 0; i < split_header.length; i++) {
            if (split_header[i].matches(options_intern.chip_atlas_column_gene_version)) {
                column_gene_version = i;
            }
            if (split_header[i].matches(options_intern.chip_atlas_column_antigen_class)) {
                column_antigen_class = i;
            }
            if (split_header[i].matches(options_intern.chip_atlas_column_antigen)) {
                column_antigen = i;
            }
            if (split_header[i].matches(options_intern.chip_atlas_column_cell_type_class)) {
                column_cell_type_class = i;
            }
            if (split_header[i].matches(options_intern.chip_atlas_column_url)) {
                column_url = i;
            }
        }

        while ((line = br.readLine()) != null) {
            String[] split = line.split(",");
            if (!split[column_gene_version].equals(options_intern.chip_atlas_genome_version)) {
                continue;
            }
            if (!split[column_antigen_class].equals("TFs and others")) {
                continue;
            }
            if (!split[column_cell_type_class].toUpperCase().matches(regex_tissue_type)) {
                continue;
            }

            String url = split[split.length - 1];
            String tf = split[column_antigen];
            if (!tf.equals("")) {
                tf_to_url.put(tf.toUpperCase(), url);
            }

        }
        br.close();

        logger.logLine("[ChIP-ATLAS] downloading available TF data and creating IGV index files");
        logger.logLine("[ChIP-ATLAS] waiting ...");

        //read in discounted cumulative gain tfs
        File f_input_dcg_result = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution +
                        File.separator + options_intern.folder_out_distribution_dcg + File.separator +
                        options_intern.file_suffix_distribution_analysis_dcg);
        ArrayList<String> ordered_tfs_dcg = new ArrayList<>();
        HashSet<String> unordered_tfs_dcg = new HashSet<>();
        HashMap<String, HashSet<String>> tf_to_splits = new HashMap<>();

        BufferedReader br_dcg = new BufferedReader(new FileReader(f_input_dcg_result));
        String line_dcg = br_dcg.readLine();
        while ((line_dcg = br_dcg.readLine()) != null) {
            String[] split = line_dcg.split("\t");

            String[] split_tf = split[1].split("\\.");

            HashSet<String> splits = new HashSet<>();
            for (String tf_name_chipatlas : split_tf) {
                if (!tf_name_chipatlas.startsWith(" ")) {
                    unordered_tfs_dcg.add(tf_name_chipatlas);
                }

                splits.add(tf_name_chipatlas);
            }

            ordered_tfs_dcg.add(split[1]);
            tf_to_splits.put(split[1].toUpperCase(), splits);


        }
        br_dcg.close();

        for (int i = 0; i < ordered_tfs_dcg.size(); i++) {
            boolean found_tf = false;

            ArrayList<String> found_tfs_to_download = new ArrayList<>();

            HashSet<String> split_of_tfs = tf_to_splits.get(ordered_tfs_dcg.get(i).toUpperCase());
            for (String split_tf : split_of_tfs) {
                if (tf_to_url.containsKey(split_tf)) {
                    found_tf = true;
                    found_tfs_to_download.add(split_tf);
                }
            }


            if (found_tf) {
                int rank = i + 1;

                File f_download_tf = new File(
                        f_output_peak_files.getAbsolutePath() + File.separator + rank + "_" + ordered_tfs_dcg.get(i));
                f_download_tf.mkdir();

                for (String tf_found_split : found_tfs_to_download) {
                    File f_download_tf_file =
                            new File(f_download_tf.getAbsolutePath() + File.separator + tf_found_split + ".bed");

                    if(f_download_tf_file.exists())
                    {
                        continue;
                    }

                    // Create a new trust manager that trust all certificates
                    TrustManager[] trustAllCerts = new TrustManager[]{new X509TrustManager() {
                        public java.security.cert.X509Certificate[] getAcceptedIssuers() {
                            return null;
                        }

                        public void checkClientTrusted(java.security.cert.X509Certificate[] certs, String authType) {
                        }

                        public void checkServerTrusted(java.security.cert.X509Certificate[] certs, String authType) {
                        }
                    }};
                    // Activate the new trust manager
                    try {
                        SSLContext sc = SSLContext.getInstance("SSL");
                        sc.init(null, trustAllCerts, new java.security.SecureRandom());
                        HttpsURLConnection.setDefaultSSLSocketFactory(sc.getSocketFactory());
                    } catch (Exception e) {
                        //e.printStackTrace();
                    }
                    // And as before now you can use URL and URLConnection
                    URL url = new URL(tf_to_url.get(tf_found_split));
                    URLConnection connection = url.openConnection();
                    InputStream is = connection.getInputStream();

                    try (OutputStream outputStream = new FileOutputStream(f_download_tf_file)) {
                        IOUtils.copy(is, outputStream);
                    } catch (FileNotFoundException e) {
                        // handle exception here
                        e.printStackTrace();
                    } catch (IOException e) {
                        // handle exception here
                        //e.printStackTrace();
                    }

                    //check if valid bed file, if not, delete
                    Path path = Paths.get(f_download_tf_file.getAbsolutePath());
                    long size = Files.size(path);
                    if (size < 300) {
                        f_download_tf_file.delete();
                        f_download_tf.delete();
                    } else {
                        DecimalFormat df = new DecimalFormat("0.00");

                        double percentage = ((i * 1.0 / ordered_tfs_dcg.size())) * 100;
                        if (percentage <= 100) {
                            logger.logLine("[ChIP-ATLAS] processed " + df.format(percentage) + "% TFs.");
                        }
                    }

                    //create index file for IGV
                    String command_create_index = options_intern.igv_path_to_igv + File.separator + "igvtools index " +
                            f_download_tf_file.getAbsolutePath();

                    Process child = Runtime.getRuntime().exec(command_create_index);
                    int code = child.waitFor();
                    switch (code) {
                        case 0:
                            break;
                        case 1:
                            String message = child.getErrorStream().toString();
                            throw new Exception(message);
                    }

                }


            }
        }


        logger.logLine("[ChIP-ATLAS] finished downloading all available ChIP-Atlas Data");

    }

    /**
     * retrieve available TF ChIP-seq from ChIP-ATLAS
     *
     * @throws Exception
     */
    public void get_chip_atlas_data_list() throws Exception {
        logger.logLine("[ChIP-ATLAS] get ChIP-Atlas Data for Discounted Cumulative Gain TFs!");
        logger.logLine("[ChIP-ATLAS] Genome version: " + options_intern.chip_atlas_genome_version);
        logger.logLine("[ChIP-ATLAS] Tissue type: " + options_intern.chip_atlas_tissue_type);

        File f_output_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_chip_atlas);
        f_output_root.mkdir();

        File f_output_list =
                new File(f_output_root.getAbsolutePath() + File.separator + options_intern.folder_out_chip_atlas_list);
        f_output_list.mkdir();

        File f_output_peak_files = new File(
                f_output_root.getAbsolutePath() + File.separator + options_intern.folder_out_chip_atlas_peak_files);
        f_output_peak_files.mkdir();

        /**
         * Download list of latest available ChIP-ATLAS files
         */
        logger.logLine("[ChIP-ATLAS] Downloading list of available TF ChIP-seq Files from " +
                options_intern.chip_atlas_url_to_list);
        logger.logLine("[ChIP-ATLAS] waiting ...");
        File f_output_list_zip = new File(
                f_output_list.getAbsolutePath() + File.separator + options_intern.file_suffix_chip_atlas_list_zipped);
        if (!f_output_list_zip.exists()) {
            // Create a new trust manager that trust all certificates
            TrustManager[] trustAllCerts = new TrustManager[]{new X509TrustManager() {
                public java.security.cert.X509Certificate[] getAcceptedIssuers() {
                    return null;
                }

                public void checkClientTrusted(java.security.cert.X509Certificate[] certs, String authType) {
                }

                public void checkServerTrusted(java.security.cert.X509Certificate[] certs, String authType) {
                }
            }};
            // Activate the new trust manager
            try {
                SSLContext sc = SSLContext.getInstance("SSL");
                sc.init(null, trustAllCerts, new java.security.SecureRandom());
                HttpsURLConnection.setDefaultSSLSocketFactory(sc.getSocketFactory());
            } catch (Exception e) {
                e.printStackTrace();
            }
            // And as before now you can use URL and URLConnection
            URL url = new URL(options_intern.chip_atlas_url_to_list);
            URLConnection connection = url.openConnection();
            InputStream is = connection.getInputStream();

            try (OutputStream outputStream = new FileOutputStream(f_output_list_zip)) {
                IOUtils.copy(is, outputStream);
            } catch (FileNotFoundException e) {
                // handle exception here
                e.printStackTrace();
            } catch (IOException e) {
                // handle exception here
                e.printStackTrace();
            }
        }

        logger.logLine("[ChIP-ATLAS] Unzipping file: " + f_output_list.getAbsolutePath());
        File f_output_list_csv = new File(
                f_output_list.getAbsolutePath() + File.separator + options_intern.file_suffix_chip_atlas_list_csv);
        if (!f_output_list_csv.exists()) {
            //unzip file
            String command_edited =
                    "unzip " + f_output_list_zip + " -d " + f_output_list_csv.getParentFile().getAbsolutePath();
            logger.logLine("[ChIP-ATLAS] executing command: " + command_edited);
            Process child = Runtime.getRuntime().exec(command_edited);
            int code = child.waitFor();
            switch (code) {
                case 0:
                    break;
                case 1:
                    String message = child.getErrorStream().toString();
                    throw new Exception(message);
            }

        }
    }

    /**
     * retrieve top k target genes of DCG found TFs
     */
    public void get_top_k_target_genes_dcg() throws IOException {
        logger.logLine("[EVALUATION-PREPROCESSING] retrieve top " + options_intern.plot_top_k_genes +
                " target genes foreach discounted cumulative gain TF");

        File f_output_results_root = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_out_target_genes_dcg);
        f_output_results_root.mkdir();

        File f_input_affinity_values_root = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_tepic_postprocessing + File.separator +
                options_intern.folder_name_tepic_postprocessing_input);

        //retrieve ordered tfs
        File f_input_dcg_result = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution +
                        File.separator + options_intern.folder_out_distribution_dcg + File.separator +
                        options_intern.file_suffix_distribution_analysis_dcg);
        ArrayList<String> ordered_tfs_dcg = new ArrayList<>();
        HashSet<String> unordered_tfs_dcg = new HashSet<>();

        BufferedReader br_dcg = new BufferedReader(new FileReader(f_input_dcg_result));
        String line_dcg = br_dcg.readLine();
        while ((line_dcg = br_dcg.readLine()) != null) {
            String[] split = line_dcg.split("\t");
            ordered_tfs_dcg.add(split[1]);
            unordered_tfs_dcg.add(split[1]);
        }
        br_dcg.close();

        HashSet<String> found_tfs = new HashSet<>();

        HashMap<String, String> ensg_gene_symbol_map = new HashMap<>();
        BufferedReader br_ensg_gene_symbol =
                new BufferedReader(new FileReader(new File(options_intern.tepic_ensg_symbol)));
        String line_ensg_gene_symbol = br_ensg_gene_symbol.readLine();
        while ((line_ensg_gene_symbol = br_ensg_gene_symbol.readLine()) != null) {
            String[] split = line_ensg_gene_symbol.split("\t");
            if (split.length > 1) {
                ensg_gene_symbol_map.put(split[0].toUpperCase(), split[1].toUpperCase());
            }
        }
        br_ensg_gene_symbol.close();


        for (File file_dir_tp : f_input_affinity_values_root.listFiles()) {
            if (file_dir_tp.isDirectory()) {
                File f_output_tp =
                        new File(f_output_results_root.getAbsolutePath() + File.separator + file_dir_tp.getName());
                f_output_tp.mkdir();

                for (File file_dir_tp_hm : file_dir_tp.listFiles()) {
                    if (file_dir_tp_hm.isDirectory()) {
                        File f_output_tp_hm =
                                new File(f_output_tp.getAbsolutePath() + File.separator + file_dir_tp_hm.getName());
                        f_output_tp_hm.mkdir();

                        HashMap<String, TF_TargetGene_DCG> tf_to_target_genes = new HashMap<>();
                        for (String key_tf : ordered_tfs_dcg) {
                            TF_TargetGene_DCG dcg = new TF_TargetGene_DCG();
                            dcg.TF = key_tf;
                            tf_to_target_genes.put(key_tf, dcg);
                        }

                        for (File file_dir_tp_hm_file : file_dir_tp_hm.listFiles()) {
                            if (file_dir_tp_hm_file.isFile()) {
                                HashMap<Integer, String> index_to_tf = new HashMap<>();
                                HashMap<String, Integer> tf_to_index = new HashMap<>();

                                BufferedReader br = new BufferedReader(new FileReader(file_dir_tp_hm_file));
                                String line = br.readLine();
                                String[] split_header = line.split("\t");
                                for (int i = 1; i < split_header.length; i++) {
                                    String name_tf = split_header[i].split("_")[0];
                                    name_tf = name_tf.replaceAll(":", "\\.");

                                    if (unordered_tfs_dcg.contains(name_tf.toUpperCase())) {
                                        index_to_tf.put(i, name_tf.toUpperCase());
                                        tf_to_index.put(name_tf.toUpperCase(), i);
                                    }
                                }

                                while ((line = br.readLine()) != null) {
                                    String[] split = line.split("\t");

                                    String gene_name = "";
                                    if (ensg_gene_symbol_map.containsKey(split[0])) {
                                        gene_name = ensg_gene_symbol_map.get(split[0]).toUpperCase();
                                    } else {
                                        continue;
                                    }

                                    for (int i = 1; i < split.length; i++) {
                                        if (index_to_tf.containsKey(i)) {
                                            String tf = index_to_tf.get(i);

                                            TF_TargetGene_DCG dcg = tf_to_target_genes.get(tf);

                                            if (dcg.target_gene_affinity_values.containsKey(gene_name)) {
                                                dcg.target_gene_affinity_values.get(gene_name).affinity_values.add(
                                                        Double.parseDouble(split[i]));
                                            } else {
                                                TargetGene_DCG tdcg = new TargetGene_DCG();
                                                tdcg.target_gene = gene_name;
                                                tdcg.affinity_values.add(Double.parseDouble(split[i]));
                                                dcg.target_gene_affinity_values.put(gene_name, tdcg);
                                            }
                                        }
                                    }

                                }
                                br.close();
                            }
                        }

                        //calculate means
                        for (String key_tf : tf_to_target_genes.keySet()) {
                            TF_TargetGene_DCG dcg = tf_to_target_genes.get(key_tf);
                            if (dcg.target_gene_affinity_values.isEmpty()) {
                                continue;
                            }

                            ArrayList<TargetGene_DCG> tdcg_ordered = dcg.get_ordered_target_gene_list();

                            if (tdcg_ordered.isEmpty()) {
                                continue;
                            }

                            File f_output_tf = new File(f_output_tp_hm.getAbsolutePath() + File.separator + key_tf);
                            f_output_tf.mkdir();
                            File f_output_tf_file =
                                    new File(f_output_tf.getAbsolutePath() + File.separator + key_tf + ".csv");
                            StringBuilder sb = new StringBuilder();
                            sb.append("TARGET_GENE\tAFFINITY_SCORE");
                            sb.append("\n");
                            for (TargetGene_DCG tdcg : tdcg_ordered) {
                                sb.append(tdcg.target_gene);
                                sb.append("\t");
                                sb.append(tdcg.mean_affinity_value);
                                sb.append("\n");
                            }

                            found_tfs.add(key_tf);

                            BufferedWriter bw = new BufferedWriter(new FileWriter(f_output_tf_file));
                            bw.write(sb.toString());
                            bw.close();
                        }
                    }
                }
            }
        }

        logger.logLine("[EVALUATION-PREPROCESSING] finished retrieving target genes.");
    }

    /**
     * calculate discounted cumulative gain rank for distribution analysis and set up html report for it
     */
    public void calculate_discounted_cumulative_gain_rank_distribution_analysis() throws Exception {
        logger.logLine("[DISTRIBUTION-ANALYSIS] Start calculate ranks overall groups.");

        //calculate ranks overall groups
        HashMap<String, HashMap<String, Integer>> group_analysis_distr_stats = new HashMap<>();
        HashMap<String, Analysis_distribution_stats_object> group_analysis_distr_stats_objects = new HashMap<>();
        HashSet<String> available_modifications = new HashSet<>();

        HashMap<String, Integer> final_ranks_only = new HashMap<>();

        File f_distr_analysis_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution);

        File f_distr_analysis_stats_root = new File(f_distr_analysis_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_stats);
        File f_distr_stats_ALL = new File(f_distr_analysis_stats_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_stats_ALL);
        File f_distr_stats_HM = new File(f_distr_analysis_stats_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_stats_HM);


        for (File fileDir_stat : f_distr_stats_HM.listFiles()) {
            Analysis_distribution_stats_object current_object = get_distribution_analysis_stats_ordered(fileDir_stat);
            group_analysis_distr_stats_objects.put(fileDir_stat.getName(), current_object);
        }

        Analysis_distribution_stats_object all_object = get_distribution_analysis_stats_ordered(f_distr_stats_ALL);
        group_analysis_distr_stats_objects.put(options_intern.distribution_analysis_all_name, all_object);

        int not_found_factor = 0;
        int group_number = 0;


        for (String key_group : group_analysis_distr_stats_objects.keySet()) {
            available_modifications.add(key_group);

            HashMap<String, Integer> tf_to_ranks = new HashMap<>();
            int rank = 0;

            ArrayList<Analysis_distribution_stats> currents_tfs_to_rank =
                    group_analysis_distr_stats_objects.get(key_group).ordered_tfs;

            for (Analysis_distribution_stats tf : currents_tfs_to_rank) {
                tf_to_ranks.put(tf.label, rank + 1);
                rank++;
            }
            group_analysis_distr_stats.put(key_group, tf_to_ranks);

            not_found_factor += currents_tfs_to_rank.size();
            group_number++;
        }

        not_found_factor /= group_number;

        ArrayList<Analysis_distribution_stats_cumulative_gain> fin_ranks = new ArrayList<>();
        HashSet<String> already_calculated_tfs = new HashSet<>();

        for (String key_group : group_analysis_distr_stats.keySet()) {
            HashMap<String, Integer> current_tfs = group_analysis_distr_stats.get(key_group);
            for (String key_tf : current_tfs.keySet()) {
                if (already_calculated_tfs.contains(key_tf)) {
                    continue;
                }

                double rank_overall = 0;

                double rank_current_group = (current_tfs.size() - current_tfs.get(key_tf)) / (current_tfs.size() * 1.0);
                rank_overall += rank_current_group;

                HashMap<String, Analysis_distribution_stats> current_group_stats = new HashMap<>();
                HashMap<String, Analysis_distribution_stats> current_group_background = new HashMap<>();

                Analysis_distribution_stats_object xx = group_analysis_distr_stats_objects.get(key_group);
                current_group_stats.put(key_group, xx.ordered_tfs.get(current_tfs.get(key_tf) - 1));
                current_group_background.put(key_group, xx.background);


                for (String key_group_clash : group_analysis_distr_stats.keySet()) {
                    if (key_group_clash.equals(key_group)) {
                        continue;
                    }
                    HashMap<String, Integer> current_tfs_clash = group_analysis_distr_stats.get(key_group_clash);
                    if (current_tfs_clash.containsKey(key_tf)) {
                        HashMap<String, Integer> current_tfs_clash_x = group_analysis_distr_stats.get(key_group_clash);

                        Analysis_distribution_stats_object xx_clash =
                                group_analysis_distr_stats_objects.get(key_group_clash);
                        current_group_stats.put(key_group_clash,
                                xx_clash.ordered_tfs.get(current_tfs_clash_x.get(key_tf) - 1));
                        current_group_background.put(key_group_clash, xx_clash.background);
                        double rank_current_group_clash =
                                (current_tfs_clash_x.size() - current_tfs_clash_x.get(key_tf)) /
                                        (current_tfs_clash_x.size() * 1.0);
                        rank_overall += rank_current_group_clash;

                    }

                }
                Analysis_distribution_stats_cumulative_gain ac = new Analysis_distribution_stats_cumulative_gain();
                ac.label = key_tf;
                ac.rank = rank_overall;
                ac.group_object = current_group_stats;
                ac.group_background = current_group_background;
                fin_ranks.add(ac);
                already_calculated_tfs.add(key_tf);
            }
        }

        Collections.sort(fin_ranks);
        //write dcg ranks
        File f_out_dcg = new File(
                f_distr_analysis_root.getAbsolutePath() + File.separator + options_intern.folder_out_distribution_dcg);
        f_out_dcg.mkdir();

        File f_out_dcg_file = new File(
                f_out_dcg.getAbsolutePath() + File.separator + options_intern.file_suffix_distribution_analysis_dcg);
        StringBuilder sb_dcg = new StringBuilder();
        sb_dcg.append("RANK\tTF\tSCORE\n");
        int i_rank = 1;
        for (Analysis_distribution_stats_cumulative_gain ac : fin_ranks) {
            sb_dcg.append(ac.toString(i_rank));
            sb_dcg.append("\n");
            i_rank++;
        }

        BufferedWriter bw_dcg = new BufferedWriter(new FileWriter(f_out_dcg_file));
        bw_dcg.write(sb_dcg.toString());
        bw_dcg.close();

        //CREATE HTML_REPORT
        File folder_website_out = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_website +
                        File.separator + options_intern.folder_out_website_htmls_distribution_analysis +
                        File.separator + options_intern.folder_out_website_htmls_distribution_analysis_cumulative_gain);
        folder_website_out.mkdir();

        String html_tail = "</body>\n" + "</html>";

        StringBuilder sb_cumulative_gain_website = new StringBuilder();
        sb_cumulative_gain_website.append(get_header_html(options_intern.html_report_levels_3_steps,
                "Distribution Analysis: " + "CUMULATIVE GAIN"));

        //get gene counts
        HashMap<String, HashMap<String, String>> tp_tf_gene_count = new HashMap<>();

        File f_input_genecounts_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_analysis_data +
                        File.separator + options_intern.folder_out_analysis_data_HM_LEVEL);
        for (File fileDir : f_input_genecounts_root.listFiles()) {
            if (fileDir.isDirectory()) {
                for (File fileDir_gc : fileDir.listFiles()) {
                    if (fileDir_gc.isFile()) {
                        BufferedReader br = new BufferedReader(new FileReader(fileDir_gc));
                        String line = br.readLine();
                        String[] split_header = line.split("\t");

                        while ((line = br.readLine()) != null) {
                            String[] split_line = line.split("\t");

                            String name_tf = split_line[0];

                            for (int i = 1; i < split_header.length; i++) {
                                HashMap<String, String> tf_genecount_for_tp;

                                if (tp_tf_gene_count.containsKey(split_header[i])) {
                                    tf_genecount_for_tp = tp_tf_gene_count.get(split_header[i]);
                                } else {
                                    tf_genecount_for_tp = new HashMap<>();
                                }

                                String x = split_line[i];
                                tf_genecount_for_tp.put(name_tf.toUpperCase(), x);

                                tp_tf_gene_count.put(split_header[i], tf_genecount_for_tp);

                            }

                        }
                        br.close();
                    }
                }
            }
        }

        DecimalFormat df = new DecimalFormat("0.00");

        for (int i = 0; i < fin_ranks.size(); i++) {
            int rank_overall = i + 1;
            Analysis_distribution_stats_cumulative_gain as = fin_ranks.get(i);
            HashMap<String, Analysis_distribution_stats> group_object = as.group_object;

            final_ranks_only.put(as.label.toUpperCase(), rank_overall);

            sb_cumulative_gain_website.append("<div id='" + as.label.toUpperCase() + "' class='w3-content'>\n");
            sb_cumulative_gain_website.append(
                    "<h2>" + rank_overall + ". <a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=" +
                            as.label.toUpperCase() + "' target='_blank'><button class='button'>" +
                            as.label.toUpperCase() + "</button></a></h2>\n");

            sb_cumulative_gain_website.append("<h5>Discounted cumulative gain details:</h5>\n");

            sb_cumulative_gain_website.append("\t\t<table style=\"width:80%;\">\n");
            sb_cumulative_gain_website.append("\t\t\t<tr>\n");
            sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
            sb_cumulative_gain_website.append("Modification");
            sb_cumulative_gain_website.append("\t\t\t\t</th>\n");
            sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
            sb_cumulative_gain_website.append("Availability");
            sb_cumulative_gain_website.append("\t\t\t\t</th>\n");
            sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
            sb_cumulative_gain_website.append("Rank");
            sb_cumulative_gain_website.append("\t\t\t\t</th>\n");
            sb_cumulative_gain_website.append("\t\t</tr>\n");

            for (String key_group : available_modifications) {
                sb_cumulative_gain_website.append("\t\t\t<tr>\n");

                if (as.group_object.containsKey(key_group)) {
                    Analysis_distribution_stats as_obj = as.group_object.get(key_group);

                    HashMap<String, Integer> actual_ranks = group_analysis_distr_stats.get(key_group);

                    sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
                    if (key_group.equals(options_intern.distribution_analysis_all_name)) {
                        sb_cumulative_gain_website.append("<a href='" + ".." + File.separator +
                                options_intern.folder_out_website_htmls_distribution_analysis_ALL + File.separator +
                                options_intern.html_report_home_regression_distribution_analysis_all + "#" +
                                as.label.toUpperCase() +
                                "' target='_blank'><button class='button' style='background-color:#67cb67;'>");

                    } else {
                        sb_cumulative_gain_website.append("<a href='" + ".." + File.separator +
                                options_intern.folder_out_website_htmls_distribution_analysis_HM + File.separator +
                                key_group + ".html#" + as.label.toUpperCase() +
                                "' target='_blank'><button class='button' style='background-color:#67cb67;'>");
                    }
                    sb_cumulative_gain_website.append(key_group);
                    sb_cumulative_gain_website.append("</button></a>\n");
                    sb_cumulative_gain_website.append("\t\t\t\t</th>\n");
                    sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
                    sb_cumulative_gain_website.append("<img src=\".." + File.separator + ".." + File.separator +
                            options_intern.folder_out_website_basics + File.separator +
                            options_intern.folder_out_website_basics_website + File.separator +
                            options_intern.folder_out_website_basics_website_images + File.separator +
                            "is_available.png" + "\" style=\"width:50px;height:50px;\"/>");
                    sb_cumulative_gain_website.append("\t\t\t\t</th>\n");
                    sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
                    sb_cumulative_gain_website.append(
                            actual_ranks.get(as.label.toUpperCase()) + "/" + actual_ranks.size());
                    sb_cumulative_gain_website.append("\t\t\t\t</th>\n");
                } else {
                    sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
                    sb_cumulative_gain_website.append(key_group);
                    sb_cumulative_gain_website.append("\t\t\t\t</th>\n");
                    sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
                    sb_cumulative_gain_website.append("<img src=\".." + File.separator + ".." + File.separator +
                            options_intern.folder_out_website_basics + File.separator +
                            options_intern.folder_out_website_basics_website + File.separator +
                            options_intern.folder_out_website_basics_website_images + File.separator +
                            "not_available.png" + "\" style=\"width:50px;height:50px;\"/>");
                    sb_cumulative_gain_website.append("\t\t\t\t</th>\n");
                    sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
                    sb_cumulative_gain_website.append("-");
                    sb_cumulative_gain_website.append("\t\t\t\t</th>\n");
                }
                sb_cumulative_gain_website.append("\t\t</tr>\n");
            }
            sb_cumulative_gain_website.append("\t\t<tr>\n");
            sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
            sb_cumulative_gain_website.append("OVERALL SCORE");
            sb_cumulative_gain_website.append("\t\t\t\t</th>\n");
            sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
            sb_cumulative_gain_website.append("");
            sb_cumulative_gain_website.append("\t\t\t\t</th>\n");
            sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
            sb_cumulative_gain_website.append(df.format(as.rank));
            sb_cumulative_gain_website.append("\t\t\t\t</th>\n");
            sb_cumulative_gain_website.append("\t\t</tr>\n");


            sb_cumulative_gain_website.append("\t\t</table>\n");

            sb_cumulative_gain_website.append("<h5>Gene Counts Table </h5>\n");

            sb_cumulative_gain_website.append("\t\t<table style=\"width:80%\">\n");
            sb_cumulative_gain_website.append("\t\t\t<tr>\n");
            sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
            sb_cumulative_gain_website.append("TF");
            sb_cumulative_gain_website.append("\t\t\t\t</th>\n");

            HashSet<String> tps = new HashSet<>();

            for (String key : tp_tf_gene_count.keySet()) {
                sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
                sb_cumulative_gain_website.append(key);
                sb_cumulative_gain_website.append("\t\t\t\t</th>\n");

                tps.add(key);
            }
            sb_cumulative_gain_website.append("\t\t\t</tr>\n");
            sb_cumulative_gain_website.append("\t\t\t<tr>\n");
            sb_cumulative_gain_website.append("\t\t\t\t<th>\n");
            sb_cumulative_gain_website.append(as.label.toUpperCase());
            sb_cumulative_gain_website.append("\t\t\t\t</th>\n");

            for (String key : tps) {
                String value = "0.0";
                HashMap<String, String> current_x = tp_tf_gene_count.get(key);
                if (current_x.containsKey(as.label.toUpperCase())) {
                    value = current_x.get(as.label.toUpperCase());
                }

                sb_cumulative_gain_website.append("\t\t\t\t<th>");
                sb_cumulative_gain_website.append(value);
                sb_cumulative_gain_website.append("\t\t\t\t</th>\n");
            }

            sb_cumulative_gain_website.append("\t\t\t</tr>\n");

            sb_cumulative_gain_website.append("\t\t</table>");
            sb_cumulative_gain_website.append("</div>");
        }

        sb_cumulative_gain_website.append(html_tail);
        BufferedWriter bw = new BufferedWriter(new FileWriter(folder_website_out.getAbsolutePath() + File.separator +
                options_intern.html_report_home_regression_distribution_analysis_all));
        bw.write(sb_cumulative_gain_website.toString());
        bw.close();

        StringBuilder sb_home_home = new StringBuilder();
        sb_home_home.append(
                get_header_html(options_intern.html_report_levels_home, options_intern.html_report_levels_home));

        //CREATE EXPANDABLE BUTTONS FOR DISTRIBUTION ANALYSIS OVERVIEW
        sb_home_home.append(
                "<button class=\"button_expandable\" id=\"button_distribution_analysis\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_distribution_analysis','table_distribution_analysis')\"> TF-TG-Score Distribution Analysis - OVERVIEW\n");
        sb_home_home.append(
                "<div style=\"display: none;background-color: white;color:black;\" id=\"table_distribution_analysis\">\n");

        ArrayList<File> stat_files = new ArrayList<>();
        stat_files.add(f_distr_stats_ALL);
        for (File fileDir_stat : f_distr_stats_HM.listFiles()) {
            stat_files.add(fileDir_stat);
        }

        sb_home_home.append(
                get_html_table_found_tfs_in_distribution_analysis(stat_files, options_intern.html_report_levels_home,
                        final_ranks_only));

        sb_home_home.append("</div>\n</button>\n");


        //CREATE EXPANDABLE BUTTONS FOR REGRESSION COEFFICIENT ANALYSIS
        sb_home_home.append(
                "<button class=\"button_expandable\" id=\"button_regression_coefficient_analysis\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_regression_coefficient_analysis','table_regression_coefficient_analysis')\"> Regression Coefficient Analysis - OVERVIEW\n");
        sb_home_home.append(
                "<div style=\"display: none;background-color: white;color:black;\" id=\"table_regression_coefficient_analysis\">\n");
        for (Double d : options_intern.plot_th_coefficient) {
            sb_home_home.append(
                    write_regression_coeffecient_analysis_found_table_html(d, options_intern.html_report_levels_home));
        }

        sb_home_home.append("</div>\n</button>\n");

        sb_home_home.append(html_tail);

        File f_output_website_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_website);
        BufferedWriter bw_home_home = new BufferedWriter(new FileWriter(
                f_output_website_root.getAbsolutePath() + File.separator + options_intern.html_report_home_home));
        bw_home_home.write(sb_home_home.toString());
        bw_home_home.close();

        logger.logLine("[DISTRIBUTION-ANALYSIS] Finish calculate ranks overall groups.");
    }

    /**
     * CREATE HTML REPORT FOR DISTRIBUTION ANALYSIS
     */
    public void create_overview_html_report_distribution() throws Exception {
        logger.logLine("[DISTRIBUTION-ANALYSIS-HTML-REPORT] Create HTML report.");

        File f_output_website_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_website);

        File f_website_distr_analysis_html_folder = new File(f_output_website_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_website_htmls_distribution_analysis);
        f_website_distr_analysis_html_folder.mkdir();
        File f_website_distr_analysis_html_folder_ALL = new File(
                f_website_distr_analysis_html_folder.getAbsolutePath() + File.separator +
                        options_intern.folder_out_website_htmls_distribution_analysis_ALL);
        f_website_distr_analysis_html_folder_ALL.mkdir();
        File f_website_distr_analysis_html_folder_HM = new File(
                f_website_distr_analysis_html_folder.getAbsolutePath() + File.separator +
                        options_intern.folder_out_website_htmls_distribution_analysis_HM);
        f_website_distr_analysis_html_folder_HM.mkdir();
        File f_website_distr_analysis_html_folder_cumulative_gain = new File(
                f_website_distr_analysis_html_folder.getAbsolutePath() + File.separator +
                        options_intern.folder_out_website_htmls_distribution_analysis_cumulative_gain);
        f_website_distr_analysis_html_folder_cumulative_gain.mkdir();

        File f_distr_analysis_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution);

        File f_distr_analysis_stats_root = new File(f_distr_analysis_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_stats);
        File f_distr_stats_ALL = new File(f_distr_analysis_stats_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_stats_ALL);
        File f_distr_stats_HM = new File(f_distr_analysis_stats_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_stats_HM);

        File f_website_distr_analysis_html_folder_ALL_FILE = new File(
                f_website_distr_analysis_html_folder_ALL.getAbsolutePath() + File.separator +
                        options_intern.html_report_home_regression_distribution_analysis_all);

        File f_distr_plots_html_report = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_website +
                        File.separator + options_intern.folder_out_website_plots_distribution_analysis);
        f_distr_plots_html_report.mkdir();

        File f_root_plots = new File(f_distr_analysis_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_plots);

        //copy plots
        String command =
                "cp -u -r " + f_root_plots.getAbsolutePath() + " " + f_distr_plots_html_report.getAbsolutePath();
        Process child = Runtime.getRuntime().exec(command);
        logger.logLine("[DISTRIBUTION-ANALYSIS-HTML-REPORT] Copy plots for transferable html report: " + command);
        int code = child.waitFor();
        switch (code) {
            case 0:
                break;
            case 1:
                String message = child.getErrorStream().toString();
                throw new Exception(message);
        }

        HashMap<String, HashMap<String, String>> tp_tf_gene_count = new HashMap<>();

        File f_input_genecounts_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_analysis_data +
                        File.separator + options_intern.folder_out_analysis_data_HM_LEVEL);
        for (File fileDir : f_input_genecounts_root.listFiles()) {
            if (fileDir.isDirectory()) {
                for (File fileDir_gc : fileDir.listFiles()) {
                    if (fileDir_gc.isFile()) {
                        BufferedReader br = new BufferedReader(new FileReader(fileDir_gc));
                        String line = br.readLine();
                        String[] split_header = line.split("\t");

                        while ((line = br.readLine()) != null) {
                            String[] split_line = line.split("\t");

                            String name_tf = split_line[0];

                            for (int i = 1; i < split_header.length; i++) {
                                HashMap<String, String> tf_genecount_for_tp;

                                if (tp_tf_gene_count.containsKey(split_header[i])) {
                                    tf_genecount_for_tp = tp_tf_gene_count.get(split_header[i]);
                                } else {
                                    tf_genecount_for_tp = new HashMap<>();
                                }

                                String x = split_line[i];
                                tf_genecount_for_tp.put(name_tf.toUpperCase(), x);

                                tp_tf_gene_count.put(split_header[i], tf_genecount_for_tp);

                            }

                        }
                        br.close();
                    }
                }
            }
        }


        String html_header = get_header_html("HOME", options_intern.analysis_types_distribution_analysis);
        String html_tail = "</body>\n" + "</html>";

        File html_home_distribution_analysis = new File(f_output_website_root.getAbsolutePath() + File.separator +
                options_intern.html_report_home_distribution_analysis);
        StringBuilder sb_home_distribution_analysis = new StringBuilder();

        sb_home_distribution_analysis.append(html_header);
        sb_home_distribution_analysis.append(
                "<div class='w3-row-padding w3-padding-64 w3-container w3-content'><a href='." + File.separator +
                        f_website_distr_analysis_html_folder_cumulative_gain.getParentFile().getName() +
                        File.separator + f_website_distr_analysis_html_folder_cumulative_gain.getName() +
                        File.separator + options_intern.html_report_home_regression_distribution_analysis_all +
                        "' target='_blank'><button class='button_expandable'>DISCOUNTED CUMULATIVE GAIN</button></a></div>");
        sb_home_distribution_analysis.append(
                "<div class='w3-row-padding w3-padding-64 w3-container w3-content'><a href='." + File.separator +
                        f_website_distr_analysis_html_folder_ALL.getParentFile().getName() + File.separator +
                        f_website_distr_analysis_html_folder_ALL.getName() + File.separator +
                        options_intern.html_report_home_regression_distribution_analysis_all +
                        "' target='_blank'><button class='button_expandable'>ALL</button></a></div>");

        //sb_home_distribution_analysis.append("<div class='w3-row-padding w3-padding-64 w3-container w3-content'><a href='"+f_website_distr_analysis_html_folder_cumulative_gain.getAbsolutePath()+File.separator+options_intern.html_report_home_regression_distribution_analysis_all+"' target='_blank'><button class='button_expandable'>DISCOUNTED CUMULATIVE GAIN</button></a></div>");
        //sb_home_distribution_analysis.append("<div class='w3-row-padding w3-padding-64 w3-container w3-content'><a href='"+f_website_distr_analysis_html_folder_ALL.getAbsolutePath()+File.separator+options_intern.html_report_home_regression_distribution_analysis_all+"' target='_blank'><button class='button_expandable'>ALL</button></a></div>");


        for (File fileDir_stat : f_distr_stats_HM.listFiles()) {
            File f_website_html_HM = new File(
                    f_website_distr_analysis_html_folder_HM.getAbsolutePath() + File.separator +
                            fileDir_stat.getName() + ".html");
            //sb_home_distribution_analysis.append("<div class='w3-row-padding w3-padding-64 w3-container w3-content'><a href='"+f_website_html_HM.getAbsolutePath()+"' target='_blank'><button class='button_expandable'>"+fileDir_stat.getName()+"</button></a></div>");
            sb_home_distribution_analysis.append(
                    "<div class='w3-row-padding w3-padding-64 w3-container w3-content'><a href='." + File.separator +
                            f_website_html_HM.getParentFile().getParentFile().getName() + File.separator +
                            f_website_html_HM.getParentFile().getName() + File.separator + f_website_html_HM.getName() +
                            "' target='_blank'><button class='button_expandable'>" + fileDir_stat.getName() +
                            "</button></a></div>");

            write_html_distribution_analysis_plots_hm_page(f_website_html_HM, fileDir_stat, "HM", tp_tf_gene_count);
        }

        write_html_distribution_analysis_plots_hm_page(f_website_distr_analysis_html_folder_ALL_FILE, f_distr_stats_ALL,
                options_intern.distribution_analysis_all_name, tp_tf_gene_count);

        sb_home_distribution_analysis.append(html_tail);

        BufferedWriter bw_html_home = new BufferedWriter(new FileWriter(html_home_distribution_analysis));
        bw_html_home.write(sb_home_distribution_analysis.toString());
        bw_html_home.close();
        logger.logLine("[DISTRIBUTION-ANALYSIS-HTML-REPORT] Finished HTML report.");
    }

    /**
     * CREATE python scripts for plots, to only accept TFs which have a higher mean TF-TG-SCORE than the background distribution
     */
    public void create_distribution_plots() throws Exception {
        logger.logLine(
                "[DISTRIBUTION-ANALYSIS-PLOTS] Start comparing background distribution to TF distribution and create plots for outstanding TFs");

        File f_analysis_distr_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution);

        File f_analysis_distr_HM_level_input = new File(f_analysis_distr_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores_background_distr + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores_background_distr_HM);

        File f_root_plots_scripts = new File(f_analysis_distr_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_plots_scripts);
        f_root_plots_scripts.mkdir();
        File f_root_plots_scripts_all = new File(f_root_plots_scripts.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_plots_script_ALL);
        f_root_plots_scripts_all.mkdir();
        File f_root_plots_scripts_hm = new File(f_root_plots_scripts.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_plots_scripts_HM);
        f_root_plots_scripts_hm.mkdir();

        File f_root_mwu_plots_scripts = new File(f_analysis_distr_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_mwu_scripts);
        f_root_mwu_plots_scripts.mkdir();
        File f_root_mwu_plots_scripts_ALL = new File(f_root_mwu_plots_scripts.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_plots_ALL);
        f_root_mwu_plots_scripts_ALL.mkdir();
        File f_root_mwu_plots_scripts_HM = new File(f_root_mwu_plots_scripts.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_plots_HM);
        f_root_mwu_plots_scripts_HM.mkdir();

        File f_root_muw_plots_plots = new File(f_analysis_distr_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_mwu_plots);
        f_root_muw_plots_plots.mkdir();
        File f_root_muw_plots_plots_ALL = new File(f_root_muw_plots_plots.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_plots_ALL);
        f_root_muw_plots_plots_ALL.mkdir();
        File f_root_muw_plots_plots_HM = new File(f_root_muw_plots_plots.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_plots_HM);
        f_root_muw_plots_plots_HM.mkdir();


        File f_root_plots = new File(f_analysis_distr_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_plots);
        f_root_plots.mkdir();
        File f_root_plots_all = new File(
                f_root_plots.getAbsolutePath() + File.separator + options_intern.folder_out_distribution_plots_ALL);
        f_root_plots_all.mkdir();
        File f_root_plots_hm = new File(
                f_root_plots.getAbsolutePath() + File.separator + options_intern.folder_out_distribution_plots_HM);
        f_root_plots_hm.mkdir();

        File f_scores_root = new File(f_analysis_distr_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores);

        File f_background_distr_root = new File(f_scores_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores_background_distr);
        File f_background_distr_tfs_ALL = new File(f_scores_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores_tf_distributions + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores_background_distr_ALL);
        File f_background_distr_tfs_HM = new File(f_scores_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores_background_distr + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores_background_distr_HM);

        File f_tf_distribution_ALL = new File(f_scores_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores_tf_distributions + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores_tf_distributions_ALL);
        File f_tf_distribution_HM = new File(f_scores_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores_tf_distributions + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores_tf_distributions_HM);

        File f_stats = new File(f_analysis_distr_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_stats);
        f_stats.mkdir();
        File f_stats_ALL =
                new File(f_stats.getAbsolutePath() + File.separator + options_intern.folder_out_distribution_stats_ALL);
        f_stats_ALL.mkdir();
        File f_stats_HM =
                new File(f_stats.getAbsolutePath() + File.separator + options_intern.folder_out_distribution_stats_HM);
        f_stats_HM.mkdir();

        ArrayList<MANN_WHITNEYU_PLOTS_FILES> mann_whitneyU_plots = new ArrayList<>();

        ArrayList<File> scripts_to_execute = new ArrayList<>();

        for (File fileDir_hm : f_analysis_distr_HM_level_input.listFiles()) {
            if (fileDir_hm.isDirectory()) {
                File f_output_script =
                        new File(f_root_plots_scripts_hm.getAbsolutePath() + File.separator + fileDir_hm.getName());
                f_output_script.mkdir();
                File f_output_script_file = new File(f_output_script.getAbsolutePath() + File.separator +
                        options_intern.file_suffix_distribution_analysis_python_script);

                File f_output_mwu_script =
                        new File(f_root_mwu_plots_scripts_HM.getAbsolutePath() + File.separator + fileDir_hm.getName());
                f_output_mwu_script.mkdir();

                File f_output_plot =
                        new File(f_root_plots_hm.getAbsolutePath() + File.separator + fileDir_hm.getName());
                f_output_plot.mkdir();

                File f_output_mwu_plot =
                        new File(f_root_muw_plots_plots_HM.getAbsolutePath() + File.separator + fileDir_hm.getName());
                f_output_mwu_plot.mkdir();

                File f_output_stats = new File(f_stats_HM.getAbsolutePath() + File.separator + fileDir_hm.getName());
                f_output_stats.mkdir();

                File f_input_background_distr = new File(
                        f_background_distr_tfs_HM.getAbsolutePath() + File.separator + fileDir_hm.getName() +
                                File.separator + options_intern.file_suffix_distribution_analysis_distributions);
                File f_input_tf_distr =
                        new File(f_tf_distribution_HM.getAbsolutePath() + File.separator + fileDir_hm.getName());

                scripts_to_execute.add(f_output_script_file);

                //File input_background_file, File input_tf_root, File output_plots, File output_script_file, File output_stats) throws IOException {

                write_python_script_distribution_analysis(f_input_background_distr, f_input_tf_distr, f_output_plot,
                        f_output_script_file, f_output_stats);


                MANN_WHITNEYU_PLOTS_FILES mwu = new MANN_WHITNEYU_PLOTS_FILES();

                mann_whitneyU_plots.add(mwu);
                mwu.background = f_input_background_distr;
                mwu.tf_data = f_input_tf_distr;
                mwu.hm = true;
                mwu.f_output_script = f_output_mwu_script;
                mwu.f_output_plot = f_output_mwu_plot;
                mwu.options_intern = options_intern;

            }
        }

        File all_background_distr = new File(f_background_distr_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores_background_distr_ALL + File.separator +
                options_intern.file_suffix_distribution_analysis_distributions);
        File f_script_out_all = new File(f_root_plots_scripts_all.getAbsolutePath() + File.separator +
                options_intern.file_suffix_distribution_analysis_python_script);

        write_python_script_distribution_analysis(all_background_distr, f_tf_distribution_ALL, f_root_plots_all,
                f_script_out_all, f_stats_ALL);
        scripts_to_execute.add(f_script_out_all);

        MANN_WHITNEYU_PLOTS_FILES mwu = new MANN_WHITNEYU_PLOTS_FILES();
        mwu.background = all_background_distr;
        mwu.tf_data = f_tf_distribution_ALL;
        mwu.hm = false;
        mwu.f_output_script = f_root_mwu_plots_scripts_ALL;
        mwu.f_output_plot = f_root_muw_plots_plots_ALL;
        mwu.options_intern = options_intern;

        mann_whitneyU_plots.add(mwu);

        for (File fileDir : scripts_to_execute) {
            String command_edited = "python3 " + fileDir.getAbsolutePath();
            logger.logLine("[DISTRIBUTION-ANALYSIS-PLOTS] Run python script: " + command_edited);


            Process child = Runtime.getRuntime().exec(command_edited);
            int code = child.waitFor();
            switch (code) {
                case 0:
                    break;
                case 1:
                    String message = child.getErrorStream().toString();
                    throw new Exception(message);
            }
        }

        //DO MANN WHITNEY U PLOTS

        for (MANN_WHITNEYU_PLOTS_FILES mwu_execute : mann_whitneyU_plots) {
            mwu_execute.write_and_execute_scripts();
        }

        logger.logLine(
                "[DISTRIBUTION-ANALYSIS-PLOTS] Finished comparing background distribution to TF distribution and create plots for outstanding TFs");
    }

    /**
     * creates data needed for distribution plots
     */
    public void perform_distribution_analysis() throws IOException {
        logger.logLine("[DISTRIBUTION-ANALYSIS] Calculate distributions for TFs");

        //search input directories and create output directories
        File f_distr_analysis = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution);
        File f_distr_analysis_analysed_tfs = new File(f_distr_analysis.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_analyzed_tfs);
        File f_distr_analysis_analysed_tfs_csv = new File(
                f_distr_analysis_analysed_tfs.getAbsolutePath() + File.separator +
                        options_intern.file_suffix_distribution_analysis_analysed_tfs);

        File f_out_distr_analysis_tf_tg_scores = new File(f_distr_analysis.getAbsolutePath() + File.separator +
                options_intern.folder_out_distribution_tf_tg_scores);
        f_out_distr_analysis_tf_tg_scores.mkdir();
        File f_out_distr_analysis_tf_tg_scores_background_distr = new File(
                f_out_distr_analysis_tf_tg_scores.getAbsolutePath() + File.separator +
                        options_intern.folder_out_distribution_tf_tg_scores_background_distr);
        f_out_distr_analysis_tf_tg_scores_background_distr.mkdir();
        File f_out_distr_analysis_tf_tg_scores_background_distr_ALL = new File(
                f_out_distr_analysis_tf_tg_scores_background_distr.getAbsolutePath() + File.separator +
                        options_intern.folder_out_distribution_tf_tg_scores_background_distr_ALL);
        f_out_distr_analysis_tf_tg_scores_background_distr_ALL.mkdir();
        File f_out_distr_analysis_tf_tg_scores_background_distr_HM = new File(
                f_out_distr_analysis_tf_tg_scores_background_distr.getAbsolutePath() + File.separator +
                        options_intern.folder_out_distribution_tf_tg_scores_background_distr_HM);
        f_out_distr_analysis_tf_tg_scores_background_distr_HM.mkdir();
        File f_out_distr_analysis_tf_tg_scores_tf_distr = new File(
                f_out_distr_analysis_tf_tg_scores.getAbsolutePath() + File.separator +
                        options_intern.folder_out_distribution_tf_tg_scores_tf_distributions);
        f_out_distr_analysis_tf_tg_scores_tf_distr.mkdir();
        File f_out_distr_analysis_tf_tg_scores_tf_distr_ALL = new File(
                f_out_distr_analysis_tf_tg_scores_tf_distr.getAbsolutePath() + File.separator +
                        options_intern.folder_out_distribution_tf_tg_scores_tf_distributions_ALL);
        f_out_distr_analysis_tf_tg_scores_tf_distr_ALL.mkdir();
        File f_out_distr_analysis_tf_tg_scores_tf_distr_HM = new File(
                f_out_distr_analysis_tf_tg_scores_tf_distr.getAbsolutePath() + File.separator +
                        options_intern.folder_out_distribution_tf_tg_scores_tf_distributions_HM);
        f_out_distr_analysis_tf_tg_scores_tf_distr_HM.mkdir();

        HashMap<String, HashSet<String>> distinct_hms_tfs = new HashMap<>();
        HashMap<String, HashSet<String>> distinct_tfs_hms = new HashMap<>();

        HashMap<String, String> ensg_symbol = new HashMap<>();
        HashMap<String, HashSet<String>> symbol_ensg = new HashMap<>();

        HashMap<String, HashMap<String, Double>> timepoint_ensg_gene_counts = new HashMap<>();
        HashMap<String, HashMap<String, Double>> group_clash_diff_gene_expression = new HashMap<>();

        HashMap<String, String> composed_tfs = new HashMap<>();
        HashMap<String, HashSet<String>> composed_tfs_tfs = new HashMap<>();

        BufferedReader br_composed_tfs = new BufferedReader(new FileReader(
                options_intern.com2pose_working_directory + File.separator +
                        options_intern.folder_name_tepic_postprocessing + File.separator +
                        options_intern.folder_name_tepic_postprocessing_tfs + File.separator +
                        options_intern.file_suffix_tepic_postprocessing_tfs_tfs));
        String line_composed_tfs = "";
        while ((line_composed_tfs = br_composed_tfs.readLine()) != null) {
            String[] split = line_composed_tfs.split("\t");

            HashSet temp_set = new HashSet();
            for (int i = 1; i < split.length; i++) {
                composed_tfs.put(split[i], split[0]);
                temp_set.add(split[i]);
            }
            composed_tfs_tfs.put(split[0], temp_set);
        }
        br_composed_tfs.close();

        BufferedReader br_distinct_hms_tfs = new BufferedReader(new FileReader(f_distr_analysis_analysed_tfs_csv));
        String line_distinct_hms_tfs = br_distinct_hms_tfs.readLine();
        while ((line_distinct_hms_tfs = br_distinct_hms_tfs.readLine()) != null) {
            String[] split = line_distinct_hms_tfs.split("\t");
            HashSet<String> current_hm_tfs;

            if (distinct_hms_tfs.containsKey(split[1])) {
                current_hm_tfs = distinct_hms_tfs.get(split[1]);
            } else {
                current_hm_tfs = new HashSet<>();

                File f_out_hms_background = new File(
                        f_out_distr_analysis_tf_tg_scores_background_distr_HM.getAbsolutePath() + File.separator +
                                split[1]);
                f_out_hms_background.mkdir();

                File f_out_hms_tf = new File(
                        f_out_distr_analysis_tf_tg_scores_tf_distr_HM.getAbsolutePath() + File.separator + split[1]);
                f_out_hms_tf.mkdir();

            }

            current_hm_tfs.add(split[0]);
            distinct_hms_tfs.put(split[1], current_hm_tfs);

            HashSet<String> current_tf_hms;
            if (distinct_tfs_hms.containsKey(split[0])) {
                current_tf_hms = distinct_tfs_hms.get(split[0]);
            } else {
                current_tf_hms = new HashSet<>();
            }

            current_tf_hms.add(split[1]);
            distinct_tfs_hms.put(split[0], current_tf_hms);
        }
        br_distinct_hms_tfs.close();

        BufferedReader br_ensg_symbol = new BufferedReader(new FileReader(new File(options_intern.tepic_ensg_symbol)));
        String line_ensg_symbol = br_ensg_symbol.readLine();
        while ((line_ensg_symbol = br_ensg_symbol.readLine()) != null) {
            String[] split = line_ensg_symbol.split("\t");
            if (split.length > 1) {
                String symbol = split[1].toUpperCase();

                if (composed_tfs.containsKey(symbol)) {
                    symbol = composed_tfs.get(symbol);
                }

                ensg_symbol.put(split[0].toUpperCase(), symbol);

                HashSet<String> current_ensgs;
                if (symbol_ensg.containsKey(symbol)) {
                    current_ensgs = symbol_ensg.get(symbol);
                } else {
                    current_ensgs = new HashSet<>();
                }
                current_ensgs.add(split[0].toUpperCase());
                symbol_ensg.put(symbol, current_ensgs);

            }
        }
        br_ensg_symbol.close();

        File f_genecounts_root = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.folder_name_deseq2_preprocessing_gene_symbols);
        for (File fileDir : f_genecounts_root.listFiles()) {
            if (fileDir.isFile()) {
                String name_tp = fileDir.getName().split("\\.")[0];

                HashMap<String, Double> ensg_genecount = new HashMap<>();

                BufferedReader br_gene_counts = new BufferedReader(new FileReader(fileDir));
                String line_gene_counts = br_gene_counts.readLine();
                while ((line_gene_counts = br_gene_counts.readLine()) != null) {
                    String[] split = line_gene_counts.split("\t");
                    ensg_genecount.put(split[1].toUpperCase(), Double.parseDouble(split[2]));
                }
                br_gene_counts.close();

                timepoint_ensg_gene_counts.put(name_tp, ensg_genecount);
            }
        }

        File f_diff_gene_expr_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_deseq2_output);
        for (File fileDir : f_diff_gene_expr_root.listFiles()) {
            if (fileDir.isFile()) {
                String[] split_name = fileDir.getName().split("\\.")[0].split("_");
                String name_group_clash = split_name[0] + "_" + split_name[1];

                HashMap<String, Double> current_diff_expr = new HashMap<>();

                BufferedReader br_diff_gene_expr = new BufferedReader(new FileReader(fileDir));
                String line_diff_gene_expr = br_diff_gene_expr.readLine();
                while ((line_diff_gene_expr = br_diff_gene_expr.readLine()) != null) {
                    String[] split = line_diff_gene_expr.split("\t");
                    current_diff_expr.put(split[0], Double.parseDouble(split[1]));
                }
                br_diff_gene_expr.close();

                group_clash_diff_gene_expression.put(name_group_clash, current_diff_expr);
            }
        }

        for (String key_tf : distinct_tfs_hms.keySet()) {
            String name_anhaengsel = "";
            if (symbol_ensg.containsKey(key_tf)) {
                HashSet<String> ensgs = symbol_ensg.get(key_tf);

                HashMap<String, Double> tp_gene_count = new HashMap<>();
                HashMap<String, Double> group_clash_diff_gene_expr = new HashMap<>();

                for (String k_tp : timepoint_ensg_gene_counts.keySet()) {
                    int gene_count_total = 0;
                    int found_gene_counts = 0;

                    HashMap<String, Double> tp_specific_genecounts = timepoint_ensg_gene_counts.get(k_tp);

                    for (String k_ensgs : ensgs) {
                        if (tp_specific_genecounts.containsKey(k_ensgs)) {
                            gene_count_total += tp_specific_genecounts.get(k_ensgs);
                            found_gene_counts++;
                        }
                    }

                    gene_count_total /= found_gene_counts;

                    tp_gene_count.put(k_tp, (double) gene_count_total);
                }

                for (String k_group_clash : group_clash_diff_gene_expression.keySet()) {

                    double diff_total = 0;
                    int diff_count = 0;

                    HashMap<String, Double> tp_specific_genecounts =
                            group_clash_diff_gene_expression.get(k_group_clash);

                    for (String k_ensgs : ensgs) {
                        if (tp_specific_genecounts.containsKey(k_ensgs)) {
                            diff_total += tp_specific_genecounts.get(k_ensgs);
                            diff_count++;
                        }
                    }

                    diff_total /= diff_count;

                    group_clash_diff_gene_expr.put(k_group_clash, diff_total);


                }

                for (String k_group_clash : group_clash_diff_gene_expr.keySet()) {
                    HashSet<String> available_hms = distinct_tfs_hms.get(key_tf);

                    for (String k_hm : available_hms) {
                        HashSet<String> available_ensgs = new HashSet<>();

                        boolean found_tgene_file = true;
                        if (!options_intern.path_tgen.equals("")) {
                            File f_tgene_input = new File(options_intern.com2pose_working_directory + File.separator +
                                    options_intern.folder_name_tgen + File.separator +
                                    options_intern.folder_name_tgen_filter_target_genes + File.separator + k_hm +
                                    File.separator + k_group_clash);

                            if (!f_tgene_input.exists()) {
                                found_tgene_file = false;
                            }
                            if (found_tgene_file) {
                                for (File fileDir : f_tgene_input.listFiles()) {
                                    if (fileDir.isFile()) {
                                        f_tgene_input = fileDir;
                                    }
                                }

                                BufferedReader br_tgene_input = new BufferedReader(new FileReader(f_tgene_input));
                                String line_tgene_input = br_tgene_input.readLine();
                                while ((line_tgene_input = br_tgene_input.readLine()) != null) {
                                    String[] split = line_tgene_input.split("\t");
                                    available_ensgs.add(split[0]);
                                }
                                br_tgene_input.close();
                            }

                        }

                        double tf_regression_coefficient = 0.0;

                        File f_hm_tf_coeff_input = new File(options_intern.com2pose_working_directory + File.separator +
                                options_intern.folder_out_put_DYNAMITE + File.separator + k_hm + File.separator +
                                k_group_clash + File.separator +
                                options_intern.file_suffix_dynamite_output_to_be_plotted);
                        if (f_hm_tf_coeff_input.exists() && f_hm_tf_coeff_input.isFile()) {
                            BufferedReader br_hm_tf_coeff = new BufferedReader(new FileReader(f_hm_tf_coeff_input));
                            String line_hm_tf_coeff = br_hm_tf_coeff.readLine();
                            while ((line_hm_tf_coeff = br_hm_tf_coeff.readLine()) != null) {
                                String[] split = line_hm_tf_coeff.split("\t");

                                String[] split_tf_long_name = split[0].split("_");

                                String tf_in_file = split_tf_long_name[0];

                                if (split_tf_long_name.length > 1) {
                                    name_anhaengsel = split_tf_long_name[1];
                                }

                                if (tf_in_file.toUpperCase().equals(key_tf.toUpperCase())) {
                                    if (split[1].equals("0")) {
                                        break;
                                    }
                                    tf_regression_coefficient = Double.parseDouble(split[1]);
                                    break;
                                }

                            }
                        } else {
                            continue;
                        }

                        if (tf_regression_coefficient == 0) {
                            continue;
                        }

                        String[] split_clash = k_group_clash.split("_");

                        double current_diff_gene_expr = group_clash_diff_gene_expr.get(k_group_clash);
                        double current_gene_counts =
                                tp_gene_count.get(split_clash[0]) + tp_gene_count.get(split_clash[1]);

                        double current_tf_score = 1;//current_gene_counts;
                        if (current_tf_score < 0) {
                            current_tf_score *= -1;
                        }

                        String[] group_clash_tps = k_group_clash.split("_");

                        File f_group1_hm_tf_target_genes_root;
                        File parent_group1;

                        if (options_intern.tepic_tpm_cutoff > 0) {
                            f_group1_hm_tf_target_genes_root = new File(
                                    options_intern.com2pose_working_directory + File.separator +
                                            options_intern.folder_name_tepic_postprocessing + File.separator +
                                            options_intern.folder_name_tepic_postprocessing_output + File.separator +
                                            k_hm + File.separator + k_group_clash + File.separator + split_clash[0]);
                            parent_group1 = new File(options_intern.com2pose_working_directory + File.separator +
                                    options_intern.folder_name_tepic_postprocessing + File.separator +
                                    options_intern.folder_name_tepic_postprocessing_output + File.separator + k_hm +
                                    File.separator + k_group_clash + File.separator + split_clash[0]);

                        } else {
                            parent_group1 = new File(options_intern.com2pose_working_directory + File.separator +
                                    options_intern.folder_name_tepic_output_raw + File.separator + group_clash_tps[0] +
                                    File.separator + k_hm + File.separator);
                            f_group1_hm_tf_target_genes_root = new File("");
                        }

                        String suffix = "";

                        if (options_intern.tepic_tpm_cutoff > 0) {
                            suffix = "_Gene_View_Filtered_TPM.txt";
                        } else {
                            suffix = "_Gene_View_Filtered.txt";
                        }


                        File f_group1_hm_tf_target_genes = new File("");

                        for (int j = 0; j < parent_group1.listFiles().length; j++) {
                            if (options_intern.tepic_tpm_cutoff > 0) {
                                j = parent_group1.listFiles().length;
                            } else {
                                f_group1_hm_tf_target_genes_root = parent_group1.listFiles()[j];
                            }

                            for (File fileDir : f_group1_hm_tf_target_genes_root.listFiles()) {
                                if (!fileDir.getName().matches(".*" + suffix + ".*")) {
                                    continue;
                                }
                                if (fileDir.isFile()) {
                                    f_group1_hm_tf_target_genes = fileDir;
                                    break;
                                }
                            }
                        }

                        File f_group2_hm_tf_target_genes_root;
                        File parent_group2;

                        if (options_intern.tepic_tpm_cutoff > 0) {
                            f_group2_hm_tf_target_genes_root = new File(
                                    options_intern.com2pose_working_directory + File.separator +
                                            options_intern.folder_name_tepic_postprocessing + File.separator +
                                            options_intern.folder_name_tepic_postprocessing_output + File.separator +
                                            k_hm + File.separator + k_group_clash + File.separator + split_clash[1]);
                            parent_group2 = new File(options_intern.com2pose_working_directory + File.separator +
                                    options_intern.folder_name_tepic_postprocessing + File.separator +
                                    options_intern.folder_name_tepic_postprocessing_output + File.separator + k_hm +
                                    File.separator + k_group_clash + File.separator + split_clash[1]);

                        } else {
                            parent_group2 = new File(options_intern.com2pose_working_directory + File.separator +
                                    options_intern.folder_name_tepic_output_raw + File.separator + group_clash_tps[1] +
                                    File.separator + k_hm + File.separator);
                            f_group2_hm_tf_target_genes_root = new File("");
                        }

                        File f_group2_hm_tf_target_genes = new File("");

                        for (int j = 0; j < parent_group2.listFiles().length; j++) {
                            if (options_intern.tepic_tpm_cutoff > 0) {

                                j = parent_group2.listFiles().length;
                            } else {
                                f_group2_hm_tf_target_genes_root = parent_group2.listFiles()[j];
                            }

                            for (File fileDir : f_group2_hm_tf_target_genes_root.listFiles()) {
                                if (!fileDir.getName().matches(".*" + suffix + ".*")) {
                                    continue;
                                }
                                if (fileDir.isFile()) {
                                    f_group2_hm_tf_target_genes = fileDir;
                                    break;
                                }
                            }
                        }

                        if (f_group1_hm_tf_target_genes.exists() && f_group1_hm_tf_target_genes.isFile() &&
                                f_group2_hm_tf_target_genes.exists() && f_group2_hm_tf_target_genes.isFile()) {

                            Map<String, Double> target_genes_scores_group1 = new HashMap<>();
                            Map<String, Double> target_genes_scores_group2 = new HashMap<>();

                            BufferedReader br_group1_hm_tf_target_genes =
                                    new BufferedReader(new FileReader(f_group1_hm_tf_target_genes));
                            String line_group1_hm_tf_target_genes = br_group1_hm_tf_target_genes.readLine();

                            String temp_key_tf = key_tf.replace(".", ":");
                            String temp_key_tf_2 = key_tf.replace(".", ":");
                            if (!name_anhaengsel.equals("")) {
                                temp_key_tf += "_" + name_anhaengsel;
                            }

                            int index_tf_group1_hm_tf_target_genes = -1;
                            String[] split_header_group1_hm_tf_target_gens = line_group1_hm_tf_target_genes.split("\t");
                            for (int i = 0; i < split_header_group1_hm_tf_target_gens.length; i++) {
                                if (split_header_group1_hm_tf_target_gens[i].toUpperCase().equals(temp_key_tf) ||
                                        split_header_group1_hm_tf_target_gens[i].toUpperCase().equals(temp_key_tf_2)) {
                                    index_tf_group1_hm_tf_target_genes = i;
                                    break;
                                }
                            }

                            while ((line_group1_hm_tf_target_genes = br_group1_hm_tf_target_genes.readLine()) != null) {
                                String[] split = line_group1_hm_tf_target_genes.split("\t");
                                target_genes_scores_group1.put(split[0],
                                        Double.parseDouble(split[index_tf_group1_hm_tf_target_genes]));
                            }
                            br_group1_hm_tf_target_genes.close();

                            BufferedReader br_group2_hm_tf_target_genes =
                                    new BufferedReader(new FileReader(f_group2_hm_tf_target_genes));
                            String line_group2_hm_tf_target_genes = br_group2_hm_tf_target_genes.readLine();

                            int index_tf_group2_hm_tf_target_genes = -1;
                            String[] split_header_group2_hm_tf_target_gens = line_group2_hm_tf_target_genes.split("\t");
                            for (int i = 0; i < split_header_group2_hm_tf_target_gens.length; i++) {
                                if (split_header_group2_hm_tf_target_gens[i].toUpperCase().equals(temp_key_tf) ||
                                        split_header_group2_hm_tf_target_gens[i].toUpperCase().equals(temp_key_tf_2)) {
                                    index_tf_group2_hm_tf_target_genes = i;
                                    break;
                                }
                            }

                            while ((line_group2_hm_tf_target_genes = br_group2_hm_tf_target_genes.readLine()) != null) {
                                String[] split = line_group2_hm_tf_target_genes.split("\t");
                                target_genes_scores_group2.put(split[0],
                                        Double.parseDouble(split[index_tf_group2_hm_tf_target_genes]));
                            }
                            br_group2_hm_tf_target_genes.close();


                            BufferedWriter bw_background_ALL;
                            File f_out_background_ALL = new File(
                                    f_out_distr_analysis_tf_tg_scores_background_distr_ALL.getAbsolutePath() +
                                            File.separator +
                                            options_intern.file_suffix_distribution_analysis_distributions);
                            if (f_out_background_ALL.exists()) {
                                bw_background_ALL = new BufferedWriter(new FileWriter(f_out_background_ALL, true));
                            } else {
                                bw_background_ALL = new BufferedWriter(new FileWriter(f_out_background_ALL));
                                bw_background_ALL.write("##HM\tALL\n");
                                bw_background_ALL.write("##TF\tALL\n");
                                bw_background_ALL.write("TARGET_GENE\tTF_TG_SCORE\tHM\tGROUPS\tTF\tTF_COEFF\n");

                            }

                            File f_out_background_HM = new File(
                                    f_out_distr_analysis_tf_tg_scores_background_distr_HM.getAbsolutePath() +
                                            File.separator + k_hm + File.separator +
                                            options_intern.file_suffix_distribution_analysis_distributions);
                            BufferedWriter bw_background_HM;
                            if (f_out_background_HM.exists()) {
                                bw_background_HM = new BufferedWriter(new FileWriter(f_out_background_HM, true));
                            } else {
                                bw_background_HM = new BufferedWriter(new FileWriter(f_out_background_HM));
                                bw_background_HM.write("##HM\t" + k_hm + "\n");
                                bw_background_HM.write("##TF\tALL\n");
                                bw_background_HM.write("TARGET_GENE\tTF_TG_SCORE\tHM\tGROUPS\tTF\tTF_COEFF\n");
                            }

                            File f_out_tf_ALL = new File(
                                    f_out_distr_analysis_tf_tg_scores_tf_distr_ALL.getAbsolutePath() + File.separator +
                                            key_tf + "_" +
                                            options_intern.file_suffix_distribution_analysis_distributions);
                            BufferedWriter bw_tf_ALL;
                            if (f_out_tf_ALL.exists()) {
                                bw_tf_ALL = new BufferedWriter(new FileWriter(f_out_tf_ALL, true));
                            } else {
                                bw_tf_ALL = new BufferedWriter(new FileWriter(f_out_tf_ALL));
                                bw_tf_ALL.write("##HM\tALL\n");
                                bw_tf_ALL.write("##TF\t" + key_tf + "\n");
                                bw_tf_ALL.write("TARGET_GENE\tTF_TG_SCORE\tHM\tGROUPS\tTF\tTF_COEFF\n");
                            }

                            BufferedWriter bw_tf_HM;
                            File f_out_tf_HM = new File(
                                    f_out_distr_analysis_tf_tg_scores_tf_distr_HM.getAbsolutePath() + File.separator +
                                            k_hm + File.separator + key_tf + "_" +
                                            options_intern.file_suffix_distribution_analysis_distributions);
                            if (f_out_tf_HM.exists()) {
                                bw_tf_HM = new BufferedWriter(new FileWriter(f_out_tf_HM, true));
                            } else {
                                bw_tf_HM = new BufferedWriter(new FileWriter(f_out_tf_HM));
                                bw_tf_HM.write("##HM\t" + k_hm + "\n");
                                bw_tf_HM.write("##TF\t" + key_tf + "\n");
                                bw_tf_HM.write("TARGET_GENE\tTF_TG_SCORE\tHM\tGROUPS\tTF\tTF_COEFF\n");
                            }

                            for (String k_target_gene : target_genes_scores_group1.keySet()) {
                                if (!options_intern.path_tgen.equals("") && found_tgene_file) {
                                    if (!available_ensgs.contains(k_target_gene)) {
                                        continue;
                                    }
                                }

                                if (target_genes_scores_group2.containsKey(k_target_gene)) {
                                    double gene_score1 = target_genes_scores_group1.get(k_target_gene);
                                    double gene_score2 = target_genes_scores_group2.get(k_target_gene);

                                    double genecount_1 = 0;
                                    double genecount_2 = 0;

                                    double ensg_diff_gene_expr = 0;
                                    if (timepoint_ensg_gene_counts.get(split_clash[0]).containsKey(k_target_gene)) {
                                        genecount_1 = timepoint_ensg_gene_counts.get(split_clash[0]).get(k_target_gene);
                                    }

                                    if (timepoint_ensg_gene_counts.get(split_clash[1]).containsKey(k_target_gene)) {
                                        genecount_2 = timepoint_ensg_gene_counts.get(split_clash[1]).get(k_target_gene);
                                    }

                                    if (group_clash_diff_gene_expr.containsKey(k_group_clash)) {
                                        ensg_diff_gene_expr = group_clash_diff_gene_expr.get(k_group_clash);
                                    }

                                    double tg_score_1 = 0;
                                    double tg_score_2 = 0;

                                    if (options_intern.plot_distribution_analysis_score_type.equals(
                                            "EXCL_GENE_COUNTS")) {
                                        tg_score_1 = ensg_diff_gene_expr * gene_score1;
                                        tg_score_2 = ensg_diff_gene_expr * gene_score2;
                                    }
                                    if (options_intern.plot_distribution_analysis_score_type.equals("GENE_COUNTS")) {
                                        tg_score_1 = genecount_1 * ensg_diff_gene_expr * gene_score1;
                                        tg_score_2 = genecount_2 * ensg_diff_gene_expr * gene_score2;
                                    }


                                    if (tg_score_1 < 0) {
                                        tg_score_1 *= -1;
                                    }
                                    if (tg_score_2 < 0) {
                                        tg_score_2 *= -1;
                                    }

                                    double tg_score_cumm = tg_score_1 + tg_score_2;

                                    double tf_tg_score = (current_tf_score * tg_score_cumm) * tf_regression_coefficient;

                                    if (tf_tg_score < 0) {
                                        tf_tg_score *= -1;
                                    }

                                    StringBuilder sb_line = new StringBuilder();
                                    sb_line.append(k_target_gene);
                                    sb_line.append("\t");
                                    sb_line.append(tf_tg_score);
                                    sb_line.append("\t");
                                    sb_line.append(k_hm);
                                    sb_line.append("\t");
                                    sb_line.append(k_group_clash);
                                    sb_line.append("\t");
                                    sb_line.append(key_tf);
                                    sb_line.append("\t");
                                    sb_line.append(tf_regression_coefficient);

                                    bw_background_HM.write(sb_line.toString());
                                    bw_background_HM.newLine();

                                    bw_background_ALL.write(sb_line.toString());
                                    bw_background_ALL.newLine();

                                    bw_tf_ALL.write(sb_line.toString());
                                    bw_tf_ALL.newLine();

                                    bw_tf_HM.write(sb_line.toString());
                                    bw_tf_HM.newLine();
                                }
                            }
                            bw_background_ALL.close();
                            bw_background_HM.close();
                            bw_tf_ALL.close();
                            bw_tf_HM.close();
                        }
                    }
                }
            } else {
                continue;
            }
        }


        logger.logLine("[DISTRIBUTION-ANALYSIS] Finished calculating distributions");
    }

    /**
     * creates the HTML report
     * it creates also one important file (available tfs in HMs and which stages) for distribution analysis
     */
    public void create_overview_html_report_before_distribution_analysis() throws Exception {
        logger.logLine("[WEBSITE] Start creating overview website.");

        File f_website_css =
                new File(options_intern.path_to_COM2POSE + File.separator + "ext" + File.separator + "WEBSITE");

        File f_move_css = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_website +
                        File.separator + options_intern.folder_out_website_basics);
        f_move_css.mkdir();

        File f_output_website = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_website);
        f_output_website.mkdir();

        File f_output_website_htmls = new File(f_output_website.getAbsolutePath() + File.separator +
                options_intern.folder_out_website_htmls_regression_coefficients);
        f_output_website_htmls.mkdir();

        String command = "cp -u -r " + f_website_css.getAbsolutePath() + " " + f_move_css.getAbsolutePath();
        Process child = Runtime.getRuntime().exec(command);
        logger.logLine("[WEBSITE] Copy CSS files: " + command);
        int code = child.waitFor();
        switch (code) {
            case 0:
                break;
            case 1:
                String message = child.getErrorStream().toString();
                throw new Exception(message);
        }

        for (int i_dont_know = 0; i_dont_know < 2; i_dont_know++) {

            HashMap<String, HashMap<String, HashSet<String>>> distinct_tf_hm_diff_same = new HashMap<>();

            String html_tail = "</body>\n" + "</html>";

            BufferedWriter bw_home = new BufferedWriter(new FileWriter(f_output_website + File.separator +
                    options_intern.html_report_home_regression_coefficient_analysis));

            //bw_home.write(sb_home_front.toString());
            bw_home.write(get_header_html(options_intern.html_report_levels_home,
                    options_intern.analysis_types_regression_coefficient_analysis));

            threshold_folders_filled = true;

            for (Double d : options_intern.plot_th_coefficient) {
                bw_home.append(write_regression_coeffecient_analysis_found_table_html(d, "HOME"));
            }

            bw_home.write(html_tail);

            bw_home.close();


            StringBuilder sb_parameter = new StringBuilder(get_header_html(options_intern.html_report_levels_home, ""));

            sb_parameter.append(" <script>\n" + " document.title = \"PARAMETERS\";\n" + " </script>\n");


            //mix_option_parameters
            {
                sb_parameter.append(
                        "<button class=\"button_expandable\" id=\"button_preprocessing_mix_options\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_preprocessing_mix_options','table_preprocessing_mix_options')\"> Preprocessing Mix Options\n");
                sb_parameter.append(
                        "<div style=\"display: none;background-color: white;color:black;\" id=\"table_preprocessing_mix_options\">\n");

                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");
                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>mix_level</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.mix_level + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: if set a mix of either the Histone Modification Level (HM_LEVEL) or the Sample Level (SAMPLE_LEVEL) will be performed, mix_option is required\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>mix_option</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.mix_option + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: set histone marks or samples will be mixed with option: UNION (all peaks of all HMs will be used), INTERSECTION (only peaks, which are in all HMs will be used)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                if (options_intern.mix_option.equals("INTERSECTION")) {
                    sb_parameter.append("\t\t\t<tr>");
                    sb_parameter.append("\t\t\t\t<th>mix_occurence_intersection</th>");
                    sb_parameter.append("\t\t\t\t<th>" + options_intern.mix_occurence_intersection + "</th>");
                    sb_parameter.append(
                            "\t\t\t\t<th>#[OPT]: minimal occurence of peaks in intersection (only applied if mix_option is set to INTERSECTION), if set to zero it means it must be in all samples /histone modifications available, default 2\n</th>");
                    sb_parameter.append("\t\t\t</tr>");
                }
                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>mix_mutually_exclusive</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.mix_mutually_exclusive + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: mutually exclusive peaks. Use only peaks, which are mutually exclusive, default: TRUE\n</th>");
                sb_parameter.append("\t\t\t</tr>");
            }

            //blacklist parameters
            {
                sb_parameter.append(
                        "<button class=\"button_expandable\" id=\"button_blacklist_parameters\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_blacklist_parameters','table_blacklist_parameters')\"> Blacklist Parameters\n");
                sb_parameter.append(
                        "<div style=\"display: none;background-color: white;color:black;\" id=\"table_blacklist_parameters\">\n");
                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>black_list_dir</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.black_list_dir + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: if blacklist filter should be used provide path to BED file here (BED files can be found at https://github.com/Boyle-Lab/Blacklist/tree/master/lists)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>black_list_signals</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.black_list_signals.toString() + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: which regions should be ignored: Low_Mappability;High_Signal_Region\n</th>");
                sb_parameter.append("\t\t\t</tr>");


                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");
            }

            //DESeq2 parameters
            {
                sb_parameter.append(
                        "<button class=\"button_expandable\" id=\"button_deseq2_parameters\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_deseq2_parameters','table_deseq2_parameters')\"> DESeq2 Parameters\n");
                sb_parameter.append(
                        "<div style=\"display: none;background-color: white;color:black;\" id=\"table_deseq2_parameters\">\n");

                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>deseq2_input_directory</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.deseq2_input_directory + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[REQ]: count files in nfcore RNA-seq format (each line is a count of one gene), directory must be ordered like: TP1 - samples_1,...,samples_n;....\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>deseq2_input_gene_id</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.deseq2_input_gene_id + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[REQ]: gene ID file from nfcore RNA-seq (each line names one gene ID (ENSG) - must be same order as deseq2_input_directory files\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>deseq2_biomart_dataset_species</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.deseq2_biomart_dataset_species + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[REQ]: biomart dataset name for the species which is used in RNA-seq and CHIP-seq data, (https://www.bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>deseq2_count_threshold</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.deseq2_count_threshold + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: minimum count over all samples of two timepoints for DESeq2, default: 0\n</th>");
                sb_parameter.append("\t\t\t</tr>");
                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");
            }

            //TEPIC parameters
            {
                sb_parameter.append(
                        "<button class=\"button_expandable\" id=\"button_tepic_parameters\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_tepic_parameters','table_tepic_parameters')\"> TEPIC Parameters\n");
                sb_parameter.append(
                        "<div style=\"display: none;background-color: white;color:black;\" id=\"table_tepic_parameters\">\n");

                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_input_directory [ORIGINAL]</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_input_original + "</th>");
                sb_parameter.append("\t\t\t\t<th></th>");
                sb_parameter.append("\t\t\t</tr>");

                if (!options_intern.tepic_input_prev.equals("")) {
                    sb_parameter.append("\t\t\t<tr>");
                    sb_parameter.append("\t\t\t\t<th>tepic_input_directory [LAST FILTER]</th>");
                    sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_input_directory + "</th>");
                    sb_parameter.append("\t\t\t\t<th></th>");
                    sb_parameter.append("\t\t\t</tr>");
                }

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_input_ref_genome</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_input_ref_genome + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[REQ]:fasta files (in RefSeq format, without \\\"chr\\\" prefix)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_path_pwms</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_path_pwms + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[REQ]: path to position specific energy matrix used for TRAP (different matrices can be found in ~/COM2POSE/ext/TEPIC/TEPIC/PWMs/2.1)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_cores</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_cores + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: number of cores\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_bed_chr_sign</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_bed_chr_sign + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: bedgraph file containing open chromatin signal, e.g. DNase1-seq, or Histone-Mark signal\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_column_bedfile</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_column_bedfile + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: column in the tepic_input_directory file containing the average per base signal within a peak. If this option is used, the tepic_bed_chr_sign option must not be used\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_gene_annot</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_gene_annot + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[REQ]: gene annotation file, required to generate the gene view, required for TPM filter\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_window_size</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_window_size + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: size of the window to be considered to generate gene view (default 50kb)]\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_onlyDNasePeaks</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_onlyDNasePeaks + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: path to annotation file and annotate only DNase peaks that are within a window specified by the tepic_window_size option around all genes contained in the gene annotation file specified by this option\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_exponential_decay</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_exponential_decay + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: indicating whether exponential decay should be used (default TRUE)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_not_norm_peak_length</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_not_norm_peak_length + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: flag to be set if affinities should not be normalised by peak length\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_not_generated</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_not_generated + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: flag to be set if peak features for peak length and peak counts should not be generated\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_original_decay</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_original_decay + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: if tepic_bed_chr_sign or tepic_column_bedfile is used together with this flag, the original (Decay-)Scaling formulation of TEPIC is used to compute gene-TF scores\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_psems_length</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_psems_length + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: path to a tab delimited file containing the length of the used PSEMs\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_entire_gene_body</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_entire_gene_body + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: flag to be set if the entire gene body should be screened for TF binding. The search window is extended by a region half of the size that is specified by the tepic_window_size option upstream of the genes 5' TSS\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_zipped</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_zipped + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: flag indicating that the output of TEPIC should be zipped\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_background_seq</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_background_seq + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: path to a set of background sequences that should be used to compute to generate a binary score for TF binding. Mutually exclusive to the tepic_2bit option\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_2bit</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_2bit + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: path to a 2bit representation of the reference genome, required to generate a binary score for TF binding. The binary score is generated in addition to the standard affinity values. Mutually exclusive to the tepic_background_seq option\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_pvalue</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_pvalue + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: p-value cut off used to determine a cut off to derive a binary score for TF binding (default 0.05)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_minutes_per_chr</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_minutes_per_chr + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: minutes that should be spend at most per chromosome to find matching random regions (default 3)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_chr_prefix</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_chr_prefix + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: flag indicating that the reference genome contains a chr prefix\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_transcript_based</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_transcript_based + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: flag indicating that the annotation should be transcript and not gene based\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_loop_list</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_loop_list + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: loop list file containing chromatin contacts\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_loop_windows</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_loop_windows + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: size of the loop window used around a genes promoter to link chromatin loops to genes (default 5000)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_only_peak_features</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_only_peak_features + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: parameter to be set if only peak features should be computed (default FALSE)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_tpm_cutoff</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_tpm_cutoff + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: set (T)ranscripts (P)er (M)illion cutoff, default: TPM filter not active\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_ensg_symbol</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_ensg_symbol + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[REQ]: IF NOT SET MAPPING WILL BE PERFORMED AUTOMATICALLY - path to input file ensg to gene symbol file\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tepic_tgene_target_genes</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tepic_tgene_target_genes + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: use only tgene linked target genes default: true\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");

            }

            //TGEN parameters
            if (!options_intern.path_tgen.equals("")) {
                sb_parameter.append(
                        "<button class=\"button_expandable\" id=\"button_tgen_parameters\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_tgen_parameters','table_tgen_parameters')\"> TGene Parameters\n");
                sb_parameter.append(
                        "<div style=\"display: none;background-color: white;color:black;\" id=\"table_tgen_parameters\">\n");

                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tgen_no_closest_locus</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tgen_no_closest_locus + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: if no locus is found within window size, the nearest locus is used, default:false, meaning locus is used\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tgen_no_closest_tss</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tgen_no_closest_tss + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: if no tss is found within window size, the nearest locus is used, default:false, meaning locus is used\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tgen_max_link_distances</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tgen_max_link_distances + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: window size of tgen, default: 500000\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tgen_pvalue</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tgen_pvalue + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: max pvalue, which is accepted, default 0.05\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>tgen_self_regulatory</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.tgen_self_regulatory + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[REQ]: should self regulatory TFs {OPT-SELF-REG} be increased? default: false\n</th>");
                sb_parameter.append("\t\t\t</tr>");
                if (options_intern.tgen_self_regulatory) {
                    sb_parameter.append("\t\t\t<tr>");
                    sb_parameter.append("\t\t\t\t<th>tgen_consensus</th>");
                    sb_parameter.append("\t\t\t\t<th>" + options_intern.tgen_consensus + "</th>");
                    sb_parameter.append(
                            "\t\t\t\t<th>#[OPT]: {OPT-SELF-REG} value for consus approach (e.g. 0.5 = 50:50 TGene:TEPIC, 0.4 = 40:60 TGene:TEPIC), default: 0.5\n</th>");
                    sb_parameter.append("\t\t\t</tr>");

                    sb_parameter.append("\t\t\t<tr>");
                    sb_parameter.append("\t\t\t\t<th>tgen_consensus_calc</th>");
                    sb_parameter.append("\t\t\t\t<th>" + options_intern.tgen_consensus_calc + "</th>");
                    sb_parameter.append(
                            "\t\t\t\t<th>#[OPT]: {OPT-SELF-REG} can be \"INCREASE_TGENE_TFS\" (increases TEPIC TF affinities for TFs found in TGENE at target gene) or \"DECREASE_NOT_TGENE_TFs (decreases TEPIC TF affinities for TFs not found in TGENE at target gene)\".\n</th>");
                    sb_parameter.append("\t\t\t</tr>");

                    sb_parameter.append("\t\t\t<tr>");
                    sb_parameter.append("\t\t\t\t<th>tgen_mt_writing</th>");
                    sb_parameter.append("\t\t\t\t<th>" + options_intern.tgen_mt_writing + "</th>");
                    sb_parameter.append(
                            "\t\t\t\t<th>#[OPT]: {OPT-SELF-REG} if TGENE consensus is used please specify the writing of the Mitochondrial DNA chromosome in Peak Files, default: MT\n</th>");
                    sb_parameter.append("\t\t\t</tr>");

                }
                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");
            }

            //DYNAMITE parameters
            {
                sb_parameter.append(
                        "<button class=\"button_expandable\" id=\"button_dynamite_parameters\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_dynamite_parameters','table_dynamite_parameters')\"> DYNAMITE parameters\n");
                sb_parameter.append(
                        "<div style=\"display: none;background-color: white;color:black;\" id=\"table_dynamite_parameters\">\n");

                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_preprocessing_integrate_data_consider_geneFile</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>" + options_intern.dynamite_preprocessing_integrate_data_consider_geneFile +
                                "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: File containing gene IDs that should be considered\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_out_var</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.dynamite_out_var + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[REQ]: Name of the response variable default: Expression (in this pipeline)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_cores</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.dynamite_cores + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: Number of the cores to use (1 as default)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_alpha</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.dynamite_alpha + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: Alpha parameter stepsize (0.1 as default)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_testsize</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.dynamite_testsize + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]:Size of test data (0.2 as default)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_Ofolds</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.dynamite_Ofolds + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: Number of outer folds for model validation (3 as default)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_Ifolds</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.dynamite_Ifolds + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: Number of inner cross validation folds (6 as default)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_balanced</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.dynamite_balanced + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: Flag indicating whether the data should be balanced through downsampling (default TRUE)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_performance</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.dynamite_performance + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: Flag indicating whether performance measures should be computed (default TRUE)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>dynamite_randomise</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.dynamite_randomise + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: Flag indicating whether a model should be learned on randomised data (default FALSE)\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");

            }

            //PLOT parameters
            {
                sb_parameter.append(
                        "<button class=\"button_expandable\" id=\"button_plot_parameters\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_plot_parameters','table_plot_parameters')\"> PLOT Parameters\n");
                sb_parameter.append(
                        "<div style=\"display: none;background-color: white;color:black;\" id=\"table_plot_parameters\">\n");

                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>plot_th_coefficient</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.plot_th_coefficient.toString() + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: thresholds for coefficent plot and overall plots, for each threshold it creates a plot for each timepoint and an overall plot for each HistoneMod\n" +
                                "#e.g. 0.1 means it creates a plot of the coefficient range ([-1;1]), it uses all TFs of [-1,-0.1] and [0.1,1]\n" +
                                "#default: 0.1;0.2;0.3;0.4;0.5;0.6</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>plot_cutoff_tps</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.plot_cutoff_tps + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: in how many TPs should a TF be found, default: 2\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>plot_cutoff_hms</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.plot_cutoff_hms + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: in how many HMs should a TF be found, default: 2\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>plot_cutoff_gcs</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.plot_cutoff_gcs + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: minimum gene counts, default: 100\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>plot_top_k_genes</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.plot_top_k_genes + "</th>");
                sb_parameter.append("\t\t\t\t<th>#[OPT]: top k target genes for TFs, default: 30\n</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");

            }

            //HTML Report PARAMETERS
            {
                sb_parameter.append(
                        "<button class=\"button_expandable\" id=\"button_html_report_parameters\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_html_report_parameters','table_html_report_parameters')\"> HTML Report Parameters\n");
                sb_parameter.append(
                        "<div style=\"display: none;background-color: white;color:black;\" id=\"table_html_report_parameters\">\n");

                sb_parameter.append("\t\t<table style=\"width:100%;font-size:15px;text-align: left;\">\n");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>Parameter</th>");
                sb_parameter.append("\t\t\t\t<th>Value</th>");
                sb_parameter.append("\t\t\t\t<th>Explanation</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t\t<tr>");
                sb_parameter.append("\t\t\t\t<th>html_report_interesting_tfs</th>");
                sb_parameter.append("\t\t\t\t<th>" + options_intern.website_interesting_tfs + "</th>");
                sb_parameter.append(
                        "\t\t\t\t<th>#[OPT]: list of TFs which you are interested in - is only used to search fast for known TFs in the results, it does not affect results\n" +
                                "#e.g. website_interesting_tfs=\"STAT3;GATA3;NFIB\"</th>");
                sb_parameter.append("\t\t\t</tr>");

                sb_parameter.append("\t\t</table>");
                sb_parameter.append("</div>\n</button>\n");
            }


            /*
            sb_parameter.append("\t\t\t<tr>");
            sb_parameter.append("\t\t\t\t<th></th>");
            sb_parameter.append("\t\t\t\t<th></th>");
            sb_parameter.append("\t\t\t\t<th></th>");
            sb_parameter.append("\t\t\t</tr>");
             */


            sb_parameter.append(html_tail);

            BufferedWriter bw_parameter =
                    new BufferedWriter(new FileWriter(new File(f_output_website + File.separator + "PARAMETERS.html")));
            bw_parameter.write(sb_parameter.toString());
            bw_parameter.close();


            HashMap<String, ArrayList<String>> tf_gene_count = new HashMap<>();

            for (File fileDir : threshold_folders) {
                HashSet<String> tfs_to_create_pages = new HashSet<>();
                HashSet<String> possible_hms = new HashSet<>();

                File f_interactive_plots_root = new File(f_output_website.getAbsolutePath() + File.separator +
                        options_intern.folder_out_website_interactive_plots + File.separator + fileDir.getName());

                StringBuilder sb_threshold = new StringBuilder(
                        get_header_html(options_intern.html_report_levels_3_steps,
                                options_intern.analysis_types_regression_coefficient_analysis));
                sb_threshold.append(
                        " <script>\n" + " document.title = \"TH: " + fileDir.getName() + "\";\n" + " </script>\n");
                sb_threshold.append(
                        "<div class=\"w3-row-padding w3-padding-64 w3-container\"><div class=\"w3-content\"><h1>Current threshold: " +
                                fileDir.getName() + "</h1></div>\n</div>\n");


                HashSet<String> total_number_tfs = new HashSet<>();

                for (File fileDir_hm : f_interactive_plots_root.listFiles()) {
                    HashSet<String> total_number_tfs_hm = new HashSet<>();

                    possible_hms.add(fileDir_hm.getName());

                    sb_threshold.append("<div class=\"w3-row-padding w3-padding-64 w3-container\">\n" +
                            "  <div class=\"w3-content\">\n" + "    <div class=\"w3-twothird\">\n" + "      <h1>");
                    sb_threshold.append(fileDir_hm.getName());
                    sb_threshold.append(
                            "</h1>\n" + "\t  \n" + "\t <h4> Different timepoints / conditions </h4> \n" + "\t  \n");
                    sb_threshold.append("</div>\n");
                    sb_threshold.append("</div>\n");
                    sb_threshold.append(
                            "<div class=\"container\" style=\"\n" + "    width: 100%;\n" + "    position: relative;\n" +
                                    "    display: block;\n" + "    overflow-x: scroll;overflow-y: scroll;\">\n" +
                                    "\t  <iframe id=\"igraph_" + fileDir_hm.getName() + "_threshold_" +
                                    fileDir.getName() + "_different_stages.html" +
                                    "\" scrolling=\"no\" style=\"border:none;\" seamless=\"seamless\" src=\"");
                    String relative_path_diff_th = ".." + File.separator + ".." + File.separator +
                            options_intern.folder_out_website_interactive_plots + File.separator + fileDir.getName() +
                            File.separator + fileDir_hm.getName() + File.separator +
                            options_intern.folder_out_website_interactive_plots_overview + File.separator +
                            fileDir_hm.getName() + "_threshold_" + fileDir.getName() + "_different_stages.html";

                    sb_threshold.append(relative_path_diff_th);
                    sb_threshold.append("\" height=\"400\" width=\"2500\" overflow=\"scroll\"></iframe>\n");
                    sb_threshold.append("</div>\n");

                    sb_threshold.append("    <div class=\"w3-twothird\">\n");
                    sb_threshold.append("  <div class=\"w3-content\">\n");

                    //gene count table for plot

                    //collapse thing
                    sb_threshold.append("<button class=\"button_expandable\" id=\"button_" + fileDir_hm.getName() +
                            "_diff\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_" +
                            fileDir_hm.getName() + "_diff','table_" + fileDir_hm.getName() + "_diff')\"> GeneCounts\n");
                    //tooltip
                    sb_threshold.append(
                            "<div class=\"tooltip\"><img src=\".." + File.separator + ".." + File.separator +
                                    options_intern.folder_out_website_basics + File.separator +
                                    options_intern.folder_out_website_basics_website + File.separator +
                                    options_intern.folder_out_website_basics_website_images + File.separator +
                                    "information.png" + "\" style=\"width:45px;height:40px;\"/>" +
                                    "  <span class=\"tooltiptext\">GeneCount is the unprocessed data counted by nfcore RNA-seq. A TF which is at least expressed with GeneCount 1000 is considered highly expressed and active.</span>\n" +
                                    "</div>");
                    sb_threshold.append(
                            "<div style=\"display: none;background-color: white;color:black;\" id=\"table_" +
                                    fileDir_hm.getName() + "_diff\">\n");

                    HashSet<String> total_numbers_tfs_hm_diff = new HashSet<>();
                    sb_threshold.append("<h4>Gene Count threshold: " + options_intern.plot_cutoff_gcs + "</h4>");
                    sb_threshold.append(
                            "<h5><i>Click on TF for detailed information - if no Button is available it means that this TF was eliminated by a filter.</i></h5>\n");
                    sb_threshold.append("\t\t<table style=\"width:100%;font-size:15px;\">\n");

                    File f_gene_counts_input_different = new File(
                            options_intern.com2pose_working_directory + File.separator +
                                    options_intern.folder_out_analysis_data + File.separator +
                                    options_intern.folder_out_analysis_data_TP_LEVEL + File.separator +
                                    fileDir_hm.getName() + File.separator + fileDir.getName() + File.separator +
                                    options_intern.file_suffix_analysis_plot_data_hm_level_different);
                    BufferedReader br_gc_different = new BufferedReader(new FileReader(f_gene_counts_input_different));
                    String line_gc_different = "";
                    while ((line_gc_different = br_gc_different.readLine()) != null) {
                        String[] split = line_gc_different.split("\t");
                        sb_threshold.append("\t\t\t<tr>\n");

                        tfs_to_create_pages.add(split[0]);

                        String tf_key = "";
                        ArrayList<String> tf_key_set = new ArrayList<>();

                        int count = 0;
                        for (String s : split) {
                            if (count == 0) {
                                tf_key = s;
                            }
                            if (count != 0) {
                                tf_key_set.add(s);
                            }
                            sb_threshold.append("\t\t\t\t<th>");
                            if (count != 0 || s.equals("TF")) {
                                sb_threshold.append(s.toUpperCase());
                            } else {
                                File f_try = new File(options_intern.com2pose_working_directory + File.separator +
                                        options_intern.folder_out_website + File.separator +
                                        options_intern.folder_out_website_htmls_regression_coefficients +
                                        File.separator + fileDir.getName() + File.separator +
                                        options_intern.folder_out_website_htmls_TFs + File.separator + s.toUpperCase() +
                                        ".html");
                                if (f_try.exists()) {
                                    sb_threshold.append("<a href='");
                                    //sb_threshold.append(f_try.getAbsolutePath());
                                    sb_threshold.append("TFs" + File.separator + f_try.getName());
                                    sb_threshold.append(
                                            "' target='_blank'><button class=\"button\">" + s.toUpperCase() +
                                                    "</button>");
                                    sb_threshold.append("</a>");

                                    total_number_tfs.add(s.toUpperCase());
                                    total_number_tfs_hm.add(s.toUpperCase());
                                    total_numbers_tfs_hm_diff.add(s.toUpperCase());

                                    HashMap<String, HashSet<String>> current_tf_hm;

                                    if (distinct_tf_hm_diff_same.containsKey(s.toUpperCase())) {
                                        current_tf_hm = distinct_tf_hm_diff_same.get(s.toUpperCase());
                                    } else {
                                        current_tf_hm = new HashMap<>();
                                    }

                                    HashSet<String> current_tf_hm_stage = new HashSet<>();

                                    if (current_tf_hm.containsKey(fileDir_hm.getName())) {
                                        current_tf_hm_stage = current_tf_hm.get(fileDir_hm.getName());
                                    } else {
                                        current_tf_hm_stage = new HashSet<>();
                                    }

                                    current_tf_hm_stage.add("DIFFERENT_TPS");
                                    current_tf_hm.put(fileDir_hm.getName(), current_tf_hm_stage);
                                    distinct_tf_hm_diff_same.put(s.toUpperCase(), current_tf_hm);

                                } else {
                                    sb_threshold.append(s.toUpperCase());
                                }

                            }
                            sb_threshold.append("\t\t\t\t</th>\n");
                            count++;
                        }

                        tf_gene_count.put(tf_key.toUpperCase(), tf_key_set);

                        sb_threshold.append("\t\t\t</tr>\n");

                    }
                    br_gc_different.close();

                    sb_threshold.append("\t\t</table>\n");

                    //collapse thing
                    sb_threshold.append("</div>\n");
                    sb_threshold.append("</button>\n");

                    sb_threshold.append("<h5> A total number of " + total_numbers_tfs_hm_diff.size() +
                            " distinct TFs are considered in different time points group </h5>\n");

                    HashSet<String> total_numbers_tfs_hm_same = new HashSet<>();
                    sb_threshold.append("\t <h4> Same timepoints / conditions </h4> \n" + "\t  \n");
                    sb_threshold.append("</div>\n");
                    sb_threshold.append("</div>\n");
                    sb_threshold.append(
                            "<div class=\"container\" style=\"\n" + "    width: 100%;\n" + "    position: relative;\n" +
                                    "    display: block;\n" + "    overflow-x: scroll;overflow-y: scroll;\">\n" +
                                    "\t  <iframe id=\"igraph_" + fileDir_hm.getName() + "_threshold_" +
                                    fileDir.getName() + "_same_stages.html" +
                                    "\" scrolling=\"no\" style=\"border:none;\" seamless=\"seamless\" src=\"");
                    String relative_path_same_th = ".." + File.separator + ".." + File.separator +
                            options_intern.folder_out_website_interactive_plots + File.separator + fileDir.getName() +
                            File.separator + fileDir_hm.getName() + File.separator +
                            options_intern.folder_out_website_interactive_plots_overview + File.separator +
                            fileDir_hm.getName() + "_threshold_" + fileDir.getName() + "_same_stages.html";
                    sb_threshold.append(relative_path_same_th);
                    sb_threshold.append("\" height=\"400\" width=\"2500\"></iframe>\n");
                    sb_threshold.append("</div>\n");

                    sb_threshold.append("    <div class=\"w3-twothird\">\n");
                    sb_threshold.append("  <div class=\"w3-content\">\n");


                    //collapse thing
                    sb_threshold.append("<button class=\"button_expandable\" id=\"button_" + fileDir_hm.getName() +
                            "_same\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_" +
                            fileDir_hm.getName() + "_diff','table_" + fileDir_hm.getName() + "_same')\"> GeneCounts\n");
                    //tooltip
                    sb_threshold.append(
                            "<div class=\"tooltip\"><img src=\".." + File.separator + ".." + File.separator +
                                    options_intern.folder_out_website_basics + File.separator +
                                    options_intern.folder_out_website_basics_website + File.separator +
                                    options_intern.folder_out_website_basics_website_images + File.separator +
                                    "information.png" + "\" style=\"width:45px;height:40px;\"/>" +
                                    "  <span class=\"tooltiptext\">GeneCount is the unprocessed data counted by nfcore RNA-seq. A TF which is at least expressed with GeneCount 1000 is considered highly expressed and active.</span>\n" +
                                    "</div>");
                    sb_threshold.append(
                            "<div style=\"display: none;background-color: white;color:black;\" id=\"table_" +
                                    fileDir_hm.getName() + "_same\">\n");

                    //gene count table for plot
                    sb_threshold.append("<h4>Gene Count threshold: " + options_intern.plot_cutoff_gcs + "</h4>");
                    sb_threshold.append(
                            "<h5><i>Click on TF for detailed information - if no Button is available it means that this TF was eliminated by a filter.</i></h5>\n");
                    sb_threshold.append("\t\t<table style=\"width:100%;font-size:15px;\">\n");

                    File f_gene_counts_input_same = new File(
                            options_intern.com2pose_working_directory + File.separator +
                                    options_intern.folder_out_analysis_data + File.separator +
                                    options_intern.folder_out_analysis_data_TP_LEVEL + File.separator +
                                    fileDir_hm.getName() + File.separator + fileDir.getName() + File.separator +
                                    options_intern.file_suffix_analysis_plot_data_hm_level_same);
                    if (f_gene_counts_input_same.exists()) {
                        BufferedReader br_gc_same = new BufferedReader(new FileReader(f_gene_counts_input_same));
                        String line_gc_same = "";
                        while ((line_gc_same = br_gc_same.readLine()) != null) {
                            String[] split = line_gc_same.split("\t");
                            sb_threshold.append("\t\t\t<tr>\n");

                            tfs_to_create_pages.add(split[0]);

                            String tf_key = "";
                            ArrayList<String> tf_key_set = new ArrayList<>();

                            int count = 0;
                            for (String s : split) {
                                if (count == 0) {
                                    tf_key = s;
                                }
                                if (count != 0) {
                                    tf_key_set.add(s);
                                }

                                sb_threshold.append("\t\t\t\t<th>");
                                if (count != 0 || s.equals("TF")) {
                                    sb_threshold.append(s.toUpperCase());
                                } else {
                                    File f_try = new File(options_intern.com2pose_working_directory + File.separator +
                                            options_intern.folder_out_website + File.separator +
                                            options_intern.folder_out_website_htmls_regression_coefficients +
                                            File.separator + fileDir.getName() + File.separator +
                                            options_intern.folder_out_website_htmls_TFs + File.separator +
                                            s.toUpperCase() + ".html");
                                    if (f_try.exists()) {

                                        sb_threshold.append("<a href='");
                                        //sb_threshold.append(f_try.getAbsolutePath());
                                        sb_threshold.append("TFs" + File.separator + f_try.getName());
                                        sb_threshold.append(
                                                "' target='_blank'><button class=\"button\">" + s.toUpperCase() +
                                                        "</button>");
                                        sb_threshold.append("</a>");

                                        total_number_tfs.add(s.toUpperCase());
                                        total_number_tfs_hm.add(s.toUpperCase());
                                        total_numbers_tfs_hm_same.add(s.toUpperCase());

                                        HashMap<String, HashSet<String>> current_tf_hm;

                                        if (distinct_tf_hm_diff_same.containsKey(s.toUpperCase())) {
                                            current_tf_hm = distinct_tf_hm_diff_same.get(s.toUpperCase());
                                        } else {
                                            current_tf_hm = new HashMap<>();
                                        }

                                        HashSet<String> current_tf_hm_stage = new HashSet<>();

                                        if (current_tf_hm.containsKey(fileDir_hm.getName())) {
                                            current_tf_hm_stage = current_tf_hm.get(fileDir_hm.getName());
                                        } else {
                                            current_tf_hm_stage = new HashSet<>();
                                        }

                                        current_tf_hm_stage.add("SAME_TPS");
                                        current_tf_hm.put(fileDir_hm.getName(), current_tf_hm_stage);
                                        distinct_tf_hm_diff_same.put(s.toUpperCase(), current_tf_hm);
                                    } else {
                                        sb_threshold.append(s.toUpperCase());
                                    }

                                }
                                sb_threshold.append("\t\t\t\t</th>\n");
                                count++;
                            }
                            tf_gene_count.put(tf_key.toUpperCase(), tf_key_set);

                            sb_threshold.append("\t\t\t</tr>\n");

                        }
                        br_gc_same.close();
                    }


                    sb_threshold.append("\t\t</table>\n");

                    //collapse thing
                    sb_threshold.append("</div>\n");
                    sb_threshold.append("</button>");

                    sb_threshold.append("<h5> A total number of " + total_numbers_tfs_hm_same.size() +
                            " distinct TFs are considered in same time points group. </h5>\n");

                    sb_threshold.append("<h5> A total number of " + total_number_tfs_hm.size() +
                            " distinct TFs are considered in Histone Modification " + fileDir_hm.getName() +
                            " group. </h5>\n");

                    //collapse thing
                    sb_threshold.append("<button class=\"button_expandable\" id=\"button_" + fileDir_hm.getName() +
                            "_distTFs\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_" +
                            fileDir_hm.getName() + "_distTFs','table_" + fileDir_hm.getName() + "_distTFs')\">" +
                            fileDir_hm.getName() + ": Distinct TFs\n");
                    //tooltip
                    sb_threshold.append(
                            "<div class=\"tooltip\"><img src=\".." + File.separator + ".." + File.separator +
                                    options_intern.folder_out_website_basics + File.separator +
                                    options_intern.folder_out_website_basics_website + File.separator +
                                    options_intern.folder_out_website_basics_website_images + File.separator +
                                    "information.png" + "\" style=\"width:45px;height:40px;\"/>" +
                                    "  <span class=\"tooltiptext\">Distinct TFs in a group are the TFs in the union of different and same conditions.</span>\n" +
                                    "</div>");
                    sb_threshold.append(
                            "<div class='container-buttons' style=\"display: none;background-color: white;color:black;width:100%;\" id=\"table_" +
                                    fileDir_hm.getName() + "_distTFs\">\n");
                    sb_threshold.append("\t\t<table style=\"width:100%;font-size:15px;\">\n");
                    for (String tf_distinct : total_number_tfs_hm) {
                        sb_threshold.append("<tr><th><a href='TFs" + File.separator + tf_distinct.toUpperCase() +
                                ".html' target='_blank'> <button class=\"button\">" + tf_distinct.toUpperCase() +
                                "</button></a></th></tr>\n");
                    }
                    sb_threshold.append("</table>");

                    //collapse thing
                    sb_threshold.append("</div>\n");
                    sb_threshold.append("</button>\n");


                    sb_threshold.append("    </div>\n" + "\n" + "    <div class=\"w3-third w3-center\">\n" +
                            "      <i class=\"fa fa-anchor w3-padding-64 w3-text-red\"></i>\n");
                    sb_threshold.append("    </div>\n" + "  </div>\n" + "</div>\n");
                }


                sb_threshold.append("<div class=\"container-buttons w3-content\"><h5> A total number of " +
                        total_number_tfs.size() + " distinct TFs are considered in Threshold " + fileDir.getName() +
                        " group. </h5>\n");

                //collapse thing
                sb_threshold.append(
                        "<button class=\"button_expandable\" id=\"button_distTFs\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_distTFs','table_distTFs')\"> Distinct TFs overall\n");
                //tooltip
                sb_threshold.append("<div class=\"tooltip\"><img src=\".." + File.separator + ".." + File.separator +
                        options_intern.folder_out_website_basics + File.separator +
                        options_intern.folder_out_website_basics_website + File.separator +
                        options_intern.folder_out_website_basics_website_images + File.separator + "information.png" +
                        "\" style=\"width:45px;height:40px;\"/>" +
                        "  <span class=\"tooltiptext\">Distinct TFs in a threshold are the TFs in the union of all timepoints and different and same conditions.</span>\n" +
                        "</div>");
                sb_threshold.append(
                        "<div class='container-buttons' style=\"display: none;background-color: white;color:black;width:100%;\" id=\"table_distTFs\">\n");
                sb_threshold.append("\t\t<table style=\"width:100%;font-size:15px;\">\n");
                for (String tf_distinct : total_number_tfs) {
                    sb_threshold.append("<tr><th><a href='TFs" + File.separator + tf_distinct.toUpperCase() +
                            ".html' target='_blank'> <button class=\"button\">" + tf_distinct.toUpperCase() +
                            "</button></a></th></tr>\n");
                }
                sb_threshold.append("</table>");

                //collapse thing
                sb_threshold.append("</div>\n");
                sb_threshold.append("</button>\n");

                sb_threshold.append(
                        write_regression_coeffecient_analysis_found_table_html(Double.parseDouble(fileDir.getName()),
                                options_intern.html_report_levels_3_steps));

                sb_threshold.append("</div>\n");
                sb_threshold.append("</div>\n");
                sb_threshold.append("</div>\n");

                sb_threshold.append(html_tail);


                BufferedWriter bw_thresholds = new BufferedWriter(new FileWriter(
                        fileDir.getAbsolutePath() + File.separator + "threshold_" + fileDir.getName() +
                                "_overview.html"));
                bw_thresholds.write(sb_threshold.toString());
                bw_thresholds.close();

                ArrayList<String> poss_hms_list = new ArrayList(possible_hms);
                ArrayList<String> poss_tfs_list = new ArrayList<>(tfs_to_create_pages);

                File f_output_tf_root = new File(
                        fileDir.getAbsolutePath() + File.separator + options_intern.folder_out_website_htmls_TFs);
                f_output_tf_root.mkdir();

                //create TF pages
                File f_root_target_genes = new File(options_intern.com2pose_working_directory + File.separator +
                        options_intern.folder_out_target_genes);

                for (String tf : tfs_to_create_pages) {
                    File f_output_tf_page =
                            new File(f_output_tf_root.getAbsolutePath() + File.separator + tf.toUpperCase() + ".html");
                    if (f_output_tf_page.exists()) {
                        //continue;
                    }

                    StringBuilder sb_tf_page = new StringBuilder();
                    sb_tf_page.append(get_header_html(options_intern.html_report_levels_4_steps,
                            options_intern.analysis_types_regression_coefficient_analysis));
                    sb_tf_page.append(" <script>\n" + " document.title = \"TF: " + tf.toUpperCase() + " TH: " +
                            fileDir.getName() + "\";\n" + " </script>\n");

                    ArrayList<File> files_to_consider = new ArrayList<>();

                    sb_tf_page.append("<div class=\"w3-row-padding w3-padding-64 w3-container\">\n" +
                            "  <div class=\"w3-content\">\n" + "    <div class=\"w3-twothird\">\n" +
                            "      <h1>Transcription Factor: ");
                    sb_tf_page.append(tf.toUpperCase());
                    sb_tf_page.append(" <i>details</i></h1>\n");
                    sb_tf_page.append("<h3>Threshold: " + fileDir.getName().toUpperCase() + "</h3>\n");
                    sb_tf_page.append("<h4><i>Click <a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=" +
                            tf.toUpperCase() + "' target='_blank'><button class=\"button\">" + tf.toUpperCase() +
                            "</button></a> to go to GeneCards</i></h4>\n");


                    sb_tf_page.append("<h4>Gene Count threshold: " + options_intern.plot_cutoff_gcs + "</h4>");
                    sb_tf_page.append("\t\t<table style=\"width:100%\">\n");

                    ArrayList<String> table_header = tf_gene_count.get("TF");
                    sb_tf_page.append("\t\t\t<tr>\n");
                    sb_tf_page.append("\t\t\t\t<th>TF</th>\n");
                    for (String s : table_header) {
                        sb_tf_page.append("\t\t\t\t<th>" + s + "</th>\n");
                    }
                    sb_tf_page.append("\t\t\t</tr>\n");

                    ArrayList<String> table_header_tf = tf_gene_count.get(tf.toUpperCase());
                    sb_tf_page.append("\t\t\t<tr>\n");
                    sb_tf_page.append("\t\t\t\t<th>" + tf.toUpperCase() + "</th>\n");
                    for (String s : table_header_tf) {

                        sb_tf_page.append("\t\t\t\t<th>" + s + "</th>\n");

                    }
                    sb_tf_page.append("\t\t\t</tr>\n");

                    sb_tf_page.append("</table>");

                    for (String hm : possible_hms) {

                        File f_interactive_plots = new File(options_intern.com2pose_working_directory + File.separator +
                                options_intern.folder_out_website + File.separator +
                                options_intern.folder_out_website_interactive_plots + File.separator +
                                fileDir.getName() + File.separator + hm + File.separator +
                                options_intern.folder_out_website_interactive_plots_tps + File.separator);

                        sb_tf_page.append("<div class=\"w3-row-padding w3-padding-64 w3-container\">\n" +
                                "  <div class=\"w3-content\">\n" + "    <div class=\"w3-twothird\">\n" + "      <h1>");
                        sb_tf_page.append(hm);
                        sb_tf_page.append("</h1>\n");

                        sb_tf_page.append("<h3> Group plots: </h3>\n");
                        sb_tf_page.append("</div>\n");
                        sb_tf_page.append(
                                "<button class=\"button_expandable\" style=\"width:1200px\" id=\"button_group_plots_" +
                                        hm +
                                        "\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_group_plots_" +
                                        hm + "','table_group_plots_" + hm + "')\">" + hm + " group plots\n");
                        sb_tf_page.append(
                                "<div class='container-buttons' style=\"display: none;background-color: white;color:black;width:100%;\" id=\"table_group_plots_" +
                                        hm + "\">\n");

                        sb_tf_page.append("\t\t<table style=\"width:100%;font-size:15px;\">\n");

                        for (File fileDir_plot : f_interactive_plots.listFiles()) {
                            String[] split_dot = fileDir_plot.getName().split("\\.");
                            String[] split_name = split_dot[0].split("_");

                            sb_tf_page.append("\t  \n" + "\t <tr><th><h4>" + split_name[1] + " VS " + split_name[2] +
                                    " </h4></th></tr> \n" + "\t  \n" + "<tr><th><div class=\"container\" style=\"\n" +
                                    "    width: 980px;\n" + "    position: relative;\n" + "    display: block;\n" +
                                    "    overflow-x: scroll;overflow-y: scroll;\">\n" + "\t  <iframe  id=\"igraph_" +
                                    fileDir_plot.getName() +
                                    "\" scrolling=\"no\" style=\"border:none;\" seamless=\"seamless\" src=\"");
                            sb_tf_page.append(".." + File.separator + ".." + File.separator + ".." + File.separator +
                                    options_intern.folder_out_website_interactive_plots + File.separator +
                                    fileDir.getName() + File.separator + hm + File.separator +
                                    options_intern.folder_out_website_interactive_plots_tps + File.separator +
                                    fileDir_plot.getName());
                            sb_tf_page.append("\" height=\"400\" width=\"2500\" overflow=\"scroll\"></iframe>\n");
                            sb_tf_page.append("</div></th></tr>\n");
                        }

                        sb_tf_page.append("</table>");


                        sb_tf_page.append("</div>\n</button>\n");

                        sb_tf_page.append("<div class=\"w3-twothird\">\n");

                        //include target genes
                        HashSet<String> already_found_tps = new HashSet<>();

                        File f_consider_different = new File(
                                f_root_target_genes.getAbsolutePath() + File.separator + hm + File.separator +
                                        fileDir.getName() + File.separator +
                                        options_intern.folder_out_target_genes_all_different);
                        for (File tps_different : f_consider_different.listFiles()) {
                            if (tps_different.exists() && !already_found_tps.contains(tps_different.getName())) {
                                files_to_consider.add(tps_different);
                                already_found_tps.add(tps_different.getName());
                            }
                        }

                        File f_consider_same = new File(
                                f_root_target_genes.getAbsolutePath() + File.separator + hm + File.separator +
                                        fileDir.getName() + File.separator +
                                        options_intern.folder_out_target_genes_same);
                        if (f_consider_same.exists()) {
                            for (File tps_different : f_consider_same.listFiles()) {
                                if (tps_different.exists() && !already_found_tps.contains(tps_different.getName())) {
                                    files_to_consider.add(tps_different);
                                    already_found_tps.add(tps_different.getName());
                                }
                            }
                            Boolean anything_to_write = false;
                            for (File f_considered_no_tf : files_to_consider) {
                                File f_considered =
                                        new File(f_considered_no_tf.getAbsolutePath() + File.separator + tf + ".csv");
                                if (f_considered.exists()) {
                                    anything_to_write = true;
                                }
                            }
                            if (anything_to_write) {
                                //tooltip

                                sb_tf_page.append("<h3> Target Genes - top " + options_intern.plot_top_k_genes +
                                        " (normalised value):\n");
                                sb_tf_page.append("<div class=\"tooltip\"><img src=\".." + File.separator + ".." +
                                        File.separator + ".." + File.separator +
                                        options_intern.folder_out_website_basics + File.separator +
                                        options_intern.folder_out_website_basics_website + File.separator +
                                        options_intern.folder_out_website_basics_website_images + File.separator +
                                        "information.png" + "\" style=\"width:35px;height:30px;\"/>" +
                                        "  <span class=\"tooltiptext\">TargetGenes are harvested from TEPIC output. The higher the normalized affinity value the higher is the ranking.</span>\n" +
                                        "</div></h3>\n");
                                sb_tf_page.append("<p><i>Click on Symbol for GeneCard</i></p>\n");
                            }

                            for (File f_considered_no_tf : files_to_consider) {
                                File f_considered =
                                        new File(f_considered_no_tf.getAbsolutePath() + File.separator + tf + ".csv");
                                if (f_considered.exists()) {
                                    sb_tf_page.append(
                                            "<button class=\"button_expandable\" style=\"width:1200px\" id=\"button_" +
                                                    hm + "_" + f_considered_no_tf.getName() +
                                                    "\" aria-expanded=\"false\" ondblclick=\"expand_collapse('button_" +
                                                    hm + "_" + f_considered_no_tf.getName() + "','table_" + hm + "_" +
                                                    f_considered_no_tf.getName() + "')\">" + hm + ": " +
                                                    f_considered_no_tf.getName() + "\n");
                                    sb_tf_page.append(
                                            "<div class='container-buttons' style=\"display: none;background-color: white;color:black;width:100%;\" id=\"table_" +
                                                    hm + "_" + f_considered_no_tf.getName() + "\">\n");


                                    sb_tf_page.append("<h4> Point: " + f_considered_no_tf.getName() + " </h4>\n");

                                    sb_tf_page.append("\t\t<table style=\"width:100%\">\n");

                                    BufferedReader br_point_target_genes =
                                            new BufferedReader(new FileReader(f_considered));
                                    String line_point_target_getnes = "";
                                    while ((line_point_target_getnes = br_point_target_genes.readLine()) != null) {
                                        String[] split = line_point_target_getnes.split("\t");

                                        sb_tf_page.append("\t\t\t<tr>\n");

                                        int i = 0;
                                        for (String xx : split) {
                                            if (xx.equals("NOT_AVAILABLE")) {
                                                sb_tf_page.append("\t\t\t\t<th>");
                                                sb_tf_page.append("-");
                                                sb_tf_page.append("\t\t\t\t</th>\n");
                                            } else {
                                                sb_tf_page.append("\t\t\t\t<th>");
                                                if (i == 1 && !xx.equals("SYMBOL")) {
                                                    sb_tf_page.append(
                                                            "<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=" +
                                                                    xx.toUpperCase() +
                                                                    "' target='_blank'><button class=\"button\">" +
                                                                    xx.toUpperCase() + "</button></a>");
                                                } else if (i == 2 && !xx.equals("AFFINITY")) {

                                                    DecimalFormat df = new DecimalFormat("0.00000000");
                                                    sb_tf_page.append(df.format(Double.parseDouble(xx)));

                                                } else {
                                                    sb_tf_page.append(xx.toUpperCase());
                                                }
                                                sb_tf_page.append("\t\t\t\t</th>\n");
                                            }
                                            i++;

                                        }
                                        sb_tf_page.append("\t\t\t</tr>\n");
                                    }
                                    br_point_target_genes.close();

                                    sb_tf_page.append("\t\t</table>\n");
                                    sb_tf_page.append("</div>\n</button>\n");
                                }
                            }
                        }

                    }


                    sb_tf_page.append(html_tail);
                    BufferedWriter bw_tf = new BufferedWriter(new FileWriter(f_output_tf_page));
                    bw_tf.write(sb_tf_page.toString());
                    bw_tf.close();
                }

                File f_output_distribution_analysis = new File(
                        options_intern.com2pose_working_directory + File.separator +
                                options_intern.folder_out_distribution);
                f_output_distribution_analysis.mkdir();

                File f_output_distribution_analysis_analyzed_tfs = new File(
                        f_output_distribution_analysis.getAbsolutePath() + File.separator +
                                options_intern.folder_out_distribution_analyzed_tfs);
                f_output_distribution_analysis_analyzed_tfs.mkdir();

                BufferedWriter bw_distinct_tf_hm_diff_same = new BufferedWriter(new FileWriter(new File(
                        f_output_distribution_analysis_analyzed_tfs.getAbsolutePath() + File.separator +
                                options_intern.file_suffix_distribution_analysis_analysed_tfs)));

                bw_distinct_tf_hm_diff_same.write("TF\tHM\tSTAGES");
                bw_distinct_tf_hm_diff_same.newLine();

                for (String key_tf : distinct_tf_hm_diff_same.keySet()) {
                    HashMap<String, HashSet<String>> hm_diff_stages = distinct_tf_hm_diff_same.get(key_tf);

                    for (String key_hm : hm_diff_stages.keySet()) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(key_tf);
                        sb.append("\t");
                        sb.append(key_hm);
                        sb.append("\t");

                        int i = 0;
                        for (String stages : hm_diff_stages.get(key_hm)) {
                            if (i == 0) {
                                sb.append(stages);
                            } else {
                                sb.append(";");
                                sb.append(stages);
                            }
                            i++;
                        }

                        bw_distinct_tf_hm_diff_same.write(sb.toString());
                        bw_distinct_tf_hm_diff_same.newLine();

                    }

                }

                bw_distinct_tf_hm_diff_same.close();
            }
        }

        logger.logLine("[WEBSITE] Finished creating overview website.");
    }

    /**
     * get target genes of tfs
     */
    public void get_top_k_target_genes_plots() throws IOException {
        logger.logLine("[PLOTS-TARGET-GENES] Start fetching top " + options_intern.plot_top_k_genes +
                " target genes for TFs under thresholds.");

        File f_input = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_data_plots);
        File f_output = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_target_genes);
        f_output.mkdir();

        File f_input_target_genes = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_tepic_postprocessing + File.separator +
                options_intern.folder_name_tepic_postprocessing_output);

        HashMap<String, String> ensg_gene_symbol_map = new HashMap<>();
        BufferedReader br_ensg_gene_symbol =
                new BufferedReader(new FileReader(new File(options_intern.tepic_ensg_symbol)));
        String line_ensg_gene_symbol = br_ensg_gene_symbol.readLine();
        while ((line_ensg_gene_symbol = br_ensg_gene_symbol.readLine()) != null) {
            String[] split = line_ensg_gene_symbol.split("\t");
            if (split.length > 1) {
                ensg_gene_symbol_map.put(split[0], split[1]);
            }
        }
        br_ensg_gene_symbol.close();

        for (File fileDir : f_input.listFiles()) {
            if (fileDir.isDirectory()) {
                String hm_name = fileDir.getName();

                File f_output_hm = new File(f_output.getAbsolutePath() + File.separator + hm_name);
                f_output_hm.mkdir();

                for (File fileDir_th : fileDir.listFiles()) {
                    if (fileDir_th.isDirectory()) {
                        String th = fileDir_th.getName();

                        File f_output_hm_th = new File(f_output_hm.getAbsolutePath() + File.separator + th);
                        f_output_hm_th.mkdir();

                        for (File fileDir_th_f : fileDir_th.listFiles()) {
                            if (fileDir_th_f.isFile()) {
                                String file_name = fileDir_th_f.getName();

                                File f_out_hm_th_file = new File(
                                        f_output_hm_th.getAbsolutePath() + File.separator + file_name.split("\\.")[0]);
                                f_out_hm_th_file.mkdir();

                                BufferedReader br_tf = new BufferedReader(new FileReader(fileDir_th_f));
                                String line_tf = br_tf.readLine();
                                String[] split_header = line_tf.split(",");
                                ArrayList<String> timepoints_in_header = new ArrayList<>();
                                for (int i = 1; i < split_header.length; i++) {
                                    String[] split_inner = split_header[i].split("VS");

                                    String[] split_left_side = split_inner[0].split(":");

                                    timepoints_in_header.add(split_left_side[1].trim());
                                    timepoints_in_header.add(split_inner[1].trim());
                                }

                                while ((line_tf = br_tf.readLine()) != null) {
                                    String[] split = line_tf.split(",", -1);

                                    ArrayList<String> timepoints_identified = new ArrayList<>();

                                    for (int i = 1; i < split.length; i++) {
                                        if (!split[i].equals("")) {
                                            timepoints_identified.add(timepoints_in_header.get(i - 1));
                                            timepoints_identified.add(timepoints_in_header.get(i));
                                        }
                                    }

                                    HashSet<String> already_done_tps = new HashSet<>();

                                    for (int i = 0; i < timepoints_identified.size(); i += 2) {
                                        if (timepoints_identified.get(i).equals(timepoints_identified.get(i + 1))) {
                                            continue;
                                        }

                                        if (!already_done_tps.contains(timepoints_identified.get(i))) {
                                            File f_input_target_genes_hm_group_clash = new File(
                                                    f_input_target_genes.getAbsolutePath() + File.separator + hm_name +
                                                            File.separator + timepoints_identified.get(i) + "_" +
                                                            timepoints_identified.get(i + 1));

                                            write_target_genes_of_tf(f_input_target_genes_hm_group_clash,
                                                    timepoints_identified.get(i), f_out_hm_th_file, split[0],
                                                    ensg_gene_symbol_map);

                                            already_done_tps.add(timepoints_identified.get(i));
                                        }

                                        if (!already_done_tps.contains(timepoints_identified.get(i + 1))) {
                                            File f_input_target_genes_hm_group_clash = new File(
                                                    f_input_target_genes.getAbsolutePath() + File.separator + hm_name +
                                                            File.separator + timepoints_identified.get(i) + "_" +
                                                            timepoints_identified.get(i + 1));

                                            write_target_genes_of_tf(f_input_target_genes_hm_group_clash,
                                                    timepoints_identified.get(i + 1), f_out_hm_th_file, split[0],
                                                    ensg_gene_symbol_map);

                                            already_done_tps.add(timepoints_identified.get(i + 1));
                                        }
                                    }
                                }
                                br_tf.close();
                            }
                        }
                    }
                }
            }
        }

        logger.logLine("[PLOTS-TARGET-GENES] Finished fetching top " + options_intern.plot_top_k_genes +
                " target genes for TFs under thresholds.");
    }

    /**
     * analyze joined dataframe data for interesting TFs
     */
    public void analyze_plots_data() throws IOException {

        logger.logLine("[PLOT-ANALYSE] Analyse plot data.");

        File input_read_counts = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.folder_name_deseq2_preprocessing_gene_symbols);
        File input_data = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_data_plots);

        File output_analysis_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_analysis_data);
        output_analysis_root.mkdir();

        File output_analysis_tp_level = new File(output_analysis_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_analysis_data_TP_LEVEL);
        output_analysis_tp_level.mkdir();
        File output_analysis_hm_level = new File(output_analysis_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_analysis_data_HM_LEVEL);
        output_analysis_hm_level.mkdir();
        File output_analysis_website_overview = new File(output_analysis_root.getAbsolutePath() + File.separator +
                options_intern.folder_out_analysis_data_WEBSITE_OVERVIEW);
        output_analysis_website_overview.mkdir();

        HashMap<String, HashMap<String, Double>> tp_gene_tpm_value = new HashMap<>();
        File f_input_tpms_root = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.folder_name_deseq2_preprocessing_tpm + File.separator +
                options_intern.folder_name_deseq2_preprocessing_tpm_results);

        for (File f_tpm_tp : f_input_tpms_root.listFiles()) {
            if (f_tpm_tp.isFile()) {
                String name_tp = f_tpm_tp.getName().split("\\.")[0].split("_")[0];

                HashMap<String, Double> gene_tpm = new HashMap<>();

                BufferedReader br = new BufferedReader(new FileReader(f_tpm_tp));
                String line = br.readLine();
                String[] split_header = line.split("\t");
                while ((line = br.readLine()) != null) {
                    String[] split = line.split("\t");
                    gene_tpm.put(split[0], Double.parseDouble(split[3]));
                }
                br.close();

                tp_gene_tpm_value.put(name_tp, gene_tpm);
            }
        }

        HashMap<String, String> symbol_ensg_map = new HashMap<>();
        BufferedReader br_ensg_gene_symbol =
                new BufferedReader(new FileReader(new File(options_intern.tepic_ensg_symbol)));
        String line_ensg_gene_symbol = br_ensg_gene_symbol.readLine();
        while ((line_ensg_gene_symbol = br_ensg_gene_symbol.readLine()) != null) {
            String[] split = line_ensg_gene_symbol.split("\t");
            if (split.length > 1) {
                symbol_ensg_map.put(split[1].toUpperCase(), split[0].toUpperCase());
            }
        }
        br_ensg_gene_symbol.close();

        HashMap<String, HashMap<String, HashMap<String, HashMap<String, Double>>>> cutoff_hm_tf_counts_DIFFERENT =
                new HashMap<>();
        HashMap<String, HashMap<String, HashMap<String, HashMap<String, Double>>>> cutoff_hm_tf_counts_SAME =
                new HashMap<>();

        HashMap<String, String> composed_tfs = new HashMap<>();
        HashMap<String, HashSet<String>> composed_tfs_tfs = new HashMap<>();
        BufferedReader br_composed_tfs = new BufferedReader(new FileReader(
                options_intern.com2pose_working_directory + File.separator +
                        options_intern.folder_name_tepic_postprocessing + File.separator +
                        options_intern.folder_name_tepic_postprocessing_tfs + File.separator +
                        options_intern.file_suffix_tepic_postprocessing_tfs_tfs));
        String line_composed_tfs = "";
        while ((line_composed_tfs = br_composed_tfs.readLine()) != null) {
            String[] split = line_composed_tfs.split("\t");

            HashSet temp_set = new HashSet();
            for (int i = 1; i < split.length; i++) {
                composed_tfs.put(split[i], split[0]);
                temp_set.add(split[i]);
            }
            composed_tfs_tfs.put(split[0], temp_set);
        }
        br_composed_tfs.close();

        HashMap<String, HashMap<String, Double>> tp_tf_gene_counts = new HashMap<>();

        for (File fileDir : input_read_counts.listFiles()) {
            String tp_name = fileDir.getName().split("\\.")[0];

            HashMap<String, Double> n_hm = new HashMap<>();
            HashMap<String, Double> composed_tfs_counts = new HashMap<>();

            BufferedReader br = new BufferedReader(new FileReader(fileDir));
            String line = br.readLine();
            while ((line = br.readLine()) != null) {
                String[] split = line.split("\t");

                if (composed_tfs.containsKey(split[0].toUpperCase())) {
                    String key_composed_tf = composed_tfs.get(split[0].toUpperCase());
                    double sum = 0;

                    if (composed_tfs_counts.containsKey(key_composed_tf)) {
                        sum += composed_tfs_counts.get(key_composed_tf);
                    }

                    sum += Double.parseDouble(split[2]);
                    composed_tfs_counts.put(key_composed_tf, sum);
                } else {
                    n_hm.put(split[0].toUpperCase(), Double.parseDouble(split[2]));
                }
            }
            br.close();

            n_hm.putAll(composed_tfs_counts);

            tp_tf_gene_counts.put(tp_name, n_hm);

        }


        HashMap<String, HashMap<String, HashMap<String, Boolean>>> hm_th_tf_found = new HashMap<>();

        for (File fileDir : input_data.listFiles()) {
            if (fileDir.isDirectory()) {
                String hm = fileDir.getName();
                File output_hm = new File(output_analysis_tp_level.getAbsolutePath() + File.separator + hm);
                output_hm.mkdir();

                File output_website_overview_hm =
                        new File(output_analysis_website_overview.getAbsolutePath() + File.separator + hm);
                output_website_overview_hm.mkdir();


                for (File fileDirHM_th : fileDir.listFiles()) {
                    if (fileDirHM_th.isDirectory()) {
                        HashMap<String, HashMap<String, Boolean>> th_tf_found;

                        HashMap<String, Boolean> tf_found;

                        if (hm_th_tf_found.containsKey(hm)) {
                            th_tf_found = hm_th_tf_found.get(hm);
                        } else {
                            th_tf_found = new HashMap<>();
                        }

                        if (th_tf_found.containsKey(fileDirHM_th.getName())) {
                            tf_found = th_tf_found.get(fileDirHM_th.getName());
                        } else {
                            tf_found = new HashMap<>();

                            for (String s : options_intern.website_interesting_tfs) {
                                tf_found.put(s, false);
                            }
                            th_tf_found.put(fileDirHM_th.getName(), tf_found);
                        }

                        hm_th_tf_found.put(hm, th_tf_found);

                        String th_name = fileDirHM_th.getName();
                        File output_th_tp_level =
                                new File(output_hm.getAbsolutePath() + File.separator + fileDirHM_th.getName());
                        output_th_tp_level.mkdir();

                        File output_website_overview_hm_th = new File(
                                output_website_overview_hm.getAbsolutePath() + File.separator + fileDirHM_th.getName());
                        output_website_overview_hm_th.mkdir();

                        for (File fileDir_HM_th_data : fileDirHM_th.listFiles()) {
                            if (fileDir_HM_th_data.isFile()) {
                                String name = fileDir_HM_th_data.getName().split("\\.")[0];

                                HashMap<String, HashMap<String, HashMap<String, Double>>> th;
                                HashMap<String, HashMap<String, Double>> hm_intern;

                                if (name.matches(".*different.*")) {
                                    if (cutoff_hm_tf_counts_DIFFERENT.containsKey(fileDirHM_th.getName())) {
                                        th = cutoff_hm_tf_counts_DIFFERENT.get(fileDirHM_th.getName());
                                        if (th.containsKey(hm)) {
                                            hm_intern = th.get(hm);
                                        } else {
                                            hm_intern = new HashMap<>();
                                        }
                                    } else {
                                        th = new HashMap<>();
                                        hm_intern = new HashMap<>();
                                    }
                                } else {
                                    if (cutoff_hm_tf_counts_SAME.containsKey(fileDirHM_th.getName())) {
                                        th = cutoff_hm_tf_counts_DIFFERENT.get(fileDirHM_th.getName());
                                        if (th.containsKey(hm)) {
                                            hm_intern = th.get(hm);
                                        } else {
                                            hm_intern = new HashMap<>();
                                        }
                                    } else {
                                        th = new HashMap<>();
                                        hm_intern = new HashMap<>();
                                    }
                                }

                                BufferedWriter bw = new BufferedWriter(new FileWriter(
                                        output_th_tp_level.getAbsolutePath() + File.separator +
                                                fileDir_HM_th_data.getName()));
                                StringBuilder sb = new StringBuilder();
                                sb.append("TF");
                                for (String k : tp_tf_gene_counts.keySet()) {
                                    sb.append("\t");
                                    sb.append(k);
                                }
                                bw.write(sb.toString());
                                bw.newLine();

                                BufferedReader br = new BufferedReader(new FileReader(fileDir_HM_th_data));
                                String line = br.readLine();
                                while ((line = br.readLine()) != null) {
                                    String[] split = line.split(",", -1);

                                    if (tf_found.containsKey(split[0].toUpperCase())) {
                                        tf_found.put(split[0].toUpperCase(), true);
                                    }

                                    int count = 0;
                                    for (String s : split) {
                                        if (!s.equals("")) {
                                            count++;
                                        }

                                    }

                                    String ensg_name = symbol_ensg_map.get(split[0]);

                                    if (count > options_intern.plot_cutoff_tps) {
                                        HashMap<String, Double> tf_gc = new HashMap<>();
                                        //tf_gc.put(split[0],tp_tf_gene_counts.get(name).get(split[0]));
                                        StringBuilder sb_intern = new StringBuilder();
                                        sb_intern.append(split[0]);
                                        int count_row = 0;

                                        for (String k : tp_tf_gene_counts.keySet()) {
                                            if (tp_tf_gene_counts.get(k).containsKey(split[0].toUpperCase())) {
                                                count_row += tp_tf_gene_counts.get(k).get(split[0].toUpperCase());

                                                tf_gc.put(k, tp_tf_gene_counts.get(k).get(split[0].toUpperCase()));
                                                sb_intern.append("\t");
                                                sb_intern.append(tp_tf_gene_counts.get(k).get(split[0].toUpperCase()));
                                            }
                                        }

                                        double cumm_tpm = 0.0;

                                        for (String key_group : tp_gene_tpm_value.keySet()) {
                                            HashMap<String, Double> lookup = tp_gene_tpm_value.get(key_group);
                                            if (lookup.containsKey(ensg_name)) {
                                                cumm_tpm += lookup.get(ensg_name);
                                            }
                                        }

                                        boolean passed_count_threshold = false;

                                        if (count_row >= options_intern.plot_cutoff_gcs) {
                                            passed_count_threshold = true;
                                        }

                                        boolean passed_tpm_threshold = false;
                                        if (cumm_tpm >= options_intern.plot_cutoff_tpms) {
                                            passed_tpm_threshold = true;
                                        }

                                        if (passed_count_threshold && passed_tpm_threshold) {
                                            hm_intern.put(split[0].toUpperCase(), tf_gc);
                                            //WRITE
                                            bw.write(sb_intern.toString());
                                            bw.newLine();
                                        }

                                    }
                                }
                                br.close();
                                bw.close();

                                th.put(hm, hm_intern);

                                //save hashmaps correspondingly
                                if (name.matches(".*different.*")) {
                                    cutoff_hm_tf_counts_DIFFERENT.put(th_name, th);
                                } else {
                                    cutoff_hm_tf_counts_SAME.put(th_name, th);
                                }
                            }
                        }

                        BufferedWriter bw_tf_key_overview = new BufferedWriter(new FileWriter(
                                output_website_overview_hm_th.getAbsolutePath() + File.separator +
                                        options_intern.file_suffix_website_analysis_tf_available));
                        bw_tf_key_overview.write("TF\tAVAILABLE");
                        bw_tf_key_overview.newLine();
                        for (String tf_key : tf_found.keySet()) {
                            bw_tf_key_overview.write(tf_key + "\t" + tf_found.get(tf_key));
                            bw_tf_key_overview.newLine();
                        }
                        bw_tf_key_overview.close();

                    }
                }
            }
        }

        for (String k : cutoff_hm_tf_counts_DIFFERENT.keySet()) {
            File f_out = new File(output_analysis_hm_level.getAbsolutePath() + File.separator + k);
            f_out.mkdir();

            HashMap<String, HashMap<String, HashMap<String, Double>>> available_hms =
                    cutoff_hm_tf_counts_DIFFERENT.get(k);

            ArrayList<HashMap<String, HashMap<String, Double>>> tf_lists = new ArrayList<>();

            for (String kk : available_hms.keySet()) {
                tf_lists.add(available_hms.get(kk));
            }

            HashMap<String, HashMap<String, Double>> results = new HashMap<>();

            for (int i = 0; i < tf_lists.size(); i++) {
                for (String tf : tf_lists.get(i).keySet()) {
                    int count_hms = 0;

                    for (int j = 0; j < tf_lists.size(); j++) {
                        HashMap<String, HashMap<String, Double>> list = tf_lists.get(j);
                        if (list.containsKey(tf)) {
                            count_hms++;
                        }
                    }

                    if (count_hms >= options_intern.plot_cutoff_hms) {
                        results.put(tf, tf_lists.get(i).get(tf));
                    }
                }
            }


            BufferedWriter bw = new BufferedWriter(new FileWriter(f_out.getAbsolutePath() + File.separator +
                    options_intern.file_suffix_analysis_plot_data_hm_level_different));
            StringBuilder sb = new StringBuilder();
            ArrayList<String> tp_order = new ArrayList<>();
            sb.append("TF");
            for (String tp : tp_tf_gene_counts.keySet()) {
                sb.append("\t");
                sb.append(tp);
                tp_order.add(tp);
            }
            bw.write(sb.toString());
            bw.newLine();

            for (String tf : results.keySet()) {
                StringBuilder sb_tf = new StringBuilder();
                sb_tf.append(tf);
                HashMap<String, Double> tf_info = results.get(tf);

                for (String order : tp_order) {
                    sb_tf.append("\t");
                    sb_tf.append(tf_info.get(order));
                }
                bw.write(sb_tf.toString());
                bw.newLine();


            }
            bw.close();
        }

        logger.logLine("[PLOT-ANALYSE] Finished analysing plot data.");
    }

    /**
     * creates python Scripts for all timepoints of all histone modifications, it also creates an overview plot of all different group clashes (e.g. P_L but not L1_L10) for a set threshold plots
     */
    public void create_tp_plots() throws Exception {

        logger.logLine("[PLOTS] start creating / running python scripts for plots");

        File folder_input = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_put_DYNAMITE);
        File folder_output =
                new File(options_intern.com2pose_working_directory + File.separator + options_intern.folder_plots);

        File data_output = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_data_plots);
        data_output.mkdir();

        File interactive_plots_output_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_website);
        interactive_plots_output_root.mkdir();

        File interactive_plots_output_interactive_plots = new File(
                interactive_plots_output_root.getAbsolutePath() + File.separator +
                        options_intern.folder_out_website_interactive_plots);
        interactive_plots_output_interactive_plots.mkdir();

        folder_output.mkdir();

        for (File fileDirHM : folder_input.listFiles()) {
            if (fileDirHM.isDirectory()) {
                File folder_outputHM = new File(folder_output.getAbsolutePath() + File.separator + fileDirHM.getName());
                folder_outputHM.mkdir();

                File folder_out_data_HM =
                        new File(data_output.getAbsolutePath() + File.separator + fileDirHM.getName());
                folder_out_data_HM.mkdir();


                for (double d : options_intern.plot_th_coefficient) {
                    File out_th = new File(folder_outputHM.getAbsolutePath() + File.separator + d);
                    out_th.mkdir();

                    File out_data_th = new File(folder_out_data_HM.getAbsolutePath() + File.separator + d);
                    out_data_th.mkdir();

                    File interactive_plots_output_coeff =
                            new File(interactive_plots_output_interactive_plots.getAbsolutePath() + File.separator + d);
                    interactive_plots_output_coeff.mkdir();

                    File interactive_plots_output_coeff_hm = new File(
                            interactive_plots_output_coeff.getAbsolutePath() + File.separator + fileDirHM.getName());
                    interactive_plots_output_coeff_hm.mkdir();

                    File interactive_plots_output_coeff_hm_tps = new File(
                            interactive_plots_output_coeff_hm.getAbsolutePath() + File.separator +
                                    options_intern.folder_out_website_interactive_plots_tps);
                    interactive_plots_output_coeff_hm_tps.mkdir();

                    File interactive_plots_output_coeff_hm_overview = new File(
                            interactive_plots_output_coeff_hm.getAbsolutePath() + File.separator +
                                    options_intern.folder_out_website_interactive_plots_overview);
                    interactive_plots_output_coeff_hm_overview.mkdir();

                    StringBuilder sb = new StringBuilder();

                    /*WEBSITE_INTERACTIVE_PLOTS*/
                    sb.append("import pip\n" + "\n" + "def import_or_install(package):\n" + "    try:\n" +
                            "        __import__(package)\n" + "    except ImportError:\n" +
                            "        pip.main(['install', package])\n\n");

                    sb.append("import io\n" + "from base64 import b64encode\n" +
                            "import_or_install(\"plotly.express\")\n" + "import plotly.express as px\n" +
                            "import_or_install(\"dash\")\n" + "import_or_install(\"dash_core_components\")\n" +
                            "import_or_install(\"dash_html_components\")\n" + "import dash_core_components as dcc\n" +
                            "import dash_html_components as html\n" + "from dash.dependencies import Input, Output\n" +
                            "import plotly.graph_objs as go\n\n");

                    sb.append("import_or_install(\"pandas\")\n");
                    sb.append("import_or_install(\"seaborn\")\n");
                    sb.append("import_or_install(\"matplotlib.pyplot\")\n");
                    /*WEBSITE_INTERACTIVE_PLOTS*/

                    sb.append("import pandas as pd\n");
                    sb.append("import seaborn as sns\n");
                    sb.append("import matplotlib.pyplot as plt\n");
                    sb.append("sns.set_context(\"notebook\")\n");
                    sb.append("color = \"#A6CEE3\"\n");
                    sb.append("sns.set_context(\"talk\")\n");
                    sb.append("sns.set_style(\"whitegrid\")\n");

                    HashSet<String> th_group_differentpoints = new HashSet<>();
                    HashSet<String> th_group_samepoints = new HashSet<>();

                    for (File fileDirHM_Group : fileDirHM.listFiles()) {

                        String[] group_split = fileDirHM_Group.getName().split("_");
                        String groupname1 = group_split[0];
                        String groupname2 = group_split[1];

                        if (groupname1.charAt(0) != groupname2.charAt(0)) {
                            th_group_differentpoints.add(fileDirHM_Group.getName());
                        } else {
                            th_group_samepoints.add(fileDirHM_Group.getName());
                        }
                        sb.append("plt.figure(figsize=(26, 20))\n");
                        sb.append("#create barplot for group: ");
                        sb.append(fileDirHM.getName());
                        sb.append("-");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("\n");
                        sb.append("#Read data\n");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("=pd.read_table('");
                        sb.append(fileDirHM_Group.getAbsolutePath() + File.separator +
                                options_intern.file_suffix_dynamite_output_to_be_plotted);
                        sb.append("').sort_values(['value'], ascending=False)\n");
                        sb.append("# Remove suffix\n");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("['TF'] = ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("['TF'].str.split('_', expand=True)[0]\n");
                        sb.append("# Sort and filter\n");
                        sb.append(fileDirHM_Group.getName() + "_temp1");
                        sb.append(" = ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("[");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("['value'] > ");
                        sb.append(d);
                        sb.append("].set_index('TF')\n");
                        sb.append(fileDirHM_Group.getName() + "_temp2");
                        sb.append(" = ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("[");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("['value'] < -");
                        sb.append(d);
                        sb.append("].set_index('TF')\n");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(" = ");
                        sb.append("pd.concat([");
                        sb.append(fileDirHM_Group.getName() + "_temp1, ");
                        sb.append(fileDirHM_Group.getName() + "_temp2],axis=0)\n");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(".columns = ['");
                        sb.append(fileDirHM.getName());
                        sb.append(": ");
                        sb.append(groupname1);
                        sb.append(" VS ");
                        sb.append(groupname2);
                        sb.append("']\n");
                        sb.append("if('Peak' in ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(".index):\n");
                        sb.append("    ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("=pd.DataFrame.drop(");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(",index='Peak')\n");
                        sb.append("if not ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(".empty:");
                        sb.append("    # Bar Plot\n");
                        sb.append("    time = '");
                        sb.append(fileDirHM.getName());
                        sb.append(": ");
                        sb.append(groupname1);
                        sb.append(" VS ");
                        sb.append(groupname2);
                        sb.append("'\n");
                        sb.append("    ax = sns.barplot(x= ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(".index, y=time, data=");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(", color=color)\n");
                        sb.append("    ax.set_title(time)\n");
                        sb.append("    ax.set_ylabel('Normalized feature value')\n");
                        sb.append("    ax.set_xlabel('')\n");
                        sb.append("    plt.xticks(rotation=90)\n");
                        sb.append("    plt.tight_layout()\n");
                        sb.append("    plt.savefig(f\"");
                        sb.append(out_th.getAbsolutePath() + File.separator + fileDirHM.getName() + "_" +
                                fileDirHM_Group.getName() + "_threshold_" + d + ".png");
                        sb.append("\")\n");
                        sb.append("    plt.clf()\n");
                        /*WEBSITE_INTERACTIVE_PLOTS*/
                        sb.append("    fig_fin = px.bar(");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(",x= ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append(".index, y=time)\n");
                        sb.append("    fig_fin.layout.xaxis.dtick=1\n");
                        sb.append("    fig_fin.write_html(f\"");
                        sb.append(interactive_plots_output_coeff_hm_tps.getAbsolutePath() + File.separator +
                                fileDirHM.getName() + "_" + fileDirHM_Group.getName() + "_threshold_" + d +
                                ".html\",include_plotlyjs=\"cdn\")\n");
                        /*WEBSITE_INTERACTIVE_PLOTS*/

                    }
                    sb.append("plt.figure(figsize=(26, 20))\n");
                    sb.append("# Heatmap different stages\n");
                    //Clear Peak values

                    for (String s : th_group_differentpoints) {
                        sb.append("if('Peak' in ");
                        sb.append(s);
                        sb.append(".index):\n");
                        sb.append("    ");
                        sb.append(s);
                        sb.append("=pd.DataFrame.drop(");
                        sb.append(s);
                        sb.append(",index='Peak')\n");
                    }

                    for (String s : th_group_differentpoints) {
                        sb.append(s);
                        sb.append("=");
                        sb.append(s);
                        sb.append("[~");
                        sb.append(s);
                        sb.append(".index.duplicated(keep='first')]\n");
                    }

                    sb.append("join_df_stages = pd.concat([");
                    int c = 0;
                    for (String s : th_group_differentpoints) {
                        if (c == 0) {
                            sb.append(s);
                        } else {
                            sb.append(", ");
                            sb.append(s);
                        }
                        c++;
                    }
                    sb.append("], axis=1)\n");

                    sb.append("join_df_stages.to_csv(r'");
                    sb.append(out_data_th.getAbsolutePath() + File.separator);
                    sb.append(options_intern.file_suffix_analysis_plot_data_hm_level_different);
                    sb.append("',index = True, header = True)");
                    sb.append("\n");

                    sb.append("if not ");
                    sb.append("join_df_stages");
                    sb.append(".empty:\n");

                    sb.append(
                            "    plot = sns.heatmap(join_df_stages.transpose(), cmap=\"Paired\",  square=True, vmin=1, vmax=1, cbar=False, linewidths=0.5, linecolor='black', xticklabels=True)\n");
                    sb.append("    plt.savefig(\"");
                    sb.append(out_th.getAbsolutePath() + File.separator + fileDirHM.getName() + "_threshold_" + d +
                            "_different_stages.png\")\n");

                    sb.append("    plt.clf()\n");

                    /*WEBSITE_INTERACTIVE_PLOTS*/
                    sb.append("    fig_fin = px.imshow(join_df_stages.transpose())\n");
                    sb.append("    fig_fin.update_layout(plot_bgcolor='rgb(136, 136, 136)')\n");
                    sb.append("    fig_fin.update_xaxes(gridcolor='rgb(136, 136, 136)')\n");
                    sb.append("    fig_fin.layout.xaxis.dtick=1\n");
                    sb.append("    fig_fin.update_yaxes(gridcolor='rgb(136, 136, 136)')\n");
                    sb.append("    fig_fin.layout.height=500\n");
                    sb.append("    fig_fin.layout.width=2500\n");
                    sb.append("    fig_fin.layout.yaxis.dtick=1\n");
                    //sb.append("    fig_fin.update_layout(legend=dict(xanchor=\"right\", x=-1))\n");
                    sb.append("    fig_fin.write_html(f\"");
                    sb.append(interactive_plots_output_coeff_hm_overview.getAbsolutePath() + File.separator +
                            fileDirHM.getName() + "_threshold_" + d +
                            "_different_stages.html\",include_plotlyjs=\"cdn\")\n");
                    /*WEBSITE_INTERACTIVE_PLOTS*/

                    sb.append("plt.figure(figsize=(26, 20))\n");

                    sb.append("# Heatmap same stages\n");
                    for (String s : th_group_samepoints) {
                        sb.append("if('Peak' in ");
                        sb.append(s);
                        sb.append(".index):\n");
                        sb.append("    ");
                        sb.append(s);
                        sb.append("=pd.DataFrame.drop(");
                        sb.append(s);
                        sb.append(",index='Peak')\n");
                    }

                    for (String s : th_group_samepoints) {
                        sb.append(s);
                        sb.append("=");
                        sb.append(s);
                        sb.append("[~");
                        sb.append(s);
                        sb.append(".index.duplicated(keep='first')]\n");
                    }

                    if (th_group_samepoints.size() != 0) {
                        sb.append("join_df_same = pd.concat([");
                        c = 0;
                        for (String s : th_group_samepoints) {
                            if (c == 0) {
                                sb.append(s);
                            } else {
                                sb.append(", ");
                                sb.append(s);
                            }
                            c++;
                        }
                        sb.append("], axis=1)\n");

                        sb.append("join_df_same.to_csv(r'");
                        sb.append(out_data_th.getAbsolutePath() + File.separator);
                        sb.append(options_intern.file_suffix_analysis_plot_data_hm_level_same);
                        sb.append("',index = True, header = True)");
                        sb.append("\n");

                        sb.append("if not ");
                        sb.append("join_df_same");
                        sb.append(".empty:\n");

                        sb.append(
                                "    plot = sns.heatmap(join_df_same.transpose(), cmap=\"Paired\",  square=True, vmin=1, vmax=1, cbar=False, linewidths=0.5, linecolor='black', xticklabels=True)\n");
                        sb.append("    plt.savefig(\"");
                        sb.append(out_th.getAbsolutePath() + File.separator + fileDirHM.getName() + "_threshold_" + d +
                                "_same_stages.png\")\n");
                        /*WEBSITE_INTERACTIVE_PLOTS*/
                        sb.append("    fig_fin = px.imshow(join_df_same.transpose())\n");
                        sb.append("    fig_fin.update_layout(plot_bgcolor='rgb(136, 136, 136)')\n");
                        sb.append("    fig_fin.update_xaxes(gridcolor='rgb(136, 136, 136)')\n");
                        sb.append("    fig_fin.layout.xaxis.dtick=1\n");
                        sb.append("    fig_fin.update_yaxes(gridcolor='rgb(136, 136, 136)')\n");
                        sb.append("    fig_fin.layout.height=500\n");
                        sb.append("    fig_fin.layout.width=2500\n");
                        //sb.append("    fig_fin.update_layout(legend=dict(yanchor=\"top\", y=0.99, xanchor=\"left\", x=0.01))\n");
                        sb.append("    fig_fin.layout.yaxis.dtick=1\n");
                        sb.append("    fig_fin.write_html(f\"");
                        sb.append(interactive_plots_output_coeff_hm_overview.getAbsolutePath() + File.separator +
                                fileDirHM.getName() + "_threshold_" + d +
                                "_same_stages.html\",include_plotlyjs=\"cdn\")\n");
                    }

                    /*WEBSITE_INTERACTIVE_PLOTS*/
                    BufferedWriter bw = new BufferedWriter(new FileWriter(new File(
                            out_th.getAbsolutePath() + File.separator + fileDirHM.getName() + "_" + d + ".py")));
                    bw.write(sb.toString());
                    bw.close();

                    String command =
                            "python3 " + out_th.getAbsolutePath() + File.separator + fileDirHM.getName() + "_" + d +
                                    ".py";

                    logger.logLine("[PLOTS] run plot: " + command);
                    Process child = Runtime.getRuntime().exec(command);
                    int code = child.waitFor();
                    switch (code) {
                        case 0:
                            break;
                        case 1:
                            String message = child.getErrorStream().toString();
                            throw new Exception(message);
                    }
                }
            }
        }
        logger.logLine("[PLOTS] end creating / running python scripts for plots");
    }

    /**
     * run DYNAMITE.R
     */
    public void run_DYNAMITE() throws Exception {
        logger.logLine("[DYNAMITE] start running DYNAMITE with parameters: ");

        String command_base = "Rscript " + options_intern.path_to_COM2POSE + File.separator +
                options_intern.directory_for_tepic_DYNAMITE + File.separator + "DYNAMITE.R";

        File folder_input;
        folder_input = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_output_preprocessing_DYNAMITE + File.separator +
                options_intern.folder_output_preprocessing_DYNAMITE_prepareClass);


        File folder_output = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_put_DYNAMITE);
        folder_output.mkdir();

        for (File fileDirHM : folder_input.listFiles()) {
            if (fileDirHM.isDirectory()) {
                File folder_outputHM = new File(folder_output.getAbsolutePath() + File.separator + fileDirHM.getName());
                //folder_outputHM.mkdir();

                for (File fileDirHM_Group : fileDirHM.listFiles()) {
                    if (fileDirHM_Group.isDirectory()) {
                        File folder_outputHM_Group = new File(
                                folder_outputHM.getAbsolutePath() + File.separator + fileDirHM_Group.getName());
                        //folder_outputHM_Group.mkdir();

                        int count_line_br = 0;

                        for (File f_x : fileDirHM_Group.listFiles()) {
                            if (f_x.isFile()) {
                                BufferedReader br = new BufferedReader(new FileReader(f_x));
                                String line = br.readLine();
                                while ((line = br.readLine()) != null) {
                                    count_line_br++;
                                }

                                br.close();
                            }
                        }

                        if (count_line_br < 2) {
                            continue;
                        }
                        folder_outputHM.mkdir();
                        folder_outputHM_Group.mkdir();

                        String input_dir = fileDirHM_Group.getAbsolutePath() + File.separator;
                        String output_dir = folder_outputHM_Group.getAbsolutePath() + File.separator;

                        String command_edited = new String(command_base);
                        command_edited += " --dataDir=" + input_dir;
                        command_edited += " --outDir=" + output_dir;
                        command_edited += " --out_var=" + options_intern.dynamite_out_var;
                        command_edited += " --Ofolds=" + options_intern.dynamite_Ofolds;
                        command_edited += " --Ifolds=" + options_intern.dynamite_Ifolds;
                        command_edited += " --performance=" + options_intern.dynamite_performance;
                        command_edited += " --alpha=" + options_intern.dynamite_alpha;
                        command_edited += " --cores=" + options_intern.dynamite_cores;
                        if (options_intern.dynamite_randomise) {
                            command_edited += " --randomise=" + options_intern.dynamite_randomise;
                        }

                        logger.logLine("[DYNAMITE] " + fileDirHM.getName() + ":" + fileDirHM_Group.getName() +
                                " DYNAMITE.R: " + command_edited);
                        logger.logLine("[DYNAMITE] ... waiting ...");
                        Process child = Runtime.getRuntime().exec(command_edited);
                        int code = child.waitFor();
                        switch (code) {
                            case 0:
                                break;
                            case 1:
                                String message = child.getErrorStream().toString();
                                throw new Exception(message);
                        }

                    }
                }
            }
        }

        logger.logLine("[DYNAMITE] finished running DYNAMITE");
    }

    /**
     * install required R packages for DYNAMITE
     */
    public void install_required_packages() throws Exception {
        logger.logLine("[DYNAMITE] installing required packages... please wait ...");
        StringBuilder sb = new StringBuilder();

        sb.append("packages <- c(\"ggplot2\", \"dplyr\", \"gplots\", \"glmnet\", \"doMC\", \"methods\")\n" + "\n" +
                "install.packages(setdiff(packages, rownames(installed.packages())))\n");

        File file_output_script = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_output_preprocessing_DYNAMITE + File.separator +
                options_intern.folder_output_preprocessing_DYNAMITE_install_required_packages);
        file_output_script.mkdir();
        File file_output_script_to_exectue = new File(file_output_script.getAbsolutePath() + File.separator +
                options_intern.file_suffix_output_preprocessing_DYNAMITE_install_required_packages_dynamite);

        BufferedWriter bw = new BufferedWriter(new FileWriter(file_output_script_to_exectue));
        bw.write(sb.toString());
        bw.close();

        String command = "Rscript " + file_output_script_to_exectue.getAbsolutePath();

        Process child = Runtime.getRuntime().exec(command);
        int code = child.waitFor();
        switch (code) {
            case 0:
                break;
            case 1:
                String message = child.getErrorStream().toString();
                throw new Exception(message);
        }
        logger.logLine("[DYNAMITE] required packages ready!");
    }

    /**
     * run integrateData.py and prepareForClassification.R
     */
    public void preprocess_dynamite() throws Exception {

        logger.logLine("[DYNAMITE]: start preprocessing data for DYNAMITE");

        File folder;
        if (options_intern.path_tgen.equals("")) {
            folder = new File(options_intern.com2pose_working_directory + File.separator +
                    options_intern.folder_name_tepic_postprocessing + File.separator +
                    options_intern.folder_name_tepic_postprocessing_output);
        } else {
            if (options_intern.tgen_self_regulatory) {
                folder = new File(
                        options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_tgen +
                                File.separator + options_intern.folder_name_tgen_integrate);
            } else {
                folder = new File(
                        options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_tgen +
                                File.separator + options_intern.folder_name_tgen_filter_target_genes);
            }
        }

        File folder_output_pre = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_output_preprocessing_DYNAMITE);
        folder_output_pre.mkdir();

        File folder_output = new File(folder_output_pre.getAbsolutePath() + File.separator +
                options_intern.folder_output_preprocessing_DYNAMITE_integrateData);
        folder_output.mkdir();

        for (File fileDirHM : folder.listFiles()) {
            if (fileDirHM.isDirectory()) {
                File folder_output_hm =
                        new File(folder_output.getAbsolutePath() + File.separator + fileDirHM.getName());
                folder_output_hm.mkdir();

                for (File fileDirGroups : fileDirHM.listFiles()) {
                    if (fileDirGroups.isDirectory()) {
                        File folder_output_group =
                                new File(folder_output_hm.getAbsolutePath() + File.separator + fileDirGroups.getName());
                        folder_output_group.mkdir();

                        String[] groups_splitted = fileDirGroups.getName().split("_");
                        String groupname1 = groups_splitted[0];
                        String groupname2 = groups_splitted[1];

                        File input_ratios;
                        if (options_intern.path_tgen.equals("")) {
                            input_ratios = new File(fileDirGroups.getAbsolutePath() + File.separator +
                                    options_intern.folder_name_tepic_postprocessing_output_ratios + File.separator +
                                    options_intern.file_suffix_tepic_postprocessing_output_ratios +
                                    fileDirGroups.getName() + ".txt");
                        } else {
                            input_ratios = new File(fileDirGroups.getAbsolutePath() + File.separator +
                                    options_intern.file_suffix_tepic_postprocessing_output_ratios + groupname1 + "_" +
                                    groupname2 + ".txt");
                        }

                        File input_diff_gene_expr = new File(
                                options_intern.com2pose_working_directory + File.separator +
                                        options_intern.folder_name_deseq2_output + File.separator +
                                        fileDirGroups.getName() + options_intern.file_suffix_deseq2_output_DYNAMITE);

                        String command = "python3 " + options_intern.path_to_COM2POSE + File.separator +
                                options_intern.directory_for_tepic_DYNAMITE + File.separator + "integrateData.py";
                        command += " " + input_ratios.getAbsolutePath();
                        command += " " + input_diff_gene_expr.getAbsolutePath();
                        command += " " + folder_output_group.getAbsolutePath() + File.separator +
                                options_intern.file_suffix_output_preprocessing_DYNAMITE_integrateData_log2coeff;
                        if (options_intern.dynamite_preprocessing_integrate_data_geneIDs != 0) {
                            command += " " + options_intern.dynamite_preprocessing_integrate_data_geneIDs;
                        }
                        if (options_intern.dynamite_preprocessing_integrate_data_log2fc != 1) {
                            command += " " + options_intern.dynamite_preprocessing_integrate_data_log2fc;
                        }
                        if (!options_intern.dynamite_preprocessing_integrate_data_consider_geneFile.equals("")) {
                            command += " " + options_intern.dynamite_preprocessing_integrate_data_consider_geneFile;
                        }

                        logger.logLine("[DYNAMITE] " + fileDirHM.getName() + ":" + fileDirGroups.getName() +
                                " preprocessing integrateData.py: " + command);
                        Process child = Runtime.getRuntime().exec(command);
                        int code = child.waitFor();
                        switch (code) {
                            case 0:
                                break;
                            case 1:
                                String message = child.getErrorStream().toString();
                                throw new Exception(message);
                        }
                    }
                }

            }
        }

        File folder_output_classification = new File(
                folder_output_pre + File.separator + options_intern.folder_output_preprocessing_DYNAMITE_prepareClass);
        folder_output_classification.mkdir();

        for (File fileDirHM : folder_output.listFiles()) {
            if (fileDirHM.isDirectory()) {
                File folder_output_classificationHM =
                        new File(folder_output_classification.getAbsolutePath() + File.separator + fileDirHM.getName());
                folder_output_classificationHM.mkdir();

                for (File fileDirGroup : fileDirHM.listFiles()) {
                    if (fileDirGroup.isDirectory()) {
                        File folder_output_classificationHM_Group = new File(
                                folder_output_classificationHM.getAbsolutePath() + File.separator +
                                        fileDirGroup.getName());
                        folder_output_classificationHM_Group.mkdir();

                        String command = "Rscript " + options_intern.path_to_COM2POSE + File.separator +
                                options_intern.directory_for_tepic_DYNAMITE + File.separator +
                                "prepareForClassificiation.R";
                        command += " " + fileDirGroup.getAbsolutePath() + File.separator +
                                options_intern.file_suffix_output_preprocessing_DYNAMITE_integrateData_log2coeff;
                        command += " " + folder_output_classificationHM_Group.getAbsolutePath() + File.separator +
                                options_intern.file_suffix_output_preprocessing_DYNAMITE_prepClass;
                        logger.logLine(
                                "[DYNAMITE] " + fileDirGroup.getName() + " preprocessing prepareForClassification.R: " +
                                        command);
                        Process child = Runtime.getRuntime().exec(command);
                        int code = child.waitFor();
                        switch (code) {
                            case 0:
                                break;
                            case 1:
                                String message = child.getErrorStream().toString();
                                throw new Exception(message);
                        }

                    }
                }

            }
        }

        logger.logLine("[DYNAMITE]: finished preprocessing data for DYNAMITE");

    }

    /**
     * integrates self-regulatory TGene TFs into TEPIC data into one affinity file so preproccesing DYNAMITE can use it.
     */
    public void integrate_self_regulatory_tgen() throws IOException {
        logger.logLine("[TGENE-SELF-REGULATORY] Integrate self-regulatoring TFs found in TGene data into TEPIC data");

        File folder_output = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_tgen +
                        File.separator + options_intern.folder_name_tgen_integrate);
        folder_output.mkdir();
        //create necessary folder structure -> TEMPLATE: postprocess TEPIC
        File folder_tepic_postprocess = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_tgen +
                        File.separator + options_intern.folder_name_tgen_filter_target_genes);

        for (File firDir : folder_tepic_postprocess.listFiles()) {
            if (firDir.isDirectory()) {
                File folder_output_hm = new File(folder_output.getAbsolutePath() + File.separator + firDir.getName());
                folder_output_hm.mkdir();
                for (File fileDirHM : firDir.listFiles()) {
                    if (fileDirHM.isDirectory()) {
                        File folder_output_HM_group =
                                new File(folder_output_hm.getAbsolutePath() + File.separator + fileDirHM.getName());
                        folder_output_HM_group.mkdir();
                    }
                }
            }
        }


        //integrate data between TEPIC and TGENE (modifying quotients)
        //based on structure build necessary files
        for (File folderDirHM : folder_output.listFiles()) {
            if (folderDirHM.isDirectory()) {
                String hm = folderDirHM.getName();
                for (File folderDirHM_Group : folderDirHM.listFiles()) {

                    if (folderDirHM_Group.isDirectory()) {
                        logger.logLine("[TGENE-SELF-REGULATORY] integrate self-regulatory TFs " + hm + ": " +
                                folderDirHM_Group.getName());

                        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(
                                folderDirHM_Group.getAbsolutePath() + File.separator +
                                        options_intern.file_suffix_tepic_postprocessing_output_ratios +
                                        folderDirHM_Group.getName() + ".txt")));

                        File input_data_TGENE = new File(options_intern.com2pose_working_directory + File.separator +
                                options_intern.folder_name_tgen + File.separator +
                                options_intern.folder_name_tgen_groups + File.separator + hm + File.separator +
                                folderDirHM_Group.getName() + File.separator +
                                options_intern.file_suffic_tgen_output_groups);
                        File input_data_TEPIC = new File(options_intern.com2pose_working_directory + File.separator +
                                options_intern.folder_name_tepic_postprocessing + File.separator +
                                options_intern.folder_name_tepic_postprocessing_output + File.separator + hm +
                                File.separator + folderDirHM_Group.getName() + File.separator +
                                options_intern.folder_name_tepic_postprocessing_output_ratios + File.separator +
                                options_intern.file_suffix_tepic_postprocessing_output_ratios +
                                folderDirHM_Group.getName() + ".txt");

                        //new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_output_preprocessing_DYNAMITE+File.separator+options_intern.folder_output_preprocessing_DYNAMITE_integrateData+File.separator+hm+File.separator+folderDirHM_Group.getName()+File.separator+options_intern.file_suffix_output_preprocessing_DYNAMITE_integrateData_log2coeff);

                        File input_map_position_ensg = new File(
                                options_intern.com2pose_working_directory + File.separator +
                                        options_intern.folder_name_tgen + File.separator +
                                        options_intern.folder_name_tgen_preprocessing + File.separator +
                                        options_intern.folder_name_tgen_preprocessing_binary_trees + File.separator +
                                        options_intern.folder_name_tgen_preprocessing_binary_trees_sorted);

                        HashMap<String, ArrayList<ENSG_ranges_binary_trees>> regions = new HashMap<>();
                        HashMap<String, Boolean> tfs_in_tgene = new HashMap<>();

                        HashMap<String, ENSG_binary_tree> chr_binary_tree = new HashMap<>();

                        for (File chr : input_map_position_ensg.listFiles()) {
                            ArrayList<ENSG_ranges_binary_trees> chr_region_list = new ArrayList<>();

                            BufferedReader br = new BufferedReader(new FileReader(chr));
                            String line = br.readLine();
                            while ((line = br.readLine()) != null) {
                                String[] split = line.split("\t");

                                ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                                iu.chromosome = chr.getName();
                                iu.number = Integer.parseInt(split[0]);
                                iu.left_border = Integer.parseInt(split[1]);
                                iu.right_border = Integer.parseInt(split[2]);
                                iu.ensgs.addAll(Arrays.asList(split[3].split(";")));

                                chr_region_list.add(iu);

                            }
                            br.close();

                            regions.put(chr.getName(), chr_region_list);
                        }

                        for (String chr : regions.keySet()) {
                            ArrayList<ENSG_ranges_binary_trees> current_regions = regions.get(chr);

                            ENSG_binary_tree_node root =
                                    new ENSG_binary_tree_node(current_regions.get(0), current_regions.get(0).number);
                            ENSG_binary_tree bin_tree = new ENSG_binary_tree(root);

                            for (int i = 1; i < current_regions.size(); i++) {
                                bin_tree.add(current_regions.get(i).number, current_regions.get(i));

                            }

                            chr_binary_tree.put(chr.split("\\.")[0], bin_tree);

                        }

                        BufferedReader br_tgene = new BufferedReader(new FileReader(input_data_TGENE));
                        String line_tgene = br_tgene.readLine();

                        HashMap<String, ENSG_ranges_binary_trees> ensgTargetGenes_TFs = new HashMap<>();

                        while ((line_tgene = br_tgene.readLine()) != null) {
                            String split[] = line_tgene.split("\t");

                            String chr = split[0];

                            ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                            iu.chromosome = chr;
                            iu.left_border = Integer.parseInt(split[1]);
                            iu.right_border = Integer.parseInt(split[2]);
                            iu.ensgs.addAll(Arrays.asList(split[3].split(";")));

                            ENSG_binary_tree tree;

                            if (chr.startsWith("M")) {
                                tree = chr_binary_tree.get("chr" + options_intern.tgen_mt_writing);
                            } else {
                                tree = chr_binary_tree.get("chr" + chr);
                            }

                            if (tree == null) {
                                continue;
                            }

                            ENSG_ranges_binary_trees ensg_match = tree.containsNode(iu);

                            if (ensg_match == null) {
                                continue;
                            }

                            for (String k : ensg_match.ensgs) {
                                if (ensgTargetGenes_TFs.containsKey(k)) {
                                    ensgTargetGenes_TFs.get(k).ensgs.addAll(iu.ensgs);
                                } else {
                                    ensgTargetGenes_TFs.put(k, iu);
                                }
                            }

                            for (String k : iu.ensgs) {
                                tfs_in_tgene.put(k.toUpperCase(), false);
                            }
                        }

                        BufferedReader br_tepic = new BufferedReader(new FileReader(input_data_TEPIC));
                        String line_tepic = br_tepic.readLine();

                        bw.write(line_tepic);
                        bw.newLine();

                        HashMap<String, Integer> tf_place = new HashMap<>();

                        ArrayList<String> header = new ArrayList<>(Arrays.asList(line_tepic.split("\t".toUpperCase())));
                        for (int i = 0; i < header.size(); i++) {
                            String key = header.get(i).split("_")[0].toUpperCase();
                            header.set(i, key);
                            String[] split_doubles = key.split("::");
                            for (int j = 0; j < split_doubles.length; j++) {
                                tf_place.put(split_doubles[j], i);
                            }
                        }

                        while ((line_tepic = br_tepic.readLine()) != null) {
                            String[] split = line_tepic.split("\t");

                            StringBuilder sb = new StringBuilder();
                            sb.append(split[0]);

                            if (ensgTargetGenes_TFs.containsKey(split[0])) {
                                HashSet<String> tgene_pref_tfs = ensgTargetGenes_TFs.get(split[0]).ensgs;
                                HashSet<Integer> pref_positions = new HashSet<>();

                                for (String k : tgene_pref_tfs) {
                                    pref_positions.add(tf_place.get(k));
                                }

                                for (int i = 1; i < split.length; i++) {
                                    if (pref_positions.contains(i)) {
                                        double score = Double.parseDouble(split[i]);
                                        if (options_intern.tgen_consensus_calc.equals("INCREASE_TGENE_TFS")) {
                                            score = score + (score * (1 - options_intern.tgen_consensus));

                                        }
                                        sb.append("\t");
                                        sb.append(score);
                                    } else {
                                        double score = Double.parseDouble(split[i]);
                                        if (options_intern.tgen_consensus_calc.equals("DECREASE_NOT_TGENE_TFs")) {
                                            score *= (1 - options_intern.tgen_consensus);
                                        }
                                        sb.append("\t");
                                        sb.append(score);
                                    }
                                }
                            } else {
                                for (int i = 1; i < split.length; i++) {
                                    double score = Double.parseDouble(split[i]);
                                    if (options_intern.tgen_consensus_calc.equals("DECREASE_NOT_TGENE_TFs")) {
                                        score *= (1 - options_intern.tgen_consensus);
                                    }
                                    sb.append("\t");
                                    sb.append(score);
                                }
                            }

                            bw.write(sb.toString());
                            bw.newLine();
                        }

                        br_tepic.close();
                        bw.close();
                    }
                }
            }
        }

        logger.logLine(
                "[TGENE-SELF-REGULATORY] Finished integrating self-regulatory TFs found in TGene data into TEPIC data");
    }

    /**
     * filter targetgenes of TEPIC using TGENE output
     */
    public void filter_target_genes_tgen() throws IOException {
        logger.logLine("[TGENE] Start filtering target genes.");

        File folder_output = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_tgen +
                        File.separator + options_intern.folder_name_tgen_filter_target_genes);
        folder_output.mkdir();
        //create necessary folder structure -> TEMPLATE: postprocess TEPIC
        File folder_tepic_postprocess = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_tepic_postprocessing + File.separator +
                options_intern.folder_name_tepic_postprocessing_output);

        HashMap<String, String> composed_tfs = new HashMap<>();
        BufferedReader br_composed_tfs = new BufferedReader(new FileReader(
                options_intern.com2pose_working_directory + File.separator +
                        options_intern.folder_name_tepic_postprocessing + File.separator +
                        options_intern.folder_name_tepic_postprocessing_tfs + File.separator +
                        options_intern.file_suffix_tepic_postprocessing_tfs_tfs));
        String line_composed_tfs = "";
        while ((line_composed_tfs = br_composed_tfs.readLine()) != null) {
            String[] split = line_composed_tfs.split("\t");

            for (int i = 1; i < split.length; i++) {
                composed_tfs.put(split[i], split[0]);
            }
        }
        br_composed_tfs.close();

        for (File firDir : folder_tepic_postprocess.listFiles()) {
            if (firDir.isDirectory()) {
                File folder_output_hm = new File(folder_output.getAbsolutePath() + File.separator + firDir.getName());
                folder_output_hm.mkdir();
                for (File fileDirHM : firDir.listFiles()) {
                    if (fileDirHM.isDirectory()) {
                        File folder_output_HM_group =
                                new File(folder_output_hm.getAbsolutePath() + File.separator + fileDirHM.getName());
                        folder_output_HM_group.mkdir();
                    }
                }
            }
        }

        HashMap<String, String> ensg_gene_symbol_map = new HashMap<>();
        BufferedReader br_ensg_gene_symbol =
                new BufferedReader(new FileReader(new File(options_intern.tepic_ensg_symbol)));
        String line_ensg_gene_symbol = br_ensg_gene_symbol.readLine();
        while ((line_ensg_gene_symbol = br_ensg_gene_symbol.readLine()) != null) {
            String[] split = line_ensg_gene_symbol.split("\t");
            if (split.length > 1) {
                if (composed_tfs.containsKey(split[0])) {
                    ensg_gene_symbol_map.put(split[1].toUpperCase(), composed_tfs.get(split[0]));
                } else {
                    ensg_gene_symbol_map.put(split[1].toUpperCase(), split[0]);
                }
            }
        }
        br_ensg_gene_symbol.close();

        for (File folderDirHM : folder_output.listFiles()) {
            if (folderDirHM.isDirectory()) {
                String hm = folderDirHM.getName();
                for (File folderDirHM_Group : folderDirHM.listFiles()) {
                    if (folderDirHM_Group.isDirectory()) {
                        String group_name = folderDirHM_Group.getName();

                        logger.logLine("[TGENE] filter target genes " + hm + ": " + group_name);

                        File input_data_TGENE = new File(options_intern.com2pose_working_directory + File.separator +
                                options_intern.folder_name_tgen + File.separator +
                                options_intern.folder_name_tgen_groups + File.separator + hm + File.separator +
                                folderDirHM_Group.getName() + File.separator +
                                options_intern.file_suffic_tgen_output_groups);
                        File input_data_TEPIC = new File(options_intern.com2pose_working_directory + File.separator +
                                options_intern.folder_name_tepic_postprocessing + File.separator +
                                options_intern.folder_name_tepic_postprocessing_output + File.separator + hm +
                                File.separator + folderDirHM_Group.getName() + File.separator +
                                options_intern.folder_name_tepic_postprocessing_output_ratios + File.separator +
                                options_intern.file_suffix_tepic_postprocessing_output_ratios +
                                folderDirHM_Group.getName() + ".txt");

                        HashSet<String> available_target_genes = new HashSet<>();

                        BufferedReader br_tgene = new BufferedReader(new FileReader(input_data_TGENE));
                        String line_tgene = br_tgene.readLine();
                        while ((line_tgene = br_tgene.readLine()) != null) {
                            String[] split = line_tgene.split("\t");

                            String[] split_genes = split[3].split(";");

                            for (String s : split_genes) {
                                if (ensg_gene_symbol_map.containsKey(s.toUpperCase())) {
                                    available_target_genes.add(ensg_gene_symbol_map.get(s.toUpperCase()));
                                }
                            }
                        }
                        br_tgene.close();

                        BufferedWriter bw_tepic = new BufferedWriter(new FileWriter(
                                folderDirHM_Group.getAbsolutePath() + File.separator + input_data_TEPIC.getName()));
                        BufferedReader br_tepic = new BufferedReader(new FileReader(input_data_TEPIC));
                        String line_tepic = br_tepic.readLine();
                        if (line_tepic == null) {
                            continue;
                        }
                        bw_tepic.write(line_tepic);
                        bw_tepic.newLine();

                        while ((line_tepic = br_tepic.readLine()) != null) {
                            String[] split = line_tepic.split("\t");
                            if (available_target_genes.contains(split[0])) {
                                bw_tepic.write(line_tepic);
                                bw_tepic.newLine();
                            }
                        }

                        br_tepic.close();
                        bw_tepic.close();
                    }
                }
            }
        }
        logger.logLine("[TGENE] Finished filtering target genes.");
    }

    /**
     * creates TGENE group clash so it can be merged into TEPIC
     */
    public void create_tgen_groups() throws IOException {
        logger.logLine("[TGENE] Create TGene group data");

        check_tepic_input_with_options();
        /*
        if(options_intern.mix_option.equals("SAMPLE_LEVEL"))
        {
            File root_mix_working_dir = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option);
            File f_sample_mix_output = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_sample_mix);
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory=f_sample_mix_output.getAbsolutePath();
        }

        if(options_intern.mix_option.equals("HM_LEVEL"))
        {
            File root_mix_working_dir = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option);
            File f_output_hm = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_hm_mix);
            f_output_hm.mkdir();

            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory=f_output_hm.getAbsolutePath();
        }
        if(!options_intern.black_list_dir.equals(""))
        {
            File output_folder = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_blacklisted_regions);
            File output_folder_new_input = new File(output_folder.getAbsolutePath()+File.separator+options_intern.folder_name_blacklisted_regions_new_input);
            output_folder_new_input.mkdir();

            //set new folder directory for tepic input and save old one
            options_intern.tepic_input_prev=options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = output_folder_new_input.getAbsolutePath();
        }
        if(options_intern.mix_mutually_exclusive)
        {
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option+File.separator+options_intern.folder_name_mix_option_mutually_exclusive+File.separator+options_intern.folder_name_mix_options_mutually_exclusive_input;
        }*/

        File folder_input = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_tgen +
                        File.separator + options_intern.folder_name_tgen_merged);
        File folder_output = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_tgen +
                        File.separator + options_intern.folder_name_tgen_groups);
        folder_output.mkdir();

        //create necessary folder structure -> TEMPLATE: postprocess TEPIC
        File folder_tepic_postprocess = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_tepic_postprocessing + File.separator +
                options_intern.folder_name_tepic_postprocessing_output);
        for (File firDir : folder_tepic_postprocess.listFiles()) {
            if (firDir.isDirectory()) {
                File folder_output_hm = new File(folder_output.getAbsolutePath() + File.separator + firDir.getName());
                folder_output_hm.mkdir();
                for (File fileDirHM : firDir.listFiles()) {
                    if (fileDirHM.isDirectory()) {
                        File folder_output_HM_group =
                                new File(folder_output_hm.getAbsolutePath() + File.separator + fileDirHM.getName());
                        folder_output_HM_group.mkdir();
                    }
                }
            }
        }

        //based on structure build necessary files
        for (File folderDirHM : folder_output.listFiles()) {
            if (folderDirHM.isDirectory()) {
                String hm = folderDirHM.getName();
                for (File folderDirHM_Group : folderDirHM.listFiles()) {

                    if (folderDirHM_Group.isDirectory()) {
                        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(
                                folderDirHM_Group.getAbsolutePath() + File.separator +
                                        options_intern.file_suffic_tgen_output_groups)));
                        bw.write("CHR\tLEFT_BORDER\tRIGHT_BORDER\tENSGS");
                        bw.newLine();

                        String[] split_name = folderDirHM_Group.getName().split("_");
                        String group1 = split_name[0];
                        String group2 = split_name[1];

                        File input_group1;
                        File input_group2;

                        if (!options_intern.mix_mutually_exclusive) {
                            input_group1 = new File(
                                    folder_input.getAbsolutePath() + File.separator + group1 + File.separator + hm +
                                            File.separator + hm + "_" + group1 + ".txt");
                            input_group2 = new File(
                                    folder_input.getAbsolutePath() + File.separator + group2 + File.separator + hm +
                                            File.separator + hm + "_" + group2 + ".txt");

                        } else {
                            input_group1 = new File(
                                    folder_input.getAbsolutePath() + File.separator + folderDirHM_Group.getName() +
                                            File.separator + hm + File.separator + hm + "_" + group1 + ".txt");
                            input_group2 = new File(
                                    folder_input.getAbsolutePath() + File.separator + folderDirHM_Group.getName() +
                                            File.separator + hm + File.separator + hm + "_" + group2 + ".txt");
                        }

                        HashMap<String, ArrayList<ENSG_ranges_binary_trees>> chr_unmerged_ius = new HashMap<>();

                        BufferedReader br_group1 = new BufferedReader(new FileReader(input_group1));
                        String line_group1 = br_group1.readLine();
                        while ((line_group1 = br_group1.readLine()) != null) {
                            String[] split_group1 = line_group1.split("\t");
                            String chr = split_group1[0];
                            int left_border = Integer.parseInt(split_group1[1]);
                            int right_border = Integer.parseInt(split_group1[2]);
                            String ensg = split_group1[3].toUpperCase();

                            ArrayList<ENSG_ranges_binary_trees> current_arr_chr;
                            if (chr_unmerged_ius.containsKey(chr)) {
                                current_arr_chr = chr_unmerged_ius.get(chr);
                            } else {
                                current_arr_chr = new ArrayList<>();
                            }

                            ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                            iu.chromosome = chr;
                            iu.ensgs.add(ensg);
                            iu.left_border = left_border;
                            iu.right_border = right_border;

                            current_arr_chr.add(iu);
                            chr_unmerged_ius.put(chr, current_arr_chr);
                        }
                        br_group1.close();

                        BufferedReader br_group2 = new BufferedReader(new FileReader(input_group2));
                        String line_group2 = br_group2.readLine();
                        while ((line_group2 = br_group2.readLine()) != null) {
                            String[] split_group2 = line_group2.split("\t");
                            String chr = split_group2[0];
                            int left_border = Integer.parseInt(split_group2[1]);
                            int right_border = Integer.parseInt(split_group2[2]);
                            String ensg = split_group2[3].toUpperCase();

                            ArrayList<ENSG_ranges_binary_trees> current_arr_chr;
                            if (chr_unmerged_ius.containsKey(chr)) {
                                current_arr_chr = chr_unmerged_ius.get(chr);
                            } else {
                                current_arr_chr = new ArrayList<>();
                            }

                            ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                            iu.chromosome = chr;
                            iu.ensgs.add(ensg);
                            iu.left_border = left_border;
                            iu.right_border = right_border;

                            current_arr_chr.add(iu);
                            chr_unmerged_ius.put(chr, current_arr_chr);
                        }
                        br_group2.close();

                        //sort all arrays
                        for (String s : chr_unmerged_ius.keySet()) {
                            ArrayList<ENSG_ranges_binary_trees> x = chr_unmerged_ius.get(s);
                            Collections.sort(x);
                            chr_unmerged_ius.put(s, x);
                        }


                        //remove duplicates or integrate ensgs if in overlap
                        for (String chr : chr_unmerged_ius.keySet()) {
                            ArrayList<ENSG_ranges_binary_trees> unmerged = chr_unmerged_ius.get(chr);

                            ArrayList<ENSG_ranges_binary_trees> unmerged_temp = new ArrayList<>();

                            for (int i = 1; i < unmerged.size(); i++) {
                                ENSG_ranges_binary_trees iu_before = unmerged.get(i - 1);
                                ENSG_ranges_binary_trees iu_after = unmerged.get(i);

                                if (iu_before.isTheSame(iu_after)) {
                                    unmerged_temp.add(iu_before);
                                    i++;
                                } else {
                                    unmerged_temp.add(iu_before);
                                }
                            }

                            unmerged = new ArrayList<>(unmerged_temp);
                            unmerged_temp.clear();

                            //check for same_ranges but not same ENSGs and merge them
                            for (int i = 1; i < unmerged.size(); i++) {
                                ENSG_ranges_binary_trees iu_before = unmerged.get(i - 1);
                                ENSG_ranges_binary_trees iu_after = unmerged.get(i);

                                if (iu_before.isSameRange(iu_after)) {
                                    for (String s : iu_after.ensgs) {
                                        iu_before.ensgs.add(s);
                                    }
                                    unmerged_temp.add(iu_before);

                                    i++;
                                } else {
                                    unmerged_temp.add(iu_before);
                                }
                            }

                            unmerged = new ArrayList<>(unmerged_temp);
                            unmerged_temp.clear();

                            /*
                            boolean something_changed = true;
                            while(something_changed)
                            {
                                something_changed=false;
                                //check for overlaps and if so join them
                                for(int i = 1; i < unmerged.size();i++)
                                {
                                    ENSG_ranges_binary_trees iu_before = unmerged.get(i - 1);
                                    ENSG_ranges_binary_trees iu_after = unmerged.get(i);

                                    if(iu_before.isOverLap(iu_after))
                                    {
                                        iu_before.right_border=iu_after.right_border;
                                        for(String s : iu_after.ensgs)
                                        {
                                            iu_before.ensgs.add(s);
                                        }
                                        something_changed=true;
                                        i++;
                                    }
                                    else
                                    {
                                        unmerged_temp.add(iu_before);
                                    }
                                }
                                unmerged=new ArrayList<>(unmerged_temp);
                                unmerged_temp.clear();
                            }*/


                            for (int i = 0; i < unmerged.size(); i++) {
                                bw.write(unmerged.get(i).toString());
                                bw.newLine();
                            }
                        }
                        bw.close();
                    }
                }
            }
        }
        logger.logLine("[TGENE] Finished creating TGene group data");
    }

    /**
     * merge different samples of one timepoint of TGENE
     */
    public void merge_tgen() throws IOException {
        logger.logLine("[TGENE] Postprocessing: Merge TGene results");

        File folder_input = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_tgen +
                        File.separator + options_intern.folder_name_tgen_output);
        File folder_output = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_tgen +
                        File.separator + options_intern.folder_name_tgen_merged);
        folder_output.mkdir();

        for (File fileDirTP : folder_input.listFiles()) {
            if (fileDirTP.isDirectory()) {
                File output_tp = new File(folder_output.getAbsolutePath() + File.separator + fileDirTP.getName());
                output_tp.mkdir();

                for (File fileDirTP_HM : fileDirTP.listFiles()) {
                    if (fileDirTP_HM.isDirectory()) {
                        File output_tp_hm =
                                new File(output_tp.getAbsolutePath() + File.separator + fileDirTP_HM.getName());
                        output_tp_hm.mkdir();

                        String name_output = fileDirTP_HM.getName() + "_" + fileDirTP.getName();


                        HashMap<String, ArrayList<ENSG_ranges_binary_trees>> tfs_to_regions = new HashMap<>();

                        //this one is needed for mutually exclusive option
                        HashMap<String, HashMap<String, ArrayList<ENSG_ranges_binary_trees>>> file_tfs_to_regions =
                                new HashMap<>();


                        for (File fileDirTP_HM_Samples : fileDirTP_HM.listFiles()) {
                            if (fileDirTP_HM_Samples.isDirectory()) {
                                //this one is needed for mutually exclusive option
                                String file_name = fileDirTP_HM_Samples.getName();
                                HashMap<String, ArrayList<ENSG_ranges_binary_trees>> current_tfs_to_regions =
                                        new HashMap<>();

                                File input_data = new File(fileDirTP_HM_Samples.getAbsolutePath() + File.separator +
                                        options_intern.file_suffix_tgen_output);

                                BufferedReader br = new BufferedReader(new FileReader(input_data));
                                String line = br.readLine();
                                while ((line = br.readLine()) != null) {
                                    if (line.startsWith("#") || line.equals("")) {
                                        continue;
                                    }

                                    String[] split = line.split("\t");

                                    ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                                    iu.ensgs.add(split[1]);

                                    String[] split_chr = split[6].split(":");
                                    String[] split_borders = split_chr[1].split("-");

                                    iu.chromosome = split_chr[0];
                                    iu.left_border = Integer.parseInt(split_borders[0]);
                                    iu.right_border = Integer.parseInt(split_borders[1]);

                                    if (tfs_to_regions.containsKey(split_chr[0])) {
                                        ArrayList<ENSG_ranges_binary_trees> x = tfs_to_regions.get(split_chr[0]);
                                        x.add(iu);
                                        tfs_to_regions.put(split_chr[0], x);

                                    } else {
                                        ArrayList<ENSG_ranges_binary_trees> x = new ArrayList<>();
                                        x.add(iu);
                                        tfs_to_regions.put(split_chr[0], x);
                                    }

                                    if (current_tfs_to_regions.containsKey(split_chr[0])) {
                                        ArrayList<ENSG_ranges_binary_trees> x =
                                                current_tfs_to_regions.get(split_chr[0]);
                                        x.add(iu);
                                        current_tfs_to_regions.put(split_chr[0], x);

                                    } else {
                                        ArrayList<ENSG_ranges_binary_trees> x = new ArrayList<>();
                                        x.add(iu);
                                        current_tfs_to_regions.put(split_chr[0], x);
                                    }
                                }

                                file_tfs_to_regions.put(file_name, current_tfs_to_regions);
                            }
                        }

                        //sort regions
                        for (String chr : tfs_to_regions.keySet()) {
                            Collections.sort(tfs_to_regions.get(chr));
                        }

                        for (String file_name : file_tfs_to_regions.keySet()) {
                            HashMap<String, ArrayList<ENSG_ranges_binary_trees>> current_tfs_to_regions =
                                    file_tfs_to_regions.get(file_name);

                            for (String chr : current_tfs_to_regions.keySet()) {
                                Collections.sort(current_tfs_to_regions.get(chr));
                            }
                        }

                        for (String chr : tfs_to_regions.keySet()) {
                            ArrayList<ENSG_ranges_binary_trees> tf_reg = tfs_to_regions.get(chr);

                            ArrayList<ENSG_ranges_binary_trees> tf_temp = new ArrayList<>();
                            for (int i = 1; i < tf_reg.size(); i++) {
                                //delete all duplicates
                                ENSG_ranges_binary_trees iu_before = tf_reg.get(i - 1);
                                ENSG_ranges_binary_trees iu_now = tf_reg.get(i);

                                if (iu_before.isTheSame(iu_now)) {
                                    int temp_i = i + 1;
                                    boolean same = true;
                                    int count_same = 0;
                                    while (temp_i < tf_reg.size() && same) {
                                        same = false;

                                        ENSG_ranges_binary_trees iu_follow = tf_reg.get(temp_i);

                                        if (iu_follow.isTheSame(iu_before)) {
                                            same = true;
                                            count_same++;
                                        }

                                        temp_i++;
                                    }
                                    i = temp_i + 1;


                                    tf_temp.add(iu_before);
                                } else {
                                    tf_temp.add(iu_before);
                                }
                            }
                            tfs_to_regions.put(chr, tf_temp);
                        }

                        for (String file_name : file_tfs_to_regions.keySet()) {
                            HashMap<String, ArrayList<ENSG_ranges_binary_trees>> current_tfs_to_regions =
                                    file_tfs_to_regions.get(file_name);

                            for (String chr : current_tfs_to_regions.keySet()) {
                                ArrayList<ENSG_ranges_binary_trees> tf_reg = current_tfs_to_regions.get(chr);

                                ArrayList<ENSG_ranges_binary_trees> tf_temp = new ArrayList<>();
                                for (int i = 1; i < tf_reg.size(); i++) {
                                    //delete all duplicates
                                    ENSG_ranges_binary_trees iu_before = tf_reg.get(i - 1);
                                    ENSG_ranges_binary_trees iu_now = tf_reg.get(i);

                                    if (iu_before.isTheSame(iu_now)) {
                                        int temp_i = i + 1;
                                        boolean same = true;
                                        int count_same = 0;
                                        while (temp_i < tf_reg.size() && same) {
                                            same = false;

                                            ENSG_ranges_binary_trees iu_follow = tf_reg.get(temp_i);

                                            if (iu_follow.isTheSame(iu_before)) {
                                                same = true;
                                                count_same++;
                                            }

                                            temp_i++;
                                        }
                                        i = temp_i + 1;


                                        tf_temp.add(iu_before);
                                    } else {
                                        tf_temp.add(iu_before);
                                    }
                                }
                                current_tfs_to_regions.put(chr, tf_temp);
                            }
                            file_tfs_to_regions.put(file_name, current_tfs_to_regions);
                        }

                        //CHECK TPM if TPM filter is set! -> non mutually exclusive
                        if (options_intern.tepic_tpm_cutoff > 0) {
                                /*command_tail += " -T " + options_intern.tepic_tpm_cutoff;
                                command_tail += " -E " + options_intern.tepic_ensg_symbol;
                                command_tail += " -A " + options_intern.deseq2_input_gene_id;
                                String n_dir = options_intern.com2pose_working_directory+File.separator+ options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_single+File.separator+dirGroup.getName()+options_intern.file_suffix_deseq2_preprocessing_meanCounts;
                                command_tail_sample += " -G " + n_dir;

                                [-G input genes count file, if set, default TPM is 1]\n
                                [-T set (T)ranscripts (P)er (M)illion cutoff must be in float form (e.g. 1.0)]
                                [-A input gene annotation desq2 file, required for TPM filter]
                                [-E input ensg to gene symbol file, required for TPM filter*/

                            logger.logLine(
                                    "[TGENE] TPM filter is set. Filtering TGENE results: " + fileDirTP_HM.getName() +
                                            " - " + fileDirTP.getName() + ".");

                            double tpm_cutoff = options_intern.tepic_tpm_cutoff;
                            File f_ensg_map_symbol = new File(options_intern.tepic_ensg_symbol);
                            File f_gene_annot_deseq2 = new File(options_intern.deseq2_input_gene_id);
                            File f_ref_genome = new File(options_intern.tepic_gene_annot);

                            ArrayList<Integer> gene_counts = new ArrayList<>();

                            if (!options_intern.mix_mutually_exclusive) {
                                File f_gene_count_group1;

                                f_gene_count_group1 = new File(
                                        options_intern.com2pose_working_directory + File.separator +
                                                options_intern.folder_name_deseq2_preprocessing + File.separator +
                                                options_intern.folder_name_deseq2_preprocessing_single +
                                                File.separator + fileDirTP.getName() +
                                                options_intern.file_suffix_deseq2_preprocessing_meanCounts);
                                BufferedReader br_gene_counts_1 =
                                        new BufferedReader(new FileReader(f_gene_count_group1));
                                String line_gene_count_1 = br_gene_counts_1.readLine();
                                while ((line_gene_count_1 = br_gene_counts_1.readLine()) != null) {
                                    gene_counts.add(Integer.parseInt(line_gene_count_1));
                                }
                                br_gene_counts_1.close();
                            } else {

                                String names[] = fileDirTP.getName().split("_");

                                File f_gene_count_group1 = new File(
                                        options_intern.com2pose_working_directory + File.separator +
                                                options_intern.folder_name_deseq2_preprocessing + File.separator +
                                                options_intern.folder_name_deseq2_preprocessing_single +
                                                File.separator + names[0] +
                                                options_intern.file_suffix_deseq2_preprocessing_meanCounts);
                                File f_gene_count_group2 = new File(
                                        options_intern.com2pose_working_directory + File.separator +
                                                options_intern.folder_name_deseq2_preprocessing + File.separator +
                                                options_intern.folder_name_deseq2_preprocessing_single +
                                                File.separator + names[1] +
                                                options_intern.file_suffix_deseq2_preprocessing_meanCounts);

                                HashMap<Integer, Integer> position_counts = new HashMap<>();

                                BufferedReader br_gene_counts_1 =
                                        new BufferedReader(new FileReader(f_gene_count_group1));
                                int counts = 0;
                                String line_gene_count_1 = br_gene_counts_1.readLine();
                                while ((line_gene_count_1 = br_gene_counts_1.readLine()) != null) {
                                    position_counts.put(counts, Integer.parseInt(line_gene_count_1));
                                    counts++;
                                }
                                br_gene_counts_1.close();

                                BufferedReader br_gene_counts_2 =
                                        new BufferedReader(new FileReader(f_gene_count_group2));
                                counts = 0;
                                String line_gene_count_2 = br_gene_counts_2.readLine();
                                while ((line_gene_count_2 = br_gene_counts_2.readLine()) != null) {
                                    int mean_count =
                                            (position_counts.get(counts) + Integer.parseInt(line_gene_count_2)) / 2;
                                    position_counts.put(counts, mean_count);
                                    counts++;
                                }
                                br_gene_counts_2.close();

                                for (int i = 0; i < position_counts.size(); i++) {
                                    gene_counts.add(position_counts.get(i));
                                }
                            }

                            ArrayList<String> ensg_numbers = new ArrayList<>();
                            BufferedReader br_ensg_numbers = new BufferedReader(new FileReader(f_gene_annot_deseq2));
                            String line_ensg_numbers = br_ensg_numbers.readLine();
                            while ((line_ensg_numbers = br_ensg_numbers.readLine()) != null) {
                                ensg_numbers.add(line_ensg_numbers.toUpperCase());
                            }
                            br_ensg_numbers.close();

                            int total_number_samples = 0;
                            HashMap<String, Integer> ensg_counts = new HashMap<>();
                            for (int i = 0; i < gene_counts.size(); i++) {
                                ensg_counts.put(ensg_numbers.get(i).toUpperCase(), gene_counts.get(i));
                                total_number_samples++;
                            }

                            HashMap<String, String> ensg_symb = new HashMap<>();

                            BufferedReader br_ensg_symb = new BufferedReader(new FileReader(f_ensg_map_symbol));
                            String line_ensg_symb = br_ensg_symb.readLine();
                            while ((line_ensg_symb = br_ensg_symb.readLine()) != null) {
                                String[] split = line_ensg_symb.split("\t");
                                if (split.length > 1) {
                                    ensg_symb.put(split[1].toUpperCase(), split[0]);
                                }
                            }
                            br_ensg_symb.close();

                            HashMap<String, Integer> gene_symbol_lengths = new HashMap<>();
                            BufferedReader br_gene_symbol_lengths = new BufferedReader(new FileReader(f_ref_genome));
                            String line_gene_symbol_lengths = "";
                            while ((line_gene_symbol_lengths = br_gene_symbol_lengths.readLine()) != null) {
                                if (line_gene_symbol_lengths.startsWith("#")) {
                                    continue;
                                }
                                String[] split = line_gene_symbol_lengths.split("\t");
                                if (split[2].equals("transcript")) {
                                    int start = Integer.parseInt(split[3]);
                                    int end = Integer.parseInt(split[4]);
                                    int diff = end - start + 1;
                                    String[] split_further = split[8].split(" ");
                                    String ensg_name = split_further[1].replace("\"", "");
                                    ensg_name = ensg_name.replace(";", "");
                                    String[] ensg_name_split = ensg_name.split("\\.");
                                    gene_symbol_lengths.put(ensg_name_split[0], diff);
                                }
                            }
                            br_gene_symbol_lengths.close();

                            double all_rpk = 0.0;

                            for (String ec : ensg_counts.keySet()) {
                                if (!ensg_counts.containsKey(ec)) {
                                    continue;
                                }

                                int ec_count = ensg_counts.get(ec);
                                int ec_lengths = 0;
                                if (!gene_symbol_lengths.containsKey(ec)) {
                                    continue;
                                }
                                ec_lengths = gene_symbol_lengths.get(ec);
                                if (ec_lengths == 0) {
                                    continue;
                                }
                                double norm_ec_lengths = ec_lengths / (1000 * 1.0);
                                double current_rpk = ec_count / (norm_ec_lengths * 1.0);
                                all_rpk += current_rpk;
                            }

                            double scaling_factor = all_rpk / (1000000 * 1.0);

                            for (String chr : tfs_to_regions.keySet()) {
                                ArrayList<ENSG_ranges_binary_trees> tf_reg = tfs_to_regions.get(chr);
                                ArrayList<ENSG_ranges_binary_trees> tf_reg_temp = new ArrayList<>();

                                for (ENSG_ranges_binary_trees iu : tf_reg) {
                                    boolean should_write = false;

                                    //check if we want that TF
                                    for (String ec : iu.ensgs) {
                                        ec = ec.toUpperCase();

                                        if (ensg_symb.containsKey(ec)) {
                                            ec = ensg_symb.get(ec);
                                        } else {
                                            continue;
                                        }

                                        if (!ensg_counts.containsKey(ec)) {
                                            continue;
                                        }

                                        int ec_count = ensg_counts.get(ec);
                                        int ec_lengths = 0;
                                        if (!gene_symbol_lengths.containsKey(ec)) {
                                            continue;
                                        }
                                        ec_lengths = gene_symbol_lengths.get(ec);
                                        if (ec_lengths == 0) {
                                            continue;
                                        }
                                        double norm_ec_lengths = ec_lengths / (1000 * 1.0);
                                        double current_rpk = ec_count / (norm_ec_lengths * 1.0);

                                        if (current_rpk >= scaling_factor) {
                                            should_write = true;
                                        }
                                    }

                                    if (should_write) {
                                        tf_reg_temp.add(iu);
                                    }
                                }

                                tfs_to_regions.put(chr, tf_reg_temp);
                            }
                        }

                        //CHECK TPM if TPM filter is set! -> mutually exclusive
                        for (String file_name : file_tfs_to_regions.keySet()) {
                            HashMap<String, ArrayList<ENSG_ranges_binary_trees>> current_tfs_to_regions =
                                    file_tfs_to_regions.get(file_name);

                            if (options_intern.tepic_tpm_cutoff > 0) {
                                /*command_tail += " -T " + options_intern.tepic_tpm_cutoff;
                                command_tail += " -E " + options_intern.tepic_ensg_symbol;
                                command_tail += " -A " + options_intern.deseq2_input_gene_id;
                                String n_dir = options_intern.com2pose_working_directory+File.separator+ options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_single+File.separator+dirGroup.getName()+options_intern.file_suffix_deseq2_preprocessing_meanCounts;
                                command_tail_sample += " -G " + n_dir;

                                [-G input genes count file, if set, default TPM is 1]\n
                                [-T set (T)ranscripts (P)er (M)illion cutoff must be in float form (e.g. 1.0)]
                                [-A input gene annotation desq2 file, required for TPM filter]
                                [-E input ensg to gene symbol file, required for TPM filter*/

                                //logger.logLine("[TGENE] TPM filter is set. Filtering TGENE results: " + fileDirTP_HM.getName() + " - " + fileDirTP.getName()+".");

                                double tpm_cutoff = options_intern.tepic_tpm_cutoff;
                                File f_ensg_map_symbol = new File(options_intern.tepic_ensg_symbol);
                                File f_gene_annot_deseq2 = new File(options_intern.deseq2_input_gene_id);
                                File f_ref_genome = new File(options_intern.tepic_gene_annot);

                                ArrayList<Integer> gene_counts = new ArrayList<>();

                                if (!options_intern.mix_mutually_exclusive) {
                                    File f_gene_count_group1;

                                    f_gene_count_group1 = new File(
                                            options_intern.com2pose_working_directory + File.separator +
                                                    options_intern.folder_name_deseq2_preprocessing + File.separator +
                                                    options_intern.folder_name_deseq2_preprocessing_single +
                                                    File.separator + fileDirTP.getName() +
                                                    options_intern.file_suffix_deseq2_preprocessing_meanCounts);
                                    BufferedReader br_gene_counts_1 =
                                            new BufferedReader(new FileReader(f_gene_count_group1));
                                    String line_gene_count_1 = br_gene_counts_1.readLine();
                                    while ((line_gene_count_1 = br_gene_counts_1.readLine()) != null) {
                                        gene_counts.add(Integer.parseInt(line_gene_count_1));
                                    }
                                    br_gene_counts_1.close();
                                } else {

                                    String names[] = fileDirTP.getName().split("_");

                                    File f_gene_count_group1 = new File(
                                            options_intern.com2pose_working_directory + File.separator +
                                                    options_intern.folder_name_deseq2_preprocessing + File.separator +
                                                    options_intern.folder_name_deseq2_preprocessing_single +
                                                    File.separator + file_name.split("_")[0] +
                                                    options_intern.file_suffix_deseq2_preprocessing_meanCounts);
                                    File f_gene_count_group2 = new File(
                                            options_intern.com2pose_working_directory + File.separator +
                                                    options_intern.folder_name_deseq2_preprocessing + File.separator +
                                                    options_intern.folder_name_deseq2_preprocessing_single +
                                                    File.separator + file_name.split("_")[0] +
                                                    options_intern.file_suffix_deseq2_preprocessing_meanCounts);

                                    HashMap<Integer, Integer> position_counts = new HashMap<>();

                                    BufferedReader br_gene_counts_1 =
                                            new BufferedReader(new FileReader(f_gene_count_group1));
                                    int counts = 0;
                                    String line_gene_count_1 = br_gene_counts_1.readLine();
                                    while ((line_gene_count_1 = br_gene_counts_1.readLine()) != null) {
                                        position_counts.put(counts, Integer.parseInt(line_gene_count_1));
                                        counts++;
                                    }
                                    br_gene_counts_1.close();

                                    BufferedReader br_gene_counts_2 =
                                            new BufferedReader(new FileReader(f_gene_count_group2));
                                    counts = 0;
                                    String line_gene_count_2 = br_gene_counts_2.readLine();
                                    while ((line_gene_count_2 = br_gene_counts_2.readLine()) != null) {
                                        int mean_count =
                                                (position_counts.get(counts) + Integer.parseInt(line_gene_count_2)) / 2;
                                        position_counts.put(counts, mean_count);
                                        counts++;
                                    }
                                    br_gene_counts_2.close();

                                    for (int i = 0; i < position_counts.size(); i++) {
                                        gene_counts.add(position_counts.get(i));
                                    }
                                }

                                ArrayList<String> ensg_numbers = new ArrayList<>();
                                BufferedReader br_ensg_numbers =
                                        new BufferedReader(new FileReader(f_gene_annot_deseq2));
                                String line_ensg_numbers = br_ensg_numbers.readLine();
                                while ((line_ensg_numbers = br_ensg_numbers.readLine()) != null) {
                                    ensg_numbers.add(line_ensg_numbers.toUpperCase());
                                }
                                br_ensg_numbers.close();

                                int total_number_samples = 0;
                                HashMap<String, Integer> ensg_counts = new HashMap<>();
                                for (int i = 0; i < gene_counts.size(); i++) {
                                    ensg_counts.put(ensg_numbers.get(i).toUpperCase(), gene_counts.get(i));
                                    total_number_samples++;
                                }

                                HashMap<String, String> ensg_symb = new HashMap<>();

                                BufferedReader br_ensg_symb = new BufferedReader(new FileReader(f_ensg_map_symbol));
                                String line_ensg_symb = br_ensg_symb.readLine();
                                while ((line_ensg_symb = br_ensg_symb.readLine()) != null) {
                                    String[] split = line_ensg_symb.split("\t");
                                    if (split.length > 1) {
                                        ensg_symb.put(split[1].toUpperCase(), split[0]);
                                    }
                                }
                                br_ensg_symb.close();

                                HashMap<String, Integer> gene_symbol_lengths = new HashMap<>();
                                BufferedReader br_gene_symbol_lengths =
                                        new BufferedReader(new FileReader(f_ref_genome));
                                String line_gene_symbol_lengths = "";
                                while ((line_gene_symbol_lengths = br_gene_symbol_lengths.readLine()) != null) {
                                    if (line_gene_symbol_lengths.startsWith("#")) {
                                        continue;
                                    }
                                    String[] split = line_gene_symbol_lengths.split("\t");
                                    if (split[2].equals("transcript")) {
                                        int start = Integer.parseInt(split[3]);
                                        int end = Integer.parseInt(split[4]);
                                        int diff = end - start + 1;
                                        String[] split_further = split[8].split(" ");
                                        String ensg_name = split_further[1].replace("\"", "");
                                        ensg_name = ensg_name.replace(";", "");
                                        String[] ensg_name_split = ensg_name.split("\\.");
                                        gene_symbol_lengths.put(ensg_name_split[0], diff);
                                    }
                                }
                                br_gene_symbol_lengths.close();

                                double all_rpk = 0.0;

                                for (String ec : ensg_counts.keySet()) {
                                    if (!ensg_counts.containsKey(ec)) {
                                        continue;
                                    }

                                    int ec_count = ensg_counts.get(ec);
                                    int ec_lengths = 0;
                                    if (!gene_symbol_lengths.containsKey(ec)) {
                                        continue;
                                    }
                                    ec_lengths = gene_symbol_lengths.get(ec);
                                    if (ec_lengths == 0) {
                                        continue;
                                    }
                                    double norm_ec_lengths = ec_lengths / (1000 * 1.0);
                                    double current_rpk = ec_count / (norm_ec_lengths * 1.0);
                                    all_rpk += current_rpk;
                                }

                                double scaling_factor = all_rpk / (1000000 * 1.0);

                                for (String chr : current_tfs_to_regions.keySet()) {
                                    ArrayList<ENSG_ranges_binary_trees> tf_reg = current_tfs_to_regions.get(chr);
                                    ArrayList<ENSG_ranges_binary_trees> tf_reg_temp = new ArrayList<>();

                                    for (ENSG_ranges_binary_trees iu : tf_reg) {
                                        boolean should_write = false;

                                        //check if we want that TF
                                        for (String ec : iu.ensgs) {
                                            ec = ec.toUpperCase();

                                            if (ensg_symb.containsKey(ec)) {
                                                ec = ensg_symb.get(ec);
                                            } else {
                                                continue;
                                            }

                                            if (!ensg_counts.containsKey(ec)) {
                                                continue;
                                            }

                                            int ec_count = ensg_counts.get(ec);
                                            int ec_lengths = 0;
                                            if (!gene_symbol_lengths.containsKey(ec)) {
                                                continue;
                                            }
                                            ec_lengths = gene_symbol_lengths.get(ec);
                                            if (ec_lengths == 0) {
                                                continue;
                                            }
                                            double norm_ec_lengths = ec_lengths / (1000 * 1.0);
                                            double current_rpk = ec_count / (norm_ec_lengths * 1.0);

                                            if (current_rpk >= scaling_factor) {
                                                should_write = true;
                                            }
                                        }

                                        if (should_write) {
                                            tf_reg_temp.add(iu);
                                        }
                                    }

                                    current_tfs_to_regions.put(chr, tf_reg_temp);
                                }
                            }

                            file_tfs_to_regions.put(file_name, current_tfs_to_regions);
                        }


                        if (!options_intern.mix_mutually_exclusive) {
                            BufferedWriter bw = new BufferedWriter(new FileWriter(
                                    new File(output_tp_hm.getAbsolutePath() + File.separator + name_output + ".txt")));
                            bw.write("CHR\tLEFT_BORDER\tRIGHT_BORDER\tENSG");
                            bw.newLine();
                            for (String chr : tfs_to_regions.keySet()) {
                                ArrayList<ENSG_ranges_binary_trees> tf_reg = tfs_to_regions.get(chr);
                                for (ENSG_ranges_binary_trees iu : tf_reg) {
                                    bw.write(iu.toString());
                                    bw.newLine();
                                }

                            }
                            bw.close();
                        } else {
                            for (String file_name : file_tfs_to_regions.keySet()) {
                                HashMap<String, ArrayList<ENSG_ranges_binary_trees>> current_tfs_to_regions =
                                        file_tfs_to_regions.get(file_name);

                                name_output = fileDirTP_HM.getName() + "_" + file_name.split("_")[0];

                                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(
                                        output_tp_hm.getAbsolutePath() + File.separator + name_output + ".txt")));
                                bw.write("CHR\tLEFT_BORDER\tRIGHT_BORDER\tENSG");
                                bw.newLine();
                                for (String chr : current_tfs_to_regions.keySet()) {
                                    ArrayList<ENSG_ranges_binary_trees> tf_reg = current_tfs_to_regions.get(chr);
                                    for (ENSG_ranges_binary_trees iu : tf_reg) {
                                        bw.write(iu.toString());
                                        bw.newLine();
                                    }

                                }
                                bw.close();


                            }

                            //create tgen groups
                            File folder_output_groups = new File(
                                    options_intern.com2pose_working_directory + File.separator +
                                            options_intern.folder_name_tgen + File.separator +
                                            options_intern.folder_name_tgen_groups);
                            folder_output_groups.mkdir();

                            File folder_output_groups_hm = new File(
                                    folder_output_groups.getAbsolutePath() + File.separator + fileDirTP_HM.getName() +
                                            File.separator + fileDirTP.getName());
                            folder_output_groups_hm.mkdirs();

                            logger.logLine("[TGENE] CREATE TGENE Group Clashes.");

                            HashMap<String, ArrayList<ENSG_ranges_binary_trees>> updated_tfs_to_region_list =
                                    new HashMap<>();

                            for (String chr : tfs_to_regions.keySet()) {
                                ArrayList<ENSG_ranges_binary_trees> old_chr_tfs_to_region_list =
                                        tfs_to_regions.get(chr);
                                ArrayList<ENSG_ranges_binary_trees> updated_chr_tfs_to_region_list = new ArrayList<>();

                                //merge same ranges
                                ENSG_ranges_binary_trees iu = old_chr_tfs_to_region_list.get(0);
                                for (int i = 1; i < old_chr_tfs_to_region_list.size(); i++) {
                                    ENSG_ranges_binary_trees current_iu = old_chr_tfs_to_region_list.get(i);

                                    if (iu.isSameRange(current_iu)) {
                                        iu.ensgs.addAll(current_iu.ensgs);
                                    } else {
                                        updated_chr_tfs_to_region_list.add(iu);
                                        iu = current_iu;
                                    }
                                }
                                updated_tfs_to_region_list.put(chr, updated_chr_tfs_to_region_list);
                            }


                            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(
                                    folder_output_groups_hm.getAbsolutePath() + File.separator +
                                            options_intern.file_suffic_tgen_output_groups)));
                            bw.write("CHR\tLEFT_BORDER\tRIGHT_BORDER\tENSGS");
                            bw.newLine();
                            for (String chr : updated_tfs_to_region_list.keySet()) {
                                ArrayList<ENSG_ranges_binary_trees> tf_reg = updated_tfs_to_region_list.get(chr);
                                for (ENSG_ranges_binary_trees iu : tf_reg) {
                                    bw.write(iu.toString());
                                    bw.newLine();
                                }

                            }
                            bw.close();
                        }

                    }
                }
            }
        }
        logger.logLine("[TGENE] Finished merging TGene results");
    }

    /**
     * run tgene for all samples
     */
    public void run_tgen() throws Exception {
        logger.logLine("[TGENE] Run TGene");

        check_tepic_input_with_options();
        /*
        if(options_intern.mix_option.equals("SAMPLE_LEVEL"))
        {
            File root_mix_working_dir = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option);
            File f_sample_mix_output = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_sample_mix);
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory=f_sample_mix_output.getAbsolutePath();
        }

        if(options_intern.mix_option.equals("HM_LEVEL"))
        {
            File root_mix_working_dir = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option);
            File f_output_hm = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_hm_mix);
            f_output_hm.mkdir();

            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory=f_output_hm.getAbsolutePath();
        }
        if(!options_intern.black_list_dir.equals(""))
        {
            File output_folder = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_blacklisted_regions);
            File output_folder_new_input = new File(output_folder.getAbsolutePath()+File.separator+options_intern.folder_name_blacklisted_regions_new_input);
            output_folder_new_input.mkdir();

            //set new folder directory for tepic input and save old one
            options_intern.tepic_input_prev=options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = output_folder_new_input.getAbsolutePath();
        }
        if(options_intern.mix_mutually_exclusive)
        {
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option+File.separator+options_intern.folder_name_mix_option_mutually_exclusive+File.separator+options_intern.folder_name_mix_options_mutually_exclusive_input;
        }*/

        logger.logLine("[TGENE] Used data: " + options_intern.tepic_input_directory);


        String command_base = options_intern.path_tgen + File.separator + "bin" + File.separator + "tgene";

        File output_tgen = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_tgen +
                        File.separator + options_intern.folder_name_tgen_output);
        output_tgen.mkdir();

        File folder = new File(options_intern.tepic_input_directory);
        for (File fileDirTP : folder.listFiles()) {
            if (fileDirTP.isDirectory()) {
                File output_tgen_tp = new File(output_tgen.getAbsolutePath() + File.separator + fileDirTP.getName());
                output_tgen_tp.mkdir();

                for (File fileDirTP_HM : fileDirTP.listFiles()) {
                    if (fileDirTP_HM.isDirectory()) {
                        File output_tgen_tp_hm =
                                new File(output_tgen_tp.getAbsolutePath() + File.separator + fileDirTP_HM.getName());
                        output_tgen_tp_hm.mkdir();

                        for (File fileDir_data : fileDirTP_HM.listFiles()) {
                            if (!fileDir_data.isDirectory()) {
                                String[] name_split = fileDir_data.getName().split("\\.");

                                File output_tgen_tp_hm_sample =
                                        new File(output_tgen_tp_hm + File.separator + name_split[0]);
                                output_tgen_tp_hm_sample.mkdir();

                                String command_execute = new String(command_base);
                                command_execute += " " + fileDir_data.getAbsolutePath();

                                File gtf_dir = new File(options_intern.com2pose_working_directory + File.separator +
                                        options_intern.folder_name_tgen + File.separator +
                                        options_intern.folder_name_tgen_preprocessing + File.separator +
                                        options_intern.folder_name_tgen_preprocessing_gtf);
                                File gtf = new File("");

                                for (File gtf_dirDir : gtf_dir.listFiles()) {
                                    if (gtf_dirDir.getName()
                                            .matches(".*" + options_intern.file_suffix_tgen_preprocess_gtf)) {
                                        gtf = gtf_dirDir;
                                        break;
                                    }
                                }

                                if (gtf.getName().equals("")) {
                                    logger.logLine("[TGENE] could not find preprocessed GTF.");
                                    System.exit(1);
                                }

                                command_execute += " " + gtf.getAbsolutePath();

                                command_execute += " -oc " + output_tgen_tp_hm_sample;

                                if (options_intern.tgen_no_closest_locus) {
                                    command_execute += " --no-closest-locus";
                                }
                                if (options_intern.tgen_no_closest_tss) {
                                    command_execute += " --no-closest-tss";
                                }

                                command_execute += " --max-link-distances " + options_intern.tgen_max_link_distances;
                                command_execute += " --max-pvalue " + options_intern.tgen_pvalue;

                                //now execute TGENE:
                                logger.logLine("[TGENE] execute TGENE with command line: " + command_execute);
                                Process child = Runtime.getRuntime().exec(command_execute);
                                int code = child.waitFor();
                                switch (code) {
                                    case 0:
                                        break;
                                    case 1:
                                        String message = child.getErrorStream().toString();
                                        logger.logLine(
                                                "[TGENE] please check chromosome names in peak file and gtf file - they need to be named the same!");
                                        throw new Exception(message);
                                }
                            }
                        }
                    }
                }
            }
        }
        logger.logLine("[TGENE] Finished TGene");
    }

    /**
     * preprocess data for tgen run
     */
    public void preprocess_tgen() throws IOException {
        logger.logLine("[TGENE] Consensus approach with TGen is used. Preprocessing ...");

        check_tepic_input_with_options();
        /*
        if(options_intern.mix_option.equals("SAMPLE_LEVEL"))
        {
            File root_mix_working_dir = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option);
            File f_sample_mix_output = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_sample_mix);
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory=f_sample_mix_output.getAbsolutePath();
        }

        if(options_intern.mix_option.equals("HM_LEVEL"))
        {
            File root_mix_working_dir = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option);
            File f_output_hm = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_hm_mix);
            f_output_hm.mkdir();

            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory=f_output_hm.getAbsolutePath();
        }
        if(options_intern.tepic_tf_binding_site_search.equals("BETWEEN"))
        {
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option + File.separator+options_intern.folder_name_mix_options_footprints_between_peaks;

        }
        if(!options_intern.black_list_dir.equals(""))
        {
            File output_folder = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_blacklisted_regions);
            File output_folder_new_input = new File(output_folder.getAbsolutePath()+File.separator+options_intern.folder_name_blacklisted_regions_new_input);
            output_folder_new_input.mkdir();

            //set new folder directory for tepic input and save old one
            options_intern.tepic_input_prev=options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = output_folder_new_input.getAbsolutePath();
        }

        if(options_intern.mix_mutually_exclusive)
        {
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option+File.separator+options_intern.folder_name_mix_option_mutually_exclusive+File.separator+options_intern.folder_name_mix_options_mutually_exclusive_input;
        }*/

        logger.logLine("[TGENE] Used data: " + options_intern.tepic_input_directory);

        //create necessary folders for preprocessing
        File f_TGEN =
                new File(options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_tgen);
        f_TGEN.mkdir();
        File f_TGEN_preprocess =
                new File(f_TGEN.getAbsolutePath() + File.separator + options_intern.folder_name_tgen_preprocessing);
        f_TGEN_preprocess.mkdir();
        File f_TGEN_preprocess_gtf = new File(f_TGEN_preprocess.getAbsolutePath() + File.separator +
                options_intern.folder_name_tgen_preprocessing_gtf);
        f_TGEN_preprocess_gtf.mkdir();
        File f_TGEN_preprocess_binary_trees = new File(f_TGEN_preprocess.getAbsolutePath() + File.separator +
                options_intern.folder_name_tgen_preprocessing_binary_trees);
        f_TGEN_preprocess_binary_trees.mkdir();
        File f_TGEN_preprocess_binary_trees_unmerged = new File(
                f_TGEN_preprocess_binary_trees.getAbsolutePath() + File.separator +
                        options_intern.folder_name_tgen_preprocessing_binary_trees_unmerged);
        f_TGEN_preprocess_binary_trees_unmerged.mkdir();
        File f_TGEN_preprocess_binary_trees_merged = new File(
                f_TGEN_preprocess_binary_trees.getAbsolutePath() + File.separator +
                        options_intern.folder_name_tgen_preprocessing_binary_trees_merged);
        f_TGEN_preprocess_binary_trees_merged.mkdir();
        File f_TGEN_preprocess_binary_trees_sorted = new File(
                f_TGEN_preprocess_binary_trees.getAbsolutePath() + File.separator +
                        options_intern.folder_name_tgen_preprocessing_binary_trees_sorted);
        f_TGEN_preprocess_binary_trees_sorted.mkdir();

        logger.logLine("[TGENE] Preprocessing GTF.");
        //restructure GTF for TGEN -> only transcript data positions are allowed

        String[] gtf_name_split_dir = options_intern.tepic_gene_annot.split(File.separator);
        String[] gtf_name_split = gtf_name_split_dir[gtf_name_split_dir.length - 1].split("\\.");
        String gtf_name = "";
        for (int i = 0; i < gtf_name_split.length - 1; i++) {
            if (i > 0) {
                gtf_name += ".";
                gtf_name += gtf_name_split[i];
            } else {
                gtf_name += gtf_name_split[i];

            }
        }

        BufferedReader br_gtf_transcripts =
                new BufferedReader(new FileReader(new File(options_intern.tepic_gene_annot)));
        BufferedWriter bw_gtf_transcripts = new BufferedWriter(new FileWriter(new File(
                f_TGEN_preprocess_gtf.getAbsolutePath() + File.separator + gtf_name +
                        options_intern.file_suffix_tgen_preprocess_gtf)));

        String line_gtf_transcripts = "";
        while ((line_gtf_transcripts = br_gtf_transcripts.readLine()) != null) {
            if (line_gtf_transcripts.startsWith("#")) {
                continue;
            }
            String[] split = line_gtf_transcripts.split("\t");
            if (split[2].equals("transcript")) {
                if (split[0].matches(".*M.*")) {
                    String line = options_intern.tgen_mt_writing;
                    for (int i = 1; i < split.length; i++) {
                        line += "\t" + split[i];
                    }
                    bw_gtf_transcripts.write(line);
                    bw_gtf_transcripts.newLine();
                } else {
                    if (split[0].matches("chr.*")) {
                        bw_gtf_transcripts.write(line_gtf_transcripts.substring(3));
                        bw_gtf_transcripts.newLine();
                    } else {
                        bw_gtf_transcripts.write(line_gtf_transcripts);
                        bw_gtf_transcripts.newLine();
                    }
                }
            }

        }
        bw_gtf_transcripts.close();
        br_gtf_transcripts.close();

        logger.logLine("[TGENE] Preprocessing Binary Trees for position to Gen Mapping.");

        BufferedReader br_gtf_bin_tree = new BufferedReader(new FileReader(new File(options_intern.tepic_gene_annot)));
        BufferedWriter bw_gtf_bin_tree = new BufferedWriter(new FileWriter(new File(
                f_TGEN_preprocess_binary_trees_unmerged.getAbsolutePath() + File.separator + "test" + ".txt")));

        String current_chr = "";
        int count = 0;

        String line_gtf_bin_tree = "";
        while ((line_gtf_bin_tree = br_gtf_bin_tree.readLine()) != null) {
            if (line_gtf_bin_tree.startsWith("#")) {
                continue;
            }
            String[] split = line_gtf_bin_tree.split("\t");
            if (split[2].equals("transcript")) {
                String[] split_line = split[8].split(";");

                String chr = split[0];

                if (!current_chr.equals(chr)) {
                    bw_gtf_bin_tree.close();
                    count = 0;
                    bw_gtf_bin_tree = new BufferedWriter(new FileWriter(new File(
                            f_TGEN_preprocess_binary_trees_unmerged.getAbsolutePath() + File.separator + chr +
                                    ".txt")));
                    bw_gtf_bin_tree.write("#\tPOS_START\tPOS_END\tENSG\n");
                    current_chr = chr;
                }

                String ensg = "";
                String start_pos = split[3];
                String end_pos = split[4];

                for (int i = 0; i < split_line.length; i++) {
                    if (split_line[i].startsWith("gene_id")) {
                        String[] split_x = split_line[i].split("\"");
                        ensg = split_x[1].substring(0, split_x[1].length()).split("\\.")[0];
                        break;
                    }
                }

                StringBuilder sb = new StringBuilder();
                sb.append(count);
                sb.append("\t");
                sb.append(start_pos);
                sb.append("\t");
                sb.append(end_pos);
                sb.append("\t");
                sb.append(ensg);

                bw_gtf_bin_tree.write(sb.toString());
                bw_gtf_bin_tree.newLine();

                count++;
            }

        }
        bw_gtf_bin_tree.close();
        br_gtf_bin_tree.close();

        //binary tree merge overlappings

        ArrayList<ENSG_ranges_binary_trees> unmerged_intervalls = new ArrayList<>();

        for (File fileDir : f_TGEN_preprocess_binary_trees_unmerged.listFiles()) {
            if (!fileDir.isDirectory() && !fileDir.getName().equals("test.txt")) {
                BufferedReader br = new BufferedReader(new FileReader(fileDir));
                String line = br.readLine();
                while ((line = br.readLine()) != null) {
                    String[] split = line.split("\t");

                    ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                    iu.number = Integer.parseInt(split[0]);
                    iu.left_border = Integer.parseInt(split[1]);
                    iu.right_border = Integer.parseInt(split[2]);
                    iu.ensgs.add(split[3]);

                    unmerged_intervalls.add(iu);
                }
                br.close();

                Collections.sort(unmerged_intervalls);

                ArrayList<ENSG_ranges_binary_trees> unmerged_intervalls_temp_list = new ArrayList<>();
                //merge intervalls with same ENSG list!
                int count_same_ENSG = 0;
                boolean last_one_merged = false;
                for (int i = 1; i < unmerged_intervalls.size(); i++) {
                    ENSG_ranges_binary_trees iu_before = unmerged_intervalls.get(i - 1);
                    ENSG_ranges_binary_trees iu_now = unmerged_intervalls.get(i);

                    HashSet<String> ensgs_before = iu_before.ensgs;
                    HashSet<String> ensgs_now = iu_now.ensgs;

                    boolean all_in = true;

                    for (String s : ensgs_before) {
                        if (!ensgs_now.contains(s)) {
                            all_in = false;
                        }
                    }
                    for (String s : ensgs_now) {
                        if (!ensgs_before.contains(s)) {
                            all_in = false;
                        }
                    }

                    if (all_in) {
                        //now merge!

                        //look for further intervalls with same ensg
                        boolean found_furthers = true;
                        int counter_further = 0;
                        int temp_i = i;
                        while (found_furthers && temp_i < unmerged_intervalls.size()) {
                            found_furthers = false;

                            HashSet<String> ensgs_further = unmerged_intervalls.get(temp_i).ensgs;


                            boolean all_in_further = true;

                            for (String s : ensgs_further) {
                                if (!ensgs_before.contains(s)) {
                                    all_in_further = false;
                                }
                            }

                            for (String s : ensgs_before) {
                                if (!ensgs_further.contains(s)) {
                                    all_in_further = false;
                                }
                            }

                            if (all_in_further) {
                                found_furthers = true;
                                iu_now = unmerged_intervalls.get(temp_i);
                            } else {
                                found_furthers = false;
                                break;
                            }

                            if (all_in_further && temp_i == unmerged_intervalls.size() - 1) {
                                last_one_merged = true;
                            }


                            counter_further++;
                            temp_i++;
                        }
                        i += counter_further + 1;


                        ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                        iu.number = count_same_ENSG;
                        iu.ensgs = new HashSet<>(ensgs_before);
                        iu.left_border = iu_before.left_border;
                        iu.right_border = iu_now.right_border;

                        unmerged_intervalls_temp_list.add(iu);

                        count_same_ENSG++;
                    } else {
                        iu_before.number = count_same_ENSG;
                        unmerged_intervalls_temp_list.add(iu_before);
                        count_same_ENSG++;
                    }
                }

                if (!last_one_merged) {
                    ENSG_ranges_binary_trees iu = unmerged_intervalls.get(unmerged_intervalls.size() - 1);
                    iu.number = unmerged_intervalls_temp_list.size();
                    unmerged_intervalls_temp_list.add(iu);
                }

                unmerged_intervalls.clear();
                unmerged_intervalls = new ArrayList<>(unmerged_intervalls_temp_list);
                unmerged_intervalls_temp_list.clear();

                //re calculate borders

                unmerged_intervalls.get(0).left_border = 0;

                for (int i = 1; i < unmerged_intervalls.size(); i++) {
                    ENSG_ranges_binary_trees iu_before = unmerged_intervalls.get(i - 1);
                    ENSG_ranges_binary_trees iu_current = unmerged_intervalls.get(i);

                    iu_current.left_border = iu_before.right_border + 1;
                }

                Collections.sort(unmerged_intervalls);

                for (int i = 0; i < unmerged_intervalls.size(); i++) {
                    unmerged_intervalls.get(i).number = i;
                }

                BufferedWriter bw = new BufferedWriter(new FileWriter(
                        new File(f_TGEN_preprocess_binary_trees_merged + File.separator + fileDir.getName())));
                bw.write("#\tPOS_START\tPOS_END\tENSG\n");
                for (int i = 1; i < unmerged_intervalls.size(); i++) {
                    ENSG_ranges_binary_trees iu = unmerged_intervalls.get(i);
                    StringBuilder sb = new StringBuilder();
                    sb.append(iu.number);
                    sb.append("\t");
                    sb.append(iu.left_border);
                    sb.append("\t");
                    sb.append(iu.right_border);
                    sb.append("\t");
                    sb.append(iu.ensgs_to_String());

                    bw.write(sb.toString());
                    bw.newLine();

                }
                bw.close();
            }
        }


        //sort for binary tree implementation

        logger.logLine("[TGENE] Preparing chromosomes for binary search.");

        for (File fileDir : f_TGEN_preprocess_binary_trees_merged.listFiles()) {
            if (!fileDir.isDirectory()) {

                ArrayList<ENSG_ranges_binary_trees> regions = new ArrayList<>();

                BufferedReader br = new BufferedReader(new FileReader(fileDir));
                String header = br.readLine();
                String line = "";

                while ((line = br.readLine()) != null) {
                    String[] split = line.split("\t");
                    ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                    iu.number = Integer.parseInt(split[0]);
                    iu.left_border = Integer.parseInt(split[1]);
                    iu.right_border = Integer.parseInt(split[2]);
                    iu.ensgs.addAll(Arrays.asList(split[3].split(";")));

                    regions.add(iu);
                }
                br.close();

                ArrayList<ENSG_ranges_binary_trees> newly_ordered = new ArrayList<>();

                ENSG_ranges_binary_trees median = regions.get(regions.size() / 2);
                newly_ordered.add(median);

                ArrayList<ENSG_ranges_binary_trees> region_left = new ArrayList<>();
                for (int i = 0; i < regions.size() / 2; i++) {
                    region_left.add(regions.get(i));
                }

                ArrayList<ENSG_ranges_binary_trees> region_right = new ArrayList<>();
                for (int i = regions.size() / 2 + 1; i < regions.size(); i++) {
                    region_right.add(regions.get(i));
                }

                newly_ordered = recursive_split(region_left, newly_ordered);
                newly_ordered = recursive_split(region_right, newly_ordered);

                BufferedWriter bw = new BufferedWriter(new FileWriter(
                        f_TGEN_preprocess_binary_trees_sorted.getAbsolutePath() + File.separator + fileDir.getName()));
                bw.write(header);
                bw.newLine();

                for (int i = 0; i < newly_ordered.size(); i++) {
                    bw.write(newly_ordered.get(i).toString_binary());
                    bw.newLine();
                }

                bw.close();
            }
        }


        logger.logLine("[TGENE] Finished preprocessing.");
    }


    /**
     * postprocesses the TEPIC output (checks for TPM filter and copies files into a structure where preprocessing of DYNAMITE can happen
     */
    public void postprocess_tepic_output() throws Exception {
        logger.logLine("Start postprocessing of TEPIC output");

        HashMap<String, HashMap<String, HashSet<String>>> groups_to_compare = checkGroupsTEPIC();

        String suffix = "";

        if (options_intern.tepic_tpm_cutoff > 0) {
            suffix = "_Gene_View_Filtered_TPM.txt";
        } else {
            suffix = "_Gene_View_Filtered.txt";
        }

        File folder_postprocessing = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_tepic_postprocessing);
        folder_postprocessing.mkdir();

        File folder_pp_input = new File(folder_postprocessing.getAbsolutePath() + File.separator +
                options_intern.folder_name_tepic_postprocessing_input);
        folder_pp_input.mkdir();

        File folder_pp_output = new File(folder_postprocessing.getAbsolutePath() + File.separator +
                options_intern.folder_name_tepic_postprocessing_output);
        folder_pp_output.mkdir();

        HashSet<String> available_hms = new HashSet<>();

        HashSet<String> check_tfs = new HashSet<>();
        String seperator = "";

        //generate input structure and copy files
        for (String s : groups_to_compare.keySet()) {
            File f_input_tp_folders = new File(folder_pp_input.getAbsolutePath() + File.separator + s);
            f_input_tp_folders.mkdir();
            HashMap<String, HashSet<String>> hm = groups_to_compare.get(s);
            for (String ss : hm.keySet()) {
                available_hms.add(ss);
                File f_input_hm_folders = new File(f_input_tp_folders.getAbsolutePath() + File.separator + ss);
                f_input_hm_folders.mkdir();

                //move coresponding sample outputs to these folders
                File f_input_samples;
                if (options_intern.tepic_randomize_tf_gene_matrix) {
                    f_input_samples = new File(options_intern.com2pose_working_directory + File.separator +
                            options_intern.folder_name_tepic_output_raw_shuffle + File.separator + File.separator + s +
                            File.separator + ss);
                } else {
                    f_input_samples = new File(options_intern.com2pose_working_directory + File.separator +
                            options_intern.folder_name_tepic_output_raw + File.separator + s + File.separator + ss);
                }

                for (File fileDir : f_input_samples.listFiles()) {
                    if (fileDir.isDirectory()) {
                        for (File fileDir2 : fileDir.listFiles()) {
                            if (!fileDir2.isDirectory()) {
                                String name = fileDir2.getName();
                                if (name.matches(".*" + suffix)) {
                                    //check tfs
                                    BufferedReader br = new BufferedReader(new FileReader(fileDir2));
                                    String line = br.readLine();
                                    String[] split = line.split("\t");
                                    for (String st : split) {
                                        if (st.matches(".*::.*") || st.matches(".*[.][.].*")) {
                                            check_tfs.add(st);
                                            if (st.matches(".*::.*")) {
                                                seperator = "::";
                                            } else {
                                                seperator = "[.][.]";
                                            }
                                        }
                                    }

                                    br.close();

                                    //COPY!!
                                    String command = "cp -u " + fileDir2.getAbsolutePath() + " " +
                                            folder_pp_input.getAbsolutePath() + File.separator + s + File.separator +
                                            ss;
                                    Process child = Runtime.getRuntime().exec(command);
                                    logger.logLine("[TEPIC] Copy files: " + command);
                                    int code = child.waitFor();
                                    switch (code) {
                                        case 0:
                                            break;
                                        case 1:
                                            String message = child.getErrorStream().toString();
                                            throw new Exception(message);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        File f_out_tfs = new File(folder_postprocessing.getAbsolutePath() + File.separator +
                options_intern.folder_name_tepic_postprocessing_tfs);
        f_out_tfs.mkdir();
        BufferedWriter bw_check_tfs = new BufferedWriter(new FileWriter(f_out_tfs.getAbsolutePath() + File.separator +
                options_intern.file_suffix_tepic_postprocessing_tfs_tfs));
        for (String key_tfs : check_tfs) {
            String[] split = key_tfs.split(seperator);

            StringBuilder sb = new StringBuilder();

            for (int i = 0; i < split.length; i++) {
                if (i < 1) {
                    sb.append(split[i].toUpperCase());
                } else {
                    sb.append("..");
                    sb.append(split[i].toUpperCase());
                }
            }

            for (int i = 0; i < split.length; i++) {
                sb.append("\t");
                sb.append(split[i].toUpperCase());
            }
            bw_check_tfs.write(sb.toString());
            bw_check_tfs.newLine();
        }
        bw_check_tfs.close();

        //generate output structure
        //HMs
        HashMap<String, File> pp_output_hms_files = new HashMap<>();
        for (String s : available_hms) {
            File f = new File(folder_pp_output.getAbsolutePath() + File.separator + s);
            f.mkdir();
            pp_output_hms_files.put(s, f);
        }

        //identify groups to compute mean ratios
        //group1_vs_group2
        HashMap<String, File> pp_output_clashedGroups = new HashMap<>();
        HashSet<String> already_checked_groups = new HashSet<>();
        for (String s : groups_to_compare.keySet()) {
            for (String ss : groups_to_compare.keySet()) {
                if (s.equals(ss)) {
                    continue;
                }
                String key1 = s + "_" + ss;
                String key2 = ss + "_" + s;
                if (already_checked_groups.contains(key1) || already_checked_groups.contains(key2)) {
                    continue;
                }
                already_checked_groups.add(key1);

                HashMap<String, HashSet<String>> group1_hms = groups_to_compare.get(s);
                HashMap<String, HashSet<String>> group2_hms = groups_to_compare.get(ss);

                for (String k : group1_hms.keySet()) {
                    if (group2_hms.containsKey(k)) {
                        if (!options_intern.mix_mutually_exclusive) {
                            File f = new File(pp_output_hms_files.get(k).getAbsolutePath() + File.separator + key1);
                            f.mkdir();
                        }
                    }
                }
            }
        }

        if (options_intern.mix_mutually_exclusive) {
            for (String key_groupClash : groups_to_compare.keySet()) {
                HashMap<String, HashSet<String>> hms = groups_to_compare.get(key_groupClash);

                for (String hm : hms.keySet()) {
                    File f = new File(
                            folder_pp_output.getAbsolutePath() + File.separator + hm + File.separator + key_groupClash);
                    f.mkdir();
                }
            }
        }

        //prepare command line for all (computeMeanRatioTFAffinities.py
        String command_base = "python3 " + options_intern.path_to_COM2POSE + File.separator +
                options_intern.directory_for_tepic_DYNAMITE + File.separator + "computeMeanRatioTFAffinities.py";

        HashSet<String> all_tfs = new HashSet<>();

        for (File fileDir : folder_pp_output.listFiles()) {
            if (fileDir.isDirectory()) {
                for (File fileDir2 : fileDir.listFiles()) {
                    if (fileDir2.isDirectory()) {
                        String current_HM = fileDir.getName();
                        String current_group_clash = fileDir2.getName();
                        logger.logLine("[TEPIC] Postprocess: " + current_HM + " - " + current_group_clash);
                        String[] group_clash_split = current_group_clash.split("_");

                        String group1 = group_clash_split[0];
                        String group2 = group_clash_split[1];

                        if (group_clash_split.length > 2) {
                            int found_non = -1;

                            for (int i = 0; i < group_clash_split.length; i++) {
                                if (group_clash_split[i].equals("non") || group_clash_split[i].equals("Non")) {
                                    found_non = i;
                                }
                            }
                            //group1= group_clash_split[0].substring(0, 1).toUpperCase() + group_clash_split[0].substring(1);
                            //group2= group_clash_split[1].substring(0, 1).toUpperCase() + group_clash_split[1].substring(1);

                            if (found_non != -1) {
                                if (found_non == 0) {
                                    group1 = group_clash_split[0] + "_" + group_clash_split[1];
                                    group2 = group_clash_split[2];
                                }
                                if (found_non == 1) {
                                    group1 = group_clash_split[0];
                                    group2 = group_clash_split[1] + "_" + group_clash_split[2];
                                }
                            }
                        }


                        String group1_input_dir = "";
                        String group2_input_dir = "";

                        File folder_group1_check_tfs;
                        File folder_group2_check_tfs;

                        if (!options_intern.mix_mutually_exclusive) {
                            folder_group1_check_tfs = new File(
                                    folder_pp_input.getAbsolutePath() + File.separator + group1 + File.separator +
                                            current_HM);
                            folder_group2_check_tfs = new File(
                                    folder_pp_input.getAbsolutePath() + File.separator + group2 + File.separator +
                                            current_HM);
                        } else {
                            folder_group1_check_tfs = new File(
                                    folder_pp_input.getAbsolutePath() + File.separator + current_group_clash +
                                            File.separator + current_HM);
                            folder_group2_check_tfs = new File(
                                    folder_pp_input.getAbsolutePath() + File.separator + current_group_clash +
                                            File.separator + current_HM);
                        }

                        File[] samples_group1_check_tfs = folder_group1_check_tfs.listFiles();
                        File[] samples_group2_check_tfs = folder_group2_check_tfs.listFiles();

                        if (!options_intern.mix_mutually_exclusive) {
                            BufferedReader br_group1_check_tfs =
                                    new BufferedReader(new FileReader((samples_group1_check_tfs[0])));
                            String header_group1_check_tfs = br_group1_check_tfs.readLine();
                            all_tfs.addAll(Arrays.asList(header_group1_check_tfs.split("\t")));
                            br_group1_check_tfs.close();

                            BufferedReader br_group2_check_tfs =
                                    new BufferedReader(new FileReader((samples_group2_check_tfs[0])));
                            String header_group2_check_tfs = br_group2_check_tfs.readLine();
                            all_tfs.addAll(Arrays.asList(header_group2_check_tfs.split("\t")));
                            br_group2_check_tfs.close();
                        } else {
                            BufferedReader br_group1_check_tfs =
                                    new BufferedReader(new FileReader((samples_group1_check_tfs[0])));
                            String header_group1_check_tfs = br_group1_check_tfs.readLine();
                            all_tfs.addAll(Arrays.asList(header_group1_check_tfs.split("\t")));
                            br_group1_check_tfs.close();

                            BufferedReader br_group2_check_tfs =
                                    new BufferedReader(new FileReader((samples_group2_check_tfs[1])));
                            String header_group2_check_tfs = br_group2_check_tfs.readLine();
                            all_tfs.addAll(Arrays.asList(header_group2_check_tfs.split("\t")));
                            br_group2_check_tfs.close();
                        }


                        //if TPM filter was used we need to postprocess the TPM files. otherwise it will not work!
                        if (options_intern.tepic_tpm_cutoff > 0) {
                            logger.logLine("[TEPIC] TPM filter > 0, start postprocessing of TPM filtered scores");
                            File output_post_group1 = new File(
                                    folder_pp_output + File.separator + current_HM + File.separator +
                                            current_group_clash + File.separator + group1);
                            output_post_group1.mkdir();
                            File output_post_group2 = new File(
                                    folder_pp_output + File.separator + current_HM + File.separator +
                                            current_group_clash + File.separator + group2);
                            output_post_group2.mkdir();

                            //build intersect and write new files with filter

                            File folder_group1;
                            File folder_group2;
                            File[] samples_group1;
                            File[] samples_group2;

                            if (!options_intern.mix_mutually_exclusive) {
                                folder_group1 = new File(
                                        folder_pp_input.getAbsolutePath() + File.separator + group1 + File.separator +
                                                current_HM);
                                folder_group2 = new File(
                                        folder_pp_input.getAbsolutePath() + File.separator + group2 + File.separator +
                                                current_HM);

                                samples_group1 = folder_group1.listFiles();
                                samples_group2 = folder_group2.listFiles();
                            } else {
                                folder_group1 = new File(
                                        folder_pp_input.getAbsolutePath() + File.separator + current_group_clash +
                                                File.separator + current_HM);
                                folder_group2 = new File(
                                        folder_pp_input.getAbsolutePath() + File.separator + current_group_clash +
                                                File.separator + current_HM);
                                samples_group1 = new File[1];
                                samples_group1[0] = folder_group1.listFiles()[0];
                                samples_group2 = new File[1];
                                samples_group2[0] = folder_group2.listFiles()[1];
                            }

                            ArrayList<Boolean> write_group1 = new ArrayList<>();
                            ArrayList<Boolean> write_group2 = new ArrayList<>();

                            ArrayList<String> header_group1;
                            ArrayList<String> header_group2;
                            HashSet<String> header_group1_set;
                            HashSet<String> header_group2_set;

                            BufferedReader br_group1_check = new BufferedReader(new FileReader(samples_group1[0]));
                            String line_group1_check = br_group1_check.readLine();
                            header_group1 = new ArrayList<>(Arrays.asList(line_group1_check.split("\t")));
                            header_group1_set = new HashSet<>(Arrays.asList(line_group1_check.split("\t")));
                            br_group1_check.close();

                            BufferedReader br_group2_check = new BufferedReader(new FileReader(samples_group2[0]));
                            String line_group2_check = br_group2_check.readLine();
                            header_group2 = new ArrayList<>(Arrays.asList(line_group2_check.split("\t")));
                            header_group2_set = new HashSet<>(Arrays.asList(line_group2_check.split("\t")));
                            br_group2_check.close();

                            for (int i = 0; i < header_group1.size(); i++) {
                                if (header_group2_set.contains(header_group1.get(i))) {
                                    write_group1.add(true);
                                } else {
                                    write_group1.add(false);
                                }
                            }
                            for (int i = 0; i < header_group2.size(); i++) {
                                if (header_group1_set.contains(header_group2.get(i))) {
                                    write_group2.add(true);
                                } else {
                                    write_group2.add(false);
                                }
                            }

                            for (File f : folder_group1.listFiles()) {
                                if (options_intern.mix_mutually_exclusive) {
                                    f = folder_group1.listFiles()[0];
                                }

                                BufferedReader br_group1 = new BufferedReader(new FileReader(f));
                                BufferedWriter bw_group1 = new BufferedWriter(new FileWriter(
                                        new File(output_post_group1.getAbsolutePath() + File.separator + f.getName())));

                                String line_group1 = "";
                                while ((line_group1 = br_group1.readLine()) != null) {
                                    String[] split_line_group1 = line_group1.split("\t");
                                    StringBuilder sb = new StringBuilder();

                                    for (int i = 0; i < split_line_group1.length; i++) {
                                        if (write_group1.get(i)) {
                                            if (i > 0) {
                                                sb.append("\t");
                                                sb.append(split_line_group1[i]);
                                            } else {
                                                sb.append(split_line_group1[i]);
                                            }
                                        }
                                    }
                                    bw_group1.write(sb.toString());
                                    bw_group1.newLine();
                                }
                                bw_group1.close();
                                br_group1.close();
                                if (options_intern.mix_mutually_exclusive) {
                                    break;
                                }
                            }


                            for (File f : folder_group2.listFiles()) {
                                if (options_intern.mix_mutually_exclusive) {
                                    f = folder_group1.listFiles()[1];
                                }

                                BufferedReader br_group2 = new BufferedReader(new FileReader(f));
                                BufferedWriter bw_group2 = new BufferedWriter(new FileWriter(
                                        new File(output_post_group2.getAbsolutePath() + File.separator + f.getName())));

                                String line_group2 = "";
                                while ((line_group2 = br_group2.readLine()) != null) {
                                    String[] split_line_group2 = line_group2.split("\t");
                                    StringBuilder sb = new StringBuilder();

                                    for (int i = 0; i < split_line_group2.length; i++) {
                                        if (write_group2.get(i)) {
                                            if (i > 0) {
                                                sb.append("\t");
                                                sb.append(split_line_group2[i]);
                                            } else {
                                                sb.append(split_line_group2[i]);
                                            }
                                        }
                                    }
                                    bw_group2.write(sb.toString());
                                    bw_group2.newLine();
                                }
                                bw_group2.close();
                                br_group2.close();
                                if (options_intern.mix_mutually_exclusive) {
                                    break;
                                }
                            }
                            logger.logLine("[TEPIC] TPM filter > 0, end postprocessing of TPM filtered scores");

                            //set to postprocessed TMP filtered data
                            group1_input_dir = output_post_group1.getAbsolutePath();
                            group2_input_dir = output_post_group2.getAbsolutePath();
                        } else {
                            if (!options_intern.mix_mutually_exclusive) {
                                group1_input_dir =
                                        folder_pp_input.getAbsolutePath() + File.separator + group1 + File.separator +
                                                current_HM;
                                group2_input_dir =
                                        folder_pp_input.getAbsolutePath() + File.separator + group2 + File.separator +
                                                current_HM;
                            } else {
                                //NOTE: THIS IS NOT TESTED AND COULD CAUSE TROUBLE WITH NO TPM + MUTUALLY EXCLUSIVE OPTION
                                if (options_intern.mix_mutually_exclusive) {
                                    File f_group1_file = new File(
                                            folder_pp_input.getAbsolutePath() + File.separator + current_group_clash +
                                                    File.separator + current_HM + File.separator +
                                                    group_clash_split[0]);
                                    f_group1_file.mkdir();

                                    File f_group2_file = new File(
                                            folder_pp_input.getAbsolutePath() + File.separator + current_group_clash +
                                                    File.separator + current_HM + File.separator +
                                                    group_clash_split[1]);
                                    f_group1_file.mkdir();

                                    File f_group1_data = new File(
                                            folder_pp_input.getAbsolutePath() + File.separator + group1 +
                                                    File.separator + current_HM).listFiles()[0];
                                    File f_group2_data = new File(
                                            folder_pp_input.getAbsolutePath() + File.separator + group1 +
                                                    File.separator + current_HM).listFiles()[1];

                                    String command_edited = "mv " + f_group1_data.getAbsolutePath() + " " +
                                            f_group1_file.getAbsolutePath() + File.separator + f_group1_data.getName();

                                    logger.logLine(
                                            "[TEPIC] Moving file due to mutually exclusive option: " + command_edited);
                                    Process child = Runtime.getRuntime().exec(command_edited);
                                    int code = child.waitFor();
                                    switch (code) {
                                        case 0:
                                            break;
                                        case 1:
                                            String message = child.getErrorStream().toString();
                                            throw new Exception(message);
                                    }

                                    String command_edited_2 = "mv " + f_group2_data.getAbsolutePath() + " " +
                                            f_group2_file.getAbsolutePath() + File.separator + f_group2_data.getName();

                                    logger.logLine("[TEPIC] Moving file due to mutually exclusive option: " +
                                            command_edited_2);
                                    Process child2 = Runtime.getRuntime().exec(command_edited_2);
                                    int code2 = child2.waitFor();
                                    switch (code2) {
                                        case 0:
                                            break;
                                        case 1:
                                            String message = child2.getErrorStream().toString();
                                            throw new Exception(message);
                                    }

                                    group1_input_dir = f_group1_file.getAbsolutePath();
                                    group2_input_dir = f_group2_file.getAbsolutePath();
                                } else {
                                    group1_input_dir = folder_pp_input.getAbsolutePath() + File.separator + group1 +
                                            File.separator + current_HM;
                                    group2_input_dir = folder_pp_input.getAbsolutePath() + File.separator + group2 +
                                            File.separator + current_HM;
                                }
                            }
                        }


                        File output_mean_affinities = new File(
                                folder_pp_output + File.separator + current_HM + File.separator + current_group_clash +
                                        File.separator +
                                        options_intern.folder_name_tepic_postprocessing_output_mean_affinities);
                        output_mean_affinities.mkdir();
                        File output_ratios = new File(
                                folder_pp_output + File.separator + current_HM + File.separator + current_group_clash +
                                        File.separator + options_intern.folder_name_tepic_postprocessing_output_ratios);
                        output_ratios.mkdir();

                        File output_mean_affinities_group1 = new File(
                                output_mean_affinities.getAbsolutePath() + File.separator +
                                        options_intern.file_suffix_tepic_postprocessing_output_mean_affinities +
                                        group1 + ".txt");
                        File output_mean_affinities_group2 = new File(
                                output_mean_affinities.getAbsolutePath() + File.separator +
                                        options_intern.file_suffix_tepic_postprocessing_output_mean_affinities +
                                        group2 + ".txt");
                        File output_ratios_group1_group2 = new File(output_ratios.getAbsolutePath() + File.separator +
                                options_intern.file_suffix_tepic_postprocessing_output_ratios + group1 + "_" + group2 +
                                ".txt");

                        HashSet<File> files_to_create = new HashSet<>();
                        files_to_create.add(output_mean_affinities_group1);
                        files_to_create.add(output_mean_affinities_group2);
                        files_to_create.add(output_ratios_group1_group2);

                        for (File f : files_to_create) {
                            BufferedWriter bw = new BufferedWriter(new FileWriter(f));
                            bw.write("");
                            bw.close();
                        }


                        String command_edited = new String(command_base);

                        command_edited += " " + group1_input_dir + File.separator;
                        command_edited += " " + group2_input_dir + File.separator;
                        command_edited += " " + output_mean_affinities_group1.getAbsolutePath();
                        command_edited += " " + output_mean_affinities_group2.getAbsolutePath();
                        command_edited += " " + output_ratios_group1_group2.getAbsolutePath();

                        String command_tail = "";

                        if (options_intern.tepic_original_decay) {
                            command_tail += " True";
                        } else {
                            command_tail += " False";
                        }
                        if (options_intern.tepic_not_generated) {
                            command_tail += " True";
                        } else {
                            command_tail += " False";
                        }
                        if (options_intern.tepic_tpm_cutoff > 0) {
                            command_tail += " True";
                        } else {
                            command_tail += " False";
                        }

                        command_edited += command_tail;

                        logger.logLine(
                                "[TEPIC] execute computeMeanRatioTFAffinities.py with command line: " + command_edited);
                        Process child = Runtime.getRuntime().exec(command_edited);
                        int code = child.waitFor();
                        switch (code) {
                            case 0:
                                break;
                            case 1:
                                String message = child.getErrorStream().toString();
                                throw new Exception(message);
                        }
                    }
                }
            }
        }

        BufferedWriter bw_all_tfs = new BufferedWriter(new FileWriter(new File(
                folder_postprocessing.getAbsolutePath() + File.separator +
                        options_intern.file_suffix_tepic_postprocessing_all_tfs)));
        for (String k : all_tfs) {
            bw_all_tfs.write(k);
            bw_all_tfs.newLine();
        }
        bw_all_tfs.close();
        logger.logLine("Finished postprocessing of TEPIC output");
    }

    /**
     * creates violine plots of lengths
     *
     * @throws IOException
     */
    public void create_open_regions_violin_plots() throws Exception {
        logger.logLine("[TEPIC-ANALYSIS] Create Violine Plots of lengths of Peaks");

        //create R dataformat
        File f_output_root = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_tepic_postprocessing + File.separator +
                options_intern.folder_name_tepic_postprocessing_open_chromatin_violins);
        f_output_root.mkdir();
        File f_output_data = new File(f_output_root.getAbsolutePath() + File.separator +
                options_intern.folder_name_tepic_postprocessing_open_chromatin_violins_data);
        f_output_data.mkdir();
        File f_output_script = new File(f_output_root.getAbsolutePath() + File.separator +
                options_intern.folder_name_tepic_postprocessing_open_chromatin_violins_script);
        f_output_script.mkdir();
        File f_output_plot = new File(f_output_root.getAbsolutePath() + File.separator +
                options_intern.folder_name_tepic_postprocessing_open_chromatin_violins_plots);
        f_output_plot.mkdir();

        File f_input_root = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_tepic_output_raw);

        HashMap<String, ArrayList<Integer>> hm_to_open_chr_lengths = new HashMap<>();

        for (File f_tp : f_input_root.listFiles()) {
            if (f_tp.isDirectory()) {
                for (File f_hm : f_tp.listFiles()) {
                    if (f_hm.isDirectory()) {
                        String hm_name = f_hm.getName();
                        ArrayList<Integer> hm_lengths = new ArrayList<>();
                        if (hm_to_open_chr_lengths.containsKey(hm_name)) {
                            hm_lengths = hm_to_open_chr_lengths.get(hm_name);
                        }
                        hm_to_open_chr_lengths.put(hm_name, hm_lengths);

                        for (File f_sample : f_hm.listFiles()) {
                            if (f_sample.isDirectory()) {
                                for (File f_regions_to_targets : f_sample.listFiles()) {
                                    if (f_regions_to_targets.isFile()) {
                                        if (f_regions_to_targets.getName()
                                                .equals(options_intern.file_suffix_tepic_output_regions_to_target_genes)) {

                                            BufferedReader br_regions_to_targets =
                                                    new BufferedReader(new FileReader(f_regions_to_targets));
                                            String line_regions_to_targets = br_regions_to_targets.readLine();
                                            String[] header_regions_to_targets = line_regions_to_targets.split("\t");
                                            while ((line_regions_to_targets = br_regions_to_targets.readLine()) !=
                                                    null) {
                                                String[] split = line_regions_to_targets.split("\t");
                                                String[] split_chr = split[0].split(":");
                                                String[] split_position = split_chr[1].split("-");

                                                int start = Integer.parseInt(split_position[0]);
                                                int end = Integer.parseInt(split_position[1]);

                                                int length = end - start;

                                                hm_lengths.add(length);

                                            }
                                            br_regions_to_targets.close();
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        File f_data_file = new File(f_output_data.getAbsolutePath() + File.separator +
                options_intern.file_suffix_tepic_postprocessing_violin_plots_data);
        BufferedWriter bw_data_file = new BufferedWriter(new FileWriter(f_data_file));
        bw_data_file.write("HM\tLENGTH\n");
        for (String hm : hm_to_open_chr_lengths.keySet()) {
            for (Integer length : hm_to_open_chr_lengths.get(hm)) {
                bw_data_file.write(hm + "\t" + length + "\n");
            }
        }

        bw_data_file.close();

        //create Rscript and execute Rscript
        File f_script = new File(f_output_script.getAbsolutePath() + File.separator +
                options_intern.file_suffix_tepic_postprocessing_violin_plots_script);
        File f_plot = new File(f_output_plot.getAbsolutePath() + File.separator +
                options_intern.file_suffix_tepic_postprocessing_violin_plots_plot);

        StringBuilder sb_script = new StringBuilder();
        sb_script.append("require(ggplot2)\n" + "require(scales)\n");
        sb_script.append("hm_to_lengths=read.csv('" + f_data_file.getAbsolutePath() + "',sep='\\t')\n");
        sb_script.append(
                "p <- ggplot(hm_to_lengths, aes(HM, LENGTH))\n" + "p<-p + geom_violin() + geom_boxplot(width=0.1)\n" +
                        "p<-p+ scale_y_log10(breaks = trans_breaks(\"log10\", function(x) 10^x),\n" +
                        "              labels = trans_format(\"log10\", math_format(10^.x)))\n" +
                        "p<- p+ theme(axis.title.x = element_blank())\n" +
                        "p <- p+ylab(\"log10 length of open chromatin\") +theme(text = element_text(size=20))\n" +
                        "#p\n");

        sb_script.append("ggsave(filename ='" + f_plot.getAbsolutePath() + "' , plot = p)\n");


        BufferedWriter bw_script = new BufferedWriter(new FileWriter(f_script));
        bw_script.write(sb_script.toString());
        bw_script.close();

        String command_execute = "Rscript " + f_script.getAbsolutePath();
        logger.logLine("[TEPIC-ANALYSIS] execute Rscript for violin plot: " + command_execute);
        Process child = Runtime.getRuntime().exec(command_execute);
        int code = child.waitFor();
        switch (code) {
            case 0:
                break;
            case 1:
                String message = child.getErrorStream().toString();
                throw new Exception(message);
        }


        logger.logLine("[TEPIC-ANALYSIS] Finished Violine Plots");
    }

    /**
     * randomizes TEPIC output for evaluation
     *
     * @throws IOException
     */
    public void randomize_tepic() throws IOException {
        logger.logLine("[TEPIC] Randomize output start.");

        File f_input_root = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_tepic_output_raw);

        File f_output_root = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_tepic_output_raw_shuffle);
        f_output_root.mkdir();

        for (File fileDir_tp : f_input_root.listFiles()) {
            File f_output_tp = new File(f_output_root.getAbsolutePath() + File.separator + fileDir_tp.getName());
            f_output_tp.mkdir();

            for (File fileDir_hm : fileDir_tp.listFiles()) {
                File f_output_tp_hm = new File(f_output_tp.getAbsolutePath() + File.separator + fileDir_hm.getName());
                f_output_tp_hm.mkdir();

                for (File fileDir_sample : fileDir_hm.listFiles()) {
                    File f_output_tp_hm_sample =
                            new File(f_output_tp_hm.getAbsolutePath() + File.separator + fileDir_sample.getName());
                    f_output_tp_hm_sample.mkdir();

                    ArrayList<File> affinity_files = new ArrayList<>();

                    for (File fileDir_affinity : fileDir_sample.listFiles()) {
                        if (fileDir_affinity.getName().matches(".*_Affinity_Gene_View_Filtered.*")) {
                            affinity_files.add(fileDir_affinity);
                        }
                    }

                    File file_to_shuffle = new File("");

                    if (affinity_files.size() == 2) {
                        for (File f : affinity_files) {
                            if (f.getName().matches(".*TPM.*")) {
                                file_to_shuffle = f;
                                break;
                            }
                        }
                    } else {
                        file_to_shuffle = affinity_files.get(0);
                    }

                    //now shuffle columns of this file and write line by line to new file
                    File f_out_shuffled = new File(
                            f_output_tp_hm_sample.getAbsolutePath() + File.separator + file_to_shuffle.getName());
                    BufferedWriter bw = new BufferedWriter(new FileWriter(f_out_shuffled));
                    BufferedReader br = new BufferedReader(new FileReader(file_to_shuffle));
                    String line = br.readLine();
                    String[] split_header = line.split("\t");

                    bw.write(line);
                    bw.newLine();

                    ArrayList<Integer> which_column_writes_which = new ArrayList<>();
                    which_column_writes_which.add(0);

                    ArrayList<Integer> not_shuffled = new ArrayList<>();
                    for (int i = 1; i < split_header.length; i++) {
                        not_shuffled.add(i);
                    }

                    Collections.shuffle(not_shuffled);
                    which_column_writes_which.addAll(not_shuffled);

                    while ((line = br.readLine()) != null) {
                        String[] split = line.split("\t");

                        StringBuilder sb = new StringBuilder();

                        int c = 0;
                        for (Integer i : which_column_writes_which) {
                            sb.append(split[i]);
                            if (c < which_column_writes_which.size() - 1) {
                                sb.append("\t");
                            }
                            c++;
                        }

                        bw.write(sb.toString());
                        bw.newLine();

                    }

                    br.close();
                    bw.close();
                }
            }
        }

        logger.logLine("[TEPIC] Randomize output finished.");

    }

    /**
     * creates the TEPIC command lines and runs all samples of all histone modification and all timepoints
     *
     * @throws Exception
     */
    public void run_tepic() throws Exception {
        logger.logLine("Start TEPIC.sh");

        check_tepic_input_with_options();

        /*
        if(options_intern.mix_level.equals("SAMPLE_LEVEL"))
        {
            File root_mix_working_dir = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option);
            File f_sample_mix_output = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_sample_mix);
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory=f_sample_mix_output.getAbsolutePath();
        }
        if(options_intern.mix_level.equals("HM_LEVEL"))
        {
            File root_mix_working_dir = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option);
            File f_output_hm = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_hm_mix);
            f_output_hm.mkdir();

            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory=f_output_hm.getAbsolutePath();
        }
        if(!options_intern.black_list_dir.equals(""))
        {
            File output_folder = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_blacklisted_regions);
            File output_folder_new_input = new File(output_folder.getAbsolutePath()+File.separator+options_intern.folder_name_blacklisted_regions_new_input);
            output_folder_new_input.mkdir();

            //set new folder directory for tepic input and save old one
            options_intern.tepic_input_prev=options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = output_folder_new_input.getAbsolutePath();
        }
        if(options_intern.mix_mutually_exclusive)
        {
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option+File.separator+options_intern.folder_name_mix_option_mutually_exclusive+File.separator+options_intern.folder_name_mix_options_mutually_exclusive_input;
        }*/
        logger.logLine("[TEPIC] Used data: " + options_intern.tepic_input_directory);


        String command = "bash";
        String tepic_path = " " + options_intern.path_to_COM2POSE + File.separator +
                options_intern.directory_for_tepic_scripts_code_tepic_sh;
        command += tepic_path;
        command += " -g " + options_intern.tepic_input_ref_genome;
        command += " -p " + options_intern.tepic_path_pwms;


        String command_tail = "";
        if (options_intern.tepic_cores > 1) {
            command_tail += " -c " + options_intern.tepic_cores;
        }
        if (!options_intern.tepic_bed_chr_sign.equals("")) {
            command_tail += " -d " + options_intern.tepic_bed_chr_sign;
        }
        if (options_intern.tepic_column_bedfile != -1) {
            command_tail += " -n " + options_intern.tepic_column_bedfile;
        }
        if (!options_intern.tepic_gene_annot.equals("")) {
            command_tail += " -a " + options_intern.tepic_gene_annot;
        }
        if (options_intern.tepic_window_size != 50000) {
            command_tail += " -w " + options_intern.tepic_window_size;

        }
        if (!options_intern.tepic_onlyDNasePeaks.equals("")) {
            command_tail += " -f " + options_intern.tepic_onlyDNasePeaks;
        }
        if (options_intern.tepic_exponential_decay) {
            command_tail += " -e TRUE";
        }
        if (options_intern.tepic_not_norm_peak_length) {
            command_tail += " -l TRUE";
        }
        if (options_intern.tepic_not_generated) {
            command_tail += " -u TRUE";
        }
        if (options_intern.tepic_original_decay) {
            command_tail += " -x TRUE";
        }
        if (!options_intern.tepic_psems_length.equals("")) {
            command_tail += " -m " + options_intern.tepic_psems_length;
        }
        if (options_intern.tepic_entire_gene_body) {
            command_tail += " -y TRUE";
        }
        if (options_intern.tepic_zipped) {
            command_tail += " -z TRUE";
        }
        if (!options_intern.tepic_2bit.equals("")) {
            command_tail += " -r " + options_intern.tepic_2bit;
        }
        if (options_intern.tepic_pvalue != 0.05) {
            command_tail += " -v " + options_intern.tepic_pvalue;

        }
        if (options_intern.tepic_minutes_per_chr != 3) {
            command_tail += " -i " + options_intern.tepic_minutes_per_chr;
        }
        if (options_intern.tepic_chr_prefix) {
            command_tail += " -j TRUE";
        }
        if (options_intern.tepic_transcript_based) {
            command_tail += " -t TRUE";
        }
        if (!options_intern.tepic_loop_list.equals("")) {
            command_tail += " -h " + options_intern.tepic_loop_list;
        }
        if (options_intern.tepic_loop_windows != 5000) {
            command_tail += " -s " + options_intern.tepic_loop_windows;
        }
        if (options_intern.tepic_only_peak_features) {
            command_tail += " -q TRUE";
        }
        if (options_intern.tepic_tpm_cutoff > 0) {
            command_tail += " -T " + options_intern.tepic_tpm_cutoff;
            command_tail += " -E " + options_intern.tepic_ensg_symbol;
            command_tail += " -A " + options_intern.deseq2_input_gene_id;
        }

        command_tail += " -B " + options_intern.tepic_tf_binding_site_search;

        File output_TEPIC = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_tepic_output_raw);
        output_TEPIC.mkdir();


        File folder = new File(options_intern.tepic_input_directory);
        for (File dirGroup : folder.listFiles()) {
            if (dirGroup.isDirectory()) {
                File output_TEPIC_group =
                        new File(output_TEPIC.getAbsolutePath() + File.separator + dirGroup.getName());
                output_TEPIC_group.mkdir();

                logger.logLine("[TEPIC] Start group " + dirGroup.getName());
                for (File dirHM : dirGroup.listFiles()) {
                    File output_TEPIC_group_hm =
                            new File(output_TEPIC_group.getAbsolutePath() + File.separator + dirHM.getName());
                    output_TEPIC_group_hm.mkdir();

                    logger.logLine("[TEPIC] Start histone modification: " + dirHM.getName());
                    for (File sample : dirHM.listFiles()) {
                        logger.logLine("[TEPIC] Start sample " + sample.getName());

                        File output_sample =
                                new File(output_TEPIC_group_hm.getAbsolutePath() + File.separator + sample.getName());
                        output_sample.mkdir();

                        String command_sample = new String(command);

                        String output_dir_combined =
                                output_sample.getAbsolutePath() + File.separator + output_sample.getName();

                        command_sample += " -b " + sample.getAbsolutePath();
                        command_sample += " -o " + output_dir_combined;

                        String command_tail_sample = new String(command_tail);
                        if (options_intern.tepic_tpm_cutoff > 0) {
                            String n_dir = "";

                            if (!options_intern.mix_mutually_exclusive) {
                                n_dir = options_intern.com2pose_working_directory + File.separator +
                                        options_intern.folder_name_deseq2_preprocessing + File.separator +
                                        options_intern.folder_name_deseq2_preprocessing_single + File.separator +
                                        dirGroup.getName() + options_intern.file_suffix_deseq2_preprocessing_meanCounts;

                            } else {
                                String[] names_sample = sample.getName().split("_");
                                String name_sample = names_sample[0];
                                n_dir = options_intern.com2pose_working_directory + File.separator +
                                        options_intern.folder_name_deseq2_preprocessing + File.separator +
                                        options_intern.folder_name_deseq2_preprocessing_single + File.separator +
                                        name_sample + options_intern.file_suffix_deseq2_preprocessing_meanCounts;
                            }

                            //n_dir=options_intern.com2pose_working_directory+File.separator+ options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_single+File.separator+dirGroup.getName()+options_intern.file_suffix_deseq2_preprocessing_meanCounts;

                            command_tail_sample += " -G " + n_dir;
                        }

                        if (!options_intern.path_tgen.equals("") && options_intern.tepic_tgene_target_genes) {
                            File f_tgene_input;

                            if (!options_intern.mix_mutually_exclusive) {
                                f_tgene_input = new File(options_intern.com2pose_working_directory + File.separator +
                                        options_intern.folder_name_tgen + File.separator +
                                        options_intern.folder_name_tgen_output + File.separator + dirGroup.getName() +
                                        File.separator + dirHM.getName() + File.separator + dirGroup.getName() + "_" +
                                        dirHM.getName() + File.separator + options_intern.file_suffix_tgen_output);
                            } else {
                                String[] names_sample = sample.getName().split("_");
                                String name_sample = names_sample[0];

                                f_tgene_input = new File(options_intern.com2pose_working_directory + File.separator +
                                        options_intern.folder_name_tgen + File.separator +
                                        options_intern.folder_name_tgen_output + File.separator + dirGroup.getName() +
                                        File.separator + dirHM.getName() + File.separator + name_sample + "_" +
                                        dirHM.getName() + File.separator + options_intern.file_suffix_tgen_output);
                            }
                            //f_tgene_input = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_output+File.separator+dirGroup.getName()+File.separator+dirHM.getName()+File.separator+dirGroup.getName()+"_"+dirHM.getName()+File.separator+options_intern.file_suffix_tgen_output);

                            command_tail_sample += " -L " + f_tgene_input.getAbsolutePath();
                            logger.logLine(
                                    "[TEPIC-TGENE] using only tgene linked target genes for TF-Gene score calculation.");
                        }

                        String command_execute = command_sample + command_tail_sample;
                        logger.logLine("[TEPIC] execute TEPIC with command line: " + command_execute);
                        Process child = Runtime.getRuntime().exec(command_execute);
                        int code = child.waitFor();
                        switch (code) {
                            case 0:
                                break;
                            case 1:
                                String message = child.getErrorStream().toString();
                                throw new Exception(message);
                        }
                    }
                }
            }
        }
        logger.logLine("[TEPIC] Finished TEPIC.sh");
    }

    /**
     * run created DESeq2 scripts and postprocess for input into DYNAMITE
     */
    public void run_and_postprocess_DESeq2() throws Exception {

        logger.logLine("Start running DESeq2 RScripts");

        File folder = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_R_scripts);

        for (File dir : folder.listFiles()) {
            if (!dir.isDirectory()) {
                String command = "Rscript " + dir.getAbsolutePath();
                logger.logLine("[DESEQ2] Running script " + dir.getName() + ": " + command);
                Process child = Runtime.getRuntime().exec(command);
                int code = child.waitFor();
                switch (code) {
                    case 0:
                        break;
                    case 1:
                        String message = child.getErrorStream().toString();
                        throw new Exception(message);
                }
            }
        }

        logger.logLine("Finished running DESeq2 RScripts");
        logger.logLine("Start postprocessing DESeq2 data for input to DYNAMITE");


        File folder_results = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_output_raw);
        File output_file = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_deseq2_output);
        output_file.mkdir();

        for (File res : folder_results.listFiles()) {
            if (!res.isDirectory()) {
                String name_res = res.getName();
                String[] split_name_res = name_res.split("_");
                String name_out =
                        split_name_res[0] + "_" + split_name_res[1] + options_intern.file_suffix_deseq2_output_DYNAMITE;

                BufferedReader br = new BufferedReader(new FileReader(res));
                BufferedWriter bw =
                        new BufferedWriter(new FileWriter(output_file.getAbsolutePath() + File.separator + name_out));

                bw.write("geneID\tlog2fc");
                bw.newLine();

                String line = br.readLine();

                while ((line = br.readLine()) != null) {
                    String[] split = line.split("\t");
                    if (!split[2].equals("NA")) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(split[0].substring(1, split[0].length() - 1));
                        sb.append("\t");
                        sb.append(split[2]);

                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                bw.close();
                br.close();


            }
        }


        logger.logLine("Finished postprocessing DESeq2 data for input to DYNAMITE");

    }

    /**
     * edit input files if TPM filter is set
     */
    public void preprocess_deseq2_input_tpm() throws IOException {
        logger.logLine("[DESEQ2-PREP-TPM] Filter genes with TPM under cutoff.");

        //read in TPMs
        HashMap<String, HashMap<String, Double>> tp_gene_tpm_value = new HashMap<>();

        File f_input_tpms_root = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.folder_name_deseq2_preprocessing_tpm + File.separator +
                options_intern.folder_name_deseq2_preprocessing_tpm_results);

        for (File f_tpm_tp : f_input_tpms_root.listFiles()) {
            if (f_tpm_tp.isFile()) {
                String name_tp = f_tpm_tp.getName().split("\\.")[0].replace("_tpms", "");

                if (name_tp.split("_").length > 1) {
                    String[] split = name_tp.split("_");
                    name_tp = split[0].substring(0, 1).toUpperCase() + split[0].substring(1) +
                            split[1].substring(0, 1).toUpperCase() + split[1].substring(1);
                }

                HashMap<String, Double> gene_tpm = new HashMap<>();

                BufferedReader br = new BufferedReader(new FileReader(f_tpm_tp));
                String line = br.readLine();
                String[] split_header = line.split("\t");
                while ((line = br.readLine()) != null) {
                    String[] split = line.split("\t");
                    if (split[3].equals("NA")) {
                        continue;
                    }
                    gene_tpm.put(split[0], Double.parseDouble(split[3]));
                }
                br.close();

                tp_gene_tpm_value.put(name_tp, gene_tpm);
            }
        }

        File f_combined_root = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.folder_name_deseq2_preprocessing_combined);

        File f_output_copy = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.folder_name_deseq2_preprocessing_combined_original);
        f_output_copy.mkdirs();

        for (File f_input : f_combined_root.listFiles()) {
            if (f_input.isDirectory()) {
                String[] groups = f_input.getName().split("_");
                String name1 = groups[0];
                String name2 = groups[1];

                if (groups.length > 2) {
                    int found_non = -1;

                    for (int i = 0; i < groups.length; i++) {
                        if (groups[i].equals("non") || groups[i].equals("Non")) {
                            found_non = i;
                        }
                    }

                    name1 = groups[0].substring(0, 1).toUpperCase() + groups[0].substring(1);
                    name2 = groups[1].substring(0, 1).toUpperCase() + groups[1].substring(1);

                    if (found_non == 0) {
                        name1 = groups[0].substring(0, 1).toUpperCase() + groups[0].substring(1) +
                                groups[1].substring(0, 1).toUpperCase() + groups[1].substring(1);
                        name2 = groups[2].substring(0, 1).toUpperCase() + groups[2].substring(1);
                    } else if (found_non == 1) {
                        name1 = groups[0].substring(0, 1).toUpperCase() + groups[0].substring(1);
                        name2 = groups[1].substring(0, 1).toUpperCase() + groups[1].substring(1) +
                                groups[2].substring(0, 1).toUpperCase() + groups[2].substring(1);
                    }

                }


                File f_output_copy_group =
                        new File(f_output_copy.getAbsolutePath() + File.separator + f_input.getName());
                f_output_copy_group.mkdirs();

                for (File f_input_read : f_input.listFiles()) {
                    if (f_input_read.isFile()) {
                        StringBuilder sb_copy = new StringBuilder();

                        StringBuilder sb_tpm_filtered = new StringBuilder();

                        BufferedReader br = new BufferedReader(new FileReader(f_input_read));
                        String line = br.readLine();

                        sb_copy.append(line);
                        sb_copy.append("\n");

                        sb_tpm_filtered.append(line);
                        sb_tpm_filtered.append("\n");

                        String[] header = line.split("\t");

                        for (int i = 1; i < header.length; i++) {
                            if (header[i].matches(name1 + ".*")) {
                                header[i] = name1;
                            }
                            if (header[i].matches(name2 + ".*")) {
                                header[i] = name2;
                            }
                        }

                        while ((line = br.readLine()) != null) {
                            sb_copy.append(line);
                            sb_copy.append("\n");

                            String[] split = line.split("\t");

                            sb_tpm_filtered.append(split[0]);

                            boolean one_changed = false;

                            StringBuilder changer = new StringBuilder();

                            for (int i = 1; i < split.length; i++) {
                                String group_name = header[i];
                                HashMap<String, Double> lookup = new HashMap<>();

                                if (tp_gene_tpm_value.containsKey(group_name)) {
                                    lookup = tp_gene_tpm_value.get(group_name);
                                } else {
                                    System.out.println("X");
                                }


                                if (!lookup.containsKey(split[0])) {
                                    changer.append("\t");
                                    changer.append("0");
                                    one_changed = true;
                                    continue;
                                }

                                double tpm_value = lookup.get(split[0]);
                                if (tpm_value < options_intern.deseq2_tpm_filter) {
                                    changer.append("\t");
                                    changer.append("0");
                                    one_changed = true;
                                } else {
                                    changer.append("\t");
                                    changer.append(split[i]);
                                }
                            }
                            changer.append("\n");

                            if (one_changed) {
                                StringBuilder sb_complete_changer = new StringBuilder();
                                for (int i = 1; i < split.length; i++) {
                                    sb_complete_changer.append("\t0");
                                }
                                sb_complete_changer.append("\n");

                                sb_tpm_filtered.append(sb_complete_changer.toString());

                            } else {
                                sb_tpm_filtered.append(changer.toString());
                            }

                        }
                        br.close();

                        BufferedWriter bw_copy = new BufferedWriter(new FileWriter(new File(
                                f_output_copy_group.getAbsolutePath() + File.separator + f_input_read.getName())));
                        bw_copy.write(sb_copy.toString());
                        bw_copy.close();

                        BufferedWriter bw_tpm_filtered = new BufferedWriter(new FileWriter(f_input_read));
                        bw_tpm_filtered.write(sb_tpm_filtered.toString());
                        bw_tpm_filtered.close();
                    }
                }
            }
        }


    }

    /**
     * created a chromosome and position mapping from biomart
     */
    public void create_gene_positions() throws Exception {

        if (!options_intern.calculcate_gene_positions) {
            return;
        }

        File f_output_positions_root = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.folder_name_deseq2_preprocessing_gene_positions);
        f_output_positions_root.mkdir();

        File f_data_input = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.file_suffix_deseq2_mapping);

        File f_script = new File(f_output_positions_root.getAbsolutePath() + File.separator +
                options_intern.file_suffix_deseq2_preprocessing_gene_positions_script);

        File f_data_prev = new File(f_output_positions_root.getAbsolutePath() + File.separator +
                options_intern.file_suffix_deseq2_preprocessing_gene_positions_data_prev);
        File f_data_version = new File(f_output_positions_root.getAbsolutePath() + File.separator +
                options_intern.file_suffix_deseq2_preprocessing_gene_positions_data_prev_version);



        logger.logLine("[PREP] Create gene positions for all RNA-seq data.");
        logger.logLine("[PREP] Get gene lengths ...");

        logger.logLine("[PREP] if not done in 5 hours:");
        logger.logLine("[PREP] 1. Please stop com2pose.");
        logger.logLine("[PREP] 2. Please run script manually in RStudio.");
        logger.logLine("[PREP] Script path: " + f_script.getAbsolutePath());
        logger.logLine(
                "[PREP] 3. Afterwards add paramter -a to com2pose command line, so this script wont be started again.");
        logger.logLine("[PREP] 4. restart com2pose");
        logger.logLine("[PREP] waiting ...");


        //get biomart gene positions

        StringBuilder sb_script = new StringBuilder();
        sb_script.append(
                "if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n" + "  install.packages(\"BiocManager\")\n" +
                        "\n" + "if (!requireNamespace(\"biomaRt\", quietly = TRUE))\n" +
                        "  BiocManager::install(\"biomaRt\")\n" + "\n" + "library('biomaRt')\n" +
                        "httr::set_config(httr::config(ssl_verifypeer = FALSE))\n\n");

        sb_script.append("input<-read.csv('" + f_data_input.getAbsolutePath() + "', sep = '\\t')\n");
        sb_script.append("input<-input$" + options_intern.deseq2_biomart_dataset_symbol_column + "\n");

        sb_script.append("not_done=TRUE\n" + "\n" + "G_list= data.frame()\n" + "\n" + "while(not_done)\n" + "{\n" +
                "  tryCatch({\n" + "    mart <- useDataset(\"" + options_intern.deseq2_biomart_dataset_species +
                "\", useMart(\"ensembl\"))\n" + "    #df$id <- NA\n" + "    G_list_intern <- getBM(filters= \"" +
                options_intern.deseq2_biomart_dataset_symbol_column + "\", attributes= c(\"" +
                options_intern.deseq2_biomart_dataset_symbol_column +
                "\",\"chromosome_name\",\"start_position\",\"end_position\",\"strand\",\"band\"),values=input,mart= mart)\n" +
                "    G_list=rbind(G_list,G_list_intern)\n" + "    not_done=FALSE\n" + "  }, warning = function(w) {\n" +
                "    print(\"WARNING SECTION\")\n" + "    print(w)\n" + "  }, error = function(e) {\n" +
                "    print(\"ERROR SECTION\")\n" + "    print(e)\n" + "  }, finally = {\n" + "  })\n" + "}\n" +
                "write.table(G_list,\"" + f_data_prev.getAbsolutePath() +
                "\", row.names = FALSE, quote = F, sep=\"\\t\")\n");

        sb_script.append("not_done=TRUE\n" + "while (not_done) {\n" + "  tryCatch({\n" +
                "    #check which genome version we have in biomart\n" + "    \n" +
                "    ensembl <- useEnsembl(biomart = \"genes\")\n" + "    datasets <- listDatasets(ensembl)\n" +
                "    used_dataset=searchDatasets(mart = ensembl, pattern = \"mmusculus_gene_ensembl\")\n" +
                "    version=used_dataset[,3]\n" + "    version=as.character(version)\n" +
                "    fileConn<-file(\""+f_data_version.getAbsolutePath()+"\")\n" +
                "    writeLines(version, fileConn)\n" + "    \n" + "    not_done=FALSE\n" +
                "  }, warning = function(w) {\n" + "    print(\"WARNING SECTION\")\n" + "    print(w)\n" +
                "  }, error = function(e) {\n" + "    print(\"ERROR SECTION\")\n" + "    print(e)\n" +
                "  }, finally = {\n" + "  })\n" + "}\n\n");

        BufferedWriter bw_script = new BufferedWriter(new FileWriter(f_script));
        bw_script.write(sb_script.toString());
        bw_script.close();

        String command = "Rscript " + f_script;
        logger.logLine("[PREP] run R script: " + command);

        Process child = Runtime.getRuntime().exec(command);
        int code = child.waitFor();
        switch (code) {
            case 0:
                break;
            case 1:
                String message = child.getErrorStream().toString();
                logger.logLine(
                        "[PREP-TPM] Script failed due to bioMart connection error. Please run script manually in RStudio.");
                logger.logLine("[PREP-TPM] Script path: " + f_script.getAbsolutePath());
                logger.logLine(
                        "[PREP-TPM] Afterwards, add paramter -b to com2pose command line, so this script wont be started again.");
                throw new Exception(message);
        }

        File f_data = new File(f_output_positions_root.getAbsolutePath() + File.separator +
                options_intern.file_suffix_deseq2_preprocessing_gene_positions_data);
        File f_uplift_script = new File(f_output_positions_root.getAbsolutePath() + File.separator +
                options_intern.file_suffix_deseq2_preprocessing_gene_positions_data_upliftpy);

        StringBuilder sb_pyuplift = new StringBuilder();
        sb_pyuplift.append("import pip\n" + "\n" + "def import_or_install(package):\n" + "    try:\n" +
                "        __import__(package)\n" + "    except ImportError:\n" +
                "        pip.main(['install', package])\n" + "\n" + "import_or_install(\"pyliftover\")\n" +
                "import_or_install(\"pandas\")\n" + "import_or_install(\"re\")\n" + "\n" + "import re\n" +
                "import pandas as pd\n" + "from pyliftover import LiftOver\n" + "\n" +
                "def convert(c, x, y, s, converter):\n" + "    first = converter.convert_coordinate(c, int(x),s)\n" +
                "    second = converter.convert_coordinate(c, int(y),s)\n" + "\n" +
                "    if(first is None or second is None):\n" + "        return None, None\n" + "\n" +
                "    if len(first) == 0 or len(second) == 0:\n" + "        return None, None\n" + "\n" +
                "    return str(first[0][1]), str(second[0][1])\n" + "\n" + "def main():\n" +
                "    path_to_X = \""+f_data_prev.getAbsolutePath()+"\"\n" +
                "    path_to_version =\""+f_data_version.getAbsolutePath()+"\"\n" +
                "    path_to_newSave = \""+f_data.getAbsolutePath()+"\"\n" +
                "\n" + "    version_of_our_dat=\""+options_intern.igv_species_ref_genome+"\"\n" + "\n" +
                "    version_of_our_dat_n= re.findall(r'\\d+', version_of_our_dat)\n" +
                "    version_of_our_dat_n=version_of_our_dat_n.__getitem__(0)\n" +
                "    name=version_of_our_dat.split(version_of_our_dat_n)\n" + "    name=name.__getitem__(0)\n" + "\n" +
                "    version = []\n" + "    with open(path_to_version) as file:\n" +
                "        version=file.readlines()\n" + "\n" + "    version = version.__getitem__(0).rstrip(\"\\n\")\n" +
                "    version_n = re.findall(r'\\d+', version)\n" + "    version_n = version_n.__getitem__(0)\n" + "\n" +
                "    version_convert_from = name+str(version_n)\n" + "\n" + "    please_convert=True\n" +
                "    if(version_of_our_dat_n == version_convert_from):\n" + "        please_convert=False\n" + "\n" +
                "\n" + "    df = pd.read_csv(path_to_X, sep=\"\\t\")\n" + "\n" + "    data_output = []\n" + "\n" +
                "    column_names=[]\n" + "\n" + "    for col in df.columns:\n" + "        column_names.append(col)\n" +
                "\n" + "    converter = LiftOver(version_convert_from, version_of_our_dat)\n" + "\n" +
                "    for index, row in df.iterrows():\n" + "        mgi_symbol = row[column_names.__getitem__(0)]\n" +
                "        chromosome_not_edited=str(row[column_names.__getitem__(1)])\n" +
                "        chromosome = \"chr\" + str(row[column_names.__getitem__(1)])\n" +
                "        start_position = row[column_names.__getitem__(2)]\n" +
                "        end_position = row[column_names.__getitem__(3)]\n" +
                "        strand = row[column_names.__getitem__(4)]\n" +
                "        band = row[column_names.__getitem__(5)]\n" + "\n" +
                "        search_string = chromosome+\":\"+str(start_position)+\"-\"+str(end_position)\n" +
                "        x_converted=start_position\n" + "        y_converted=end_position\n" + "\n" +
                "        if please_convert:\n" +
                "            x_converted, y_converted = convert(chromosome, start_position, end_position, strand, converter)\n" +
                "            if(x_converted is None or y_converted is None):\n" + "                continue\n" + "\n" +
                "        start_position=x_converted\n" + "        end_position=y_converted\n" + "\n" +
                "        data_output.append([mgi_symbol,chromosome_not_edited,start_position,end_position,strand,band])\n" +
                "\n" + "        #string_converted = chromosome+\":\"+x_converted+\"-\"+y_converted\n" +
                "        #print(\"X\")\n" + "\n" +
                "    df_converted = pd.DataFrame(data_output, columns=column_names)\n" +
                "    df_converted.to_csv(path_to_newSave,sep=\"\\t\",index=False)\n\n" +
                "main()\n\n");

        BufferedWriter bw_pyuplift = new BufferedWriter(new FileWriter(f_uplift_script));
        bw_pyuplift.write(sb_pyuplift.toString());
        bw_pyuplift.close();

        logger.logLine("[PREP] Uplift positions to correct genome version");

        String command_pyuplift = "python3 "+ f_uplift_script.getAbsolutePath();

        logger.logLine("[PREP] executing command: " + command_pyuplift);

        Process child_pyuplift = Runtime.getRuntime().exec(command_pyuplift);
        int code_pyuplift = child_pyuplift.waitFor();
        switch (code_pyuplift) {
            case 0:
                break;
            case 1:
                String message = child_pyuplift.getErrorStream().toString();
                throw new Exception(message);
        }



        logger.logLine("[PREP] Finished creating gene positions for all RNA-seq data.");


    }

    /**
     * creates a mapping for TPM values for all RNA-seq genes
     */
    public void create_TPM_mappings() throws Exception {
        logger.logLine("[PREP] Create TPM values for all RNA-seq data.");

        File f_output_ot = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.folder_name_deseq2_preprocessing_tpm);
        f_output_ot.mkdirs();

        File f_output_script_lengths = new File(f_output_ot.getAbsolutePath() + File.separator +
                options_intern.file_suffix_deseq2_preprocessing_tpm_mapping_get_geneLengths_script);
        File f_output_lengths = new File(f_output_ot.getAbsolutePath() + File.separator +
                options_intern.file_suffix_deseq2_preprocessing_tpm_mapping_geneLengths_file);

        logger.logLine("[PREP] Get gene lengths ...");

        logger.logLine("[PREP] if not done in 5 hours:");
        logger.logLine("[PREP] 1. Please stop com2pose.");
        logger.logLine("[PREP] 2. Please run script manually in RStudio.");
        logger.logLine("[PREP] Script path: " + f_output_script_lengths.getAbsolutePath());
        logger.logLine(
                "[PREP] 3. Afterwards add paramter -a to com2pose command line, so this script wont be started again.");
        logger.logLine("[PREP] 4. restart com2pose");
        logger.logLine("[PREP] waiting ...");


        if (options_intern.calculate_tpm_lengths) {
            StringBuilder sb_lengths = new StringBuilder();
            sb_lengths.append("if(!\"EDASeq\" %in% rownames(installed.packages()))\n" + "{\n" +
                    "  if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n" + "  {\n" +
                    "    install.packages(\"BiocManager\")\n" + "    \n" + "  }\n" + "  \n" +
                    "  BiocManager::install(\"EDASeq\")\n" + "}\n" + "\n" +
                    "if(!\"biomaRt\" %in% rownames(installed.packages()))\n" + "{\n" +
                    "  if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n" + "  {\n" +
                    "    install.packages(\"BiocManager\")\n" + "  }\n" + "  BiocManager::install(\"biomaRt\")\n" +
                    "}\n" + "require(\"dplyr\")\n" + "require(\"biomaRt\")\n" + "require(\"EDASeq\")\n" +
                    "require(\"data.table\")\n" + "\n" + "httr::set_config(httr::config(ssl_verifypeer = FALSE))\n" +
                    "ensg_names<-read.csv('" + options_intern.deseq2_input_gene_id + "',sep='\\t')\n" +
                    "ensg_names=dplyr::filter(ensg_names, !grepl(\"_\",Geneid))\n" +
                    "list_of_something<-split(ensg_names, (0:nrow(ensg_names) %/% 5000))  # modulo division\n" + "\n" +
                    "ensembl_to_lenth<- data.frame()\n" + "\n" + "counter=1\n" +
                    "counter_end=length(list_of_something)\n" + "\n" + "while(counter<=counter_end)\n" + "{\n" +
                    "  print(\"starting\")\n" + "  print(counter)\n" + "  df<-unlist(list_of_something[counter])\n" +
                    "  df<-as.data.frame(df)\n" + "  df_names<-df$df\n" + "  result = tryCatch({\n" +
                    "    ensembl_to_lenth_slice=getGeneLengthAndGCContent(df_names, \"" +
                    options_intern.deseq2_biomart_dataset_species + "\")\n" +
                    "    ensembl_to_lenth<-rbind(ensembl_to_lenth,ensembl_to_lenth_slice)\n" +
                    "    print(\"Success: I made it!!\")\n" + "    print(\"current length of output csv\")\n" +
                    "    print(length(ensembl_to_lenth))\n" + "    counter=counter+1\n" +
                    "  }, warning = function(w) {\n" + "    print(\"WARNING SECTION\")\n" + "    print(w)\n" +
                    "  }, error = function(e) {\n" + "    print(\"ERROR SECTION\")\n" + "    print(e)\n" +
                    "  }, finally = {\n" + "  })\n" + "}\n\n");
            sb_lengths.append("ensembl_to_lenth=setDT(ensembl_to_lenth, keep.rownames = TRUE)[]\n" +
                    "colnames(ensembl_to_lenth)<-c(\"ENSG\",\"length\",\"gc\")\n" +
                    "write.table(ensembl_to_lenth, file=\"" + f_output_lengths.getAbsolutePath() +
                    "\", sep='\\t', quote=FALSE, row.names=FALSE)\n");

            /*
            StringBuilder sb_lengths = new StringBuilder();
            sb_lengths.append("if(!\"EDASeq\" %in% rownames(installed.packages()))\n" +
                    "{\n" +
                    "  if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n" +
                    "  {\n" +
                    "    install.packages(\"BiocManager\")\n" +
                    "  \n" +
                    "  }\n" +
                    "  \n" +
                    "  BiocManager::install(\"EDASeq\")\n" +
                    "}\n" +
                    "\n" +
                    "if(!\"biomaRt\" %in% rownames(installed.packages()))\n" +
                    "{\n" +
                    "  if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n" +
                    "  {\n" +
                    "    install.packages(\"BiocManager\")\n" +
                    "  }\n" +
                    "  BiocManager::install(\"biomaRt\")\n" +
                    "}\n");

            sb_lengths.append("require(\"dplyr\")\n");
            sb_lengths.append("require(\"biomaRt\")\n");
            sb_lengths.append("require(\"EDASeq\")\n");
            sb_lengths.append("require(\"data.table\")\n");
            sb_lengths.append("httr::set_config(httr::config(ssl_verifypeer = FALSE))\n");

            sb_lengths.append("ensg_names<-read.csv('"+options_intern.deseq2_input_gene_id+"',sep='\\t')\n");
            sb_lengths.append("ensembl_list<-ensg_names$Geneid\n");
            sb_lengths.append("ensembl_to_lenth=getGeneLengthAndGCContent(ensembl_list, \""+options_intern.deseq2_biomart_dataset_species+"\")\n");
            sb_lengths.append("ensembl_to_lenth=as.data.frame(ensembl_to_lenth)\n");
            sb_lengths.append("ensembl_to_lenth=setDT(ensembl_to_lenth, keep.rownames = TRUE)[]\n");
            sb_lengths.append("colnames(ensembl_to_lenth)<-c(\"ENSG\",\"length\",\"gc\")\n");
            sb_lengths.append("write.table(ensembl_to_lenth, file=\""+f_output_lengths.getAbsolutePath()+"\", sep='\\t', quote=FALSE, row.names=FALSE)\n");*/

            BufferedWriter bw_lengths_script = new BufferedWriter(new FileWriter(f_output_script_lengths));
            bw_lengths_script.append(sb_lengths.toString());
            bw_lengths_script.close();

            String command = "Rscript " + f_output_script_lengths;
            logger.logLine("[PREP-TPM] run R script: " + command);

            Process child = Runtime.getRuntime().exec(command);
            int code = child.waitFor();
            switch (code) {
                case 0:
                    break;
                case 1:
                    String message = child.getErrorStream().toString();
                    logger.logLine(
                            "[PREP-TPM] Script failed due to bioMart connection error. Please run script manually in RStudio.");
                    logger.logLine("[PREP-TPM] Script path: " + f_output_script_lengths.getAbsolutePath());
                    logger.logLine(
                            "[PREP-TPM] Afterwards, add paramter -a to com2pose command line, so this script wont be started again.");
                    throw new Exception(message);
            }
        } else {
            logger.logLine(
                    "[PREP-TPM] -a is set, get gene_lengths script is not executed, as it was executed manually.");
        }

        //now do this for all RNA-seq data

        File f_output_scripts_rna_seq = new File(f_output_ot.getAbsolutePath() + File.separator +
                options_intern.folder_name_deseq2_preprocessing_tpm_scripts);
        f_output_scripts_rna_seq.mkdirs();

        File f_output_data_rna_seq = new File(f_output_ot.getAbsolutePath() + File.separator +
                options_intern.folder_name_deseq2_preprocessing_tpm_results);
        f_output_data_rna_seq.mkdirs();

        File f_input_rna_seq_root = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.folder_name_deseq2_preprocessing_gene_symbols);

        HashSet<File> files_to_execute = new HashSet<>();

        for (File f_rnaseq : f_input_rna_seq_root.listFiles()) {
            if (f_rnaseq.isFile()) {
                File f_output_script = new File(f_output_scripts_rna_seq.getAbsolutePath() + File.separator +
                        f_rnaseq.getName().split("\\.")[0] +
                        options_intern.file_suffix_deseq2_preprocessing_tpm_mapping_get_tpm_mappings_script);
                files_to_execute.add(f_output_script);

                File f_output_data = new File(
                        f_output_data_rna_seq.getAbsolutePath() + File.separator + f_rnaseq.getName().split("\\.")[0] +
                                options_intern.file_suffix_deseq2_preprocessing_tpm_mapping_get_tpm_mappings_data);

                StringBuilder sb_rscript_calc_tpm = new StringBuilder();
                sb_rscript_calc_tpm.append("require(\"dplyr\")\n");
                sb_rscript_calc_tpm.append("require(\"tidyr\")\n");
                sb_rscript_calc_tpm.append("require(\"data.table\")\n");
                sb_rscript_calc_tpm.append("require(\"tibble\")\n");

                sb_rscript_calc_tpm.append("tpm3 <- function(counts,len) {\n" + "  x <- counts/len\n" +
                        "  return(t(t(x)*1e6/colSums(x)))\n" + "}\n");

                sb_rscript_calc_tpm.append(
                        "expression_matrix=read.csv(\"" + f_rnaseq.getAbsolutePath() + "\",sep=\"\\t\")\n");
                sb_rscript_calc_tpm.append(
                        "gene_lengths=read.csv(\"" + f_output_lengths.getAbsolutePath() + "\",sep=\"\\t\")\n");
                sb_rscript_calc_tpm.append("##CHECK NAs\n" + "expression_matrix=na.omit(expression_matrix)\n" +
                        "gene_lengths=na.omit(gene_lengths)\n");
                sb_rscript_calc_tpm.append(
                        "##CHECK DISTINCT VALUES\n" + "gene_lengths=gene_lengths%>%dplyr::select(ENSG,length)\n" +
                                "gene_lengths=distinct(gene_lengths)\n" + "\n" +
                                "gene_counts=dplyr::select(expression_matrix, c(\"ENSG\",\"MEAN_COUNT\"))\n" +
                                "gene_counts=distinct(gene_counts)\n" + "\n" +
                                "##CHECK VALUES INSIDE gene_counts and gene_lengths\n" +
                                "available_ensg_counts=gene_counts$ENSG\n" +
                                "available_ensg_lenghts=gene_lengths$ENSG\n" +
                                "gene_counts=filter(gene_counts, ENSG %in% available_ensg_lenghts)\n" +
                                "gene_lengths=filter(gene_lengths, ENSG %in% available_ensg_counts)\n" + "\n" +
                                "gene_counts$ENSG=sort(gene_counts$ENSG)\n" +
                                "gene_lengths$ENSG=sort(gene_lengths$ENSG)\n" + "\n" +
                                "gene_lengths=column_to_rownames(gene_lengths,\"ENSG\")\n" +
                                "gene_counts=column_to_rownames(gene_counts,\"ENSG\")\n" + "\n" +
                                "tpm_value=tpm3(gene_counts,gene_lengths)\n" + "tpm_value=as.data.frame(tpm_value)\n" +
                                "\n" + "tpm_value=setDT(tpm_value, keep.rownames = TRUE)[]\n" +
                                "colnames(tpm_value)=c(\"ENSG\",\"TPM\")\n" + "\n" +
                                "gene_lengths=setDT(gene_lengths, keep.rownames = TRUE)[]\n" +
                                "colnames(gene_lengths)=c(\"ENSG\",\"length\")\n" + "\n" + "\n" +
                                "merged_output=merge(expression_matrix,tpm_value)\n" +
                                "merged_output=merge(merged_output,gene_lengths)\n" + "\n" +
                                "write.table(merged_output, file=\"" + f_output_data.getAbsolutePath() +
                                "\", sep=\"\\t\", row.names = FALSE, quote = FALSE)\n");

                BufferedWriter bw = new BufferedWriter(new FileWriter(f_output_script));
                bw.write(sb_rscript_calc_tpm.toString());
                bw.close();

                String command_intern = "Rscript " + f_output_script;
                logger.logLine("[PREP-TPM] run R script: " + command_intern);

                Process child_intern = Runtime.getRuntime().exec(command_intern);
                int code_intern = child_intern.waitFor();
                switch (code_intern) {
                    case 0:
                        break;
                    case 1:
                        String message = child_intern.getErrorStream().toString();
                        throw new Exception(message);
                }
            }
        }


        logger.logLine("[PREP] TPM values finished.");
    }

    /**
     * create DESeq2 scripts based on input directory for DESeq2 - each group against each group, save intermediate steps and R Scripts
     */
    public void create_DESeq2_scripts() throws IOException {

        logger.logLine("Start preprocessing nfcore RNA-seq for DESeq2 input");
        HashMap<Integer, String> row_ensg_name = new HashMap<>();

        File ensg_names = new File(options_intern.deseq2_input_gene_id);
        BufferedReader br_ensg_per_line = new BufferedReader(new FileReader(ensg_names));
        String line_ensg_per_line = "";
        line_ensg_per_line = br_ensg_per_line.readLine();
        int count_ensg_lines = 0;
        while ((line_ensg_per_line = br_ensg_per_line.readLine()) != null) {
            row_ensg_name.put(count_ensg_lines, line_ensg_per_line);
            count_ensg_lines++;
        }
        br_ensg_per_line.close();

        File output_intermediate_steps = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing);
        if (output_intermediate_steps.exists()) {
            //logger.logLine("Working directory was already used - please us another one or empty this one completely");
            //TODO: after debugging use system exit !!!
            //System.exit(1);
        }
        output_intermediate_steps.mkdir();
        File output_inter_steps_combined = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.folder_name_deseq2_preprocessing_combined);
        output_inter_steps_combined.mkdir();
        File output_inter_steps_single = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing + File.separator +
                options_intern.folder_name_deseq2_preprocessing_single);
        output_inter_steps_single.mkdir();
        File output_inter_steps_symbols_ensg_mean_counts = new File(
                options_intern.com2pose_working_directory + File.separator +
                        options_intern.folder_name_deseq2_preprocessing + File.separator +
                        options_intern.folder_name_deseq2_preprocessing_gene_symbols);
        output_inter_steps_symbols_ensg_mean_counts.mkdir();

        File folder = new File(options_intern.deseq2_input_directory);

        HashMap<String, String> ensg_symbol = new HashMap<>();
        BufferedReader br_ensg_symbol = new BufferedReader(new FileReader(new File(options_intern.tepic_ensg_symbol)));
        String line_ensg_symbol = br_ensg_symbol.readLine();
        while ((line_ensg_symbol = br_ensg_symbol.readLine()) != null) {
            String[] split = line_ensg_symbol.split("\t");
            if (split.length > 1) {
                ensg_symbol.put(split[0], split[1]);
            }
        }

        br_ensg_symbol.close();

        HashMap<String, HashSet<String>> timepoints_samples = new HashMap<>();

        //CREATE SINGLES
        for (File fileDir : folder.listFiles()) {
            if (fileDir.isDirectory()) {
                HashMap<Integer, Integer> mean_line_counts = new HashMap<>();
                int count_samples = 0;
                File group = new File(output_inter_steps_single.getAbsolutePath() + File.separator + fileDir.getName());
                group.mkdir();

                HashSet<String> x = new HashSet<>();
                timepoints_samples.put(group.getName(), x);
                for (File sample : fileDir.listFiles()) {
                    String name = sample.getName();
                    timepoints_samples.get(group.getName()).add(name);

                    BufferedReader br = new BufferedReader(new FileReader(sample));
                    BufferedWriter bw =
                            new BufferedWriter(new FileWriter(group.getAbsolutePath() + File.separator + name));
                    String line = br.readLine();
                    bw.write(line);
                    bw.newLine();
                    int count = 0;
                    while ((line = br.readLine()) != null) {
                        bw.write(row_ensg_name.get(count) + "\t" + line);
                        bw.newLine();
                        int count_line = Integer.parseInt(line);
                        if (mean_line_counts.containsKey(count)) {
                            int z = mean_line_counts.get(count);
                            z += count_line;
                            mean_line_counts.put(count, z);
                        } else {
                            mean_line_counts.put(count, count_line);
                        }
                        count++;
                    }
                    if (count != count_ensg_lines) {
                        logger.logLine("Error in nfcore RNA-seq data: File " + sample.getName() +
                                " has not the same number of rows as in File " + ensg_names.getName());
                        System.exit(1);
                    }
                    br.close();
                    bw.close();

                    count_samples++;
                }

                BufferedWriter bw_means = new BufferedWriter(new FileWriter(new File(
                        output_inter_steps_single.getAbsolutePath() + File.separator + group.getName() +
                                options_intern.file_suffix_deseq2_preprocessing_meanCounts)));
                bw_means.write(group.getName() + "_MEANS");
                bw_means.newLine();
                for (int i = 0; i < count_ensg_lines; i++) {
                    int mean_count = mean_line_counts.get(i) / count_samples;
                    bw_means.write("" + mean_count);
                    bw_means.newLine();

                }
                bw_means.close();

                BufferedReader br_ensg =
                        new BufferedReader(new FileReader(new File(options_intern.deseq2_input_gene_id)));
                BufferedWriter bw_symbol = new BufferedWriter(new FileWriter(
                        output_inter_steps_symbols_ensg_mean_counts.getAbsolutePath() + File.separator +
                                group.getName() + ".csv"));
                bw_symbol.write("SYMBOL\tENSG\tMEAN_COUNT");
                bw_symbol.newLine();

                String line = br_ensg.readLine();

                int i = 0;
                while ((line = br_ensg.readLine()) != null) {
                    StringBuilder sb = new StringBuilder();

                    int mean_count = mean_line_counts.get(i) / count_samples;
                    String gene_symbol_name = "NO_SYMBOL";

                    if (ensg_symbol.containsKey(line)) {
                        gene_symbol_name = ensg_symbol.get(line);
                    }

                    sb.append(gene_symbol_name);
                    sb.append("\t");
                    sb.append(line);
                    sb.append("\t");
                    sb.append(mean_count);


                    bw_symbol.write(sb.toString());
                    bw_symbol.newLine();
                    i++;
                }
                br_ensg.close();
                bw_symbol.close();


            }
        }
        //CREATE ALL COMBINED FILES
        HashSet<String> already_combined = new HashSet<>();
        for (String k : timepoints_samples.keySet()) {
            for (String kk : timepoints_samples.keySet()) {
                String key1 = k + "_" + kk;
                String key2 = kk + "_" + k;

                if (already_combined.contains(key1) || already_combined.contains(key2) || k.equals(kk)) {
                    continue;
                }

                already_combined.add(key1);

                HashSet<String> samples_group1 = timepoints_samples.get(k);
                HashSet<String> samples_group2 = timepoints_samples.get(kk);

                File combined_file = new File(output_inter_steps_combined.getAbsolutePath() + File.separator + key1);
                combined_file.mkdir();

                StringBuilder sb_header = new StringBuilder();
                sb_header.append("geneID");

                ArrayList<BufferedReader> bufferedReaders = new ArrayList<>();
                for (String s : samples_group1) {
                    sb_header.append("\t");
                    sb_header.append(s);
                    BufferedReader br = new BufferedReader(new FileReader(
                            new File(output_inter_steps_single + File.separator + k + File.separator + s)));
                    br.readLine();
                    bufferedReaders.add(br);
                }

                for (String s : samples_group2) {
                    sb_header.append("\t");
                    sb_header.append(s);
                    BufferedReader br = new BufferedReader(new FileReader(
                            new File(output_inter_steps_single + File.separator + kk + File.separator + s)));
                    br.readLine();
                    bufferedReaders.add(br);
                }

                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(
                        output_inter_steps_combined + File.separator + key1 + File.separator + key1 + ".csv")));
                bw.write(sb_header.toString());
                bw.newLine();

                BufferedReader br_ensg =
                        new BufferedReader(new FileReader(new File(options_intern.deseq2_input_gene_id)));
                String line = "";
                br_ensg.readLine();
                while ((line = br_ensg.readLine()) != null) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(line);
                    for (BufferedReader br : bufferedReaders) {
                        String line_intern = br.readLine();

                        String[] split = line_intern.split("\t");
                        sb.append("\t");
                        sb.append(split[1]);
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br_ensg.close();
                bw.close();
            }
        }

        logger.logLine("Finished preprocessing of nfcore RNA-seq data - no errors detected");
        logger.logLine("Started creating RScripts for running DESeq2");

        //CREATE RScripts for every folder in here
        File r_scripts = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_R_scripts);
        r_scripts.mkdir();
        for (File dir : output_inter_steps_combined.listFiles()) {
            if (dir.isDirectory()) {
                String name_rscript = dir.getName() + ".R";

                String group1 = dir.getName().split("_")[0];
                String group2 = dir.getName().split("_")[1];


                BufferedWriter bw = new BufferedWriter(
                        new FileWriter(new File(r_scripts.getAbsolutePath() + File.separator + name_rscript)));

                StringBuilder sb = new StringBuilder();
                sb.append("if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n");
                sb.append("  install.packages(\"BiocManager\")\n");
                sb.append("if (!requireNamespace(\"DESeq2\", quietly = TRUE))\n");
                sb.append("  BiocManager::install(\"DESeq2\")\n");
                sb.append("library(DESeq2)\n");
                sb.append("#" + group1 + " VS " + group2 + "\n");
                logger.logLine("[DESEQ2] Group " + group1 + " VS " + group2);

                BufferedReader br_header = new BufferedReader(
                        new FileReader(dir.getAbsolutePath() + File.separator + dir.getName() + ".csv"));
                String line_header = br_header.readLine();
                br_header.close();

                ArrayList<String> groups = new ArrayList<>();
                ArrayList<String> samples = new ArrayList<>();

                String[] split_header = line_header.split("\t");
                for (int i = 1; i < split_header.length; i++) {
                    String[] split_samples = split_header[i].split("-");
                    String[] split_groups = split_samples[0].split("_");

                    samples.add(split_samples[0]);
                    groups.add(split_groups[0]);
                }

                sb.append("metadata_df<-data.frame(sample_id = c(");
                int count_s = 0;
                for (String s : samples) {
                    if (count_s > 0) {
                        sb.append(" ,\"");
                    } else {
                        sb.append("\"");
                    }
                    sb.append(s);
                    sb.append("\"");
                    count_s++;
                }
                sb.append("), group = c(");
                count_s = 0;
                for (String s : groups) {
                    if (count_s > 0) {
                        sb.append(" ,\"");
                    } else {
                        sb.append("\"");
                    }
                    sb.append(s);
                    sb.append("\"");
                    count_s++;
                }
                sb.append("))\n");

                sb.append("input_groups = \"");
                sb.append(dir.getName());
                sb.append("\"\n");

                sb.append("group_one = strsplit(input_groups, \"_\")[[1]][1]\n");
                sb.append("group_two = strsplit(input_groups, \"_\")[[1]][2]\n");

                sb.append("rownames(metadata_df) <- metadata_df$sample_id\n");
                sb.append("metadata_df$sample_id <- NULL\n");

                sb.append("count_path = \"");
                sb.append(dir.getAbsolutePath());
                sb.append(File.separator + dir.getName());
                sb.append(".csv\"\n");
                sb.append("count_df = read.csv(count_path, sep = \"\\t\", header = T, row.names = NULL)\n");
                sb.append("count_df=count_df[!duplicated(count_df$geneID), ]\n" +
                        "row.names(count_df) <- count_df[,1]\n" + "count_df$geneID<-NULL\n");

                sb.append("dds <- DESeqDataSetFromMatrix(countData=count_df, \n");
                sb.append("                              colData=metadata_df, \n");
                sb.append("                              design=~group)\n");

                if (options_intern.deseq2_count_threshold > 0) {
                    sb.append("threshold = ");
                    sb.append(options_intern.deseq2_count_threshold);
                    sb.append("\n");
                    sb.append("keep <- rowSums(counts(dds)) >= threshold\n");
                    sb.append("dds <- dds[keep,]\n");
                }

                File output_deseq2 = new File(options_intern.com2pose_working_directory + File.separator +
                        options_intern.folder_name_deseq2_output_raw);
                output_deseq2.mkdir();

                sb.append("output_path = \"");
                sb.append(output_deseq2.getAbsolutePath());
                sb.append(File.separator);
                sb.append(dir.getName());
                sb.append("_raw_output_deseq2.tsv");
                sb.append("\"\n");

                sb.append("dds <- DESeq(dds)\n");
                sb.append("res <- results(dds)\n");
                sb.append("summary(res)\n");
                sb.append("write.table(res,file=output_path,sep=\"\\t\")\n");

                bw.write(sb.toString());
                bw.close();
            }
        }

        logger.logLine("Finished creating RScripts for running DESeq2");


    }

    public void get_ensg_symbol_mapping() throws Exception {

        File script_output_dir = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_deseq2_preprocessing);
        script_output_dir.mkdir();
        File script_output = new File(script_output_dir.getAbsolutePath() + File.separator + "ENSG_SYMBOL_MAP.R");

        logger.logLine("[DESEQ2] Start mapping ENSG to GENE SYMBOLS.");
        logger.logLine("[DESEQ2] if not done in 2 hours:");
        logger.logLine("[DESEQ2] 1. Please stop com2pose.");
        logger.logLine("[DESEQ2] 2. Please run script manually in RStudio.");
        logger.logLine("[DESEQ2] Script path: " + script_output.getAbsolutePath());
        logger.logLine(
                "[DESEQ2] 3. Afterwards add paramter -m to com2pose command line, so this script wont be started again.");
        logger.logLine("[DESEQ2] 4. restart com2pose");
        logger.logLine("[DESEQ2] waiting ...");


        File results = new File(
                script_output_dir.getAbsolutePath() + File.separator + options_intern.file_suffix_deseq2_mapping);

        StringBuilder sb = new StringBuilder();
        sb.append(
                "if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n" + "  install.packages(\"BiocManager\")\n" +
                        "\n" + "if (!requireNamespace(\"biomaRt\", quietly = TRUE))\n" +
                        "  BiocManager::install(\"biomaRt\")\n" + "\n" + "library('biomaRt')\n" +
                        "httr::set_config(httr::config(ssl_verifypeer = FALSE))\n" + "df <- read.csv('" +
                        options_intern.deseq2_input_gene_id + "')\n" + "\n" + "not_done=TRUE\n" + "\n" +
                        "G_list= data.frame()\n" + "\n" + "while(not_done)\n" + "{\n" + "  tryCatch({\n" +
                        "    mart <- useDataset(\"" + options_intern.deseq2_biomart_dataset_species +
                        "\", useMart(\"ensembl\"))\n" + "    df$id <- NA\n" +
                        "    G_list_intern <- getBM(filters= \"ensembl_gene_id\", attributes= c(\"ensembl_gene_id\",\"" +
                        options_intern.deseq2_biomart_dataset_symbol_column + "\"),values=df$Geneid,mart= mart)\n" +
                        "    G_list=rbind(G_list,G_list_intern)\n" + "    not_done=FALSE\n" +
                        "  }, warning = function(w) {\n" + "    print(\"WARNING SECTION\")\n" + "    print(w)\n" +
                        "  }, error = function(e) {\n" + "    print(\"ERROR SECTION\")\n" + "    print(e)\n" +
                        "  }, finally = {\n" + "  })\n" + "}\n" + "\n" + "write.table(G_list,\"" +
                        results.getAbsolutePath() + "\", row.names = FALSE, quote = F, sep=\"\\t\")\n");

        /*
        StringBuilder sb = new StringBuilder();
        sb.append("if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n" +
                "  install.packages(\"BiocManager\")\n" +
                "\n" +
                "if (!requireNamespace(\"biomaRt\", quietly = TRUE))\n" +
                "  BiocManager::install(\"biomaRt\")\n" +
                "\n" +
                "library('biomaRt')\n");

        sb.append("httr::set_config(httr::config(ssl_verifypeer = FALSE))\n");

        sb.append("df <- read.csv('"+options_intern.deseq2_input_gene_id+"')\n");

        sb.append("mart <- useDataset(\""+options_intern.deseq2_biomart_dataset_species+"\", useMart(\"ensembl\"))");
        sb.append("\n");

        sb.append("df$id <- NA\n" +
                "G_list <- getBM(filters= \"ensembl_gene_id\", attributes= c(\"ensembl_gene_id\",\""+options_intern.deseq2_biomart_dataset_symbol_column+"\"),values=df$Geneid,mart= mart)\n" +
                "write.table(G_list,\""+results.getAbsolutePath()+"\", row.names = FALSE, quote = F, sep=\"\\t\")\n");*/


        BufferedWriter bw = new BufferedWriter(new FileWriter(script_output));
        bw.write(sb.toString());
        bw.close();

        String command = "Rscript " + script_output.getAbsolutePath();
        Process child = Runtime.getRuntime().exec(command);
        logger.logLine("[DESEQ2] Running script " + script_output.getName() + ": " + command);
        int code = child.waitFor();
        switch (code) {
            case 0:
                break;
            case 1:
                String message = child.getErrorStream().toString();
                logger.logLine(
                        "[DESEQ2] Script failed due to bioMart connection error. Please run script manually in RStudio.");
                logger.logLine("[DESEQ2] Script path: " + script_output.getAbsolutePath());
                logger.logLine(
                        "[DESEQ2] Afterwards add paramter -m to com2pose command line, so this script wont be started again.");
                System.exit(1);
        }

        options_intern.tepic_ensg_symbol = results.getAbsolutePath();


        logger.logLine("[DESEQ2] Finished mapping ENSG to GENE SYMBOLS.");

    }

    public void mix_mutually_exclusive_peaks() throws IOException {
        logger.logLine("[MUTUALLY-EXCLUSIVE-PEAKS] Start mutually exclusive peaks calculation.");
        logger.logLine("[MUTUALLY-EXCLUSIVE-PEAKS] Preprocessing mutually exclusive peaks for binary tree comparison.");

        logger.logLine("[MUTUALLY-EXCLUSIVE-PEAKS] Used data: " + options_intern.tepic_input_directory);

        options_intern.tepic_input_prev = options_intern.tepic_input_directory;
        File f_annotation_check = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option +
                        File.separator + options_intern.folder_name_mix_option_preprocessing_check_chr);
        options_intern.tepic_input_directory = f_annotation_check.getAbsolutePath();

        if (options_intern.mix_level.equals("SAMPLE_LEVEL")) {
            File root_mix_working_dir = new File(
                    options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option);
            File f_sample_mix_output = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                    options_intern.folder_name_mix_option_sample_mix);
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = f_sample_mix_output.getAbsolutePath();
        }

        if (options_intern.mix_level.equals("HM_LEVEL")) {
            File root_mix_working_dir = new File(
                    options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option);
            File f_output_hm = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                    options_intern.folder_name_mix_option_hm_mix);
            f_output_hm.mkdir();

            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = f_output_hm.getAbsolutePath();
        }

        if (options_intern.tepic_tf_binding_site_search.equals("BETWEEN")) {
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory =
                    options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option +
                            File.separator + options_intern.folder_name_mix_options_footprints_between_peaks;

        }

        if (!options_intern.black_list_dir.equals("")) {
            File output_folder = new File(options_intern.com2pose_working_directory + File.separator +
                    options_intern.folder_name_blacklisted_regions);
            File output_folder_new_input = new File(output_folder.getAbsolutePath() + File.separator +
                    options_intern.folder_name_blacklisted_regions_new_input);
            output_folder_new_input.mkdir();

            //set new folder directory for tepic input and save old one
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = output_folder_new_input.getAbsolutePath();
        }


        File f_output_mix_option = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option);
        f_output_mix_option.mkdir();

        File f_output_mix_option_mutually_exclusive = new File(f_output_mix_option.getAbsolutePath() + File.separator +
                options_intern.folder_name_mix_option_mutually_exclusive);
        f_output_mix_option_mutually_exclusive.mkdir();

        File f_output_mix_option_mutually_exclusive_preprocessing = new File(
                f_output_mix_option_mutually_exclusive.getAbsolutePath() + File.separator +
                        options_intern.folder_name_mix_option_mutually_exclusive_preprocessing);
        f_output_mix_option_mutually_exclusive_preprocessing.mkdir();

        File f_output_mix_option_mutually_exclusive_new_input = new File(
                f_output_mix_option_mutually_exclusive.getAbsolutePath() + File.separator +
                        options_intern.folder_name_mix_options_mutually_exclusive_input);
        f_output_mix_option_mutually_exclusive_new_input.mkdir();

        HashMap<String, HashSet<String>> tps_to_hms = new HashMap<>();

        File f_input_peaks = new File(options_intern.tepic_input_directory);

        for (File fileDir_tp : f_input_peaks.listFiles()) {
            if (fileDir_tp.isDirectory()) {
                HashSet<String> hms;
                if (tps_to_hms.containsKey(fileDir_tp.getName())) {
                    hms = tps_to_hms.get(fileDir_tp.getName());
                } else {
                    hms = new HashSet<>();
                }

                File output_tp = new File(
                        f_output_mix_option_mutually_exclusive_preprocessing.getAbsolutePath() + File.separator +
                                fileDir_tp.getName());
                output_tp.mkdir();

                for (File fileDir_hm : fileDir_tp.listFiles()) {
                    if (fileDir_hm.isDirectory()) {
                        hms.add(fileDir_hm.getName());

                        File output_hm = new File(output_tp.getAbsolutePath() + File.separator + fileDir_hm.getName());
                        output_hm.mkdir();

                        //we can only use this option with one input file
                        if (fileDir_hm.listFiles().length > 1) {
                            logger.logLine(
                                    "[MUTUALLY-EXCLUSIVE-PEAKS] A UNION OVER SAMPLES MUST BE MADE! Start MIX-OPTION now");
                            logger.logLine(
                                    "[MUTUALLY-EXCLUSIVE-PEAKS] Setting parameters to \"UNION, SAMPLE_LEVEL\" now.");

                            options_intern.mix_option = "UNION";
                            options_intern.mix_level = "SAMPLE_LEVEL";

                            mix_option();
                            if (!options_intern.black_list_dir.equals("")) {
                                logger.logLine("[MUTUALLY-EXCLUSIVE-PEAKS] Redo blacklist on unioned samples");
                                preprocess_blacklist();
                                filter_blacklist();
                            }
                            mix_mutually_exclusive_peaks();
                        }

                        BufferedReader br_peaks = new BufferedReader(new FileReader(fileDir_hm.listFiles()[0]));
                        String line_peaks = "";
                        String current_chr = "-1";
                        int position = 0;
                        ArrayList<BL_ranges_binary_tree> current_chr_ranges = new ArrayList<>();
                        while ((line_peaks = br_peaks.readLine()) != null) {
                            String[] split = line_peaks.split("\t");

                            if (!current_chr.equals(split[0])) {
                                if (!current_chr.equals("-1")) {
                                    ArrayList<BL_ranges_binary_tree> newly_ordered = new ArrayList<>();
                                    //recursive_split
                                    current_chr_ranges = recursive_split_BL(current_chr_ranges, newly_ordered);
                                    //write_file
                                    BufferedWriter bw_peaks = new BufferedWriter(new FileWriter(new File(
                                            output_hm.getAbsolutePath() + File.separator + current_chr + ".txt")));
                                    bw_peaks.append("#\tCHR\tLEFT_BORDER\tRIGHT_BORDER\tPEAK_SCORE");
                                    bw_peaks.newLine();
                                    for (BL_ranges_binary_tree bl : newly_ordered) {
                                        bw_peaks.write(bl.toString());
                                        bw_peaks.newLine();
                                    }
                                    bw_peaks.close();


                                    current_chr = split[0];
                                    current_chr_ranges = new ArrayList<>();
                                    position = 0;
                                    continue;
                                } else {
                                    current_chr = split[0];
                                }
                            }

                            BL_ranges_binary_tree range = new BL_ranges_binary_tree();
                            range.left_border = Integer.parseInt(split[1]);
                            range.right_border = Integer.parseInt(split[2]);
                            range.chr = current_chr;
                            //range.signal=split[0];
                            range.number = position;
                            range.peak_score = Double.parseDouble(split[4]);

                            current_chr_ranges.add(range);

                            position++;


                        }
                        br_peaks.close();
                        ArrayList<BL_ranges_binary_tree> newly_ordered = new ArrayList<>();
                        //recursive_split
                        current_chr_ranges = recursive_split_BL(current_chr_ranges, newly_ordered);
                        //write_file
                        BufferedWriter bw_peaks = new BufferedWriter(new FileWriter(
                                new File(output_hm.getAbsolutePath() + File.separator + current_chr + ".txt")));
                        bw_peaks.append("#\tCHR\tLEFT_BORDER\tRIGHT_BORDER");
                        bw_peaks.newLine();
                        for (BL_ranges_binary_tree bl : newly_ordered) {
                            bw_peaks.write(bl.toString());
                            bw_peaks.newLine();
                        }
                        bw_peaks.close();
                    }
                }
                tps_to_hms.put(fileDir_tp.getName(), hms);
            }
        }
        logger.logLine("[MUTUALLY-EXCLUSIVE-PEAKS] Finished preprocessing. Starting binary tree filtering now.");

        HashSet<String> already_worked_tp_groups = new HashSet<>();

        for (String key_tp_1 : tps_to_hms.keySet()) {
            for (String key_tp_2 : tps_to_hms.keySet()) {
                if (key_tp_1.equals(key_tp_2)) {
                    continue;
                }

                String key_clash1 = key_tp_1 + "_" + key_tp_2;
                String key_clash2 = key_tp_2 + "_" + key_tp_1;

                if (already_worked_tp_groups.contains(key_clash1) || already_worked_tp_groups.contains(key_clash2)) {
                    continue;
                } else {
                    already_worked_tp_groups.add(key_clash1);
                }

                File clash_output = new File(
                        f_output_mix_option_mutually_exclusive_new_input.getAbsolutePath() + File.separator +
                                key_clash1);
                clash_output.mkdir();

                HashSet<String> hms_tp1 = tps_to_hms.get(key_tp_1);
                HashSet<String> hms_tp2 = tps_to_hms.get(key_tp_2);

                for (String key_hm : hms_tp1) {
                    if (hms_tp2.contains(key_hm)) {
                        File clash_output_hm = new File(clash_output.getAbsolutePath() + File.separator + key_hm);
                        clash_output_hm.mkdir();

                        //output input files
                        File f_input_tp1_dir = new File(
                                f_input_peaks.getAbsolutePath() + File.separator + key_tp_1 + File.separator + key_hm);
                        File f_input_tp2_dir = new File(
                                f_input_peaks.getAbsolutePath() + File.separator + key_tp_2 + File.separator + key_hm);
                        File f_input_tp1 = f_input_tp1_dir.listFiles()[0];
                        File f_input_tp2 = f_input_tp2_dir.listFiles()[0];


                        //build binary tree for tp1
                        HashMap<String, BL_binary_tree> chr_tree_tp1 = new HashMap<>();
                        File input_folder_tp1 = new File(
                                f_output_mix_option_mutually_exclusive_preprocessing.getAbsolutePath() +
                                        File.separator + key_tp_1 + File.separator + key_hm);
                        for (File fileChr : input_folder_tp1.listFiles()) {
                            if (fileChr.isFile()) {
                                String name = fileChr.getName().split("\\.")[0];

                                ArrayList<BL_ranges_binary_tree> region = new ArrayList<>();

                                BufferedReader br_chr = new BufferedReader(new FileReader(fileChr));
                                String line_chr = br_chr.readLine();
                                while ((line_chr = br_chr.readLine()) != null) {
                                    String[] split = line_chr.split("\t");

                                    BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
                                    iu.number = Integer.parseInt(split[0]);
                                    iu.chr = split[1];
                                    iu.left_border = Integer.parseInt(split[2]);
                                    iu.right_border = Integer.parseInt(split[3]);
                                    //iu.signal = split[4];
                                    iu.peak_score = Double.parseDouble(split[4]);

                                    region.add(iu);
                                }
                                br_chr.close();

                                BL_binary_tree_node root = new BL_binary_tree_node(region.get(0), region.get(0).number);
                                BL_binary_tree tree = new BL_binary_tree(root);

                                for (int i = 1; i < region.size(); i++) {
                                    tree.add(region.get(i).number, region.get(i));
                                    tree.peak_signals.add(region.get(i).peak_score);
                                }

                                tree.calcualte_average_peak_score();

                                chr_tree_tp1.put(name, tree);

                            }
                        }

                        //filter tp2
                        BufferedWriter bw_tp2 = new BufferedWriter(new FileWriter(
                                clash_output_hm.getAbsolutePath() + File.separator + f_input_tp2.getName()));
                        BufferedReader br_tp2 = new BufferedReader(new FileReader(f_input_tp2));
                        String line_tp2 = "";
                        int count_line_tp_2 = 0;
                        while ((line_tp2 = br_tp2.readLine()) != null) {
                            String[] split = line_tp2.split("\t");

                            String chr = split[0];

                            BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
                            iu.chr = chr;
                            iu.left_border = Integer.parseInt(split[1]);
                            iu.right_border = Integer.parseInt(split[2]);
                            iu.peak_score = Double.parseDouble(split[4]);

                            if (!chr_tree_tp1.containsKey(chr)) {
                                count_line_tp_2++;
                                continue;
                            }
                            BL_binary_tree tree = chr_tree_tp1.get(chr);

                            if (options_intern.mix_mutually_exclusive_diff_peak_signals) {
                                BL_ranges_binary_tree current_match = tree.containsNode(iu);
                                if (current_match == null) {
                                    //this peak is mutually exclusive
                                    bw_tp2.write(line_tp2);
                                    bw_tp2.newLine();
                                    continue;
                                }

                                double peak_difference = current_match.peak_score - iu.peak_score;

                                if (peak_difference < 0) {
                                    peak_difference *= -1;
                                }

                                if (peak_difference > tree.average_peak_score) {
                                    //this peak is mutually exclusive in peak score
                                    bw_tp2.write(line_tp2);
                                    bw_tp2.newLine();
                                }
                            } else {
                                if (tree.containsNode(iu) == null) {
                                    //this peak is mutually exclusive
                                    bw_tp2.write(line_tp2);
                                    bw_tp2.newLine();
                                }
                            }


                        }
                        br_tp2.close();
                        bw_tp2.close();

                        //build binary tree fpr tp2
                        HashMap<String, BL_binary_tree> chr_tree_tp2 = new HashMap<>();
                        File input_folder_tp2 = new File(
                                f_output_mix_option_mutually_exclusive_preprocessing.getAbsolutePath() +
                                        File.separator + key_tp_2 + File.separator + key_hm);
                        for (File fileChr : input_folder_tp2.listFiles()) {
                            if (fileChr.isFile()) {
                                String name = fileChr.getName().split("\\.")[0];

                                ArrayList<BL_ranges_binary_tree> region = new ArrayList<>();

                                BufferedReader br_chr = new BufferedReader(new FileReader(fileChr));
                                String line_chr = br_chr.readLine();
                                while ((line_chr = br_chr.readLine()) != null) {
                                    String[] split = line_chr.split("\t");

                                    BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
                                    iu.number = Integer.parseInt(split[0]);
                                    iu.chr = split[1];
                                    iu.left_border = Integer.parseInt(split[2]);
                                    iu.right_border = Integer.parseInt(split[3]);
                                    //iu.signal = split[4];
                                    iu.peak_score = Double.parseDouble(split[4]);

                                    region.add(iu);
                                }
                                br_chr.close();

                                BL_binary_tree_node root = new BL_binary_tree_node(region.get(0), region.get(0).number);
                                BL_binary_tree tree = new BL_binary_tree(root);

                                for (int i = 1; i < region.size(); i++) {
                                    tree.add(region.get(i).number, region.get(i));
                                    tree.peak_signals.add(region.get(i).peak_score);
                                }
                                tree.calcualte_average_peak_score();

                                chr_tree_tp2.put(name, tree);

                            }
                        }


                        //filter tp1
                        BufferedWriter bw_tp1 = new BufferedWriter(new FileWriter(
                                clash_output_hm.getAbsolutePath() + File.separator + f_input_tp1.getName()));
                        BufferedReader br_tp1 = new BufferedReader(new FileReader(f_input_tp1));
                        String line_tp1 = "";
                        int count_line_tp_1 = 0;
                        while ((line_tp1 = br_tp1.readLine()) != null) {
                            String[] split = line_tp1.split("\t");

                            String chr = split[0];

                            BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
                            iu.chr = chr;
                            iu.left_border = Integer.parseInt(split[1]);
                            iu.right_border = Integer.parseInt(split[2]);
                            iu.peak_score = Double.parseDouble(split[4]);

                            if (!chr_tree_tp1.containsKey(chr)) {
                                count_line_tp_2++;
                                continue;
                            }
                            BL_binary_tree tree = chr_tree_tp2.get(chr);

                            if (options_intern.mix_mutually_exclusive_diff_peak_signals) {
                                BL_ranges_binary_tree current_match = tree.containsNode(iu);

                                if (current_match == null) {
                                    //this peak is mutually exclusive
                                    bw_tp1.write(line_tp1);
                                    bw_tp1.newLine();
                                    continue;
                                }

                                double peak_difference = current_match.peak_score - iu.peak_score;

                                if (peak_difference < 0) {
                                    peak_difference *= -1;
                                }

                                if (peak_difference > tree.average_peak_score) {
                                    //this peak is mutually exclusive
                                    bw_tp1.write(line_tp1);
                                    bw_tp1.newLine();
                                }
                            } else {
                                if (tree.containsNode(iu) == null) {
                                    //this peak is mutually exclusive
                                    bw_tp1.write(line_tp1);
                                    bw_tp1.newLine();
                                }
                            }
                        }
                        br_tp1.close();
                        bw_tp1.close();

                    }
                }
            }
        }

        logger.logLine("[MUTUALLY-EXCLUSIVE-PEAKS] Finished binary tree filtering.");
    }

    /**
     * Filter TEPIC input files for blacklisted regions
     */
    public void filter_blacklist() throws IOException {

        options_intern.tepic_input_prev = options_intern.tepic_input_directory;
        File f_annotation_check = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option +
                        File.separator + options_intern.folder_name_mix_option_preprocessing_check_chr);
        options_intern.tepic_input_directory = f_annotation_check.getAbsolutePath();

        if (options_intern.tepic_tf_binding_site_search.equals("BETWEEN")) {
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory =
                    options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option +
                            File.separator + options_intern.folder_name_mix_options_footprints_between_peaks;

        }
        File folder_input = new File(options_intern.tepic_input_directory);
        File output_folder = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_blacklisted_regions);
        File output_folder_new_input = new File(output_folder.getAbsolutePath() + File.separator +
                options_intern.folder_name_blacklisted_regions_new_input);
        output_folder_new_input.mkdir();

        logger.logLine("[BLACKLIST] Used input: " + options_intern.tepic_input_directory);

        //set new folder directory for tepic input and save old one
        options_intern.tepic_input_prev = options_intern.tepic_input_directory;
        options_intern.tepic_input_directory = output_folder_new_input.getAbsolutePath();


        logger.logLine("[BLACKLIST] Create chromosome binary trees.");
        //CREATE BINARY TREES
        HashMap<String, BL_binary_tree> chr_tree = new HashMap<>();

        File input_folder_chr = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_blacklisted_regions + File.separator +
                options_intern.folder_name_blacklisted_regions_preprocessing + File.separator +
                options_intern.folder_name_blacklisted_regions_preprocessing_sorted);
        for (File fileChr : input_folder_chr.listFiles()) {
            if (fileChr.isFile()) {
                String name = fileChr.getName().split("\\.")[0];

                ArrayList<BL_ranges_binary_tree> region = new ArrayList<>();

                BufferedReader br_chr = new BufferedReader(new FileReader(fileChr));
                String line_chr = br_chr.readLine();
                while ((line_chr = br_chr.readLine()) != null) {
                    String[] split = line_chr.split("\t");

                    BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
                    iu.number = Integer.parseInt(split[0]);
                    iu.chr = split[1];
                    iu.left_border = Integer.parseInt(split[2]);
                    iu.right_border = Integer.parseInt(split[3]);
                    iu.signal = split[4];

                    region.add(iu);
                }
                br_chr.close();

                BL_binary_tree_node root = new BL_binary_tree_node(region.get(0), region.get(0).number);
                BL_binary_tree tree = new BL_binary_tree(root);

                for (int i = 1; i < region.size(); i++) {
                    tree.add(region.get(i).number, region.get(i));
                }

                chr_tree.put(name, tree);

            }
        }

        logger.logLine("[BLACKLIST] Filter input files for blacklisted regions.");
        //now filter all files
        int count_matched = 0;

        for (File fileDirTP : folder_input.listFiles()) {
            if (fileDirTP.isDirectory()) {
                File output_folder_new_input_TP =
                        new File(output_folder_new_input.getAbsolutePath() + File.separator + fileDirTP.getName());
                output_folder_new_input_TP.mkdir();

                for (File fileDirTP_HM : fileDirTP.listFiles()) {
                    if (fileDirTP_HM.isDirectory()) {
                        File output_folder_new_input_TP_HM = new File(
                                output_folder_new_input_TP.getAbsolutePath() + File.separator + fileDirTP_HM.getName());
                        output_folder_new_input_TP_HM.mkdir();

                        for (File filrDirTP_HM_sample : fileDirTP_HM.listFiles()) {
                            if (filrDirTP_HM_sample.isFile()) {
                                logger.logLine(
                                        "[BLACKLIST] Filter: " + fileDirTP.getName() + ": " + fileDirTP_HM.getName() +
                                                " - " + filrDirTP_HM_sample.getName());
                                //FILTER HERE WITH BINARY TREE

                                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(
                                        output_folder_new_input_TP_HM.getAbsolutePath() + File.separator +
                                                filrDirTP_HM_sample.getName())));
                                BufferedReader br = new BufferedReader(new FileReader(filrDirTP_HM_sample));
                                String line = "";
                                int count_line = 0;
                                while ((line = br.readLine()) != null) {
                                    String[] split = line.split("\t");

                                    String chr = split[0];
                                    if (!chr.matches(".*chr.*")) {
                                        chr = "chr" + chr;
                                    }

                                    BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
                                    iu.chr = chr;
                                    iu.left_border = Integer.parseInt(split[1]);
                                    iu.right_border = Integer.parseInt(split[2]);

                                    if (!chr_tree.containsKey(chr)) {
                                        count_line++;
                                        continue;
                                    }

                                    BL_binary_tree tree = chr_tree.get(chr);
                                    if (tree.containsNode(iu) == null) {
                                        bw.write(line);
                                        bw.newLine();
                                    } else {
                                        count_matched++;
                                    }
                                    count_line++;
                                }
                                br.close();
                                bw.close();
                            }
                        }
                    }
                }
            }
        }
        logger.logLine("[BLACKLIST] Finished filtering for blacklisted regions.");
    }

    /**
     * preprocesses the blacklist file for binary search
     */
    public void preprocess_blacklist() throws IOException {

        logger.logLine("[BLACKLIST] start preprocessing blacklist");

        File f_blacklist = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_blacklisted_regions);
        f_blacklist.mkdir();
        File f_blacklist_pre = new File(f_blacklist.getAbsolutePath() + File.separator +
                options_intern.folder_name_blacklisted_regions_preprocessing);
        f_blacklist_pre.mkdir();
        File f_blacklist_pre_chr = new File(f_blacklist_pre.getAbsolutePath() + File.separator +
                options_intern.folder_name_blacklisted_regions_preprocessing_perChr);
        f_blacklist_pre_chr.mkdir();

        if (options_intern.mix_option.equals("SAMPLE_LEVEL")) {
            File root_mix_working_dir = new File(
                    options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option);
            File f_sample_mix_output = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                    options_intern.folder_name_mix_option_sample_mix);
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = f_sample_mix_output.getAbsolutePath();
        }

        if (options_intern.mix_option.equals("HM_LEVEL")) {
            File root_mix_working_dir = new File(
                    options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option);
            File f_output_hm = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                    options_intern.folder_name_mix_option_hm_mix);
            f_output_hm.mkdir();

            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = f_output_hm.getAbsolutePath();
        }

        BufferedReader br_chr = new BufferedReader(new FileReader(new File(options_intern.black_list_dir)));
        ArrayList<BL_ranges_binary_tree> chr_ius = new ArrayList<>();

        String line_chr = "";
        String current_chr = "";
        BufferedWriter bw_chr = new BufferedWriter(
                new FileWriter(new File(f_blacklist_pre_chr.getAbsolutePath() + File.separator + "test.txt")));
        while ((line_chr = br_chr.readLine()) != null) {
            String[] split = line_chr.split("\t");

            String chr = split[0];
            if (!chr.equals(current_chr)) {
                Collections.sort(chr_ius);

                int i = 0;

                for (BL_ranges_binary_tree iu : chr_ius) {
                    iu.number = i;
                    bw_chr.write(iu.toString());
                    bw_chr.newLine();
                    i++;
                }

                chr_ius.clear();

                if (!chr.matches(".*chr.*")) {
                    chr = "chr" + chr;
                }

                bw_chr.close();
                bw_chr = new BufferedWriter(new FileWriter(
                        new File(f_blacklist_pre_chr.getAbsolutePath() + File.separator + chr + ".txt")));
                bw_chr.write("#\tCHR\tLEFT_BORDER\tRIGHT_BORDER\tSIGNAL");
                bw_chr.newLine();
                current_chr = split[0];
            }

            BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
            iu.chr = chr;
            iu.left_border = Integer.parseInt(split[1]);
            iu.right_border = Integer.parseInt(split[2]);
            iu.signal = split[3].replace(' ', '_').toUpperCase();

            if (options_intern.black_list_signals.contains(iu.signal)) {
                chr_ius.add(iu);
            }
        }

        Collections.sort(chr_ius);

        int i = 0;

        for (BL_ranges_binary_tree iu : chr_ius) {
            iu.number = i;
            bw_chr.write(iu.toString());
            bw_chr.newLine();
            i++;
        }

        bw_chr.close();
        br_chr.close();

        logger.logLine("[BLACKLIST] prepare chromosomes for binary tree");

        File f_blacklist_pre_sorted = new File(f_blacklist_pre.getAbsolutePath() + File.separator +
                options_intern.folder_name_blacklisted_regions_preprocessing_sorted);
        f_blacklist_pre_sorted.mkdir();

        for (File fileDir : f_blacklist_pre_chr.listFiles()) {
            if (fileDir.isFile() && !fileDir.getName().equals("test.txt")) {
                ArrayList<BL_ranges_binary_tree> current_ius = new ArrayList<>();
                String header;

                BufferedReader br_sort = new BufferedReader(new FileReader(fileDir));
                header = br_sort.readLine();
                String line_sort = "";
                while ((line_sort = br_sort.readLine()) != null) {
                    String[] split = line_sort.split("\t");

                    BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
                    iu.number = Integer.parseInt(split[0]);
                    iu.chr = split[1];
                    iu.left_border = Integer.parseInt(split[2]);
                    iu.right_border = Integer.parseInt(split[3]);
                    iu.signal = split[4];

                    current_ius.add(iu);
                }
                br_sort.close();

                ArrayList<BL_ranges_binary_tree> newly_ordered = new ArrayList<>();

                recursive_split_BL(current_ius, newly_ordered);

                BufferedWriter bw_sort = new BufferedWriter(
                        new FileWriter(f_blacklist_pre_sorted.getAbsolutePath() + File.separator + fileDir.getName()));
                bw_sort.write(header);
                bw_sort.newLine();

                for (int j = 0; j < newly_ordered.size(); j++) {
                    bw_sort.write(newly_ordered.get(j).toString());
                    bw_sort.newLine();
                }
                bw_sort.close();
            }

        }


        logger.logLine("[BLACKLIST] finished preprocessing blacklist");


    }

    /**
     * CREATES FOOTPRINTS OF GAPS BETWEEN PEAKS
     *
     * @throws IOException
     */
    public void create_footprints_between_peaks() throws IOException {
        logger.logLine(
                "[FOOTPRINTS] OPTION tepic_tf_binding_site_search=\"" + options_intern.tepic_tf_binding_site_search +
                        "\" was set!");
        logger.logLine("[FOOTPRINTS] creating footprints");


        File root_mix_working_dir = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option);
        root_mix_working_dir.mkdir();

        File f_sample_mix_preprocess = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                options_intern.folder_name_mix_option_sample_mix_preprocessing);
        f_sample_mix_preprocess.mkdir();

        File f_sample_mix_output = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                options_intern.folder_name_mix_option_sample_mix);
        f_sample_mix_output.mkdir();

        if (options_intern.mix_level.equals("SAMPLE_LEVEL")) {
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = f_sample_mix_output.getAbsolutePath();
        }

        logger.logLine("[FOOTPRINTS] Used data: " + options_intern.tepic_input_directory);

        options_intern.tepic_input_prev = options_intern.tepic_input_directory;
        File f_annotation_check = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option +
                        File.separator + options_intern.folder_name_mix_option_preprocessing_check_chr);
        options_intern.tepic_input_directory = f_annotation_check.getAbsolutePath();

        File folder_input = new File(options_intern.tepic_input_directory);

        File output_folder = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option);
        File output_folder_new_input = new File(output_folder.getAbsolutePath() + File.separator +
                options_intern.folder_name_mix_options_footprints_between_peaks);
        output_folder_new_input.mkdirs();

        //set new folder directory for tepic input and save old one
        options_intern.tepic_input_prev = options_intern.tepic_input_directory;
        options_intern.tepic_input_directory = output_folder_new_input.getAbsolutePath();


        for (File f_tp : folder_input.listFiles()) {
            if (f_tp.isDirectory()) {
                String key_tp = f_tp.getName();
                File f_output_tp = new File(output_folder_new_input.getAbsolutePath() + File.separator + key_tp);
                f_output_tp.mkdir();

                for (File f_hm : f_tp.listFiles()) {
                    if (f_hm.isDirectory()) {
                        String key_hm = f_hm.getName();
                        File f_output_tp_hm = new File(f_output_tp.getAbsolutePath() + File.separator + key_hm);
                        f_output_tp_hm.mkdir();

                        for (File f_sample : f_hm.listFiles()) {
                            if (f_sample.isFile()) {
                                File f_output_peak_file = new File(
                                        f_output_tp_hm.getAbsolutePath() + File.separator + f_sample.getName());

                                //read in all regions in array list per chromosome
                                HashMap<String, ArrayList<Footprint_Interval>> chr_regions = new HashMap<>();

                                BufferedReader br_input_regions_chr = new BufferedReader(new FileReader(f_sample));
                                String line_input_regions_chr = "";
                                while ((line_input_regions_chr = br_input_regions_chr.readLine()) != null) {
                                    String[] split = line_input_regions_chr.split("\t");
                                    if (chr_regions.containsKey(split[0])) {
                                        continue;
                                    } else {
                                        ArrayList<Footprint_Interval> fp_al = new ArrayList<>();
                                        chr_regions.put(split[0], fp_al);
                                    }

                                }
                                br_input_regions_chr.close();

                                BufferedReader br_input_regions = new BufferedReader(new FileReader(f_sample));
                                String line_input_regions = "";
                                while ((line_input_regions = br_input_regions.readLine()) != null) {
                                    String[] split = line_input_regions.split("\t");

                                    String chromosome = split[0];
                                    int start = Integer.parseInt(split[1]);
                                    int end = Integer.parseInt(split[2]);
                                    String name = split[3];
                                    int score = Integer.parseInt(split[4]);
                                    String strand = split[5];
                                    double signalValue = Double.parseDouble(split[6]);
                                    double pValue = Double.parseDouble(split[7]);
                                    double qValue = Double.parseDouble(split[8]);
                                    Footprint_Interval fp =
                                            new Footprint_Interval(chromosome, start, end, name, score, strand,
                                                    signalValue, pValue, qValue, line_input_regions);

                                    ArrayList<Footprint_Interval> fp_al = chr_regions.get(chromosome);
                                    fp_al.add(fp);
                                }
                                br_input_regions.close();

                                //sort region array lists per chromosome
                                for (String key_chr : chr_regions.keySet()) {
                                    Collections.sort(chr_regions.get(key_chr));
                                }

                                //determine chromosome order
                                ArrayList<String> chromosome_ordered = new ArrayList<>();
                                for (int i = 0; i < 100; i++) {
                                    String search_chr = "" + i;
                                    if (chr_regions.containsKey(search_chr)) {
                                        chromosome_ordered.add(search_chr);
                                    }
                                }
                                ArrayList<String> character_chromosomes = new ArrayList<>();
                                for (String key_chr : chr_regions.keySet()) {
                                    if (!chromosome_ordered.contains(key_chr)) {
                                        character_chromosomes.add(key_chr);
                                    }
                                }
                                Collections.sort(character_chromosomes);
                                if (character_chromosomes.get(0).matches(".*M.*")) {
                                    String mt_chr = character_chromosomes.get(0);
                                    character_chromosomes.remove(0);
                                    character_chromosomes.add(mt_chr);
                                }
                                chromosome_ordered.addAll(character_chromosomes);


                                //see if gaps between Footprint Intervals are < tepic_between_max_bps if so connect to one region for TEPIC

                                if (options_intern.tepic_tf_binding_site_search.equals("BETWEEN")) {
                                    for (String key_chr : chr_regions.keySet()) {
                                        ArrayList<Footprint_Interval> intervals = chr_regions.get(key_chr);

                                        ArrayList<Footprint_Interval> intervals_combined = new ArrayList<>();

                                        for (int i = 0; i < intervals.size(); i++) {
                                            int start_i = i;
                                            int end_i = i;

                                            int start_position = intervals.get(i).start;
                                            int end_position = intervals.get(i).end;

                                            String name = intervals.get(i).name;
                                            int score = intervals.get(i).score;
                                            String strand = intervals.get(i).strand;
                                            double signalValue = intervals.get(i).signalValue;
                                            double pValue = intervals.get(i).pValue;
                                            double qValue = intervals.get(i).qValue;


                                            boolean include_next_peak = true;

                                            while (include_next_peak) {
                                                if (end_i + 1 >= intervals.size()) {
                                                    break;
                                                }

                                                int start_position_next = intervals.get(end_i + 1).start;
                                                int end_position_next = intervals.get(end_i + 1).end;

                                                int distance = start_position_next - end_position;

                                                if (distance < options_intern.tepic_between_max_bps) {
                                                    include_next_peak = true;
                                                    end_position = end_position_next;
                                                    name += ";" + intervals.get(end_i + 1).name;
                                                    strand += ";" + intervals.get(end_i + 1).strand;


                                                    score += intervals.get(end_i + 1).score;
                                                    signalValue += intervals.get(end_i + 1).signalValue;
                                                    pValue += intervals.get(end_i + 1).pValue;
                                                    qValue += intervals.get(end_i + 1).qValue;


                                                    end_i += 1;
                                                } else {
                                                    include_next_peak = false;
                                                }
                                            }

                                            if (start_i == end_i) {
                                                intervals_combined.add(intervals.get(i));
                                            } else {
                                                i = end_i + 1;

                                                int quotient = end_i - start_i + 1;

                                                score /= quotient;
                                                signalValue /= quotient;
                                                pValue /= quotient;
                                                qValue /= quotient;


                                                Footprint_Interval combined =
                                                        new Footprint_Interval(key_chr, start_position, end_position,
                                                                name, score, strand, signalValue, pValue, qValue);
                                                combined.make_line();

                                                intervals_combined.add(combined);
                                            }
                                        }
                                        chr_regions.put(key_chr, intervals_combined);
                                    }
                                }

                                if (options_intern.tepic_tf_binding_site_search.equals("EXCL_BETWEEN")) {
                                    for (String key_chr : chr_regions.keySet()) {
                                        ArrayList<Footprint_Interval> intervals = chr_regions.get(key_chr);

                                        ArrayList<Footprint_Interval> intervals_combined = new ArrayList<>();

                                        for (int i = 0; i < intervals.size(); i++) {
                                            int start_i = i;
                                            int end_i = i;

                                            int start_position = intervals.get(i).start;
                                            int end_position = intervals.get(i).end;

                                            String name = intervals.get(i).name;
                                            int score = intervals.get(i).score;
                                            String strand = intervals.get(i).strand;
                                            double signalValue = intervals.get(i).signalValue;
                                            double pValue = intervals.get(i).pValue;
                                            double qValue = intervals.get(i).qValue;

                                            if (end_i + 1 >= intervals.size()) {
                                                intervals_combined.add(intervals.get(i));
                                                break;
                                            } else {
                                                int start_position_next = intervals.get(end_i + 1).start;
                                                int end_position_next = intervals.get(end_i + 1).end;

                                                int distance = start_position_next - end_position;

                                                if (distance < options_intern.tepic_between_max_bps) {
                                                    start_position = end_position;
                                                    end_position = start_position_next;

                                                    name += ";" + intervals.get(end_i + 1).name;
                                                    strand += ";" + intervals.get(end_i + 1).strand;


                                                    score += intervals.get(end_i + 1).score;
                                                    signalValue += intervals.get(end_i + 1).signalValue;
                                                    pValue += intervals.get(end_i + 1).pValue;
                                                    qValue += intervals.get(end_i + 1).qValue;


                                                    end_i += 1;

                                                }
                                            }

                                            if (start_i == end_i) {
                                                intervals_combined.add(intervals.get(i));
                                            } else {
                                                i = end_i + 1;

                                                int quotient = end_i - start_i + 1;

                                                score /= quotient;
                                                signalValue /= quotient;
                                                pValue /= quotient;
                                                qValue /= quotient;


                                                Footprint_Interval combined =
                                                        new Footprint_Interval(key_chr, start_position, end_position,
                                                                name, score, strand, signalValue, pValue, qValue);
                                                combined.make_line();

                                                intervals_combined.add(combined);
                                            }
                                        }
                                        chr_regions.put(key_chr, intervals_combined);
                                    }
                                }


                                //write output file

                                BufferedWriter bw_out = new BufferedWriter(new FileWriter(f_output_peak_file));

                                for (String key_chr : chromosome_ordered) {
                                    ArrayList<Footprint_Interval> footprints = chr_regions.get(key_chr);

                                    for (Footprint_Interval fpi : footprints) {
                                        bw_out.write(fpi.line);
                                        bw_out.newLine();
                                    }

                                }

                                bw_out.close();

                            }
                        }
                    }
                }
            }
        }

        logger.logLine("[FOOTPRINTS] finished footprints");
    }

    /**
     * preprocess mix histones, search for same peaks and use either the union or the intersection of all
     */
    public void mix_option() throws IOException {

        options_intern.tepic_input_prev = options_intern.tepic_input_directory;
        File f_annotation_check = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option +
                        File.separator + options_intern.folder_name_mix_option_preprocessing_check_chr);
        options_intern.tepic_input_directory = f_annotation_check.getAbsolutePath();

        File file_root_input = new File(options_intern.tepic_input_directory);
        File root_mix_working_dir = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option);
        root_mix_working_dir.mkdir();

        File f_sample_mix_preprocess = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                options_intern.folder_name_mix_option_sample_mix_preprocessing);
        f_sample_mix_preprocess.mkdir();

        File f_sample_mix_output = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                options_intern.folder_name_mix_option_sample_mix);
        f_sample_mix_output.mkdir();

        logger.logLine("[MIX] Used data: " + options_intern.tepic_input_directory);

        if (options_intern.mix_level.equals("SAMPLE_LEVEL")) {
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = f_sample_mix_output.getAbsolutePath();
        }

        logger.logLine("[MIX] Preprocess input data for sample mix - split chromosomes.");

        for (File fileDir : file_root_input.listFiles()) {
            if (fileDir.isDirectory()) {
                String timepoint = fileDir.getName();
                File file_output_tp = new File(f_sample_mix_preprocess + File.separator + timepoint);
                file_output_tp.mkdir();

                for (File fileDirHM : fileDir.listFiles()) {
                    if (fileDirHM.isDirectory()) {
                        File file_output_tp_hm =
                                new File(file_output_tp.getAbsolutePath() + File.separator + fileDirHM.getName());
                        file_output_tp_hm.mkdir();

                        for (File fileDirHM_sample : fileDirHM.listFiles()) {
                            if (fileDirHM_sample.isFile()) {
                                BufferedWriter bw = new BufferedWriter(new FileWriter(
                                        new File(file_output_tp_hm.getAbsolutePath() + File.separator + "test.txt")));
                                BufferedReader br = new BufferedReader(new FileReader(fileDirHM_sample));
                                String line = "";
                                String currentChr = "";
                                while ((line = br.readLine()) != null) {
                                    String[] split = line.split("\t");
                                    if (!split[0].equals(currentChr)) {
                                        bw.close();
                                        currentChr = split[0];
                                        File f_output_chr = new File(
                                                file_output_tp_hm.getAbsolutePath() + File.separator + currentChr);
                                        f_output_chr.mkdir();

                                        bw = new BufferedWriter(new FileWriter(new File(
                                                f_output_chr.getAbsolutePath() + File.separator +
                                                        fileDirHM_sample.getName())));
                                    }
                                    bw.write(line);
                                    bw.newLine();
                                }
                                bw.close();
                            }
                        }
                    }
                }
            }
        }

        logger.logLine("[MIX] Create " + options_intern.mix_option + " of samples");

        for (File fileDir : f_sample_mix_preprocess.listFiles()) {
            if (fileDir.isDirectory()) {
                String timepoint = fileDir.getName();
                File f_output_union_samples_tp =
                        new File(f_sample_mix_output.getAbsolutePath() + File.separator + timepoint);
                f_output_union_samples_tp.mkdir();

                for (File fileDirHM : fileDir.listFiles()) {
                    if (fileDirHM.isDirectory()) {
                        String hm = fileDirHM.getName();
                        File f_output_union_samples_tp_hm =
                                new File(f_output_union_samples_tp.getAbsolutePath() + File.separator + hm);
                        f_output_union_samples_tp_hm.mkdir();

                        String file_ending = "";

                        HashMap<String, ArrayList<MIX_Interval>> all_chromosomes = new HashMap();
                        ArrayList<Integer> chr_alpha = new ArrayList<>();
                        ArrayList<String> chr_str = new ArrayList<>();

                        for (File fileDirHM_Chr : fileDirHM.listFiles()) {
                            if (fileDirHM_Chr.isDirectory()) {
                                String chr = fileDirHM_Chr.getName();

                                ArrayList<MIX_Interval> all_intervals = new ArrayList<>();

                                int sample_number = 0;

                                for (File fileDirHM_Chr_sample : fileDirHM_Chr.listFiles()) {
                                    if (fileDirHM_Chr_sample.isFile()) {
                                        String[] f_ending = fileDirHM_Chr_sample.getName().split("\\.");
                                        file_ending = f_ending[f_ending.length - 1];
                                        //build array of a union of all peaks
                                        BufferedReader br = new BufferedReader(new FileReader(fileDirHM_Chr_sample));
                                        String line = "";
                                        while ((line = br.readLine()) != null) {
                                            String[] split = line.split("\t");

                                            MIX_Interval_Object mio =
                                                    new MIX_Interval_Object(split[0], Integer.parseInt(split[1]),
                                                            Integer.parseInt(split[2]), split[3],
                                                            Integer.parseInt(split[4]), split[5],
                                                            Double.parseDouble(split[6]), Double.parseDouble(split[7]),
                                                            Double.parseDouble(split[8]));
                                            MIX_Interval mi = new MIX_Interval(Integer.parseInt(split[1]),
                                                    Integer.parseInt(split[2]));
                                            mi.merged_intervals.add(mio);

                                            all_intervals.add(mi);

                                        }
                                        br.close();
                                        sample_number++;

                                        /*
                                        if(options_intern.mix_option.equals("INTERSECTION"))
                                        {
                                            ArrayList<MIX_Interval> intersec_al = new ArrayList<>();

                                            BufferedReader br_intersec = new BufferedReader(new FileReader(fileDirHM_Chr_sample));
                                            String line_intersec = "";
                                            while((line_intersec=br_intersec.readLine())!=null)
                                            {
                                                String[] split = line_intersec.split("\t");

                                                MIX_Interval_Object mio = new MIX_Interval_Object(split[0],Integer.parseInt(split[1]),Integer.parseInt(split[2]),split[3],Integer.parseInt(split[4]),split[5],Double.parseDouble(split[6]),Double.parseDouble(split[7]),Double.parseDouble(split[8]));
                                                MIX_Interval mi = new MIX_Interval(Integer.parseInt(split[1]),Integer.parseInt(split[2]));
                                                mi.merged_intervals.add(mio);

                                                intersec_al.add(mi);

                                            }
                                            br.close();


                                            intersec_vals.put(fileDirHM_Chr_sample.getName(),intersec_al);
                                        }*/

                                    }
                                }

                                //intersections sorting!
                                /*if(options_intern.mix_option.equals("INTERSECTION"))
                                {
                                    for(String k: intersec_vals.keySet())
                                    {
                                        Collections.sort(intersec_vals.get(k));
                                    }
                                }*/

                                //unions sorting
                                Collections.sort(all_intervals);

                                Stack<MIX_Interval> stack_union = mergeIntervals(all_intervals);

                                ArrayList<MIX_Interval> chr_unions = new ArrayList<>();

                                int min_occurence = 0;
                                if (options_intern.mix_occurence_intersection == -1) {
                                    min_occurence = sample_number;
                                } else {
                                    if (sample_number < options_intern.mix_occurence_intersection) {
                                        min_occurence = sample_number;

                                    } else {
                                        min_occurence = options_intern.mix_occurence_intersection;
                                    }
                                }

                                while (!stack_union.isEmpty()) {
                                    MIX_Interval t = stack_union.pop();
                                    t.calculate_mean("SAMPLE_LEVEL");

                                    if (options_intern.mix_option.equals("INTERSECTION")) {
                                        if (t.merged_intervals.size() >= min_occurence) {
                                            chr_unions.add(t);
                                        }
                                    } else {
                                        chr_unions.add(t);
                                    }
                                }

                                Collections.sort(chr_unions);

                                all_chromosomes.put(chr, chr_unions);

                                try {
                                    chr_alpha.add(Integer.parseInt(chr));
                                } catch (Exception e) {
                                    chr_str.add(chr);
                                }
                            }
                        }

                        Collections.sort(chr_alpha);
                        Collections.sort(chr_str);


                        //print Sample_unions => can be used if not HM_Level is used
                        BufferedWriter bw = new BufferedWriter(new FileWriter(
                                f_output_union_samples_tp_hm.getAbsolutePath() + File.separator + timepoint + "_" + hm +
                                        "." + file_ending));

                        int peak_counter = 1;

                        for (int chr : chr_alpha) {
                            ArrayList<MIX_Interval> x = all_chromosomes.get("" + chr);
                            for (int i = 0; i < x.size(); i++) {
                                bw.write(x.get(i).meanToString(peak_counter));
                                bw.newLine();
                                peak_counter++;
                            }
                        }

                        for (String chr : chr_str) {
                            ArrayList<MIX_Interval> x = all_chromosomes.get(chr);
                            for (int i = 0; i < x.size(); i++) {
                                bw.write(x.get(i).meanToString(peak_counter));
                                bw.newLine();
                                peak_counter++;
                            }
                        }
                        bw.close();


                    }
                }
            }
        }

        if (options_intern.mix_level.equals("HM_LEVEL")) {

            logger.logLine("[MIX] Preprocess sample unions for HM mix");

            File f_output_preprocessing_hm = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                    options_intern.folder_name_mix_option_preprocess_hm_mix);
            f_output_preprocessing_hm.mkdir();

            File f_output_hm = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                    options_intern.folder_name_mix_option_hm_mix);
            f_output_hm.mkdir();

            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = f_output_hm.getAbsolutePath();

            logger.logLine("[MIX] Identify possible Timepoints with same Histone Modifications");

            HashMap<String, ArrayList<String>> timepoints_histone_modifications = new HashMap<>();
            HashMap<String, ArrayList<String>> deleted_tps = new HashMap<>();
            HashSet<String> available_hms = new HashSet<>();

            //identify possible timepoints
            for (File fileDir : f_sample_mix_output.listFiles()) {
                if (fileDir.isDirectory()) {
                    String timepoint = fileDir.getName();
                    ArrayList hm = new ArrayList();

                    for (File fileDirHM : fileDir.listFiles()) {
                        if (fileDirHM.isDirectory()) {
                            hm.add(fileDirHM.getName());
                            available_hms.add(fileDirHM.getName());
                        }
                    }
                    timepoints_histone_modifications.put(timepoint, hm);
                }
            }


            for (String tp : timepoints_histone_modifications.keySet()) {
                if (timepoints_histone_modifications.get(tp).size() < available_hms.size()) {
                    deleted_tps.put(tp, timepoints_histone_modifications.get(tp));
                    timepoints_histone_modifications.remove(tp);
                    continue;
                }

                boolean all_in = true;
                for (String hm : available_hms) {
                    if (!timepoints_histone_modifications.get(tp).contains(hm)) {
                        all_in = false;
                    }
                }

                if (!all_in) {
                    deleted_tps.put(tp, timepoints_histone_modifications.get(tp));
                    timepoints_histone_modifications.remove(tp);
                }
            }

            StringBuilder sb_found_mixing_tps = new StringBuilder();
            sb_found_mixing_tps.append("[MIX] Can perform complete mix for HMs (");
            for (String hm : available_hms) {
                sb_found_mixing_tps.append(hm);
                sb_found_mixing_tps.append(" ");
            }
            sb_found_mixing_tps.append(") in timepoints (");
            for (String tp : timepoints_histone_modifications.keySet()) {
                sb_found_mixing_tps.append(tp);
                sb_found_mixing_tps.append(" ");
            }
            sb_found_mixing_tps.append("). Can perform part-mix or no-mix for timepoints (");
            for (String tp : deleted_tps.keySet()) {
                sb_found_mixing_tps.append(tp);
                sb_found_mixing_tps.append(" ");
            }
            sb_found_mixing_tps.append(").");
            logger.logLine(sb_found_mixing_tps.toString());


            for (File fileDir : f_sample_mix_output.listFiles()) {
                if (fileDir.isDirectory()) {
                    String timepoint = fileDir.getName();
                    File f_output_hm_prepro =
                            new File(f_output_preprocessing_hm.getAbsolutePath() + File.separator + timepoint);
                    f_output_hm_prepro.mkdir();

                    for (File fileDirHM : fileDir.listFiles()) {
                        if (fileDirHM.isDirectory()) {
                            String hm = fileDirHM.getName();
                            File f_output_hm_prepro_hm =
                                    new File(f_output_hm_prepro.getAbsolutePath() + File.separator + "MIX");
                            f_output_hm_prepro_hm.mkdir();

                            for (File fileDirHM_samples : fileDirHM.listFiles()) {
                                if (fileDirHM_samples.isFile()) {
                                    BufferedWriter bw = new BufferedWriter(new FileWriter(new File(
                                            f_output_hm_prepro_hm.getAbsolutePath() + File.separator + "test.txt")));
                                    BufferedReader br = new BufferedReader(new FileReader(fileDirHM_samples));
                                    String line = "";
                                    String currentChr = "";
                                    while ((line = br.readLine()) != null) {
                                        String[] split = line.split("\t");
                                        if (!split[0].equals(currentChr)) {
                                            bw.close();
                                            currentChr = split[0];
                                            File f_output_chr = new File(
                                                    f_output_hm_prepro_hm.getAbsolutePath() + File.separator +
                                                            currentChr);
                                            f_output_chr.mkdir();

                                            bw = new BufferedWriter(new FileWriter(new File(
                                                    f_output_chr.getAbsolutePath() + File.separator +
                                                            fileDirHM_samples.getName())));
                                        }
                                        bw.write(line);
                                        bw.newLine();
                                    }
                                    bw.close();
                                }
                            }
                        }
                    }
                }
            }

            logger.logLine("[MIX] Create " + options_intern.mix_option + " of HMs");

            for (File fileDir : f_output_preprocessing_hm.listFiles()) {
                if (fileDir.isDirectory()) {
                    String timepoint = fileDir.getName();
                    File f_output_union_samples_tp =
                            new File(f_output_hm.getAbsolutePath() + File.separator + timepoint);
                    f_output_union_samples_tp.mkdir();

                    for (File fileDirHM : fileDir.listFiles()) {
                        if (fileDirHM.isDirectory()) {
                            String hm = fileDirHM.getName();
                            File f_output_union_samples_tp_hm =
                                    new File(f_output_union_samples_tp.getAbsolutePath() + File.separator + hm);
                            f_output_union_samples_tp_hm.mkdir();

                            String file_ending = "";

                            HashMap<String, ArrayList<MIX_Interval>> all_chromosomes = new HashMap();
                            ArrayList<Integer> chr_alpha = new ArrayList<>();
                            ArrayList<String> chr_str = new ArrayList<>();

                            for (File fileDirHM_Chr : fileDirHM.listFiles()) {
                                if (fileDirHM_Chr.isDirectory()) {
                                    String chr = fileDirHM_Chr.getName();

                                    ArrayList<MIX_Interval> all_intervals = new ArrayList<>();

                                    int sample_number = 0;

                                    for (File fileDirHM_Chr_sample : fileDirHM_Chr.listFiles()) {
                                        if (fileDirHM_Chr_sample.isFile()) {
                                            String[] f_ending = fileDirHM_Chr_sample.getName().split("\\.");
                                            file_ending = f_ending[f_ending.length - 1];
                                            //build array of a union of all peaks
                                            BufferedReader br =
                                                    new BufferedReader(new FileReader(fileDirHM_Chr_sample));
                                            String line = "";
                                            while ((line = br.readLine()) != null) {
                                                String[] split = line.split("\t");

                                                MIX_Interval_Object mio =
                                                        new MIX_Interval_Object(split[0], Integer.parseInt(split[1]),
                                                                Integer.parseInt(split[2]), split[3],
                                                                Integer.parseInt(split[4]), split[5],
                                                                Double.parseDouble(split[6]),
                                                                Double.parseDouble(split[7]),
                                                                Double.parseDouble(split[8]));
                                                MIX_Interval mi = new MIX_Interval(Integer.parseInt(split[1]),
                                                        Integer.parseInt(split[2]));
                                                mi.merged_intervals.add(mio);

                                                all_intervals.add(mi);

                                            }
                                            br.close();
                                            sample_number++;

                                        /*
                                        if(options_intern.mix_option.equals("INTERSECTION"))
                                        {
                                            ArrayList<MIX_Interval> intersec_al = new ArrayList<>();

                                            BufferedReader br_intersec = new BufferedReader(new FileReader(fileDirHM_Chr_sample));
                                            String line_intersec = "";
                                            while((line_intersec=br_intersec.readLine())!=null)
                                            {
                                                String[] split = line_intersec.split("\t");

                                                MIX_Interval_Object mio = new MIX_Interval_Object(split[0],Integer.parseInt(split[1]),Integer.parseInt(split[2]),split[3],Integer.parseInt(split[4]),split[5],Double.parseDouble(split[6]),Double.parseDouble(split[7]),Double.parseDouble(split[8]));
                                                MIX_Interval mi = new MIX_Interval(Integer.parseInt(split[1]),Integer.parseInt(split[2]));
                                                mi.merged_intervals.add(mio);

                                                intersec_al.add(mi);

                                            }
                                            br.close();


                                            intersec_vals.put(fileDirHM_Chr_sample.getName(),intersec_al);
                                        }*/

                                        }
                                    }

                                    //intersections sorting!
                                /*if(options_intern.mix_option.equals("INTERSECTION"))
                                {
                                    for(String k: intersec_vals.keySet())
                                    {
                                        Collections.sort(intersec_vals.get(k));
                                    }
                                }*/

                                    //unions sorting
                                    Collections.sort(all_intervals);

                                    Stack<MIX_Interval> stack_union = mergeIntervals(all_intervals);

                                    ArrayList<MIX_Interval> chr_unions = new ArrayList<>();

                                    int min_occurence = 0;
                                    if (options_intern.mix_occurence_intersection == -1) {
                                        min_occurence = sample_number;
                                    } else {
                                        if (sample_number < options_intern.mix_occurence_intersection) {
                                            min_occurence = sample_number;

                                        } else {
                                            min_occurence = options_intern.mix_occurence_intersection;
                                        }
                                    }

                                    while (!stack_union.isEmpty()) {
                                        MIX_Interval t = stack_union.pop();
                                        t.calculate_mean("HM_LEVEL");

                                        if (options_intern.mix_option.equals("INTERSECTION")) {
                                            if (t.merged_intervals.size() >= min_occurence) {
                                                chr_unions.add(t);
                                            }
                                        } else {
                                            chr_unions.add(t);
                                        }
                                    }

                                    Collections.sort(chr_unions);

                                    all_chromosomes.put(chr, chr_unions);

                                    try {
                                        chr_alpha.add(Integer.parseInt(chr));
                                    } catch (Exception e) {
                                        chr_str.add(chr);
                                    }
                                }
                            }

                            Collections.sort(chr_alpha);
                            Collections.sort(chr_str);


                            //print Sample_unions => can be used if not HM_Level is used
                            BufferedWriter bw = new BufferedWriter(new FileWriter(
                                    f_output_union_samples_tp_hm.getAbsolutePath() + File.separator + timepoint + "_" +
                                            hm + "." + file_ending));

                            int peak_counter = 1;

                            for (int chr : chr_alpha) {
                                ArrayList<MIX_Interval> x = all_chromosomes.get("" + chr);
                                for (int i = 0; i < x.size(); i++) {
                                    bw.write(x.get(i).meanToString(peak_counter));
                                    bw.newLine();
                                    peak_counter++;
                                }
                            }

                            for (String chr : chr_str) {
                                ArrayList<MIX_Interval> x = all_chromosomes.get(chr);
                                for (int i = 0; i < x.size(); i++) {
                                    bw.write(x.get(i).meanToString(peak_counter));
                                    bw.newLine();
                                    peak_counter++;
                                }
                            }
                            bw.close();
                        }
                    }
                }
            }
        }
    }

    /**
     * check if chromosome are annotated without a chr prefix
     *
     * @throws IOException
     */
    public void check_chromosomes() throws IOException {
        logger.logLine("[CHECK CHROMOSOMES] Check chromosome for naming convention and alter if necessary!");

        File file_root_input = new File(options_intern.tepic_input_directory);

        File root_mix_working_dir = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option);
        root_mix_working_dir.mkdir();

        File chr_annotation_output = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                options_intern.folder_name_mix_option_preprocessing_check_chr);
        chr_annotation_output.mkdir();

        options_intern.tepic_input_prev = options_intern.tepic_input_directory;
        options_intern.tepic_input_directory = chr_annotation_output.getAbsolutePath();

        for (File fileDir : file_root_input.listFiles()) {
            if (fileDir.isDirectory()) {
                String timepoint = fileDir.getName();
                File file_output_tp = new File(chr_annotation_output.getAbsolutePath() + File.separator + timepoint);
                file_output_tp.mkdir();

                for (File fileDirHM : fileDir.listFiles()) {
                    if (fileDirHM.isDirectory()) {
                        File file_output_tp_hm =
                                new File(file_output_tp.getAbsolutePath() + File.separator + fileDirHM.getName());
                        file_output_tp_hm.mkdir();

                        for (File fileDirHM_sample : fileDirHM.listFiles()) {
                            if (fileDirHM_sample.isFile()) {
                                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(
                                        file_output_tp_hm.getAbsolutePath() + File.separator +
                                                fileDirHM_sample.getName())));
                                BufferedReader br = new BufferedReader(new FileReader(fileDirHM_sample));
                                String line = "";
                                while ((line = br.readLine()) != null) {
                                    String[] split = line.split("\t");

                                    if (split[0].startsWith("chr")) {
                                        String[] split_chr = split[0].split("chr");
                                        String chr_name = split_chr[1];

                                        StringBuilder sb_line = new StringBuilder();
                                        sb_line.append(chr_name);

                                        for (int i = 1; i < split.length; i++) {
                                            sb_line.append("\t");
                                            sb_line.append(split[i]);
                                        }

                                        bw.write(sb_line.toString());
                                        bw.newLine();

                                    } else {
                                        bw.write(line);
                                        bw.newLine();
                                    }
                                }
                                bw.close();
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * read config file
     *
     * @param check_options should options be checked for validity? Should be true if pipeline is run, should be false if analyse programms are run
     */
    public void read_config_file(boolean check_options) throws IOException {

        logger.logLine("Start reading config file at " + options_intern.config_data_path);

        File f_check_config = new File(options_intern.config_data_path);
        if (!f_check_config.exists()) {
            logger.logLine("[CHECK] config file does not exist. Please check path!");
            logger.logLine("[CHECK] " + options_intern.config_data_path);
            System.exit(1);
        }

        File f_check_working_dir = new File(options_intern.com2pose_working_directory);
        if (!f_check_working_dir.exists()) {
            logger.logLine("[CHECK] working directory does not exist. Please check path!");
            logger.logLine("[CHECK] " + options_intern.com2pose_working_directory);
            System.exit(1);
        }

        BufferedReader br = new BufferedReader(new FileReader(new File(options_intern.config_data_path)));
        String line = "";
        while ((line = br.readLine()) != null) {
            if (line.startsWith("#")) {
                continue;
            }

            String[] split = line.split("=");

            switch (split[0]) {
                case "mix_level":
                    options_intern.mix_level = split[1].substring(1, split[1].length() - 1);
                    break;
                case "mix_option":
                    options_intern.mix_option = split[1].substring(1, split[1].length() - 1);
                    break;
                case "mix_occurence_intersection":
                    options_intern.mix_occurence_intersection = Integer.parseInt(split[1]);
                    break;
                case "mix_mutually_exclusive":
                    options_intern.mix_mutually_exclusive = Boolean.parseBoolean(split[1]);
                    break;
                case "mix_mutually_exclusive_diff_peak_signals":
                    options_intern.mix_mutually_exclusive_diff_peak_signals = Boolean.parseBoolean(split[1]);
                    break;
                case "black_list_dir":
                    options_intern.black_list_dir = split[1].substring(1, split[1].length() - 1);
                    break;
                case "black_list_signals":
                    options_intern.black_list_signals = new HashSet<>(
                            Arrays.asList(split[1].substring(1, split[1].length() - 1).toUpperCase().split(";")));
                    break;
                case "deseq2_input_directory":
                    options_intern.deseq2_input_directory = split[1].substring(1, split[1].length() - 1);
                    break;
                case "deseq2_input_gene_id":
                    options_intern.deseq2_input_gene_id = split[1].substring(1, split[1].length() - 1);
                    break;
                case "deseq2_biomart_dataset_species":
                    options_intern.deseq2_biomart_dataset_species = split[1].substring(1, split[1].length() - 1);
                    break;
                case "deseq2_biomart_dataset_symbol_column":
                    options_intern.deseq2_biomart_dataset_symbol_column = split[1].substring(1, split[1].length() - 1);
                    break;
                case "deseq2_count_threshold":
                    options_intern.deseq2_count_threshold = Integer.parseInt(split[1]);
                    break;
                case "deseq2_tpm_filter":
                    options_intern.deseq2_tpm_filter = Double.parseDouble(split[1]);
                    break;
                case "tepic_input_directory":
                    options_intern.tepic_input_directory = split[1].substring(1, split[1].length() - 1);
                    options_intern.tepic_input_original = split[1].substring(1, split[1].length() - 1);
                    break;
                case "tepic_input_ref_genome":
                    options_intern.tepic_input_ref_genome = split[1].substring(1, split[1].length() - 1);
                    break;
                case "tepic_path_pwms":
                    options_intern.tepic_path_pwms = split[1].substring(1, split[1].length() - 1);
                    break;
                case "tepic_cores":
                    options_intern.tepic_cores = Integer.parseInt(split[1]);
                    break;
                case "tepic_bed_chr_sign":
                    options_intern.tepic_bed_chr_sign = split[1].substring(1, split[1].length() - 1);
                    break;
                case "tepic_column_bedfile":
                    options_intern.tepic_column_bedfile = Integer.parseInt(split[1]);
                    break;
                case "tepic_gene_annot":
                    options_intern.tepic_gene_annot = split[1].substring(1, split[1].length() - 1);
                    break;
                case "tepic_window_size":
                    options_intern.tepic_window_size = Integer.parseInt(split[1]);
                    break;
                case "tepic_onlyDNasePeaks":
                    options_intern.tepic_onlyDNasePeaks = split[1].substring(1, split[1].length() - 1);
                    break;
                case "tepic_exponential_decay":
                    options_intern.tepic_exponential_decay = Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_not_norm_peak_length":
                    options_intern.tepic_not_norm_peak_length = Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_not_generated":
                    options_intern.tepic_not_generated = Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_original_decay":
                    options_intern.tepic_original_decay = Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_psems_length":
                    options_intern.tepic_psems_length = split[1].substring(1, split[1].length() - 1);
                    break;
                case "tepic_entire_gene_body":
                    options_intern.tepic_entire_gene_body = Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_zipped":
                    options_intern.tepic_zipped = Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_background_seq":
                    options_intern.tepic_background_seq = split[1].substring(1, split[1].length() - 1);
                    break;
                case "tepic_2bit":
                    options_intern.tepic_2bit = split[1].substring(1, split[1].length() - 1);
                    break;
                case "tepic_pvalue":
                    options_intern.tepic_pvalue = Double.parseDouble(split[1]);
                    break;
                case "tepic_minutes_per_chr":
                    options_intern.tepic_minutes_per_chr = Integer.parseInt(split[1]);
                    break;
                case "tepic_chr_prefix":
                    options_intern.tepic_chr_prefix = Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_transcript_based":
                    options_intern.tepic_transcript_based = Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_loop_list":
                    options_intern.tepic_loop_list = split[1].substring(1, split[1].length() - 1);
                    break;
                case "tepic_loop_windows":
                    options_intern.tepic_loop_windows = Integer.parseInt(split[1]);
                    break;
                case "tepic_only_peak_features":
                    options_intern.tepic_only_peak_features = Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_tpm_cutoff":
                    options_intern.tepic_tpm_cutoff = Double.parseDouble(split[1]);
                    break;
                case "tepic_ensg_symbol":
                    options_intern.tepic_ensg_symbol = split[1].substring(1, split[1].length() - 1);
                    break;
                case "tepic_tgene_target_genes":
                    options_intern.tepic_tgene_target_genes = Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_randomize_tf_gene_matrix":
                    options_intern.tepic_randomize_tf_gene_matrix = Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_tf_binding_site_search":
                    options_intern.tepic_tf_binding_site_search = split[1].substring(1, split[1].length() - 1);
                    break;
                case "tepic_between_max_bps":
                    options_intern.tepic_between_max_bps = Integer.parseInt(split[1]);
                    break;
                case "tgen_consensus":
                    options_intern.tgen_consensus = Double.parseDouble(split[1]);
                    break;
                case "tgen_consensus_calc":
                    options_intern.tgen_consensus_calc = split[1].substring(1, split[1].length() - 1);
                    break;
                case "tgen_no_closest_locus":
                    options_intern.tgen_no_closest_locus = Boolean.parseBoolean(split[1]);
                    break;
                case "tgen_self_regulatory":
                    options_intern.tgen_self_regulatory = Boolean.parseBoolean(split[1]);
                    break;
                case "tgen_no_closest_tss":
                    options_intern.tgen_no_closest_tss = Boolean.parseBoolean(split[1]);
                    break;
                case "tgen_max_link_distances":
                    options_intern.tgen_max_link_distances = Integer.parseInt(split[1]);
                    break;
                case "tgen_pvalue":
                    options_intern.tgen_pvalue = Double.parseDouble(split[1]);
                    break;
                case "tgen_mt_writing":
                    options_intern.tgen_mt_writing = split[1].substring(1, split[1].length() - 1);
                    break;
                case "dynamite_preprocessing_integrate_data_geneIDs":
                    options_intern.dynamite_preprocessing_integrate_data_geneIDs = Integer.parseInt(split[1]);
                    break;
                case "dynamite_preprocessing_integrate_data_log2fc":
                    options_intern.dynamite_preprocessing_integrate_data_log2fc = Integer.parseInt(split[1]);
                    break;
                case "dynamite_preprocessing_integrate_data_consider_geneFile":
                    options_intern.dynamite_preprocessing_integrate_data_consider_geneFile =
                            split[1].substring(1, split[1].length() - 1);
                    break;
                case "dynamite_out_var":
                    options_intern.dynamite_out_var = split[1].substring(1, split[1].length() - 1);
                    break;
                case "dynamite_cores":
                    options_intern.dynamite_cores = Integer.parseInt(split[1]);
                    break;
                case "dynamite_alpha":
                    options_intern.dynamite_alpha = Double.parseDouble(split[1]);
                    break;
                case "dynamite_testsize":
                    options_intern.dynamite_testsize = Double.parseDouble(split[1]);
                    break;
                case "dynamite_Ofolds":
                    options_intern.dynamite_Ofolds = Integer.parseInt(split[1]);
                    break;
                case "dynamite_Ifolds":
                    options_intern.dynamite_Ifolds = Integer.parseInt(split[1]);
                    break;
                case "dynamite_balanced":
                    options_intern.dynamite_balanced = Boolean.parseBoolean(split[1]);
                    break;
                case "dynamite_performance":
                    options_intern.dynamite_performance = Boolean.parseBoolean(split[1]);
                    break;
                case "dynamite_randomise":
                    options_intern.dynamite_randomise = Boolean.parseBoolean(split[1]);
                    break;
                case "plot_th_coefficient":
                    String[] split_coefficient_ths = split[1].split(";");
                    options_intern.plot_th_coefficient.clear();
                    for (String s : split_coefficient_ths) {
                        options_intern.plot_th_coefficient.add(Double.parseDouble(s));
                    }
                    break;
                case "plot_cutoff_tps":
                    options_intern.plot_cutoff_tps = Integer.parseInt(split[1]);
                    break;
                case "plot_cutoff_hms":
                    options_intern.plot_cutoff_hms = Integer.parseInt(split[1]);
                    break;
                case "plot_cutoff_gcs":
                    options_intern.plot_cutoff_gcs = Integer.parseInt(split[1]);
                    break;
                case "plot_cutoff_tpms":
                    options_intern.plot_cutoff_tpms = Double.parseDouble(split[1]);
                    break;
                case "plot_top_k_genes":
                    options_intern.plot_top_k_genes = Integer.parseInt(split[1]);
                    break;
                case "plot_distribution_analysis_score_type":
                    options_intern.plot_distribution_analysis_score_type = split[1].substring(1, split[1].length() - 1);
                    break;
                case "website_interesting_tfs":
                    String[] split_interesting_tfs = split[1].substring(1, split[1].length() - 1).split(";");
                    options_intern.website_interesting_tfs.addAll(Arrays.asList(split_interesting_tfs));
                    break;
                case "html_report_interesting_tfs":
                    String[] split_interesting_tfs_2 = split[1].substring(1, split[1].length() - 1).split(";");
                    options_intern.website_interesting_tfs.addAll(Arrays.asList(split_interesting_tfs_2));
                    break;
                case "plot_mann_whitneyU_pvalue_cutoff":
                    options_intern.plot_mann_whitneyU_pvalue_cutoff = Double.parseDouble(split[1]);
                    break;
                case "chip_atlas_genome_version":
                    options_intern.chip_atlas_genome_version = split[1].substring(1, split[1].length() - 1);
                    break;
                case "chip_atlas_tissue_type":
                    options_intern.chip_atlas_tissue_type = split[1].substring(1, split[1].length() - 1);
                    break;
                case "igv_path_to_igv":
                    options_intern.igv_path_to_igv = split[1].substring(1, split[1].length() - 1);
                    break;
                case "igv_path_to_tf_chip_seq":
                    options_intern.igv_path_to_tf_chip_seq = split[1].substring(1, split[1].length() - 1);
                    break;
                case "igv_include_prediction_data":
                    String predicition_data_list = split[1].substring(1, split[1].length() - 1);
                    String[] split_prediction_data_list = predicition_data_list.split(";");
                    options_intern.igv_include_prediction_data.addAll(Arrays.asList(split_prediction_data_list));
                    break;
                case "igv_important_locus_all_prio_tf":
                    String loci_of_interest = split[1].substring(1, split[1].length() - 1);
                    String[] split_loci_of_interest = loci_of_interest.split(";");
                    options_intern.igv_important_locus_all_prio_tf.addAll(Arrays.asList(split_loci_of_interest));
                    break;
                case "igv_port_number":
                    options_intern.igv_port_number = Integer.parseInt(split[1]);
                    break;
                case "igv_species_ref_genome":
                    options_intern.igv_species_ref_genome = split[1].substring(1, split[1].length() - 1);
                    break;
                default:
                    logger.logLine(
                            "Misformed cfg file - please use template of: /COM2POSE/config_templates/com2pose_template.cfg");
                    logger.logLine("Do not delete unused parameters in config data!");
                    System.exit(1);
            }
        }
        br.close();

        if (check_options) {
            boolean all_set = checkOptions();
            logger.logLine("Check config file parameters for validity");
            if (!all_set) {
                logger.logLine("Not all [REQ]uired options set. Please set them in config file");
                logger.logLine("Aborting COM2POSE");
                System.exit(1);
            } else {
                logger.logLine("Parameters in config file valid");
            }
        }


        logger.logLine("Reading config file finished - no errors detected");


    }

    private boolean checkOptions() throws IOException {

        boolean all_set = true;

        /**
         * mix histone options
         */

        if (!options_intern.mix_level.equals("")) {
            if (!options_intern.mix_level.equals("HM_LEVEL") && !options_intern.mix_level.equals("SAMPLE_LEVEL")) {
                logger.logLine("[MIX]: mix_level parameter must be either HM_LEVEL or SAMPLE_LEVEL");
                all_set = false;
            }
            if (!options_intern.mix_option.equals("UNION") && !options_intern.mix_option.equals("INTERSECTION")) {
                logger.logLine("[MIX]: mix_option parameter must be either UNION or INTERSECTION");
                all_set = false;
            }
        }

        /**
         * black list options
         */

        if (!options_intern.black_list_dir.equals("")) {
            File f = new File(options_intern.black_list_dir);
            if (!f.exists()) {
                logger.logLine("[BLACKLIST] blacklist path does not exist!");
                all_set = false;
            }
            if (f.isDirectory()) {
                logger.logLine("[BLACKLIST] blacklist path is a directory but must be a file");
                all_set = false;
            }

            if (options_intern.black_list_signals.isEmpty()) {
                logger.logLine("[BLACKLIST] signals must not be empty");
                all_set = false;
            }
        }

        /**
         * DESEQ2 options
         */

        if (options_intern.deseq2_input_directory.equals("")) {
            logger.logLine("[DESEQ2] input directory is not given");
            all_set = false;
        } else {
            File f = new File(options_intern.deseq2_input_directory);
            if (!f.exists()) {
                logger.logLine("[DESEQ2] input path does not exist!");
                all_set = false;
            }
        }

        if (options_intern.deseq2_input_gene_id.equals("")) {
            logger.logLine("[DESEQ2] gene ID file from nfcore RNA-seq is not given");
            all_set = false;
        } else {
            File f = new File(options_intern.deseq2_input_gene_id);
            if (!f.exists()) {
                logger.logLine("[DESEQ2] gene ID file from nfcore RNA-seq path does not exist!");
                all_set = false;
            }
        }

        if (options_intern.deseq2_biomart_dataset_species.equals("") && options_intern.tepic_ensg_symbol.equals("")) {
            logger.logLine("[DESEQ2] deseq2_biomart_dataset_species must be filled if tepic_ensg_symbol is empty!");
            all_set = false;
        }

        if (!options_intern.deseq2_biomart_dataset_species.equals("")) {
            options_intern.tepic_ensg_symbol = options_intern.com2pose_working_directory + File.separator +
                    options_intern.folder_name_deseq2_preprocessing + File.separator +
                    options_intern.file_suffix_deseq2_mapping;
        }


        /**
         * TEPIC options
         */
        if (options_intern.tepic_input_ref_genome.equals("")) {
            logger.logLine("[TEPIC] Reference genome path is not given");
            all_set = false;
        } else {
            File f = new File(options_intern.tepic_input_ref_genome);
            if (!f.exists()) {
                logger.logLine("[TEPIC] Reference genome path does not exist!");
                all_set = false;
            }
        }

        if (options_intern.tepic_gene_annot.equals("")) {
            logger.logLine("[TEPIC] Gene annotation path is not given");
            all_set = false;
        } else {
            File f = new File(options_intern.tepic_gene_annot);
            if (!f.exists()) {
                logger.logLine("[TEPIC] Gene annotation path does not exist!");
                all_set = false;
            }
        }

        if (options_intern.tepic_input_directory.equals("")) {
            logger.logLine("[TEPIC] nfcore ChIP-seq data directory is not given");
            all_set = false;
        } else {
            File f = new File(options_intern.tepic_input_directory);
            if (!f.exists()) {
                logger.logLine("[TEPIC] nfcore ChIP-seq data path does not exist!");
                all_set = false;
            }
        }

        if (options_intern.tepic_path_pwms.equals("")) {
            logger.logLine("[TEPIC] position specific energy matrix directory is not given");
            all_set = false;
        } else {
            File f = new File(options_intern.tepic_path_pwms);
            if (!f.exists()) {
                logger.logLine("[TEPIC] position specific energy matrix path does not exist!");
                all_set = false;
            }
        }

        if (options_intern.tepic_tpm_cutoff > 0) {
            //check for gene annotation
            if (options_intern.tepic_gene_annot.equals("")) {
                logger.logLine("[TEPIC] TPM cutoff set, but no annotation file is given");
                all_set = false;
            } else {
                File f = new File(options_intern.tepic_gene_annot);
                if (!f.exists()) {
                    logger.logLine("[TEPIC] TPM cutoff set and gene annotation file path does not exist!");
                    all_set = false;
                }
            }
        }

        if (options_intern.tepic_tf_binding_site_search.equals("") ||
                !(options_intern.tepic_tf_binding_site_search.equals("INSIDE") ||
                        options_intern.tepic_tf_binding_site_search.equals("BETWEEN") ||
                        options_intern.tepic_tf_binding_site_search.equals("EXCL_BETWEEN"))) {
            //check tf_binding_site_search option
            logger.logLine("[TEPIC] tepic_tf_binding_site_search must be either INSIDE or BETWEEN or EXCL_BETWEEN");
            all_set = false;
        }
        if (options_intern.tepic_tf_binding_site_search.equals("BETWEEN") &&
                options_intern.tepic_between_max_bps <= 0) {
            logger.logLine(
                    "[TEPIC] tepic_tf_binding_site_search is set to 'BETWEEN', tepic_between_max_bps must be set > 0");
            all_set = false;
        }


        /**
         *check for map ensg symbol
         */
        if (options_intern.tepic_ensg_symbol.equals("") && options_intern.deseq2_biomart_dataset_species.equals("")) {
            logger.logLine("[TEPIC] No map of ENSG to Gene Symbol is given!");
        } else if (options_intern.deseq2_biomart_dataset_species.equals("")) {
            File f = new File(options_intern.tepic_ensg_symbol);
            if (!f.exists()) {
                logger.logLine("[TEPIC] Map of ENSG to Gene Symbol file path does not exist!");
                all_set = false;
            }
        }

        if (!options_intern.tepic_background_seq.equals("") && !options_intern.tepic_2bit.equals("")) {
            logger.logLine("[TEPIC] parameters: tepic_background_seq and tepic_2bit are mutually exclusive");
            all_set = false;
        }

        if (options_intern.tepic_column_bedfile != -1 && !options_intern.tepic_bed_chr_sign.equals("")) {
            logger.logLine("[TEPIC] parameters: tepic_column_bedfile and tepic_bed_chr_sign are mutually exclusive");
            all_set = false;
        }

        if (!options_intern.tepic_bed_chr_sign.equals("")) {
            File f = new File(options_intern.tepic_bed_chr_sign);
            if (!f.exists()) {
                logger.logLine("[TEPIC] tepic_bed_chr_sign file path does not exists!");
                all_set = false;
            }
        }

        if (!options_intern.tepic_psems_length.equals("")) {
            File f = new File(options_intern.tepic_psems_length);
            if (!f.exists()) {
                logger.logLine("[TEPIC] tepic_psems_length file path does not exists!");
                all_set = false;
            }
        }

        if (!options_intern.tepic_loop_list.equals("")) {
            File f = new File(options_intern.tepic_loop_list);
            if (!f.exists()) {
                logger.logLine("[TEPIC] tepic_loop_list file path does not exists!");
                all_set = false;
            }
        }

        /**
         * TGEN options
         */
        if (!options_intern.path_tgen.equals("")) {
            File file_tgen = new File(options_intern.path_tgen);
            if (!file_tgen.exists() || !file_tgen.isDirectory()) {
                logger.logLine("[TGENE] TGene file directory does not exist or is not a directory!");
                all_set = false;
            }

            File tgene_dir = new File(options_intern.path_tgen + File.separator + "bin");

            if (!tgene_dir.exists()) {
                logger.logLine("[TGENE] TGene binary directory cannot be found: " + tgene_dir.getAbsolutePath());
                all_set = false;
            }

            if (options_intern.tgen_mt_writing.equals("")) {
                logger.logLine("[TGENE] Please specify spelling of Mitochondrial DNA, e.g. M or MT (default: MT)");
                all_set = false;
            }

            if (options_intern.tgen_consensus == 0.0) {
                logger.logLine(
                        "[TGENE] tgen_consensus must be in range ]0.0,1.0], it cannot be 0.0, if you do not want to use consensus set path_tgen=\"\"");
                all_set = false;
            }

            if (options_intern.tgen_self_regulatory) {
                if (options_intern.tgen_consensus_calc.equals("")) {
                    logger.logLine("[TGENE] tgen_consensus_calc must be set.");
                    all_set = false;
                }
            }
        }

        /**
         * DYNAMITE OPTIONS
         */
        if (!options_intern.dynamite_preprocessing_integrate_data_consider_geneFile.equals("")) {
            File f = new File(options_intern.dynamite_preprocessing_integrate_data_consider_geneFile);
            if (!f.exists()) {
                logger.logLine("[DYNAMITE] preprocessing consider genes file for integrateData.py does not exist!");
                all_set = false;
            }
        }

        /**
         * PLOT OPTIONS
         */
        if (options_intern.plot_th_coefficient.isEmpty()) {
            logger.logLine("[PLOTS] plot th coefficients is empty, please use at least one coefficient.");
            all_set = false;
        }
        if (options_intern.plot_cutoff_tps < 1) {
            logger.logLine("[PLOTS] plot_cutoff_tps must be >= 1");
            all_set = false;
        }
        if (options_intern.plot_cutoff_hms < 1) {
            logger.logLine("[PLOTS] plot_cutoff_hms must be >= 1");
            all_set = false;
        }
        if (options_intern.plot_cutoff_gcs < 0) {
            logger.logLine("[PLOTS] plot_cutoff_gcs must be >= 0");
            all_set = false;
        }
        if (options_intern.plot_cutoff_tpms < 0) {
            logger.logLine("[PLOTS] plot_cutoff_tpms must be >= 0.0");
            all_set = false;
        }
        if (options_intern.plot_top_k_genes < 1) {
            logger.logLine("[PLOTS] plot_top_k_genes must be >= 1");
            all_set = false;
        }
        if (options_intern.plot_mann_whitneyU_pvalue_cutoff <= 0) {
            logger.logLine("[PLOTS] Mann WhitneyU pvalue cutoff must be > 0");
            all_set = false;
        }
        if (options_intern.plot_distribution_analysis_score_type.equals("")) {
            logger.logLine(
                    "[PLOTS] plot_distribution_analysis_score_type must be either set to GENE_COUNTS or EXCL_GENE_COUNTS");
            all_set = false;
        }
        if (!options_intern.plot_distribution_analysis_score_type.equals("GENE_COUNTS") &&
                !options_intern.plot_distribution_analysis_score_type.equals("EXCL_GENE_COUNTS")) {
            logger.logLine(
                    "[PLOTS] plot_distribution_analysis_score_type must be either set to GENE_COUNTS or EXCL_GENE_COUNTS");
            all_set = false;
        }

        /**
         * ChIP ATLAS OPTIONS
         */
        if (!options_intern.chip_atlas_genome_version.equals("") && !options_intern.chip_atlas_tissue_type.equals("")) {
            options_intern.chip_atlas_activated_chip_atlas = true;
        }
        if (options_intern.chip_atlas_tissue_type.equals("") && !options_intern.chip_atlas_genome_version.equals("")) {
            logger.logLine(
                    "[ChIP-ATLAS] Genome version specified, but no tissue type specified. (chip_atlas_tissue_type).");
            all_set = false;
        }
        if (options_intern.chip_atlas_tissue_type.equals("") && !options_intern.chip_atlas_genome_version.equals("")) {
            logger.logLine(
                    "[ChIP-ATLAS] Tissue type specified, but no genome version specified. (chip_atlas_genome_version).");
            all_set = false;
        }
        if (options_intern.chip_atlas_activated_chip_atlas) {
            if (options_intern.igv_path_to_igv.equals("")) {
                logger.logLine("[ChIP-ATLAS-IGV] ChIP-ATLAS activated, but no IGV path set.");
                all_set = false;
            }
            if (options_intern.igv_species_ref_genome.equals("")) {
                logger.logLine("[ChIP-ATLAS-IGV] ChIP-ATLAS activated, but no IGV ref genome set set.");
                all_set = false;
            }

            try {
                Socket socket = new Socket("127.0.0.1", options_intern.igv_port_number);
                PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
                BufferedReader in = new BufferedReader(new InputStreamReader(socket.getInputStream()));
                out.println("genome " + options_intern.igv_species_ref_genome);
                String response = in.readLine();
                if (!response.equals("OK")) {
                    logger.logLine("[IGV] igv_species_ref_genome not OK!");
                    all_set = false;
                }
                socket.close();
            } catch (Exception e) {
                logger.logLine(
                        "[IGV] IGV option is used, but IGV is not started. Please start IGV, so COM2POSE can listen to its port: " +
                                options_intern.igv_port_number);
                System.exit(1);
            }

        }

        /**
         * IGV OPTIONS
         */
        if (!options_intern.igv_path_to_igv.equals("")) {
            if (!options_intern.igv_path_to_tf_chip_seq.equals("")) {
                File f = new File(options_intern.igv_path_to_tf_chip_seq);
                if (!f.exists()) {
                    logger.logLine("[IGV] path to TF ChIP-seq data does not exist!");
                    all_set = false;
                }
            }

            if (options_intern.igv_path_to_tf_chip_seq.equals("") && !options_intern.chip_atlas_activated_chip_atlas) {
                logger.logLine("[IGV] path to IGV set, but no TF ChIP-seq data provided.");
                all_set = false;
            }

            if (options_intern.igv_species_ref_genome.equals("")) {
                logger.logLine("[IGV] reference genome not provided!");
                all_set = false;
            }

            File f_igv_tools = new File(options_intern.igv_path_to_igv + File.separator + "igvtools");
            if (!f_igv_tools.exists()) {
                all_set = false;
                logger.logLine(
                        "[IGV] igvtools does not exists. Please make sure that you download the full installation of IGV");
            }

            Socket socket = new Socket("127.0.0.1", options_intern.igv_port_number);
            PrintWriter out = new PrintWriter(socket.getOutputStream(), true);
            BufferedReader in = new BufferedReader(new InputStreamReader(socket.getInputStream()));
            out.println("genome " + options_intern.igv_species_ref_genome);
            String response = in.readLine();
            if (!response.equals("OK")) {
                logger.logLine("[IGV] igv_species_ref_genome not OK!");
                all_set = false;
            }
            socket.close();
        }

        if(options_intern.igv_include_prediction_data.size()!=0)
        {
            HashSet<String> available_hms = new HashSet<>();

            //get all histone marks
            File f_input = new File(options_intern.tepic_input_directory);
            for(File f_tp : f_input.listFiles())
            {
                if(f_tp.isDirectory())
                {
                    for(File f_hm : f_tp.listFiles())
                    {
                        available_hms.add(f_hm.getName());
                    }
                }
            }


            for(String s : options_intern.igv_include_prediction_data)
            {
                if(!available_hms.contains(s))
                {
                    all_set=false;
                    logger.logLine("[IGV] ERROR in igv_include_prediction_data");
                    logger.logLine("[IGV] Histone Modification " + s + " does not exist, please check in config file");
                }

            }

        }

        /**
         * CHECK NAMING CONVENTIONS
         */

        if (!all_set) {
            return all_set;
        }

        all_set = check_naming_conventions();


        return all_set;
    }

    /**
     * CHECK naming conventions between RNA-seq and ChIP-seq/ATAC-seq data
     */
    private boolean check_naming_conventions() throws IOException {
        boolean all_set = true;

        HashSet<String> chipseq_tp_names = new HashSet<>();
        HashSet<String> rnaseq_tp_names = new HashSet<>();
        HashSet<String> chipseq_hm_names = new HashSet<>();
        HashSet<String> rna_seq_hm_names = new HashSet<>();


        File f_chipseq_root = new File(options_intern.tepic_input_directory);
        for (File f_tp : f_chipseq_root.listFiles()) {
            if (f_tp.isDirectory()) {
                String tp_name = f_tp.getName();
                chipseq_tp_names.add(tp_name);

                for (File f_hm : f_tp.listFiles()) {
                    if (f_hm.isDirectory()) {
                        String hm_name = f_hm.getName();
                        chipseq_hm_names.add(hm_name);
                    }
                }
            }
        }

        File f_rnaseq_root = new File(options_intern.deseq2_input_directory);

        for (File f_tp : f_rnaseq_root.listFiles()) {
            if (f_tp.isDirectory()) {
                String tp_name = f_tp.getName();
                rnaseq_tp_names.add(tp_name);

                for (File f_hm : f_tp.listFiles()) {
                    if (f_hm.isDirectory()) {
                        String hm_name = f_hm.getName();
                        rna_seq_hm_names.add(hm_name);
                    }
                }
            }
        }

        /*
        for(String control_hm: chipseq_hm_names)
        {
            if(!rna_seq_hm_names.contains(control_hm))
            {
                all_set=false;
                logger.logLine("[NAMING-CONVENTIONS] ERROR: there is a naming mistake for: " + control_hm);
                logger.logLine("[NAMING-CONVENTIONS] ERROR: RNA-seq namings do not match ChIP-seq namings. (Timepoints or Groups)");
                break;
            }
        }

        for(String control_hm: rna_seq_hm_names)
        {
            if(!chipseq_hm_names.contains(control_hm))
            {
                all_set=false;
                logger.logLine("[NAMING-CONVENTIONS] ERROR: there is a naming mistake for: " + control_hm);
                logger.logLine("[NAMING-CONVENTIONS] ERROR: RNA-seq namings do not match ChIP-seq namings. (Histone Modification or ATAC-seq)");
                break;
            }
        }*/

        for (String control_tp : chipseq_tp_names) {
            if (!rnaseq_tp_names.contains(control_tp)) {
                all_set = false;
                logger.logLine("[NAMING-CONVENTIONS] ERROR: there is a naming mistake for: " + control_tp);
                logger.logLine(
                        "[NAMING-CONVENTIONS] ERROR: RNA-seq namings do not match ChIP-seq namings. (Timepoints or Groups)");
                break;
            }
        }

        for (String control_tp : rnaseq_tp_names) {
            if (!chipseq_tp_names.contains(control_tp)) {
                all_set = false;
                logger.logLine("[NAMING-CONVENTIONS] ERROR: there is a naming mistake for: " + control_tp);
                logger.logLine(
                        "[NAMING-CONVENTIONS] ERROR: RNA-seq namings do not match ChIP-seq namings. (Timepoints or Groups)");
                break;
            }
        }


        return all_set;
    }

    private HashMap<String, HashMap<String, HashSet<String>>> checkGroupsTEPIC() {
        HashMap<String, HashMap<String, HashSet<String>>> groups = new HashMap<>();

        File folder = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_tepic_output_raw);

        for (File fileDir : folder.listFiles()) {
            if (fileDir.isDirectory()) {
                HashMap<String, HashSet<String>> n_timepoint = new HashMap<>();
                String n_tp_name = fileDir.getName();
                for (File fileDir2 : fileDir.listFiles()) {
                    String n_hm_name = "";
                    HashSet<String> n_hm = new HashSet<>();
                    n_hm_name = fileDir2.getName();

                    for (File fileDir3 : fileDir2.listFiles()) {
                        //samples
                        if (fileDir3.isDirectory()) {
                            n_hm.add(fileDir3.getName());
                        }
                    }
                    n_timepoint.put(n_hm_name, n_hm);

                }
                groups.put(n_tp_name, n_timepoint);
            }
        }
        return groups;
    }

    /**
     * BLACKLIST preprocessing - creates binary tree friendly files of chromosomes
     */
    private ArrayList<BL_ranges_binary_tree> recursive_split_BL(ArrayList<BL_ranges_binary_tree> region,
                                                                ArrayList<BL_ranges_binary_tree> newly_ordered) {
        if (region.size() > 0) {
            BL_ranges_binary_tree median = region.get(region.size() / 2);
            newly_ordered.add(median);

            ArrayList<BL_ranges_binary_tree> region_left = new ArrayList<>();

            for (int i = 0; i < region.size() / 2; i++) {
                region_left.add(region.get(i));
            }
            ArrayList<BL_ranges_binary_tree> region_right = new ArrayList<>();
            for (int i = region.size() / 2 + 1; i < region.size(); i++) {
                region_right.add(region.get(i));
            }

            recursive_split_BL(region_left, newly_ordered);
            recursive_split_BL(region_right, newly_ordered);

        }
        return newly_ordered;
    }

    /**
     * TGENE preprocessing - creates binary tree friendly files of chromosomes
     */
    private ArrayList<ENSG_ranges_binary_trees> recursive_split(ArrayList<ENSG_ranges_binary_trees> region,
                                                                ArrayList<ENSG_ranges_binary_trees> newly_ordered) {
        if (region.size() > 0) {
            ENSG_ranges_binary_trees median = region.get(region.size() / 2);
            newly_ordered.add(median);

            ArrayList<ENSG_ranges_binary_trees> region_left = new ArrayList<>();
            for (int i = 0; i < region.size() / 2; i++) {
                region_left.add(region.get(i));
            }
            ArrayList<ENSG_ranges_binary_trees> region_right = new ArrayList<>();
            for (int i = region.size() / 2 + 1; i < region.size(); i++) {
                region_right.add(region.get(i));
            }

            recursive_split(region_left, newly_ordered);
            recursive_split(region_right, newly_ordered);

        }
        return newly_ordered;
    }

    /**
     * mixes the samples of one folder into one file, based on mix_option (UNION or INTERSECTION)
     */
    private static Stack<MIX_Interval> mergeIntervals(ArrayList<MIX_Interval> interval) {
        Stack<MIX_Interval> stack = new Stack<>();

        if (stack.empty()) {
            stack.push(interval.get(0));
        }

        for (int i = 1; i < interval.size(); i++) {
            MIX_Interval top = stack.peek();

            if (top.end < interval.get(i).start) {
                stack.push(interval.get(i));
            } else if (top.end < interval.get(i).end) {
                top.end = interval.get(i).end;
                top.merged_intervals.addAll(interval.get(i).merged_intervals);

                stack.pop();
                stack.push(top);
            } else {
                top.merged_intervals.addAll(interval.get(i).merged_intervals);

                stack.pop();
                stack.push(top);
            }
        }
        return stack;
    }

    private void write_target_genes_of_tf(File f_input_target_genes_hm_group_clash, String timepoint,
                                          File f_out_hm_th_file, String tf,
                                          HashMap<String, String> ensg_gene_symbol_map) throws IOException {

        File f_output = new File(f_out_hm_th_file + File.separator + timepoint);
        f_output.mkdir();

        if (!f_input_target_genes_hm_group_clash.exists()) {
            String[] split = f_input_target_genes_hm_group_clash.getName().split("_");

            String[] split_slashes = f_input_target_genes_hm_group_clash.getAbsolutePath().split(File.separator);
            String path = "";
            for (int i = 0; i < split_slashes.length - 1; i++) {
                path += File.separator + split_slashes[i];
            }
            path += File.separator + split[1] + "_" + split[0];

            f_input_target_genes_hm_group_clash = new File(path);
        }

        File f_input;
        File parent;

        if (options_intern.tepic_tpm_cutoff > 0) {
            f_input = new File(f_input_target_genes_hm_group_clash.getAbsolutePath() + File.separator + timepoint);
            parent = new File(f_input_target_genes_hm_group_clash.getAbsolutePath() + File.separator + timepoint);
        } else {
            parent = new File(options_intern.com2pose_working_directory + File.separator +
                    options_intern.folder_name_tepic_output_raw + File.separator + timepoint + File.separator +
                    f_input_target_genes_hm_group_clash.getParentFile().getName() + File.separator);
            f_input = new File("");
        }

        String suffix = "";

        if (options_intern.tepic_tpm_cutoff > 0) {
            suffix = "_Gene_View_Filtered_TPM.txt";
        } else {
            suffix = "_Gene_View_Filtered.txt";
        }

        for (int j = 0; j < parent.listFiles().length; j++) {
            if (options_intern.tepic_tpm_cutoff > 0) {
                j = parent.listFiles().length;
            } else {
                f_input = parent.listFiles()[j];
            }

            for (File fileDir : f_input.listFiles()) {
                if (!fileDir.getName().matches(".*" + suffix + ".*")) {
                    continue;
                }
                if (fileDir.isFile()) {
                    BufferedReader br = new BufferedReader(new FileReader(fileDir));
                    String line = br.readLine();

                    int interesting_column = -1;

                    String[] header = line.split("\t");
                    for (int i = 0; i < header.length; i++) {
                        if (header[i].toUpperCase().matches(".*" + tf.toUpperCase() + ".*")) {
                            interesting_column = i;
                        }
                    }

                    if (interesting_column == -1) {
                        return;
                    }

                    ArrayList<Gene_Affinity_Value> all_affinities = new ArrayList<>();

                    while ((line = br.readLine()) != null) {
                        String[] split = line.split("\t");

                        Gene_Affinity_Value gav = new Gene_Affinity_Value();
                        gav.gene_name = split[0];
                        if (ensg_gene_symbol_map.containsKey(gav.gene_name)) {
                            gav.gene_symbol = ensg_gene_symbol_map.get(gav.gene_name);
                        }
                        gav.affinity_value = Double.parseDouble(split[interesting_column]);
                        all_affinities.add(gav);
                    }
                    br.close();

                    Collections.sort(all_affinities);

                    BufferedWriter bw = new BufferedWriter(
                            new FileWriter(new File(f_output.getAbsolutePath() + File.separator + tf + ".csv")));

                    bw.write("ENSG\tSYMBOL\tAFFINITY");
                    bw.newLine();

                    for (int i = 0; i < options_intern.plot_top_k_genes; i++) {
                        if (all_affinities.size() > i) {
                            bw.write(all_affinities.get(i).toString());
                            bw.newLine();
                        }

                    }

                    bw.close();
                }
            }
        }


    }

    private String write_regression_coeffecient_analysis_found_table_html(Double d, String level) throws IOException {

        File input_dir_root = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_analysis_data +
                        File.separator + options_intern.folder_out_analysis_data_WEBSITE_OVERVIEW);
        //HashMap<String,File> hm_file = new HashMap<>();
        HashMap<String, HashMap<String, Boolean>> hm_tf_found = new HashMap<>();

        String hm = "";

        for (File fileDir : input_dir_root.listFiles()) {
            HashMap<String, Boolean> tf_found = new HashMap<>();

            File f = new File(fileDir.getAbsolutePath() + File.separator + d + File.separator +
                    options_intern.file_suffix_website_analysis_tf_available);
            //hm_file.put(fileDir.getName(),f);
            BufferedReader br = new BufferedReader(new FileReader(f));
            String line = br.readLine();
            while ((line = br.readLine()) != null) {
                String[] split = line.split("\t");
                tf_found.put(split[0], Boolean.parseBoolean(split[1]));
            }
            hm_tf_found.put(fileDir.getName(), tf_found);
            hm = fileDir.getName();


            br.close();
        }

        HashMap<String, Boolean> first = hm_tf_found.get(hm);

        StringBuilder sb = new StringBuilder();
        sb.append("<div class=\"w3-content\">");
        sb.append("<h2> Threshold:" + d + "</h2>\n");
        sb.append("\t\t<table style=\"width:80%\">\n");
        sb.append("\t\t\t<tr>\n");
        sb.append("\t\t\t\t<th>");
        sb.append("TF");
        sb.append("\t\t\t\t</th>\n");

        for (String hms_key : hm_tf_found.keySet()) {
            sb.append("\t\t\t\t<th>");
            sb.append(hms_key);
            sb.append("\t\t\t\t</th>\n");

        }
        sb.append("\t\t\t</tr>\n");

        for (String tf : first.keySet()) {
            sb.append("\t\t\t<tr>\n");

            sb.append("\t\t\t\t<th>");
            sb.append(tf);
            sb.append("\t\t\t\t</th>\n");

            for (String key_hm : hm_tf_found.keySet()) {
                HashMap<String, Boolean> current_tf_list = hm_tf_found.get(key_hm);
                if (current_tf_list.get(tf)) {
                    sb.append("\t\t\t\t<th>");
                    if (level.equals(options_intern.html_report_levels_home)) {
                        sb.append("<img src=\"" + options_intern.folder_out_website_basics + File.separator +
                                options_intern.folder_out_website_basics_website + File.separator +
                                options_intern.folder_out_website_basics_website_images + File.separator +
                                "is_available.png" + "\" style=\"width:50px;height:50px;\"/>");

                    }
                    if (level.equals(options_intern.html_report_levels_3_steps)) {
                        sb.append("<img src=\".." + File.separator + ".." + File.separator +
                                options_intern.folder_out_website_basics + File.separator +
                                options_intern.folder_out_website_basics_website + File.separator +
                                options_intern.folder_out_website_basics_website_images + File.separator +
                                "is_available.png" + "\" style=\"width:50px;height:50px;\"/>");
                    }
                    sb.append("\t\t\t\t</th>\n");
                } else {
                    sb.append("\t\t\t\t<th>");
                    if (level.equals(options_intern.html_report_levels_home)) {
                        sb.append("<img src=\"" + options_intern.folder_out_website_basics + File.separator +
                                options_intern.folder_out_website_basics_website + File.separator +
                                options_intern.folder_out_website_basics_website_images + File.separator +
                                "not_available.png" + "\" style=\"width:50px;height:50px;\"/>");

                    }
                    if (level.equals(options_intern.html_report_levels_3_steps)) {
                        sb.append("<img src=\".." + File.separator + ".." + File.separator +
                                options_intern.folder_out_website_basics + File.separator +
                                options_intern.folder_out_website_basics_website + File.separator +
                                options_intern.folder_out_website_basics_website_images + File.separator +
                                "not_available.png" + "\" style=\"width:50px;height:50px;\"/>");
                    }
                    sb.append("\t\t\t\t</th>\n");
                }
            }
            sb.append("\t\t\t</tr>\n");
        }

        sb.append("</table>");
        sb.append("<div>");

        return sb.toString();
    }

    private String get_header_html(String level, String which_analysis) {

        File f_website_css = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_website +
                        File.separator + options_intern.folder_out_website_basics + File.separator +
                        options_intern.folder_out_website_basics_website + File.separator +
                        options_intern.folder_out_website_basics_website_css);

        File f_output_website = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_website);
        File f_output_website_htmls = new File(f_output_website.getAbsolutePath() + File.separator +
                options_intern.folder_out_website_htmls_regression_coefficients);

        StringBuilder sb_home_front = new StringBuilder();
        sb_home_front.append(
                "<!DOCTYPE html>\n" + "<html lang=\"en\">\n" + "<title>COM2POSE: " + which_analysis + "</title>\n" +
                        "<meta charset=\"UTF-8\">\n" +
                        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n");
        sb_home_front.append("<head>\n<script>\n" + "  function expand_collapse(id, div_id) {\n" +
                "  var x = document.getElementById(id).getAttribute(\"aria-expanded\"); \n" +
                "  if (x == \"true\") \n" + "  {\n" + "  x = \"false\"\n" + "  } else {\n" + "  x = \"true\"\n" +
                "  }\n" + "  document.getElementById(id).setAttribute(\"aria-expanded\", x);\n" +

                "  var y = document.getElementById(div_id); \n" + "  if (y.style.display === \"none\") \n" + "  {\n" +
                "  y.style.display = \"block\"\n" + "  } else {\n" + "  y.style.display = \"none\"\n" + "  }\n" +
                "  }\n" + "</script>\n</head>\n");

        for (File fileDir : f_website_css.listFiles()) {
            if (fileDir.isFile()) {
                sb_home_front.append("<link rel=\"stylesheet\" href=\"");
                if (level.equals(options_intern.html_report_levels_home)) {
                    sb_home_front.append(options_intern.folder_out_website_basics + File.separator +
                            options_intern.folder_out_website_basics_website + File.separator +
                            options_intern.folder_out_website_basics_website_css + File.separator + fileDir.getName());
                }
                if (level.equals(options_intern.html_report_levels_2_steps)) {
                    sb_home_front.append(
                            ".." + File.separator + options_intern.folder_out_website_basics + File.separator +
                                    options_intern.folder_out_website_basics_website + File.separator +
                                    options_intern.folder_out_website_basics_website_css + File.separator +
                                    fileDir.getName());

                }
                if (level.equals(options_intern.html_report_levels_3_steps)) {
                    sb_home_front.append(
                            ".." + File.separator + ".." + File.separator + options_intern.folder_out_website_basics +
                                    File.separator + options_intern.folder_out_website_basics_website + File.separator +
                                    options_intern.folder_out_website_basics_website_css + File.separator +
                                    fileDir.getName());

                }
                if (level.equals(options_intern.html_report_levels_4_steps)) {
                    sb_home_front.append(".." + File.separator + ".." + File.separator + ".." + File.separator +
                            options_intern.folder_out_website_basics + File.separator +
                            options_intern.folder_out_website_basics_website + File.separator +
                            options_intern.folder_out_website_basics_website_css + File.separator + fileDir.getName());
                }
                sb_home_front.append("\">\n");
            }
        }

        String rel_path = "";
        if (level.equals(options_intern.html_report_levels_home)) {
            rel_path = options_intern.html_report_home_regression_coefficient_analysis;
        }
        if (level.equals(options_intern.html_report_levels_2_steps)) {
            rel_path = ".." + File.separator + options_intern.html_report_home_regression_coefficient_analysis;
        }
        if (level.equals(options_intern.html_report_levels_3_steps)) {
            rel_path = ".." + File.separator + ".." + File.separator +
                    options_intern.html_report_home_regression_coefficient_analysis;
        }
        if (level.equals(options_intern.html_report_levels_4_steps)) {
            rel_path = ".." + File.separator + ".." + File.separator + ".." + File.separator +
                    options_intern.html_report_home_regression_coefficient_analysis;
        }


        String rel_path_distribution_analysis = "";
        if (level.equals(options_intern.html_report_levels_home)) {
            rel_path_distribution_analysis = options_intern.html_report_home_distribution_analysis;
        }
        if (level.equals(options_intern.html_report_levels_2_steps)) {
            rel_path_distribution_analysis =
                    ".." + File.separator + options_intern.html_report_home_distribution_analysis;
        }
        if (level.equals(options_intern.html_report_levels_3_steps)) {
            rel_path_distribution_analysis = ".." + File.separator + ".." + File.separator +
                    options_intern.html_report_home_distribution_analysis;
        }
        if (level.equals(options_intern.html_report_levels_4_steps)) {
            rel_path_distribution_analysis = ".." + File.separator + ".." + File.separator + ".." + File.separator +
                    options_intern.html_report_home_distribution_analysis;
        }


        sb_home_front.append("<style>\n" + "body,h1,h2,h3,h4,h5,h6 {font-family: \"Lato\"; width: \"max-content\";}\n" +
                ".w3-bar,h1\n" +
                //,button {font-family: "Montserrat", sans-serif}
                ".fa-anchor,.fa-coffee {font-size:200px}\n" + ".button {\n" +
                "  background-color: #19631c; /* Green */\n" + "  border: none;\n" + "  color: white;\n" +
                "  padding: 14px 31px;\n" + "  text-align: center;\n" + "  text-decoration: none;\n" +
                "  display: inline-block;\n" + "  font-size: 15px;\n" + "  margin: 4px 2px;\n" +
                "  cursor: pointer;\n" + "  width: 180px;\n" + "}\n" + "</style>\n" + "<body>\n" + "  \n" +
                "<!-- Header -->\n" +
                "<header class=\"w3-container w3-red w3-center\" style=\"padding:128px 16px\">\n" + " <a href='" +
                rel_path + "' target='_blank'><button class=\"button\">HOME coefficient analysis</button></a>\n" +
                " <a href='" + rel_path_distribution_analysis +
                "' target='_blank'><button class=\"button\">HOME distribution analysis</button></a>\n");

        String rel_path_parameter = "";
        if (level.equals(options_intern.html_report_levels_home)) {
            rel_path_parameter = "PARAMETERS.html";
        }
        if (level.equals(options_intern.html_report_levels_2_steps)) {
            rel_path_parameter = ".." + File.separator + "PARAMETERS.html";
        }
        if (level.equals(options_intern.html_report_levels_3_steps)) {
            rel_path_parameter = ".." + File.separator + ".." + File.separator + "PARAMETERS.html";
        }
        if (level.equals(options_intern.html_report_levels_4_steps)) {
            rel_path_parameter =
                    ".." + File.separator + ".." + File.separator + ".." + File.separator + "PARAMETERS.html";
        }

        sb_home_front.append("<a href='" + rel_path_parameter +
                "' target='_blank'><button class=\"button\" id=\"button_parameters\" > Parameters</button></a>\n");

        sb_home_front.append("  <h1 class=\"w3-margin w3-jumbo\">COM2POSE: <i>results overview</i></h1>\n" + "  \n");

        if (which_analysis.equals(options_intern.analysis_types_distribution_analysis)) {
            sb_home_front.append(
                    "  <h2 class=\"w3-center\" style='margin-left:27%'><i>TF-TG Score Distribution Analysis</i></h2>\n");
        }
        if (!which_analysis.equals(options_intern.analysis_types_distribution_analysis) &&
                !which_analysis.equals(options_intern.analysis_types_regression_coefficient_analysis)) {
            sb_home_front.append(
                    "  <h2 class=\"w3-center\" style='margin-left:27%'><i>" + which_analysis + "</i></h2>\n");
        }

        if (which_analysis.equals(options_intern.analysis_types_regression_coefficient_analysis)) {
            sb_home_front.append(
                    "  <h2 class=\"w3-center\" style='margin-left:27%'><i>Regression Coefficient Analysis</i></h2>\n");

            sb_home_front.append("    <div class=\"dropdown\">\n" + "  <h2>Please choose a threshold</h2>\n" +
                    "  <button class=\"dropbtn\">Threshold</button>\n" + "  <div class=\"dropdown-content\">\n");


            for (Double d : options_intern.plot_th_coefficient) {
                File f_output_website_htmls_th =
                        new File(f_output_website_htmls.getAbsolutePath() + File.separator + d);
                f_output_website_htmls_th.mkdir();

                if (!threshold_folders_filled) {
                    threshold_folders.add(f_output_website_htmls_th);
                }

                String rel_path_2 = "";
                if (level.equals(options_intern.html_report_levels_home)) {
                    rel_path_2 = options_intern.folder_out_website_htmls_regression_coefficients + File.separator + d +
                            File.separator + "threshold_" + d + "_overview.html";
                }
                if (level.equals(options_intern.html_report_levels_2_steps)) {
                    rel_path_2 = d + File.separator + "threshold_" + d + "_overview.html";
                }
                if (level.equals(options_intern.html_report_levels_3_steps)) {
                    rel_path_2 = ".." + File.separator + d + File.separator + "threshold_" + d + "_overview.html";
                }
                if (level.equals(options_intern.html_report_levels_4_steps)) {
                    rel_path_2 = ".." + File.separator + ".." + File.separator + d + File.separator + "threshold_" + d +
                            "_overview.html";
                }


                sb_home_front.append("<a href=\"");
                sb_home_front.append(rel_path_2 + "\" target='_blank'>Coefficient " + d);
                sb_home_front.append("</a>\n");
            }
            sb_home_front.append("  </div>\n");

        }


        sb_home_front.append("</div>\n" + "</header>\n");

        return sb_home_front.toString();
    }

    private void write_python_script_distribution_analysis(File input_background_file, File input_tf_root,
                                                           File output_plots, File output_script_file,
                                                           File output_stats) throws IOException {
        HashMap<String, String> composed_tfs = new HashMap<>();
        HashMap<String, HashSet<String>> composed_tfs_tfs = new HashMap<>();

        BufferedReader br_composed_tfs = new BufferedReader(new FileReader(
                options_intern.com2pose_working_directory + File.separator +
                        options_intern.folder_name_tepic_postprocessing + File.separator +
                        options_intern.folder_name_tepic_postprocessing_tfs + File.separator +
                        options_intern.file_suffix_tepic_postprocessing_tfs_tfs));
        String line_composed_tfs = "";
        while ((line_composed_tfs = br_composed_tfs.readLine()) != null) {
            String[] split = line_composed_tfs.split("\t");

            HashSet temp_set = new HashSet();
            for (int i = 1; i < split.length; i++) {
                composed_tfs.put(split[i], split[0]);
                temp_set.add(split[i]);
            }
            composed_tfs_tfs.put(split[0], temp_set);
        }
        br_composed_tfs.close();

        String imports = "import pip\n" + "\n" + "def import_or_install(package):\n" + "    try:\n" +
                "        __import__(package)\n" + "    except ImportError:\n" +
                "        pip.main(['install', package])\n" + "\n" + "import io\n" + "from base64 import b64encode\n" +
                "import_or_install(\"plotly.express\")\n" + "import plotly.express as px\n" +
                "import_or_install(\"dash\")\n" + "import_or_install(\"dash_core_components\")\n" +
                "import_or_install(\"dash_html_components\")\n" + "import dash_core_components as dcc\n" +
                "import dash_html_components as html\n" + "from dash.dependencies import Input, Output\n" +
                "import plotly.graph_objs as go\n" + "\n" + "import_or_install(\"pandas\")\n" +
                "import_or_install(\"seaborn\")\n" + "import_or_install(\"matplotlib.pyplot\")\n" +
                "import pandas as pd\n" + "import seaborn as sns\n" + "import matplotlib.pyplot as plt\n" +
                "sns.set_context(\"notebook\")\n" + "color = \"#A6CEE3\"\n" + "sns.set_context(\"talk\")\n" +
                "sns.set_style(\"whitegrid\")\n" + "\n" + "import_or_install(\"numpy\")\n" + "import numpy as np\n" +
                "import_or_install(\"sts\")\n" + "import statistics as sts\n" + "import scipy.stats as scp\n" +
                "plt.figure(figsize=(20, 17))\n\n\n" +
                "df_interesting_stats=pd.DataFrame(columns=['label','sum_all_values','number_target_genes','mean','median','95_quantile','99_quantile'])\n" +
                "row_counter=0\n\n" + "";

        StringBuilder sb_all = new StringBuilder(imports);

        sb_all.append("background=pd.read_table('");
        sb_all.append(input_background_file);
        sb_all.append("', comment=\"#\", usecols=['TF_TG_SCORE']).sort_values(['TF_TG_SCORE'], ascending=False)\n");
        sb_all.append("background[\"label\"] = \"background\"\n" +
                "background.dropna(subset = [\"TF_TG_SCORE\"], inplace=True)\n" + "\n" +
                "background_sum = sum(background[\"TF_TG_SCORE\"])\n" + "background_length = len(background)\n" +
                "background_mean = background_sum/background_length\n" +
                "background_median=sts.median(background['TF_TG_SCORE'])\n" +
                "background_quantile=np.percentile(background[\"TF_TG_SCORE\"],95)\n" +
                "background_quantile_99=np.percentile(background[\"TF_TG_SCORE\"],99)\n" +
                "df_interesting_stats.loc[0]=['background',background_sum,background_length,background_mean,background_median,background_quantile,background_quantile_99]\n" +
                "row_counter=1" + "\n\n\n");

        for (File fileDir : input_tf_root.listFiles()) {
            String name_tf = fileDir.getName().split("\\.")[0].split("_")[0];
            String name_composed = "";

            if (composed_tfs.containsKey(name_tf)) {
                name_composed = composed_tfs.get(name_tf);
            }

            sb_all.append(name_tf);
            sb_all.append("=pd.read_table('");
            sb_all.append(fileDir.getAbsolutePath());
            sb_all.append(
                    "', comment=\"#\", usecols=['TF_TG_SCORE','TF']).sort_values(['TF_TG_SCORE'], ascending=False)\n");
            sb_all.append(name_tf);
            sb_all.append(".columns=['TF_TG_SCORE','label']\n");
            sb_all.append(name_tf + ".dropna(subset = [\"TF_TG_SCORE\"], inplace=True)\n");
            sb_all.append(
                    name_tf + "_sum = sum(" + name_tf + "['TF_TG_SCORE'])\n" + name_tf + "_length = len(" + name_tf +
                            ")\n" + name_tf + "_mean=0\n" + name_tf + "_quantile=0\n" + name_tf + "_quantile_95=0\n" +
                            name_tf + "_median=0\n" + name_tf +
                            "_mannwhitneyU=scp.mannwhitneyu(background['TF_TG_SCORE']," + name_tf +
                            "['TF_TG_SCORE'])\n" + "if(" + name_tf + "_length>0):\n" + "    " + name_tf + "_mean=" +
                            name_tf + "_sum/" + name_tf + "_length\n");
            sb_all.append("    " + name_tf + "_quantile=np.percentile(" + name_tf + "['TF_TG_SCORE'], 99)\n");
            sb_all.append("    " + name_tf + "_quantile_95=np.percentile(" + name_tf + "['TF_TG_SCORE'], 95)\n");
            sb_all.append("    " + name_tf + "_median=sts.median(" + name_tf + "['TF_TG_SCORE'])\n");
            sb_all.append("if(" + name_tf + "_median > background_median and " + name_tf + "_mannwhitneyU[1]<" +
                    options_intern.plot_mann_whitneyU_pvalue_cutoff + "):\n");
            sb_all.append(
                    "    background_" + name_tf + " = pd.concat([background," + name_tf + "],axis=0)\n" + "    ax_" +
                            name_tf + " = sns.boxplot(x=\"label\", y=\"TF_TG_SCORE\",data=background_" + name_tf +
                            ",palette=\"Set3\")\n" + "    ax_" + name_tf + ".set_yscale(\"log\")\n");
            if (name_composed.equals("")) {
                sb_all.append(
                        "    plt.savefig(f'" + output_plots.getAbsolutePath() + File.separator + name_tf + ".png')\n");
            } else {
                sb_all.append("    plt.savefig(f'" + output_plots.getAbsolutePath() + File.separator + name_composed +
                        ".png')\n");
            }
            sb_all.append("    del background_" + name_tf + "\n" + "    plt.clf()\n");
            if (name_composed.equals("")) {
                sb_all.append(
                        "    df_interesting_stats.loc[row_counter]=['" + name_tf + "'," + name_tf + "_sum," + name_tf +
                                "_length," + name_tf + "_mean," + name_tf + "_median," + name_tf + "_quantile_95," +
                                name_tf + "_quantile]\n");
            } else {
                sb_all.append(
                        "    df_interesting_stats.loc[row_counter]=['" + name_composed + "'," + name_tf + "_sum," +
                                name_tf + "_length," + name_tf + "_mean," + name_tf + "_median," + name_tf +
                                "_quantile_95," + name_tf + "_quantile]\n");
            }
            sb_all.append("    row_counter=row_counter+1\n");
            sb_all.append("del " + name_tf + "_sum\n");
            sb_all.append("del " + name_tf + "_length\n");
            sb_all.append("del " + name_tf + "_mean\n");
            sb_all.append("del " + name_tf + "_median\n");
            sb_all.append("del " + name_tf + "_quantile_95\n");
            sb_all.append("del " + name_tf + "_quantile\n");
            sb_all.append("del " + name_tf + "_mannwhitneyU\n");
            sb_all.append("plt.clf()\n");
            sb_all.append("del " + name_tf + "\n" + "#plt.figure(figsize=(20, 17))\n\n\n");

        }

        File f_stats_all_print = new File(output_stats.getAbsolutePath() + File.separator +
                options_intern.file_suffix_distribution_analysis_plot_stats);

        sb_all.append("df_interesting_stats.to_csv('" + f_stats_all_print.getAbsolutePath() + "',sep='\t')\n");

        BufferedWriter bw_all = new BufferedWriter(new FileWriter(output_script_file.getAbsolutePath()));
        bw_all.write(sb_all.toString());
        bw_all.close();
    }

    private void write_html_distribution_analysis_plots_hm_page(File f_website_html_hm, File fileDir_stat, String level,
                                                                HashMap<String, HashMap<String, String>> tp_tf_gene_count)
            throws IOException {

        Analysis_distribution_stats_object stats_object = get_distribution_analysis_stats_ordered(fileDir_stat);

        ArrayList<Analysis_distribution_stats> all_considered_tfs = stats_object.ordered_tfs;
        Analysis_distribution_stats background = stats_object.background;

        /*
        //read in stats
        BufferedReader br = new BufferedReader(new FileReader(fileDir_stat+File.separator+options_intern.file_suffix_distribution_analysis_plot_stats));
        String line = br.readLine();
        line = br.readLine();

        ArrayList<Analysis_distribution_stats> all_considered_tfs = new ArrayList<>();

        String[] split_background = line.split("\t");
        Analysis_distribution_stats background = new Analysis_distribution_stats();
        background.label = split_background[1];
        background.sum_all_values = Double.parseDouble(split_background[2]);
        background.number_target_genes = Double.parseDouble(split_background[3]);
        background.mean=Double.parseDouble(split_background[4]);
        background.median=Double.parseDouble(split_background[5]);
        background.quantile_95=Double.parseDouble(split_background[6]);
        background.quantile_99=Double.parseDouble(split_background[7]);

        while((line=br.readLine())!=null)
        {
            String[] split = line.split("\t");

            Analysis_distribution_stats tf = new Analysis_distribution_stats();
            tf.label = split[1];
            tf.sum_all_values = Double.parseDouble(split[2]);
            tf.number_target_genes = Double.parseDouble(split[3]);
            tf.mean=Double.parseDouble(split[4]);
            tf.median=Double.parseDouble(split[5]);
            tf.quantile_95=Double.parseDouble(split[6]);
            tf.quantile_99=Double.parseDouble(split[7]);
            all_considered_tfs.add(tf);
        }

        //sort stats - median and then quantile - mean?

        Collections.sort(all_considered_tfs);*/

        //setup webpage basics
        String html_tail = "</body>\n" + "</html>";

        StringBuilder sb = new StringBuilder(get_header_html(options_intern.html_report_levels_3_steps,
                "Distribution Analysis: " + fileDir_stat.getName()));

        sb.append(" <script>\n" + " document.title = \"" + fileDir_stat.getName() + "\";\n" + " </script>\n");
        //write stats table inc. comparison to background

        int rank = 1;
        for (Analysis_distribution_stats as : all_considered_tfs) {
            sb.append("<div id='" + as.label.toUpperCase() + "' class='w3-content'>\n");
            sb.append("<h2>" + rank + ". <a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=" +
                    as.label.toUpperCase() + "' target='_blank'><button class='button'>" + as.label.toUpperCase() +
                    "</button></a></h2>\n");

            sb.append("<h5>Distribution Analysis Details:</h5>\n");
            sb.append(as.to_html_String(background));

            ArrayList<String> tps = new ArrayList<>();

            sb.append("<h5>Gene Counts Table </h5>\n");

            sb.append("\t\t<table style=\"width:80%\">\n");
            sb.append("\t\t\t<tr>\n");
            sb.append("\t\t\t\t<th>\n");
            sb.append("TF");
            sb.append("\t\t\t\t</th>\n");

            for (String key : tp_tf_gene_count.keySet()) {
                sb.append("\t\t\t\t<th>\n");
                sb.append(key);
                sb.append("\t\t\t\t</th>\n");

                tps.add(key);
            }
            sb.append("\t\t\t</tr>\n");
            sb.append("\t\t\t<tr>\n");
            sb.append("\t\t\t\t<th>\n");
            sb.append(as.label.toUpperCase());
            sb.append("\t\t\t\t</th>\n");

            for (String key : tps) {
                String value = "0.0";
                HashMap<String, String> current_x = tp_tf_gene_count.get(key);
                if (current_x.containsKey(as.label.toUpperCase())) {
                    value = current_x.get(as.label.toUpperCase());
                }

                sb.append("\t\t\t\t<th>");
                sb.append(value);
                sb.append("\t\t\t\t</th>\n");
            }

            sb.append("\t\t\t</tr>\n");

            sb.append("\t\t</table>");

            String rel_path = ".." + File.separator + ".." + File.separator +
                    options_intern.folder_out_website_plots_distribution_analysis + File.separator +
                    options_intern.folder_out_distribution_plots;
            if (level.equals("HM")) {
                rel_path += File.separator + options_intern.folder_out_distribution_plots_HM + File.separator +
                        fileDir_stat.getName() + File.separator + as.label.toUpperCase() + ".png";
            } else {
                rel_path += File.separator + options_intern.folder_out_distribution_plots_ALL + File.separator +
                        as.label.toUpperCase() + ".png";

            }
            sb.append("<img src=\"" + rel_path + "\" style=\"width:900px;height:800px;\"/>\n");
            sb.append("</div>\n");

            rank++;
        }


        //include img

        sb.append(html_tail);

        BufferedWriter bw = new BufferedWriter(new FileWriter(f_website_html_hm));
        bw.write(sb.toString());
        bw.close();


    }

    private String get_html_table_found_tfs_in_distribution_analysis(ArrayList<File> files_stat, String level,
                                                                     HashMap<String, Integer> fin_rank)
            throws Exception {
        StringBuilder sb_table = new StringBuilder();

        HashMap<String, ArrayList<Analysis_distribution_stats>> hm_distribution_stats = new HashMap<>();
        ArrayList<String> hm_names_order = new ArrayList<>();

        for (File fileDir_hm : files_stat) {
            String name = fileDir_hm.getName();

            //File stats_file = new File(fileDir_hm.getAbsolutePath()+File.separator+options_intern.file_suffix_distribution_analysis_plot_stats);

            Analysis_distribution_stats_object stats_object = get_distribution_analysis_stats_ordered(fileDir_hm);

            ArrayList<Analysis_distribution_stats> all_considered_tfs = stats_object.ordered_tfs;

            hm_distribution_stats.put(name, all_considered_tfs);
            hm_names_order.add(name);
        }
        sb_table.append("<div class='w3-content'>\n");
        sb_table.append("\t\t<table style=\"width:80%;\">\n");
        sb_table.append("\t\t\t<tr>\n");
        sb_table.append("\t\t\t\t<th>\n");
        sb_table.append("TF");
        sb_table.append("\t\t\t\t</th>\n");
        for (String hms : hm_names_order) {
            sb_table.append("\t\t\t\t<th>\n");
            sb_table.append(hms);
            sb_table.append("\t\t\t\t</th>\n");
            sb_table.append("\t\t\t\t<th>\n");
            sb_table.append("RANK");
            sb_table.append("\t\t\t\t</th>\n");
        }
        sb_table.append("\t\t\t\t<th>\n");
        sb_table.append("DISCOUNTED CUMULATIVE GAIN");
        sb_table.append("\t\t\t\t</th>\n");
        sb_table.append("\t\t\t\t<th>\n");
        sb_table.append("RANK");
        sb_table.append("\t\t\t\t</th>\n");
        sb_table.append("\t\t</tr>\n");

        String rel_path = "";
        if (level.equals(options_intern.html_report_levels_home)) {
            rel_path += options_intern.folder_out_website_basics + File.separator +
                    options_intern.folder_out_website_basics_website + File.separator +
                    options_intern.folder_out_website_basics_website_images;
        }
        if (level.equals(options_intern.html_report_levels_2_steps)) {
            rel_path += ".." + File.separator + options_intern.folder_out_website_basics + File.separator +
                    options_intern.folder_out_website_basics_website + File.separator +
                    options_intern.folder_out_website_basics_website_images;
        }
        if (level.equals(options_intern.html_report_levels_3_steps)) {
            rel_path += ".." + File.separator + ".." + File.separator + options_intern.folder_out_website_basics +
                    File.separator + options_intern.folder_out_website_basics_website + File.separator +
                    options_intern.folder_out_website_basics_website_images;
        }
        if (level.equals(options_intern.html_report_levels_4_steps)) {
            rel_path += ".." + File.separator + ".." + File.separator + ".." + File.separator +
                    options_intern.folder_out_website_basics + File.separator +
                    options_intern.folder_out_website_basics_website + File.separator +
                    options_intern.folder_out_website_basics_website_images;
        }

        HashMap<String, HashSet<String>> hm_drawn_white_balls = new HashMap<>();
        HashMap<String, Integer> hm_sample_size = new HashMap<>();

        hm_sample_size.put("DISCOUNTED CUMULATIVE GAIN", fin_rank.size());

        for (String k_tf : options_intern.website_interesting_tfs) {
            sb_table.append("\t\t\t<tr>\n");
            sb_table.append("\t\t\t\t<th>\n");
            sb_table.append(k_tf);
            sb_table.append("\t\t\t\t</th>\n");

            for (String hms : hm_names_order) {
                ArrayList<Analysis_distribution_stats> temp_stats = hm_distribution_stats.get(hms);

                int rank_intern = -1;
                int rank_max = 1;

                for (int i = 0; i < temp_stats.size(); i++) {
                    if (temp_stats.get(i).label.equals(k_tf)) {
                        rank_intern = i + 1;
                    }
                    rank_max++;
                }

                hm_sample_size.put(hms, rank_max);


                if (rank_intern != -1) {
                    sb_table.append("\t\t\t\t<th>\n");
                    sb_table.append("<img src=\"" + rel_path + File.separator +
                            "is_available.png\" style=\"width:50px;height:50px;\"/>");
                    sb_table.append("\t\t\t\t</th>\n");
                    sb_table.append("\t\t\t\t<th>\n");
                    sb_table.append(rank_intern + "/" + rank_max);
                    sb_table.append("\t\t\t\t</th>\n");

                    HashSet<String> drawn_white_balls;
                    if (hm_drawn_white_balls.containsKey(hms)) {
                        drawn_white_balls = hm_drawn_white_balls.get(hms);
                    } else {
                        drawn_white_balls = new HashSet<>();
                    }
                    drawn_white_balls.add(k_tf);
                    hm_drawn_white_balls.put(hms, drawn_white_balls);
                } else {
                    sb_table.append("\t\t\t\t<th>\n");
                    sb_table.append("<img src=\"" + rel_path + File.separator +
                            "not_available.png\" style=\"width:50px;height:50px;\"/>");
                    sb_table.append("\t\t\t\t</th>\n");
                    sb_table.append("\t\t\t\t<th>\n");
                    sb_table.append("-");
                    sb_table.append("\t\t\t\t</th>\n");
                }

            }
            if (fin_rank.containsKey(k_tf)) {
                sb_table.append("\t\t\t\t<th>\n");
                sb_table.append("<img src=\"" + rel_path + File.separator +
                        "is_available.png\" style=\"width:50px;height:50px;\"/>");
                sb_table.append("\t\t\t\t</th>\n");
                sb_table.append("\t\t\t\t<th>\n");
                sb_table.append(fin_rank.get(k_tf) + "/" + fin_rank.size());
                sb_table.append("\t\t\t\t</th>\n");

                HashSet<String> drawn_white_balls;
                if (hm_drawn_white_balls.containsKey("DISCOUNTED CUMULATIVE GAIN")) {
                    drawn_white_balls = hm_drawn_white_balls.get("DISCOUNTED CUMULATIVE GAIN");
                } else {
                    drawn_white_balls = new HashSet<>();
                }
                drawn_white_balls.add(k_tf);
                hm_drawn_white_balls.put("DISCOUNTED CUMULATIVE GAIN", drawn_white_balls);
            } else {
                sb_table.append("\t\t\t\t<th>\n");
                sb_table.append("<img src=\"" + rel_path + File.separator +
                        "not_available.png\" style=\"width:50px;height:50px;\"/>");
                sb_table.append("\t\t\t\t</th>\n");
                sb_table.append("\t\t\t\t<th>\n");
                sb_table.append("-");
                sb_table.append("\t\t\t\t</th>\n");
            }


            sb_table.append("\t\t\t</tr>\n");

        }

        //last row with hypergeometric test results for each group
        int size_white_balls = options_intern.website_interesting_tfs.size();

        //read in all_tfs.csv
        HashSet<String> all_tfs = new HashSet<>();
        File f_input_all_tfs = new File(options_intern.com2pose_working_directory + File.separator +
                options_intern.folder_name_tepic_postprocessing + File.separator +
                options_intern.file_suffix_tepic_postprocessing_all_tfs);
        BufferedReader br_all_tfs = new BufferedReader(new FileReader(f_input_all_tfs));
        String line_all_tfs = "";
        while ((line_all_tfs = br_all_tfs.readLine()) != null) {
            all_tfs.add(line_all_tfs);
        }
        br_all_tfs.close();

        int size_all_balls = all_tfs.size();

        //WRITE AND EXECUTE R_SCRIPT FOR HYPERGEOMETRIC TEST
        File f_out_rscript_hypergeometric_test_directory = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_out_distribution +
                        File.separator + options_intern.folder_out_distribution_hypergeometric_test);
        f_out_rscript_hypergeometric_test_directory.mkdir();

        File f_out_hypergeometric_test_results = new File(
                f_out_rscript_hypergeometric_test_directory.getAbsolutePath() + File.separator +
                        options_intern.file_suffix_distribution_analysis_hypergeometric_test_output);

        StringBuilder sb_rscript_hypergeometric_test = new StringBuilder();

        sb_rscript_hypergeometric_test.append("results<-data.frame(Name=character(),p_value=double())\n");

        for (String key_hm : hm_sample_size.keySet()) {
            if (!hm_drawn_white_balls.containsKey(key_hm)) {
                continue;
            }
            String name = key_hm.replace(" ", "_");
            if (name.startsWith("01")) {
                name = name.substring(3);
            }
            sb_rscript_hypergeometric_test.append(name + "_name<-");
            sb_rscript_hypergeometric_test.append("\"" + key_hm + "\"\n");
            sb_rscript_hypergeometric_test.append(name);
            sb_rscript_hypergeometric_test.append("<-1-phyper(");
            sb_rscript_hypergeometric_test.append(hm_drawn_white_balls.get(key_hm).size());
            sb_rscript_hypergeometric_test.append(",");
            sb_rscript_hypergeometric_test.append(size_white_balls);
            sb_rscript_hypergeometric_test.append(",");
            sb_rscript_hypergeometric_test.append(size_all_balls - size_white_balls);
            sb_rscript_hypergeometric_test.append(",");
            sb_rscript_hypergeometric_test.append(hm_sample_size.get(key_hm));
            sb_rscript_hypergeometric_test.append(")\n");
            sb_rscript_hypergeometric_test.append("results<-rbind(results,c(");
            sb_rscript_hypergeometric_test.append(name + "_name," + name);
            sb_rscript_hypergeometric_test.append("))\n");
        }
        sb_rscript_hypergeometric_test.append(
                "write.csv(results,\"" + f_out_hypergeometric_test_results.getAbsolutePath() + "\", sep=\"\\t\")");


        File f_out_rscript_hypergeometric_test = new File(
                f_out_rscript_hypergeometric_test_directory.getAbsolutePath() + File.separator +
                        options_intern.file_suffix_distribution_analysis_hypergeometric_test_rscript);
        BufferedWriter bw_rscript_hypergeometric_test =
                new BufferedWriter(new FileWriter(f_out_rscript_hypergeometric_test));
        bw_rscript_hypergeometric_test.write(sb_rscript_hypergeometric_test.toString());
        bw_rscript_hypergeometric_test.close();

        String command = "Rscript " + f_out_rscript_hypergeometric_test;
        logger.logLine("[DISTRIBUTION ANALYSIS] Perform hypergeometric test: " + command);
        Process child = Runtime.getRuntime().exec(command);
        int code = child.waitFor();
        switch (code) {
            case 0:
                break;
            case 1:
                String message = child.getErrorStream().toString();
                throw new Exception(message);
        }

        HashMap<String, Double> results_hypergeometric_test = new HashMap<>();
        BufferedReader br_results_hyp_geo_test = new BufferedReader(new FileReader(f_out_hypergeometric_test_results));
        String line_results_hyp_geo_test = br_results_hyp_geo_test.readLine();
        while ((line_results_hyp_geo_test = br_results_hyp_geo_test.readLine()) != null) {
            line_results_hyp_geo_test = line_results_hyp_geo_test.replace("\"", "");
            String[] split = line_results_hyp_geo_test.split(",");
            results_hypergeometric_test.put(split[1], Double.parseDouble(split[2]));
        }
        br_results_hyp_geo_test.close();

        DecimalFormat df = new DecimalFormat("0.00000");

        sb_table.append("\t\t\t<tr>\n");
        sb_table.append("\t\t\t\t<th>\n");
        sb_table.append("HYP.GEO.TEST p_value");
        sb_table.append("\t\t\t\t</th>\n");
        for (String key_hm : hm_names_order) {
            if (!results_hypergeometric_test.containsKey(key_hm)) {
                sb_table.append("\t\t\t\t<th>\n");
                sb_table.append("-");
                sb_table.append("\t\t\t\t</th>\n");
                sb_table.append("\t\t\t\t<th>\n");
                sb_table.append("\t\t\t\t</th>\n");
                continue;
            }
            sb_table.append("\t\t\t\t<th>\n");
            if (results_hypergeometric_test.containsKey(key_hm)) {
                sb_table.append(df.format(results_hypergeometric_test.get(key_hm)));
            } else {
                sb_table.append("NOT PERFORMED");
            }
            sb_table.append("\t\t\t\t</th>\n");
            sb_table.append("\t\t\t\t<th>\n");
            sb_table.append("\t\t\t\t</th>\n");
        }
        sb_table.append("\t\t\t\t<th>\n");
        String key_hypergeometric_test = "DISCOUNTED CUMULATIVE GAIN";
        if (results_hypergeometric_test.containsKey(key_hypergeometric_test)) {
            sb_table.append(df.format(results_hypergeometric_test.get(key_hypergeometric_test)));
        } else {
            sb_table.append("NOT PERFORMED");
        }
        sb_table.append("\t\t\t\t</th>\n");
        sb_table.append("\t\t\t\t<th>\n");
        sb_table.append("\t\t\t\t</th>\n");


        sb_table.append("\t\t\t</tr>\n");


        sb_table.append("</table>\n");
        sb_table.append("</div>\n");

        return sb_table.toString();
    }

    private Analysis_distribution_stats_object get_distribution_analysis_stats_ordered(File fileDir_stat)
            throws IOException {

        Analysis_distribution_stats background = new Analysis_distribution_stats();

        //read in stats
        BufferedReader br = new BufferedReader(new FileReader(
                fileDir_stat + File.separator + options_intern.file_suffix_distribution_analysis_plot_stats));
        String line = br.readLine();
        line = br.readLine();

        ArrayList<Analysis_distribution_stats> all_considered_tfs = new ArrayList<>();

        String[] split_background = line.split("\t");
        background.label = split_background[1];
        background.sum_all_values = Double.parseDouble(split_background[2]);
        background.number_target_genes = Double.parseDouble(split_background[3]);
        background.mean = Double.parseDouble(split_background[4]);
        background.median = Double.parseDouble(split_background[5]);
        background.quantile_95 = Double.parseDouble(split_background[6]);
        background.quantile_99 = Double.parseDouble(split_background[7]);

        while ((line = br.readLine()) != null) {
            String[] split = line.split("\t");

            Analysis_distribution_stats tf = new Analysis_distribution_stats();
            tf.label = split[1];
            tf.sum_all_values = Double.parseDouble(split[2]);
            tf.number_target_genes = Double.parseDouble(split[3]);
            tf.mean = Double.parseDouble(split[4]);
            tf.median = Double.parseDouble(split[5]);
            tf.quantile_95 = Double.parseDouble(split[6]);
            tf.quantile_99 = Double.parseDouble(split[7]);
            all_considered_tfs.add(tf);
        }

        //sort stats - median and then quantile - mean?

        Collections.sort(all_considered_tfs);


        Analysis_distribution_stats_object ret_obj = new Analysis_distribution_stats_object();
        ret_obj.ordered_tfs = all_considered_tfs;
        ret_obj.background = background;

        return ret_obj;
    }

    private void check_tepic_input_with_options() {
        options_intern.tepic_input_prev = options_intern.tepic_input_directory;
        File f_annotation_check = new File(
                options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option +
                        File.separator + options_intern.folder_name_mix_option_preprocessing_check_chr);
        options_intern.tepic_input_directory = f_annotation_check.getAbsolutePath();

        if (options_intern.mix_level.equals("SAMPLE_LEVEL")) {
            File root_mix_working_dir = new File(
                    options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option);
            File f_sample_mix_output = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                    options_intern.folder_name_mix_option_sample_mix);
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = f_sample_mix_output.getAbsolutePath();
        }

        if (options_intern.mix_level.equals("HM_LEVEL")) {
            File root_mix_working_dir = new File(
                    options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option);
            File f_output_hm = new File(root_mix_working_dir.getAbsolutePath() + File.separator +
                    options_intern.folder_name_mix_option_hm_mix);
            f_output_hm.mkdir();

            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = f_output_hm.getAbsolutePath();
        }
        if (options_intern.tepic_tf_binding_site_search.equals("BETWEEN")) {
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory =
                    options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option +
                            File.separator + options_intern.folder_name_mix_options_footprints_between_peaks;

        }
        if (!options_intern.black_list_dir.equals("")) {
            File output_folder = new File(options_intern.com2pose_working_directory + File.separator +
                    options_intern.folder_name_blacklisted_regions);
            File output_folder_new_input = new File(output_folder.getAbsolutePath() + File.separator +
                    options_intern.folder_name_blacklisted_regions_new_input);
            output_folder_new_input.mkdir();

            //set new folder directory for tepic input and save old one
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory = output_folder_new_input.getAbsolutePath();
        }

        if (options_intern.mix_mutually_exclusive) {
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory =
                    options_intern.com2pose_working_directory + File.separator + options_intern.folder_name_mix_option +
                            File.separator + options_intern.folder_name_mix_option_mutually_exclusive + File.separator +
                            options_intern.folder_name_mix_options_mutually_exclusive_input;
        }
    }

    private void igv_save_sessions(File f_output_file, String load_command) throws IOException {
        //save sessions

        File f_session_tp_name = f_output_file;

        Socket socket_session = new Socket("127.0.0.1", options_intern.igv_port_number);
        PrintWriter out_session = new PrintWriter(socket_session.getOutputStream(), true);
        BufferedReader in_session = new BufferedReader(new InputStreamReader(socket_session.getInputStream()));

        out_session.println("genome " + options_intern.igv_species_ref_genome);
        String response_session = in_session.readLine();

        //logger.logLine("[IGV] " + load_tf_chip_seq);
        out_session.println(load_command);
        String response_load_sessoin = in_session.readLine();
        //logger.logLine("[IGV] " + response_load);

        //logger.logLine("[IGV] "+ "genome "+options_intern.igv_species_ref_genome);

        //logger.logLine("[IGV] "+ response);

        //out_session.println(load_command);
        //response_session = in_session.readLine();
        //out_session.println("expand");
        //response_session = in_session.readLine();

        out_session.println("saveSession " + f_session_tp_name.getAbsolutePath());
        response_session = in_session.readLine();
        out_session.println("new");
        response_session = in_session.readLine();

        socket_session.close();
    }
}
