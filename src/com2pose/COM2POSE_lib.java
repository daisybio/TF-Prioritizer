package com2pose;
import util.*;

import javax.swing.*;
import java.io.*;
import java.lang.reflect.Array;
import java.util.*;

public class COM2POSE_lib
{
    public Options_intern options_intern;
    public Logger logger;

    /**
     * constructor for analysis programms, you cannot run the pipeline here
     * @param options_intern
     * @param logger
     */
    public COM2POSE_lib(Options_intern options_intern, Logger logger)
    {
        this.options_intern = options_intern;
        this.logger=logger;
    }

    /**
     * constructor for the pipeline
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
        logger.logLine("COM2POSE path set to: "+ options_intern.path_to_COM2POSE);
    }

    /**
     * analyze joined dataframe data for interesting TFs
     */
    public void analyze_plots_data() throws IOException {

        logger.logLine("[PLOT-ANALYSE] Analyse plot data.");

        File input_read_counts = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_gene_symbols);
        File input_data = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_data_plots);

        File output_analysis_root = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_analysis_data);
        output_analysis_root.mkdir();

        File output_analysis_tp_level = new File(output_analysis_root.getAbsolutePath()+File.separator+options_intern.folder_out_analysis_data_TP_LEVEL);
        output_analysis_tp_level.mkdir();
        File output_analysis_hm_level = new File(output_analysis_root.getAbsolutePath()+File.separator+options_intern.folder_out_analysis_data_HM_LEVEL);
        output_analysis_hm_level.mkdir();

        HashMap<String,HashMap<String,HashMap<String,HashMap<String,Double>>>> cutoff_hm_tf_counts_DIFFERENT = new HashMap<>();
        HashMap<String,HashMap<String,HashMap<String,HashMap<String,Double>>>> cutoff_hm_tf_counts_SAME = new HashMap<>();


        HashMap<String,HashMap<String,Double>> tp_tf_gene_counts = new HashMap<>();

        for(File fileDir:input_read_counts.listFiles())
        {
            String tp_name = fileDir.getName().split("\\.")[0];
            HashMap<String,Double> n_hm = new HashMap<>();

            BufferedReader br = new BufferedReader(new FileReader(fileDir));
            String line = br.readLine();
            while((line=br.readLine())!=null)
            {
                String[] split = line.split("\t");
                n_hm.put(split[0].toUpperCase(),Double.parseDouble(split[2]));

            }
            br.close();

            tp_tf_gene_counts.put(tp_name,n_hm);

        }

        for(File fileDir: input_data.listFiles())
        {
            if(fileDir.isDirectory())
            {
                String hm = fileDir.getName();
                File output_hm = new File(output_analysis_tp_level.getAbsolutePath()+File.separator+hm);
                output_hm.mkdir();

                for(File fileDirHM_th :fileDir.listFiles())
                {
                    if(fileDirHM_th.isDirectory())
                    {
                        String th_name = fileDirHM_th.getName();
                        File output_th_tp_level = new File(output_hm.getAbsolutePath()+File.separator+fileDirHM_th.getName());
                        output_th_tp_level.mkdir();

                        for(File fileDir_HM_th_data : fileDirHM_th.listFiles())
                        {
                            if(fileDir_HM_th_data.isFile())
                            {
                                String name = fileDir_HM_th_data.getName().split("\\.")[0];

                                HashMap<String,HashMap<String,HashMap<String,Double>>> th;
                                HashMap<String,HashMap<String,Double>> hm_intern;

                                if(name.matches(".*different.*"))
                                {
                                    if(cutoff_hm_tf_counts_DIFFERENT.containsKey(fileDirHM_th.getName()))
                                    {
                                        th=cutoff_hm_tf_counts_DIFFERENT.get(fileDirHM_th.getName());
                                        if(th.containsKey(hm))
                                        {
                                            hm_intern=th.get(hm);
                                        }
                                        else
                                        {
                                            hm_intern=new HashMap<>();
                                        }
                                    }
                                    else
                                    {
                                        th=new HashMap<>();
                                        hm_intern=new HashMap<>();
                                    }
                                }
                                else
                                {
                                    if(cutoff_hm_tf_counts_SAME.containsKey(fileDirHM_th.getName()))
                                    {
                                        th=cutoff_hm_tf_counts_DIFFERENT.get(fileDirHM_th.getName());
                                        if(th.containsKey(hm))
                                        {
                                            hm_intern=th.get(hm);
                                        }
                                        else
                                        {
                                            hm_intern=new HashMap<>();
                                        }
                                    }
                                    else
                                    {
                                        th=new HashMap<>();
                                        hm_intern=new HashMap<>();
                                    }
                                }

                                BufferedWriter bw = new BufferedWriter(new FileWriter(output_th_tp_level.getAbsolutePath()+File.separator+fileDir_HM_th_data.getName()));
                                StringBuilder sb = new StringBuilder();
                                sb.append("TF");
                                for(String k: tp_tf_gene_counts.keySet())
                                {
                                    sb.append("\t");
                                    sb.append(k);
                                }
                                bw.write(sb.toString());
                                bw.newLine();

                                BufferedReader br = new BufferedReader(new FileReader(fileDir_HM_th_data));
                                String line = br.readLine();
                                while((line=br.readLine())!=null)
                                {
                                    String[] split = line.split(",");

                                    int count = 0;
                                    for(String s:split)
                                    {
                                        if(!s.equals(""))
                                        {
                                            count++;
                                        }

                                    }

                                    if(count>options_intern.plot_cutoff_tps)
                                    {
                                        HashMap<String,Double> tf_gc = new HashMap<>();
                                        //tf_gc.put(split[0],tp_tf_gene_counts.get(name).get(split[0]));
                                        StringBuilder sb_intern = new StringBuilder();
                                        sb_intern.append(split[0]);
                                        int count_row = 0;

                                        for(String k: tp_tf_gene_counts.keySet())
                                        {
                                            if(tp_tf_gene_counts.get(k).containsKey(split[0].toUpperCase()))
                                            {
                                                count_row+=tp_tf_gene_counts.get(k).get(split[0].toUpperCase());

                                                tf_gc.put(k,tp_tf_gene_counts.get(k).get(split[0].toUpperCase()));
                                                sb_intern.append("\t");
                                                sb_intern.append(tp_tf_gene_counts.get(k).get(split[0].toUpperCase()));
                                            }
                                        }

                                        if(count_row>=options_intern.plot_cutoff_gcs)
                                        {
                                            hm_intern.put(split[0].toUpperCase(),tf_gc);
                                            //WRITE
                                            bw.write(sb_intern.toString());
                                            bw.newLine();
                                        }
                                    }
                                }
                                br.close();
                                bw.close();

                                th.put(hm,hm_intern);

                                //save hashmaps correspondingly
                                if(name.matches(".*different.*"))
                                {
                                    cutoff_hm_tf_counts_DIFFERENT.put(th_name,th);
                                }
                                else
                                {
                                    cutoff_hm_tf_counts_SAME.put(th_name,th);
                                }
                            }
                        }
                    }
                }
            }
        }

        for(String k : cutoff_hm_tf_counts_DIFFERENT.keySet())
        {
            File f_out = new File(output_analysis_hm_level.getAbsolutePath()+File.separator+k);
            f_out.mkdir();

            HashMap<String,HashMap<String,HashMap<String,Double>>> available_hms = cutoff_hm_tf_counts_DIFFERENT.get(k);

            ArrayList<HashMap<String,HashMap<String,Double>>> tf_lists = new ArrayList<>();

            for(String kk:available_hms.keySet())
            {
                tf_lists.add(available_hms.get(kk));
            }

            HashMap<String,HashMap<String,Double>> results = new HashMap<>();

            for(int i = 0; i < tf_lists.size();i++)
            {
                for(String tf:tf_lists.get(i).keySet())
                {
                    int count_hms = 0;

                    for(int j = 0; j < tf_lists.size(); j++)
                    {
                        HashMap<String,HashMap<String,Double>> list = tf_lists.get(j);
                        if(list.containsKey(tf))
                        {
                            count_hms++;
                        }
                    }

                    if(count_hms>=options_intern.plot_cutoff_hms)
                    {
                        results.put(tf,tf_lists.get(i).get(tf));
                    }
                }
            }


            BufferedWriter bw = new BufferedWriter(new FileWriter(f_out.getAbsolutePath()+File.separator+options_intern.file_suffix_analysis_plot_data_hm_level_different));
            StringBuilder sb = new StringBuilder();
            ArrayList<String> tp_order = new ArrayList<>();
            sb.append("TF");
            for(String tp : tp_tf_gene_counts.keySet())
            {
                sb.append("\t");
                sb.append(tp);
                tp_order.add(tp);
            }
            bw.write(sb.toString());
            bw.newLine();

            for(String tf : results.keySet())
            {
                StringBuilder sb_tf = new StringBuilder();
                sb_tf.append(tf);
                HashMap<String,Double> tf_info = results.get(tf);

                for(String order: tp_order)
                {
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

        File folder_input = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_put_DYNAMITE);
        File folder_output = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_plots);

        File data_output = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_data_plots);
        data_output.mkdir();

        folder_output.mkdir();

        for(File fileDirHM:folder_input.listFiles())
        {
            if(fileDirHM.isDirectory())
            {
                File folder_outputHM = new File(folder_output.getAbsolutePath()+File.separator+fileDirHM.getName());
                folder_outputHM.mkdir();

                File folder_out_data_HM = new File(data_output.getAbsolutePath()+File.separator+fileDirHM.getName());
                folder_out_data_HM.mkdir();

                for(double d: options_intern.plot_th_coefficient)
                {
                    File out_th = new File(folder_outputHM.getAbsolutePath()+File.separator+d);
                    out_th.mkdir();

                    File out_data_th = new File(folder_out_data_HM.getAbsolutePath()+File.separator+d);
                    out_data_th.mkdir();

                    StringBuilder sb = new StringBuilder();
                    sb.append("import pandas as pd\n");
                    sb.append("import seaborn as sns\n");
                    sb.append("import matplotlib.pyplot as plt\n");
                    sb.append("sns.set_context(\"notebook\")\n");
                    sb.append("color = \"#A6CEE3\"\n");
                    sb.append("sns.set_context(\"talk\")\n");
                    sb.append("sns.set_style(\"whitegrid\")\n");

                    HashSet<String> th_group_differentpoints = new HashSet<>();
                    HashSet<String> th_group_samepoints = new HashSet<>();

                    for(File fileDirHM_Group : fileDirHM.listFiles())
                    {

                        String[] group_split = fileDirHM_Group.getName().split("_");
                        String groupname1 = group_split[0];
                        String groupname2 = group_split[1];

                        if(groupname1.charAt(0)!=groupname2.charAt(0))
                        {
                            th_group_differentpoints.add(fileDirHM_Group.getName());
                        }
                        else
                        {
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
                        sb.append(fileDirHM_Group.getAbsolutePath()+File.separator+options_intern.file_suffix_dynamite_output_to_be_plotted);
                        sb.append("').sort_values(['value'], ascending=False)\n");
                        sb.append("# Remove suffix\n");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("['TF'] = ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("['TF'].str.split('_', expand=True)[0]\n");
                        sb.append("# Sort and filter\n");
                        sb.append(fileDirHM_Group.getName()+"_temp1");
                        sb.append(" = ");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("[");
                        sb.append(fileDirHM_Group.getName());
                        sb.append("['value'] > ");
                        sb.append(d);
                        sb.append("].set_index('TF')\n");
                        sb.append(fileDirHM_Group.getName()+"_temp2");
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
                        sb.append(fileDirHM_Group.getName()+"_temp1, ");
                        sb.append(fileDirHM_Group.getName()+"_temp2],axis=0)\n");
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
                        sb.append(out_th.getAbsolutePath()+File.separator+fileDirHM.getName()+"_"+fileDirHM_Group.getName()+"_threshold_"+d+".png");
                        sb.append("\")\n");
                        sb.append("    plt.clf()\n");
                    }
                    sb.append("plt.figure(figsize=(26, 20))\n");
                    sb.append("# Heatmap different stages\n");
                    //Clear Peak values

                    for(String s: th_group_differentpoints)
                    {
                        sb.append("if('Peak' in ");
                        sb.append(s);
                        sb.append(".index):\n");
                        sb.append("    ");
                        sb.append(s);
                        sb.append("=pd.DataFrame.drop(");
                        sb.append(s);
                        sb.append(",index='Peak')\n");
                    }

                    for(String s: th_group_differentpoints)
                    {
                        sb.append(s);
                        sb.append("=");
                        sb.append(s);
                        sb.append("[~");
                        sb.append(s);
                        sb.append(".index.duplicated(keep='first')]\n");
                    }

                    sb.append("join_df_stages = pd.concat([");
                    int c = 0;
                    for(String s : th_group_differentpoints)
                    {
                        if(c==0)
                        {
                            sb.append(s);
                        }
                        else
                        {
                            sb.append(", ");
                            sb.append(s);
                        }
                        c++;
                    }
                    sb.append("], axis=1)\n");

                    sb.append("join_df_stages.to_csv(r'");
                    sb.append(out_data_th.getAbsolutePath()+File.separator);
                    sb.append("all_data_different_stages.csv");
                    sb.append("',index = True, header = True)");
                    sb.append("\n");

                    sb.append("if not ");
                    sb.append("join_df_stages");
                    sb.append(".empty:\n");

                    sb.append("    plot = sns.heatmap(join_df_stages.transpose(), cmap=\"Paired\",  square=True, vmin=1, vmax=1, cbar=False, linewidths=0.5, linecolor='black', xticklabels=True)\n");
                    sb.append("    plt.savefig(\"");
                    sb.append(out_th.getAbsolutePath()+File.separator+fileDirHM.getName()+"_threshold_"+d+"_different_stages.png\")\n");

                    sb.append("    plt.clf()\n");
                    sb.append("plt.figure(figsize=(26, 20))\n");

                    sb.append("# Heatmap same stages\n");
                    for(String s: th_group_samepoints)
                    {
                        sb.append("if('Peak' in ");
                        sb.append(s);
                        sb.append(".index):\n");
                        sb.append("    ");
                        sb.append(s);
                        sb.append("=pd.DataFrame.drop(");
                        sb.append(s);
                        sb.append(",index='Peak')\n");
                    }

                    for(String s: th_group_samepoints)
                    {
                        sb.append(s);
                        sb.append("=");
                        sb.append(s);
                        sb.append("[~");
                        sb.append(s);
                        sb.append(".index.duplicated(keep='first')]\n");
                    }

                    sb.append("join_df_same = pd.concat([");
                    c = 0;
                    for(String s : th_group_samepoints)
                    {
                        if(c==0)
                        {
                            sb.append(s);
                        }
                        else
                        {
                            sb.append(", ");
                            sb.append(s);
                        }
                        c++;
                    }
                    sb.append("], axis=1)\n");

                    sb.append("join_df_same.to_csv(r'");
                    sb.append(out_data_th.getAbsolutePath()+File.separator);
                    sb.append("all_data_same_stages.csv");
                    sb.append("',index = True, header = True)");
                    sb.append("\n");

                    sb.append("if not ");
                    sb.append("join_df_same");
                    sb.append(".empty:\n");

                    sb.append("    plot = sns.heatmap(join_df_same.transpose(), cmap=\"Paired\",  square=True, vmin=1, vmax=1, cbar=False, linewidths=0.5, linecolor='black', xticklabels=True)\n");
                    sb.append("    plt.savefig(\"");
                    sb.append(out_th.getAbsolutePath()+File.separator+fileDirHM.getName()+"_threshold_"+d+"_same_stages.png\")\n");


                    BufferedWriter bw = new BufferedWriter(new FileWriter(new File(out_th.getAbsolutePath()+File.separator+fileDirHM.getName()+"_"+d+".py")));
                    bw.write(sb.toString());
                    bw.close();

                    String command = "python3 " + out_th.getAbsolutePath()+File.separator+fileDirHM.getName()+"_"+d+".py";

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

        String command_base = "Rscript " + options_intern.path_to_COM2POSE+File.separator+options_intern.directory_for_tepic_DYNAMITE+File.separator+"DYNAMITE.R";

        File folder_input;
        folder_input=new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_output_preprocessing_DYNAMITE+File.separator+options_intern.folder_output_preprocessing_DYNAMITE_prepareClass);


        File folder_output = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_out_put_DYNAMITE);
        folder_output.mkdir();

        for(File fileDirHM : folder_input.listFiles())
        {
            if(fileDirHM.isDirectory())
            {
                File folder_outputHM = new File(folder_output.getAbsolutePath()+File.separator+fileDirHM.getName());
                folder_outputHM.mkdir();

                for(File fileDirHM_Group : fileDirHM.listFiles())
                {
                    if(fileDirHM_Group.isDirectory())
                    {
                        File folder_outputHM_Group = new File(folder_outputHM.getAbsolutePath()+File.separator+fileDirHM_Group.getName());
                        folder_outputHM_Group.mkdir();

                        String input_dir = fileDirHM_Group.getAbsolutePath()+File.separator;
                        String output_dir = folder_outputHM_Group.getAbsolutePath()+File.separator;

                        String command_edited = new String(command_base);
                        command_edited += " --dataDir="+input_dir;
                        command_edited += " --outDir="+output_dir;
                        command_edited += " --out_var="+options_intern.dynamite_out_var;
                        command_edited += " --Ofolds="+options_intern.dynamite_Ofolds;
                        command_edited += " --Ifolds="+options_intern.dynamite_Ifolds;
                        command_edited += " --performance="+options_intern.dynamite_performance;
                        command_edited += " --alpha="+options_intern.dynamite_alpha;
                        command_edited += " --cores="+options_intern.dynamite_cores;
                        if(options_intern.dynamite_randomise)
                        {
                            command_edited += " --randomise="+options_intern.dynamite_randomise;
                        }

                        logger.logLine("[DYNAMITE] "+ fileDirHM.getName()+":"+ fileDirHM_Group.getName()+" DYNAMITE.R: " + command_edited);
                        logger.logLine("[DYNAMITE] ... waiting ...");
                        Process child = Runtime.getRuntime().exec(command_edited);
                        int code = child.waitFor();
                        switch (code){
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
     * run integrateData.py and prepareForClassification.R
     */
    public void preprocess_dynamite() throws Exception {

        logger.logLine("[DYNAMITE]: start preprocessing data for DYNAMITE");

        File folder;
        if(options_intern.path_tgen.equals(""))
        {
            folder= new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing+File.separator+options_intern.folder_name_tepic_postprocessing_output);
        }
        else
        {
            folder= new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_integrate);
        }

        File folder_output_pre = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_output_preprocessing_DYNAMITE);
        folder_output_pre.mkdir();

        File folder_output = new File(folder_output_pre.getAbsolutePath()+File.separator+options_intern.folder_output_preprocessing_DYNAMITE_integrateData);
        folder_output.mkdir();

        for(File fileDirHM : folder.listFiles())
        {
            if(fileDirHM.isDirectory())
            {
                File folder_output_hm = new File(folder_output.getAbsolutePath()+File.separator+fileDirHM.getName());
                folder_output_hm.mkdir();

                for(File fileDirGroups : fileDirHM.listFiles())
                {
                    if(fileDirGroups.isDirectory())
                    {
                        File folder_output_group = new File(folder_output_hm.getAbsolutePath()+ File.separator + fileDirGroups.getName());
                        folder_output_group.mkdir();

                        String[] groups_splitted = fileDirGroups.getName().split("_");
                        String groupname1 = groups_splitted[0];
                        String groupname2 = groups_splitted[1];

                        File input_ratios;
                        if(options_intern.path_tgen.equals(""))
                        {
                            input_ratios = new File(fileDirGroups.getAbsolutePath()+File.separator+options_intern.folder_name_tepic_postprocessing_output_ratios+File.separator+options_intern.file_suffix_tepic_postprocessing_output_ratios+fileDirGroups.getName()+".txt");
                        }
                        else
                        {
                            input_ratios = new File(fileDirGroups.getAbsolutePath()+File.separator+options_intern.file_suffix_tepic_postprocessing_output_ratios+groupname1+"_"+groupname2+".txt");
                        }

                        File input_diff_gene_expr = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_output+File.separator+fileDirGroups.getName()+options_intern.file_suffix_deseq2_output_DYNAMITE);

                        String command = "python " + options_intern.path_to_COM2POSE+File.separator+options_intern.directory_for_tepic_DYNAMITE+File.separator+"integrateData.py";
                        command += " "+ input_ratios.getAbsolutePath();
                        command += " " + input_diff_gene_expr.getAbsolutePath();
                        command += " " + folder_output_group.getAbsolutePath()+File.separator+options_intern.file_suffix_output_preprocessing_DYNAMITE_integrateData_log2coeff;
                        if(options_intern.dynamite_preprocessing_integrate_data_geneIDs!=0)
                        {
                            command += " " + options_intern.dynamite_preprocessing_integrate_data_geneIDs;
                        }
                        if(options_intern.dynamite_preprocessing_integrate_data_log2fc!=1)
                        {
                            command += " " + options_intern.dynamite_preprocessing_integrate_data_log2fc;
                        }
                        if(!options_intern.dynamite_preprocessing_integrate_data_consider_geneFile.equals(""))
                        {
                            command += " " + options_intern.dynamite_preprocessing_integrate_data_consider_geneFile;
                        }

                        logger.logLine("[DYNAMITE] " +fileDirHM.getName()+":"+ fileDirGroups.getName()+" preprocessing integrateData.py: " + command);
                        Process child = Runtime.getRuntime().exec(command);
                        int code = child.waitFor();
                        switch (code){
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

        File folder_output_classification = new File(folder_output_pre+File.separator+options_intern.folder_output_preprocessing_DYNAMITE_prepareClass);
        folder_output_classification.mkdir();

        for(File fileDirHM:folder_output.listFiles())
        {
            if(fileDirHM.isDirectory())
            {
                File folder_output_classificationHM = new File(folder_output_classification.getAbsolutePath()+File.separator+fileDirHM.getName());
                folder_output_classificationHM.mkdir();

                for(File fileDirGroup:fileDirHM.listFiles())
                {
                    if(fileDirGroup.isDirectory())
                    {
                        File folder_output_classificationHM_Group= new File(folder_output_classificationHM.getAbsolutePath()+File.separator+fileDirGroup.getName());
                        folder_output_classificationHM_Group.mkdir();

                        String command = "Rscript "+ options_intern.path_to_COM2POSE+File.separator+options_intern.directory_for_tepic_DYNAMITE+File.separator+"prepareForClassificiation.R";
                        command += " " + fileDirGroup.getAbsolutePath()+File.separator+options_intern.file_suffix_output_preprocessing_DYNAMITE_integrateData_log2coeff;
                        command += " " + folder_output_classificationHM_Group.getAbsolutePath()+File.separator+options_intern.file_suffix_output_preprocessing_DYNAMITE_prepClass;
                        logger.logLine("[DYNAMITE] " + fileDirGroup.getName()+" preprocessing prepareForClassification.R: " + command);
                        Process child = Runtime.getRuntime().exec(command);
                        int code = child.waitFor();
                        switch (code){
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
     * integrates TGene and TEPIC data into one affinity file so preproccesing DYNAMITE can use it.
     */
    public void integrate_tgen() throws IOException {
        logger.logLine("[TGENE] Integrate TGene data into TEPIC data");

        File folder_output = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_integrate);
        folder_output.mkdir();
        //create necessary folder structure -> TEMPLATE: postprocess TEPIC
        File folder_tepic_postprocess = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing+File.separator+options_intern.folder_name_tepic_postprocessing_output);
        for(File firDir:folder_tepic_postprocess.listFiles())
        {
            if(firDir.isDirectory())
            {
                File folder_output_hm = new File(folder_output.getAbsolutePath()+File.separator+firDir.getName());
                folder_output_hm.mkdir();
                for(File fileDirHM: firDir.listFiles())
                {
                    if(fileDirHM.isDirectory())
                    {
                        File folder_output_HM_group= new File(folder_output_hm.getAbsolutePath()+File.separator+fileDirHM.getName());
                        folder_output_HM_group.mkdir();
                    }
                }
            }
        }


        //integrate data between TEPIC and TGENE (modifying quotients)
        //based on structure build necessary files
        for(File folderDirHM:folder_output.listFiles())
        {
            if (folderDirHM.isDirectory())
            {
                String hm = folderDirHM.getName();
                for (File folderDirHM_Group : folderDirHM.listFiles())
                {

                    if (folderDirHM_Group.isDirectory())
                    {
                        logger.logLine("[TGENE] integrate " + hm +": " + folderDirHM_Group.getName());

                        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(folderDirHM_Group.getAbsolutePath() + File.separator + options_intern.file_suffix_tepic_postprocessing_output_ratios+folderDirHM_Group.getName()+".txt")));

                        File input_data_TGENE = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_groups+File.separator+hm+File.separator+folderDirHM_Group.getName()+File.separator+options_intern.file_suffic_tgen_output_groups);
                        File input_data_TEPIC = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing+File.separator+options_intern.folder_name_tepic_postprocessing_output+File.separator+hm+File.separator+folderDirHM_Group.getName()+File.separator+options_intern.folder_name_tepic_postprocessing_output_ratios+File.separator+options_intern.file_suffix_tepic_postprocessing_output_ratios+folderDirHM_Group.getName()+".txt");

                        //new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_output_preprocessing_DYNAMITE+File.separator+options_intern.folder_output_preprocessing_DYNAMITE_integrateData+File.separator+hm+File.separator+folderDirHM_Group.getName()+File.separator+options_intern.file_suffix_output_preprocessing_DYNAMITE_integrateData_log2coeff);

                        File input_map_position_ensg = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_preprocessing+File.separator+options_intern.folder_name_tgen_preprocessing_binary_trees+File.separator+options_intern.folder_name_tgen_preprocessing_binary_trees_sorted);

                        HashMap<String,ArrayList<ENSG_ranges_binary_trees>> regions = new HashMap<>();
                        HashMap<String,Boolean> tfs_in_tgene = new HashMap<>();

                        HashMap<String, ENSG_binary_tree> chr_binary_tree = new HashMap<>();

                        for(File chr : input_map_position_ensg.listFiles())
                        {
                            ArrayList<ENSG_ranges_binary_trees> chr_region_list = new ArrayList<>();

                            BufferedReader br = new BufferedReader(new FileReader(chr));
                            String line = br.readLine();
                            while((line=br.readLine())!=null)
                            {
                                String[] split = line.split("\t");

                                ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                                iu.chromosome=chr.getName();
                                iu.number = Integer.parseInt(split[0]);
                                iu.left_border=Integer.parseInt(split[1]);
                                iu.right_border=Integer.parseInt(split[2]);
                                iu.ensgs.addAll(Arrays.asList(split[3].split(";")));

                                chr_region_list.add(iu);

                            }
                            br.close();

                            regions.put(chr.getName(),chr_region_list);
                        }

                        for(String chr : regions.keySet())
                        {
                            ArrayList<ENSG_ranges_binary_trees> current_regions = regions.get(chr);

                            ENSG_binary_tree_node root = new ENSG_binary_tree_node(current_regions.get(0),current_regions.get(0).number);
                            ENSG_binary_tree bin_tree = new ENSG_binary_tree(root);

                            for(int i = 1; i < current_regions.size();i++)
                            {
                                bin_tree.add(current_regions.get(i).number,current_regions.get(i));

                            }

                            chr_binary_tree.put(chr.split("\\.")[0],bin_tree);

                        }

                        BufferedReader br_tgene = new BufferedReader(new FileReader(input_data_TGENE));
                        String line_tgene = br_tgene.readLine();

                        HashMap<String,ENSG_ranges_binary_trees> ensgTargetGenes_TFs = new HashMap<>();

                        while((line_tgene= br_tgene.readLine())!=null)
                        {
                            String split[] = line_tgene.split("\t");

                            String chr = split[0];

                            ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                            iu.chromosome=chr;
                            iu.left_border= Integer.parseInt(split[1]);
                            iu.right_border=Integer.parseInt(split[2]);
                            iu.ensgs.addAll(Arrays.asList(split[3].split(";")));

                            ENSG_binary_tree tree;

                            if(chr.startsWith("M"))
                            {
                                tree = chr_binary_tree.get("chr"+options_intern.tgen_mt_writing);
                            }
                            else
                            {
                                tree = chr_binary_tree.get("chr"+chr);
                            }

                            if(tree == null)
                            {
                                continue;
                            }

                            ENSG_ranges_binary_trees ensg_match = tree.containsNode(iu);

                            if(ensg_match == null)
                            {
                                continue;
                            }

                            for(String k: ensg_match.ensgs)
                            {
                                if(ensgTargetGenes_TFs.containsKey(k))
                                {
                                    ensgTargetGenes_TFs.get(k).ensgs.addAll(iu.ensgs);
                                }
                                else
                                {
                                    ensgTargetGenes_TFs.put(k,iu);
                                }
                            }

                            for(String k: iu.ensgs)
                            {
                                tfs_in_tgene.put(k.toUpperCase(),false);
                            }
                        }

                        BufferedReader br_tepic = new BufferedReader(new FileReader(input_data_TEPIC));
                        String line_tepic = br_tepic.readLine();

                        bw.write(line_tepic);
                        bw.newLine();

                        HashMap<String,Integer> tf_place = new HashMap<>();

                        ArrayList<String> header = new ArrayList<>(Arrays.asList(line_tepic.split("\t".toUpperCase())));
                        for(int i = 0; i < header.size();i++)
                        {
                            String key = header.get(i).split("_")[0].toUpperCase();
                            header.set(i,key);
                            String[] split_doubles = key.split("::");
                            for(int j = 0; j < split_doubles.length; j++)
                            {
                                tf_place.put(split_doubles[j],i);
                            }
                        }

                        while((line_tepic= br_tepic.readLine())!=null)
                        {
                            String[] split = line_tepic.split("\t");

                            StringBuilder sb = new StringBuilder();
                            sb.append(split[0]);

                            if(ensgTargetGenes_TFs.containsKey(split[0]))
                            {
                                HashSet<String> tgene_pref_tfs =  ensgTargetGenes_TFs.get(split[0]).ensgs;
                                HashSet<Integer> pref_positions = new HashSet<>();

                                for(String k: tgene_pref_tfs)
                                {
                                    pref_positions.add(tf_place.get(k));
                                }

                                for(int i = 1; i < split.length;i++)
                                {
                                    if(pref_positions.contains(i))
                                    {
                                        double score = Double.parseDouble(split[i]);
                                        if(options_intern.tgen_consensus_calc.equals("INCREASE_TGENE_TFS"))
                                        {
                                            score = score + (score * (1-options_intern.tgen_consensus));

                                        }
                                        sb.append("\t");
                                        sb.append(score);
                                    }
                                    else
                                    {
                                        double score = Double.parseDouble(split[i]);
                                        if(options_intern.tgen_consensus_calc.equals("DECREASE_NOT_TGENE_TFs"))
                                        {
                                            score *= (1-options_intern.tgen_consensus);
                                        }
                                        sb.append("\t");
                                        sb.append(score);
                                    }
                                }
                            }
                            else
                            {
                                for(int i = 1; i < split.length;i++)
                                {
                                    double score = Double.parseDouble(split[i]);
                                    if(options_intern.tgen_consensus_calc.equals("DECREASE_NOT_TGENE_TFs"))
                                    {
                                        score *= (1-options_intern.tgen_consensus);
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

        logger.logLine("[TGENE] Finished integrating TGene data into TEPIC data");
    }

    /**
     * creates TGENE group clash so it can be merged into TEPIC
     */
    public void create_tgen_groups() throws IOException {
        logger.logLine("[TGENE] Create TGene groupe data");

        File folder_input = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_merged);
        File folder_output = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_groups);
        folder_output.mkdir();

        //create necessary folder structure -> TEMPLATE: postprocess TEPIC
        File folder_tepic_postprocess = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing+File.separator+options_intern.folder_name_tepic_postprocessing_output);
        for(File firDir:folder_tepic_postprocess.listFiles())
        {
            if(firDir.isDirectory())
            {
                File folder_output_hm = new File(folder_output.getAbsolutePath()+File.separator+firDir.getName());
                folder_output_hm.mkdir();
                for(File fileDirHM: firDir.listFiles())
                {
                    if(fileDirHM.isDirectory())
                    {
                        File folder_output_HM_group= new File(folder_output_hm.getAbsolutePath()+File.separator+fileDirHM.getName());
                        folder_output_HM_group.mkdir();
                    }
                }
            }
        }

        //based on structure build necessary files
        for(File folderDirHM:folder_output.listFiles())
        {
            if(folderDirHM.isDirectory())
            {
                String hm = folderDirHM.getName();
                for(File folderDirHM_Group : folderDirHM.listFiles())
                {

                    if(folderDirHM_Group.isDirectory())
                    {
                        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(folderDirHM_Group.getAbsolutePath()+File.separator+options_intern.file_suffic_tgen_output_groups)));
                        bw.write("CHR\tLEFT_BORDER\tRIGHT_BORDER\tENSGS");
                        bw.newLine();

                        String[] split_name = folderDirHM_Group.getName().split("_");
                        String group1 = split_name[0];
                        String group2 = split_name[1];

                        File input_group1 = new File(folder_input.getAbsolutePath()+File.separator+group1+File.separator+hm+File.separator+hm+"_"+group1+".txt");
                        File input_group2 = new File(folder_input.getAbsolutePath()+File.separator+group2+File.separator+hm+File.separator+hm+"_"+group2+".txt");

                        HashMap<String,ArrayList<ENSG_ranges_binary_trees>> chr_unmerged_ius = new HashMap<>();

                        BufferedReader br_group1 = new BufferedReader(new FileReader(input_group1));
                        String line_group1 = br_group1.readLine();
                        while((line_group1=br_group1.readLine())!=null)
                        {
                            String[] split_group1 = line_group1.split("\t");
                            String chr = split_group1[0];
                            int left_border = Integer.parseInt(split_group1[1]);
                            int right_border = Integer.parseInt(split_group1[2]);
                            String ensg = split_group1[3].toUpperCase();

                            ArrayList<ENSG_ranges_binary_trees> current_arr_chr;
                            if(chr_unmerged_ius.containsKey(chr))
                            {
                                current_arr_chr = chr_unmerged_ius.get(chr);
                            }
                            else
                            {
                                current_arr_chr = new ArrayList<>();
                            }

                            ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                            iu.chromosome=chr;
                            iu.ensgs.add(ensg);
                            iu.left_border=left_border;
                            iu.right_border=right_border;

                            current_arr_chr.add(iu);
                            chr_unmerged_ius.put(chr,current_arr_chr);
                        }
                        br_group1.close();

                        BufferedReader br_group2 = new BufferedReader(new FileReader(input_group2));
                        String line_group2 = br_group2.readLine();
                        while((line_group2=br_group2.readLine())!=null)
                        {
                            String[] split_group2 = line_group2.split("\t");
                            String chr = split_group2[0];
                            int left_border = Integer.parseInt(split_group2[1]);
                            int right_border = Integer.parseInt(split_group2[2]);
                            String ensg = split_group2[3].toUpperCase();

                            ArrayList<ENSG_ranges_binary_trees> current_arr_chr;
                            if(chr_unmerged_ius.containsKey(chr))
                            {
                                current_arr_chr = chr_unmerged_ius.get(chr);
                            }
                            else
                            {
                                current_arr_chr = new ArrayList<>();
                            }

                            ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                            iu.chromosome=chr;
                            iu.ensgs.add(ensg);
                            iu.left_border=left_border;
                            iu.right_border=right_border;

                            current_arr_chr.add(iu);
                            chr_unmerged_ius.put(chr,current_arr_chr);
                        }
                        br_group2.close();

                        //sort all arrays
                        for(String s : chr_unmerged_ius.keySet())
                        {
                            ArrayList<ENSG_ranges_binary_trees> x = chr_unmerged_ius.get(s);
                            Collections.sort(x);
                            chr_unmerged_ius.put(s,x);
                        }


                        //remove duplicates or integrate ensgs if in overlap
                        for(String chr:chr_unmerged_ius.keySet())
                        {
                            ArrayList<ENSG_ranges_binary_trees> unmerged = chr_unmerged_ius.get(chr);

                            ArrayList<ENSG_ranges_binary_trees> unmerged_temp = new ArrayList<>();

                            for(int i = 1; i < unmerged.size();i++)
                            {
                                ENSG_ranges_binary_trees iu_before = unmerged.get(i-1);
                                ENSG_ranges_binary_trees iu_after = unmerged.get(i);

                                if(iu_before.isTheSame(iu_after))
                                {
                                    unmerged_temp.add(iu_before);
                                    i++;
                                }
                                else
                                {
                                    unmerged_temp.add(iu_before);
                                }
                            }

                            unmerged = new ArrayList<>(unmerged_temp);
                            unmerged_temp.clear();

                            //check for same_ranges but not same ENSGs and merge them
                            for(int i = 1; i < unmerged.size();i++)
                            {
                                ENSG_ranges_binary_trees iu_before = unmerged.get(i-1);
                                ENSG_ranges_binary_trees iu_after = unmerged.get(i);

                                if(iu_before.isSameRange(iu_after))
                                {
                                    for(String s : iu_after.ensgs)
                                    {
                                        iu_before.ensgs.add(s);
                                    }
                                    unmerged_temp.add(iu_before);

                                    i++;
                                }
                                else
                                {
                                    unmerged_temp.add(iu_before);
                                }
                            }

                            unmerged=new ArrayList<>(unmerged_temp);
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


                            for(int i = 0; i < unmerged.size(); i++)
                            {
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

        File folder_input = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_output);
        File folder_output = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_merged);
        folder_output.mkdir();

        for(File fileDirTP : folder_input.listFiles())
        {
            if(fileDirTP.isDirectory())
            {
                File output_tp = new File(folder_output.getAbsolutePath()+File.separator+fileDirTP.getName());
                output_tp.mkdir();

                for(File fileDirTP_HM : fileDirTP.listFiles())
                {
                    if(fileDirTP_HM.isDirectory())
                    {
                        File output_tp_hm = new File(output_tp.getAbsolutePath()+File.separator+fileDirTP_HM.getName());
                        output_tp_hm.mkdir();

                        String name_output = fileDirTP_HM.getName() + "_"+ fileDirTP.getName();


                        HashMap<String,ArrayList<ENSG_ranges_binary_trees>> tfs_to_regions = new HashMap<>();

                        for(File fileDirTP_HM_Samples:fileDirTP_HM.listFiles())
                        {
                            if(fileDirTP_HM_Samples.isDirectory())
                            {
                                File input_data = new File(fileDirTP_HM_Samples.getAbsolutePath()+File.separator+options_intern.file_suffix_tgen_output);

                                BufferedReader br = new BufferedReader(new FileReader(input_data));
                                String line = br.readLine();
                                while((line=br.readLine())!=null)
                                {
                                    if(line.startsWith("#") || line.equals(""))
                                    {
                                        continue;
                                    }

                                    String[] split = line.split("\t");

                                    ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                                    iu.ensgs.add(split[1]);

                                    String[] split_chr = split[6].split(":");
                                    String[] split_borders = split_chr[1].split("-");

                                    iu.chromosome = split_chr[0];
                                    iu.left_border=Integer.parseInt(split_borders[0]);
                                    iu.right_border=Integer.parseInt(split_borders[1]);

                                    if(tfs_to_regions.containsKey(split_chr[0]))
                                    {
                                        ArrayList<ENSG_ranges_binary_trees> x = tfs_to_regions.get(split_chr[0]);
                                        x.add(iu);
                                        tfs_to_regions.put(split_chr[0],x);

                                    }
                                    else
                                    {
                                        ArrayList<ENSG_ranges_binary_trees> x = new ArrayList<>();
                                        x.add(iu);
                                        tfs_to_regions.put(split_chr[0],x);
                                    }
                                }
                            }
                        }

                        for(String chr: tfs_to_regions.keySet())
                        {
                            Collections.sort(tfs_to_regions.get(chr));
                        }

                        for(String chr: tfs_to_regions.keySet())
                        {
                            ArrayList<ENSG_ranges_binary_trees> tf_reg = tfs_to_regions.get(chr);

                            ArrayList<ENSG_ranges_binary_trees> tf_temp = new ArrayList<>();
                            for(int i = 1; i < tf_reg.size(); i++)
                            {
                                //delete all duplicates
                                ENSG_ranges_binary_trees iu_before = tf_reg.get(i-1);
                                ENSG_ranges_binary_trees iu_now = tf_reg.get(i);

                                if(iu_before.isTheSame(iu_now))
                                {
                                    int temp_i = i+1;
                                    boolean same = true;
                                    int count_same = 0;
                                    while(temp_i < tf_reg.size()&&same)
                                    {
                                        same=false;

                                        ENSG_ranges_binary_trees iu_follow = tf_reg.get(temp_i);

                                        if(iu_follow.isTheSame(iu_before))
                                        {
                                            same=true;
                                            count_same++;
                                        }

                                        temp_i++;
                                    }
                                    i=temp_i+1;


                                    tf_temp.add(iu_before);
                                }
                                else
                                {
                                    tf_temp.add(iu_before);
                                }
                            }
                            tfs_to_regions.put(chr,tf_temp);
                        }

                        //CHECK TPM if TPM filter is set!
                        if(options_intern.tepic_tpm_cutoff>0)
                        {
                                /*command_tail += " -T " + options_intern.tepic_tpm_cutoff;
                                command_tail += " -E " + options_intern.tepic_ensg_symbol;
                                command_tail += " -A " + options_intern.deseq2_input_gene_id;
                                String n_dir = options_intern.com2pose_working_directory+File.separator+ options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_single+File.separator+dirGroup.getName()+options_intern.file_suffix_deseq2_preprocessing_meanCounts;
                                command_tail_sample += " -G " + n_dir;

                                [-G input genes count file, if set, default TPM is 1]\n
                                [-T set (T)ranscripts (P)er (M)illion cutoff must be in float form (e.g. 1.0)]
                                [-A input gene annotation desq2 file, required for TPM filter]
                                [-E input ensg to gene symbol file, required for TPM filter*/

                            logger.logLine("[TGENE] TPM filter is set. Filtering TGENE results: " + fileDirTP_HM.getName() + " - " + fileDirTP.getName()+".");

                            double tpm_cutoff = options_intern.tepic_tpm_cutoff;
                            File f_ensg_map_symbol = new File(options_intern.tepic_ensg_symbol);
                            File f_gene_annot_deseq2 = new File(options_intern.deseq2_input_gene_id);
                            File f_gene_count_group1 = new File(options_intern.com2pose_working_directory+File.separator+ options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_single+File.separator+fileDirTP.getName()+options_intern.file_suffix_deseq2_preprocessing_meanCounts);
                            File f_ref_genome = new File(options_intern.tepic_gene_annot);


                            ArrayList<Integer> gene_counts = new ArrayList<>();
                            BufferedReader br_gene_counts_1 = new BufferedReader(new FileReader(f_gene_count_group1));
                            String line_gene_count_1 = br_gene_counts_1.readLine();
                            while ((line_gene_count_1= br_gene_counts_1.readLine())!=null)
                            {
                                gene_counts.add(Integer.parseInt(line_gene_count_1));
                            }
                            br_gene_counts_1.close();


                            ArrayList<String> ensg_numbers = new ArrayList<>();
                            BufferedReader br_ensg_numbers = new BufferedReader(new FileReader(f_gene_annot_deseq2));
                            String line_ensg_numbers = br_ensg_numbers.readLine();
                            while((line_ensg_numbers= br_ensg_numbers.readLine())!=null)
                            {
                                ensg_numbers.add(line_ensg_numbers.toUpperCase());
                            }
                            br_ensg_numbers.close();

                            int total_number_samples=0;
                            HashMap<String,Integer> ensg_counts = new HashMap<>();
                            for(int i = 0; i < gene_counts.size();i++)
                            {
                                ensg_counts.put(ensg_numbers.get(i).toUpperCase(),gene_counts.get(i));
                                total_number_samples++;
                            }

                            HashMap<String,String> ensg_symb = new HashMap<>();

                            BufferedReader br_ensg_symb = new BufferedReader(new FileReader(f_ensg_map_symbol));
                            String line_ensg_symb = br_ensg_symb.readLine();
                            while((line_ensg_symb=br_ensg_symb.readLine())!=null)
                            {
                                String[] split = line_ensg_symb.split("\t");
                                if(split.length>1)
                                {
                                    ensg_symb.put(split[1].toUpperCase(),split[0]);
                                }
                            }
                            br_ensg_symb.close();

                            HashMap<String,Integer> gene_symbol_lengths = new HashMap<>();
                            BufferedReader br_gene_symbol_lengths = new BufferedReader(new FileReader(f_ref_genome));
                            String line_gene_symbol_lengths = "";
                            while((line_gene_symbol_lengths=br_gene_symbol_lengths.readLine())!=null)
                            {
                                if(line_gene_symbol_lengths.startsWith("#"))
                                {
                                    continue;
                                }
                                String[] split = line_gene_symbol_lengths.split("\t");
                                if(split[2].equals("transcript"))
                                {
                                    int start=Integer.parseInt(split[3]);
                                    int end=Integer.parseInt(split[4]);
                                    int diff=end-start+1;
                                    String[] split_further = split[8].split(" ");
                                    String ensg_name=split_further[1].replace("\"","");
                                    ensg_name=ensg_name.replace(";","");
                                    String[] ensg_name_split=ensg_name.split("\\.");
                                    gene_symbol_lengths.put(ensg_name_split[0],diff);
                                }
                            }
                            br_gene_symbol_lengths.close();

                            double all_rpk = 0.0;

                            for(String ec : ensg_counts.keySet())
                            {
                                if(!ensg_counts.containsKey(ec))
                                {
                                    continue;
                                }

                                int ec_count = ensg_counts.get(ec);
                                int ec_lengths = 0;
                                if(!gene_symbol_lengths.containsKey(ec))
                                {
                                    continue;
                                }
                                ec_lengths=gene_symbol_lengths.get(ec);
                                if(ec_lengths==0)
                                {
                                    continue;
                                }
                                double norm_ec_lengths= ec_lengths/(1000*1.0);
                                double current_rpk = ec_count/(norm_ec_lengths*1.0);
                                all_rpk+=current_rpk;
                            }

                            double 	scaling_factor=all_rpk/(1000000*1.0);

                            for(String chr: tfs_to_regions.keySet())
                            {
                                ArrayList<ENSG_ranges_binary_trees> tf_reg = tfs_to_regions.get(chr);
                                ArrayList<ENSG_ranges_binary_trees> tf_reg_temp = new ArrayList<>();

                                for(ENSG_ranges_binary_trees iu : tf_reg)
                                {
                                    boolean should_write = false;

                                    //check if we want that TF
                                    for(String ec: iu.ensgs)
                                    {
                                        ec = ec.toUpperCase();

                                        if(ensg_symb.containsKey(ec))
                                        {
                                            ec = ensg_symb.get(ec);
                                        }
                                        else
                                        {
                                            continue;
                                        }

                                        if(!ensg_counts.containsKey(ec))
                                        {
                                            continue;
                                        }

                                        int ec_count = ensg_counts.get(ec);
                                        int ec_lengths = 0;
                                        if(!gene_symbol_lengths.containsKey(ec))
                                        {
                                            continue;
                                        }
                                        ec_lengths=gene_symbol_lengths.get(ec);
                                        if(ec_lengths==0)
                                        {
                                            continue;
                                        }
                                        double norm_ec_lengths= ec_lengths/(1000*1.0);
                                        double current_rpk = ec_count/(norm_ec_lengths*1.0);

                                        if(current_rpk>=scaling_factor)
                                        {
                                            should_write=true;
                                        }
                                    }

                                    if(should_write)
                                    {
                                        tf_reg_temp.add(iu);
                                    }
                                }

                                tfs_to_regions.put(chr,tf_reg_temp);
                            }
                        }


                        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(output_tp_hm.getAbsolutePath()+File.separator+name_output+".txt")));
                        bw.write("CHR\tLEFT_BORDER\tRIGHT_BORDER\tENSG");
                        bw.newLine();
                        for(String chr: tfs_to_regions.keySet())
                        {
                            ArrayList<ENSG_ranges_binary_trees> tf_reg = tfs_to_regions.get(chr);
                            for(ENSG_ranges_binary_trees iu : tf_reg)
                            {
                                bw.write(iu.toString());
                                bw.newLine();
                            }

                        }
                        bw.close();
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

        String command_base =options_intern.path_tgen+File.separator +"bin"+File.separator+"tgene";

        File output_tgen = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_output);
        output_tgen.mkdir();

        File folder = new File(options_intern.tepic_input_directory);
        for(File fileDirTP: folder.listFiles())
        {
            if(fileDirTP.isDirectory())
            {
                File output_tgen_tp = new File(output_tgen.getAbsolutePath()+File.separator+fileDirTP.getName());
                output_tgen_tp.mkdir();

                for(File fileDirTP_HM : fileDirTP.listFiles())
                {
                    if(fileDirTP_HM.isDirectory())
                    {
                        File output_tgen_tp_hm = new File(output_tgen_tp.getAbsolutePath()+File.separator+fileDirTP_HM.getName());
                        output_tgen_tp_hm.mkdir();

                        for(File fileDir_data : fileDirTP_HM.listFiles())
                        {
                            if(!fileDir_data.isDirectory())
                            {
                                String[] name_split = fileDir_data.getName().split("\\.");

                                File output_tgen_tp_hm_sample = new File(output_tgen_tp_hm+File.separator+name_split[0]);
                                output_tgen_tp_hm_sample.mkdir();

                                String command_execute = new String(command_base);
                                command_execute += " " + fileDir_data.getAbsolutePath();

                                File gtf_dir = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen+File.separator+options_intern.folder_name_tgen_preprocessing+File.separator+options_intern.folder_name_tgen_preprocessing_gtf);
                                File gtf = new File("");

                                for(File gtf_dirDir : gtf_dir.listFiles())
                                {
                                    if(gtf_dirDir.getName().matches(".*"+options_intern.file_suffix_tgen_preprocess_gtf))
                                    {
                                        gtf=gtf_dirDir;
                                        break;
                                    }
                                }

                                if(gtf.getName().equals(""))
                                {
                                    logger.logLine("[TGENE] could not find preprocessed GTF.");
                                    System.exit(1);
                                }

                                command_execute += " " + gtf.getAbsolutePath();

                                command_execute += " -oc " + output_tgen_tp_hm_sample;

                                if(options_intern.tgen_no_closest_locus)
                                {
                                    command_execute += " --no-closest-locus";
                                }
                                if(options_intern.tgen_no_closest_tss)
                                {
                                    command_execute += " --no-closest-tss";
                                }

                                command_execute+= " --max-link-distances " + options_intern.tgen_max_link_distances;
                                command_execute+= " --max-pvalue " + options_intern.tgen_pvalue;

                                //now execute TGENE:
                                logger.logLine("[TGENE] execute TGENE with command line: " + command_execute);
                                Process child = Runtime.getRuntime().exec(command_execute);
                                int code = child.waitFor();
                                switch (code){
                                    case 0:
                                        break;
                                    case 1:
                                        String message = child.getErrorStream().toString();
                                        logger.logLine("[TGENE] please check chromosome names in peak file and gtf file - they need to be named the same!");
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

        //create necessary folders for preprocessing
        File f_TGEN = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tgen);
        f_TGEN.mkdir();
        File f_TGEN_preprocess = new File(f_TGEN.getAbsolutePath()+File.separator+options_intern.folder_name_tgen_preprocessing);
        f_TGEN_preprocess.mkdir();
        File f_TGEN_preprocess_gtf = new File(f_TGEN_preprocess.getAbsolutePath()+File.separator+options_intern.folder_name_tgen_preprocessing_gtf);
        f_TGEN_preprocess_gtf.mkdir();
        File f_TGEN_preprocess_binary_trees = new File(f_TGEN_preprocess.getAbsolutePath()+File.separator+options_intern.folder_name_tgen_preprocessing_binary_trees);
        f_TGEN_preprocess_binary_trees.mkdir();
        File f_TGEN_preprocess_binary_trees_unmerged = new File(f_TGEN_preprocess_binary_trees.getAbsolutePath()+File.separator+options_intern.folder_name_tgen_preprocessing_binary_trees_unmerged);
        f_TGEN_preprocess_binary_trees_unmerged.mkdir();
        File f_TGEN_preprocess_binary_trees_merged = new File(f_TGEN_preprocess_binary_trees.getAbsolutePath()+File.separator+options_intern.folder_name_tgen_preprocessing_binary_trees_merged);
        f_TGEN_preprocess_binary_trees_merged.mkdir();
        File f_TGEN_preprocess_binary_trees_sorted = new File(f_TGEN_preprocess_binary_trees.getAbsolutePath()+File.separator+options_intern.folder_name_tgen_preprocessing_binary_trees_sorted);
        f_TGEN_preprocess_binary_trees_sorted.mkdir();

        logger.logLine("[TGENE] Preprocessing GTF.");
        //restructure GTF for TGEN -> only transcript data positions are allowed

        String[] gtf_name_split_dir= options_intern.tepic_gene_annot.split(File.separator);
        String[] gtf_name_split = gtf_name_split_dir[gtf_name_split_dir.length-1].split("\\.");
        String gtf_name = "";
        for(int i = 0; i < gtf_name_split.length-1; i++)
        {
            if(i>0)
            {
                gtf_name+=".";
                gtf_name+=gtf_name_split[i];
            }
            else
            {
                gtf_name+=gtf_name_split[i];

            }
        }

        BufferedReader br_gtf_transcripts = new BufferedReader(new FileReader(new File(options_intern.tepic_gene_annot)));
        BufferedWriter bw_gtf_transcripts = new BufferedWriter(new FileWriter(new File(f_TGEN_preprocess_gtf.getAbsolutePath()+File.separator+gtf_name+options_intern.file_suffix_tgen_preprocess_gtf)));

        String line_gtf_transcripts = "";
        while((line_gtf_transcripts= br_gtf_transcripts.readLine())!=null)
        {
            if(line_gtf_transcripts.startsWith("#"))
            {
                continue;
            }
            String[] split = line_gtf_transcripts.split("\t");
            if(split[2].equals("transcript"))
            {
                if(split[0].matches(".*M.*"))
                {
                    String line = options_intern.tgen_mt_writing;
                    for(int i = 1; i < split.length;i++)
                    {
                        line+="\t"+split[i];
                    }
                    bw_gtf_transcripts.write(line);
                    bw_gtf_transcripts.newLine();
                }
                else
                {
                    if(split[0].matches("chr.*"))
                    {
                        bw_gtf_transcripts.write(line_gtf_transcripts.substring(3));
                        bw_gtf_transcripts.newLine();
                    }
                    else
                    {
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
        BufferedWriter bw_gtf_bin_tree=new BufferedWriter(new FileWriter(new File(f_TGEN_preprocess_binary_trees_unmerged.getAbsolutePath()+File.separator+"test"+".txt")));

        String current_chr="";
        int count = 0;

        String line_gtf_bin_tree = "";
        while((line_gtf_bin_tree= br_gtf_bin_tree.readLine())!=null)
        {
            if(line_gtf_bin_tree.startsWith("#"))
            {
                continue;
            }
            String[] split = line_gtf_bin_tree.split("\t");
            if(split[2].equals("transcript"))
            {
                String[] split_line = split[8].split(";");

                String chr = split[0];

                if(!current_chr.equals(chr))
                {
                    bw_gtf_bin_tree.close();
                    count=0;
                    bw_gtf_bin_tree=new BufferedWriter(new FileWriter(new File(f_TGEN_preprocess_binary_trees_unmerged.getAbsolutePath()+File.separator+chr+".txt")));
                    bw_gtf_bin_tree.write("#\tPOS_START\tPOS_END\tENSG\n");
                    current_chr=chr;
                }

                String ensg = "";
                String start_pos = split[3];
                String end_pos = split[4];

                for(int i = 0; i < split_line.length; i++)
                {
                    if(split_line[i].startsWith("gene_id"))
                    {
                        String[] split_x = split_line[i].split("\"");
                        ensg=split_x[1].substring(0,split_x[1].length()).split("\\.")[0];
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

        for(File fileDir: f_TGEN_preprocess_binary_trees_unmerged.listFiles())
        {
            if(!fileDir.isDirectory()&&!fileDir.getName().equals("test.txt"))
            {
                BufferedReader br = new BufferedReader(new FileReader(fileDir));
                String line = br.readLine();
                while((line=br.readLine())!=null)
                {
                    String[] split = line.split("\t");

                    ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                    iu.number=Integer.parseInt(split[0]);
                    iu.left_border=Integer.parseInt(split[1]);
                    iu.right_border=Integer.parseInt(split[2]);
                    iu.ensgs.add(split[3]);

                    unmerged_intervalls.add(iu);
                }
                br.close();

                Collections.sort(unmerged_intervalls);

                ArrayList<ENSG_ranges_binary_trees> unmerged_intervalls_temp_list = new ArrayList<>();
                //merge intervalls with same ENSG list!
                int count_same_ENSG = 0;
                boolean last_one_merged=false;
                for(int i = 1; i < unmerged_intervalls.size();i++)
                {
                    ENSG_ranges_binary_trees iu_before = unmerged_intervalls.get(i-1);
                    ENSG_ranges_binary_trees iu_now = unmerged_intervalls.get(i);

                    HashSet<String> ensgs_before = iu_before.ensgs;
                    HashSet<String> ensgs_now = iu_now.ensgs;

                    boolean all_in = true;

                    for(String s: ensgs_before)
                    {
                        if(!ensgs_now.contains(s))
                        {
                            all_in=false;
                        }
                    }
                    for(String s: ensgs_now)
                    {
                        if(!ensgs_before.contains(s))
                        {
                            all_in=false;
                        }
                    }

                    if(all_in)
                    {
                        //now merge!

                        //look for further intervalls with same ensg
                        boolean found_furthers = true;
                        int counter_further = 0;
                        int temp_i = i;
                        while(found_furthers && temp_i < unmerged_intervalls.size())
                        {
                            found_furthers=false;

                            HashSet<String> ensgs_further = unmerged_intervalls.get(temp_i).ensgs;


                            boolean all_in_further = true;

                            for(String s: ensgs_further)
                            {
                                if(!ensgs_before.contains(s))
                                {
                                    all_in_further=false;
                                }
                            }

                            for(String s: ensgs_before)
                            {
                                if(!ensgs_further.contains(s))
                                {
                                    all_in_further=false;
                                }
                            }

                            if(all_in_further)
                            {
                                found_furthers=true;
                                iu_now=unmerged_intervalls.get(temp_i);
                            }
                            else
                            {
                                found_furthers=false;
                                break;
                            }

                            if(all_in_further&&temp_i==unmerged_intervalls.size()-1)
                            {
                                last_one_merged=true;
                            }


                            counter_further++;
                            temp_i++;
                        }
                        i+=counter_further+1;


                        ENSG_ranges_binary_trees iu = new ENSG_ranges_binary_trees();
                        iu.number=count_same_ENSG;
                        iu.ensgs = new HashSet<>(ensgs_before);
                        iu.left_border=iu_before.left_border;
                        iu.right_border=iu_now.right_border;

                        unmerged_intervalls_temp_list.add(iu);

                        count_same_ENSG++;
                    }
                    else
                    {
                        iu_before.number=count_same_ENSG;
                        unmerged_intervalls_temp_list.add(iu_before);
                        count_same_ENSG++;
                    }
                }

                if(!last_one_merged)
                {
                    ENSG_ranges_binary_trees iu = unmerged_intervalls.get(unmerged_intervalls.size()-1);
                    iu.number=unmerged_intervalls_temp_list.size();
                    unmerged_intervalls_temp_list.add(iu);
                }

                unmerged_intervalls.clear();
                unmerged_intervalls = new ArrayList<>(unmerged_intervalls_temp_list);
                unmerged_intervalls_temp_list.clear();

                //re calculate borders

                unmerged_intervalls.get(0).left_border=0;

                for(int i = 1; i < unmerged_intervalls.size(); i++)
                {
                    ENSG_ranges_binary_trees iu_before = unmerged_intervalls.get(i-1);
                    ENSG_ranges_binary_trees iu_current = unmerged_intervalls.get(i);

                    iu_current.left_border=iu_before.right_border+1;
                }

                Collections.sort(unmerged_intervalls);

                for(int i = 0; i < unmerged_intervalls.size();i++)
                {
                    unmerged_intervalls.get(i).number = i;
                }

                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(f_TGEN_preprocess_binary_trees_merged+File.separator+fileDir.getName())));
                bw.write("#\tPOS_START\tPOS_END\tENSG\n");
                for(int i = 1; i < unmerged_intervalls.size(); i++)
                {
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

        for(File fileDir: f_TGEN_preprocess_binary_trees_merged.listFiles())
        {
            if(!fileDir.isDirectory())
            {

                ArrayList<ENSG_ranges_binary_trees> regions = new ArrayList<>();

                BufferedReader br = new BufferedReader(new FileReader(fileDir));
                String header = br.readLine();
                String line = "";

                while((line= br.readLine())!=null)
                {
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

                ENSG_ranges_binary_trees median = regions.get(regions.size()/2);
                newly_ordered.add(median);

                ArrayList<ENSG_ranges_binary_trees> region_left=new ArrayList<>();
                for(int i = 0; i < regions.size()/2; i++)
                {
                    region_left.add(regions.get(i));
                }

                ArrayList<ENSG_ranges_binary_trees> region_right=new ArrayList<>();
                for(int i = regions.size()/2+1; i < regions.size(); i++)
                {
                    region_right.add(regions.get(i));
                }

                newly_ordered = recursive_split(region_left,newly_ordered);
                newly_ordered = recursive_split(region_right,newly_ordered);

                BufferedWriter bw = new BufferedWriter(new FileWriter(f_TGEN_preprocess_binary_trees_sorted.getAbsolutePath()+File.separator+fileDir.getName()));
                bw.write(header);
                bw.newLine();

                for(int i = 0; i < newly_ordered.size();i++)
                {
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

        HashMap<String,HashMap<String,HashSet<String>>> groups_to_compare = checkGroupsTEPIC();

        String suffix="";

        if(options_intern.tepic_tpm_cutoff>0)
        {
            suffix="_Gene_View_Filtered_TPM.txt";
        }
        else
        {
            suffix="_Gene_View_Filtered.txt";
        }

        File folder_postprocessing = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_postprocessing);
        folder_postprocessing.mkdir();

        File folder_pp_input = new File(folder_postprocessing.getAbsolutePath()+File.separator+options_intern.folder_name_tepic_postprocessing_input);
        folder_pp_input.mkdir();

        File folder_pp_output = new File(folder_postprocessing.getAbsolutePath()+File.separator+options_intern.folder_name_tepic_postprocessing_output);
        folder_pp_output.mkdir();

        HashSet<String> available_hms = new HashSet<>();

        //generate input structure and copy files
        for(String s : groups_to_compare.keySet())
        {
            File f_input_tp_folders = new File(folder_pp_input.getAbsolutePath()+File.separator+s);
            f_input_tp_folders.mkdir();
            HashMap<String,HashSet<String>> hm = groups_to_compare.get(s);
            for(String ss: hm.keySet())
            {
                available_hms.add(ss);
                File f_input_hm_folders = new File(f_input_tp_folders.getAbsolutePath()+File.separator+ss);
                f_input_hm_folders.mkdir();

                //move coresponding sample outputs to these folders
                File f_input_samples = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_output_raw+File.separator+s+File.separator+ss);

                for(File fileDir : f_input_samples.listFiles())
                {
                    if (fileDir.isDirectory()) {
                        for (File fileDir2 : fileDir.listFiles())
                        {
                            if (!fileDir2.isDirectory()) {
                                String name = fileDir2.getName();
                                if (name.matches(".*" + suffix))
                                {
                                    //COPY!!
                                    String command = "cp -u " + fileDir2.getAbsolutePath() + " " + folder_pp_input.getAbsolutePath() + File.separator + s + File.separator + ss;
                                    Process child = Runtime.getRuntime().exec(command);
                                    logger.logLine("[TEPIC] Copy files: " + command);
                                    int code = child.waitFor();
                                    switch (code){
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

        //generate output structure
        //HMs
        HashMap<String,File> pp_output_hms_files = new HashMap<>();
        for(String s: available_hms)
        {
            File f = new File(folder_pp_output.getAbsolutePath()+File.separator+s);
            f.mkdir();
            pp_output_hms_files.put(s,f);
        }

        //identify groups to compute mean ratios
        //group1_vs_group2
        HashMap<String,File> pp_output_clashedGroups = new HashMap<>();
        HashSet<String> already_checked_groups = new HashSet<>();
        for(String s: groups_to_compare.keySet())
        {
            for(String ss: groups_to_compare.keySet())
            {
                if(s.equals(ss))
                {
                    continue;
                }
                String key1 = s+"_"+ss;
                String key2 = ss+"_"+s;
                if(already_checked_groups.contains(key1) || already_checked_groups.contains(key2))
                {
                    continue;
                }
                already_checked_groups.add(key1);

                HashMap<String,HashSet<String>> group1_hms = groups_to_compare.get(s);
                HashMap<String,HashSet<String>> group2_hms = groups_to_compare.get(ss);

                for(String k: group1_hms.keySet())
                {
                    if(group2_hms.containsKey(k))
                    {
                        File f = new File(pp_output_hms_files.get(k).getAbsolutePath()+File.separator+key1);
                        f.mkdir();
                    }
                }
            }
        }

        //prepare command line for all (computeMeanRatioTFAffinities.py
        String command_base = "python " + options_intern.path_to_COM2POSE+File.separator+options_intern.directory_for_tepic_DYNAMITE+File.separator+"computeMeanRatioTFAffinities.py";

        for(File fileDir : folder_pp_output.listFiles())
        {
            if(fileDir.isDirectory())
            {
                for(File fileDir2 : fileDir.listFiles())
                {
                    if(fileDir2.isDirectory())
                    {
                        String current_HM = fileDir.getName();
                        String current_group_clash = fileDir2.getName();
                        String[] group_clash_split = current_group_clash.split("_");

                        String group1 = group_clash_split[0];
                        String group2 = group_clash_split[1];

                        String group1_input_dir ="";
                        String group2_input_dir="";
                        //if TPM filter was used we need to postprocess the TPM files. otherwise it will not work!
                        if(options_intern.tepic_tpm_cutoff>0)
                        {
                            logger.logLine("[TEPIC] TPM filter > 0, start postprocessing of TPM filtered scores");
                            File output_post_group1 = new File(folder_pp_output+File.separator+current_HM+File.separator+current_group_clash+File.separator+group1);
                            output_post_group1.mkdir();
                            File output_post_group2 = new File(folder_pp_output+File.separator+current_HM+File.separator+current_group_clash+File.separator+group2);
                            output_post_group2.mkdir();

                            //build intersect and write new files with filter

                            File folder_group1 = new File(folder_pp_input.getAbsolutePath()+File.separator+group1+File.separator+current_HM);
                            File folder_group2 = new File(folder_pp_input.getAbsolutePath()+File.separator+group2+File.separator+current_HM);

                            File[] samples_group1 = folder_group1.listFiles();
                            File[] samples_group2 = folder_group2.listFiles();

                            ArrayList<Boolean> write_group1 = new ArrayList<>();
                            ArrayList<Boolean> write_group2 = new ArrayList<>();

                            ArrayList<String> header_group1;
                            ArrayList<String> header_group2;
                            HashSet<String> header_group1_set;
                            HashSet<String> header_group2_set;

                            BufferedReader br_group1_check = new BufferedReader(new FileReader(samples_group1[0]));
                            String line_group1_check =br_group1_check.readLine();
                            header_group1=new ArrayList<>(Arrays.asList(line_group1_check.split("\t")));
                            header_group1_set=new HashSet<>(Arrays.asList(line_group1_check.split("\t")));
                            br_group1_check.close();

                            BufferedReader br_group2_check = new BufferedReader(new FileReader(samples_group2[0]));
                            String line_group2_check =br_group2_check.readLine();
                            header_group2=new ArrayList<>(Arrays.asList(line_group2_check.split("\t")));
                            header_group2_set=new HashSet<>(Arrays.asList(line_group2_check.split("\t")));
                            br_group2_check.close();

                            for(int i = 0; i < header_group1.size();i++)
                            {
                                if(header_group2_set.contains(header_group1.get(i)))
                                {
                                    write_group1.add(true);
                                }
                                else
                                {
                                    write_group1.add(false);
                                }
                            }
                            for(int i = 0; i < header_group2.size(); i++)
                            {
                                if(header_group1_set.contains(header_group2.get(i)))
                                {
                                    write_group2.add(true);
                                }
                                else
                                {
                                    write_group2.add(false);
                                }
                            }

                            for(File f:folder_group1.listFiles())
                            {
                                BufferedReader br_group1 = new BufferedReader(new FileReader(f));
                                BufferedWriter bw_group1 = new BufferedWriter(new FileWriter(new File(output_post_group1.getAbsolutePath()+File.separator+f.getName())));

                                String line_group1 = "";
                                while((line_group1=br_group1.readLine())!=null)
                                {
                                    String[] split_line_group1 = line_group1.split("\t");
                                    StringBuilder sb = new StringBuilder();

                                    for(int i = 0; i < split_line_group1.length;i++)
                                    {
                                        if(write_group1.get(i))
                                        {
                                            if(i>0)
                                            {
                                                sb.append("\t");
                                                sb.append(split_line_group1[i]);
                                            }
                                            else
                                            {
                                                sb.append(split_line_group1[i]);
                                            }
                                        }
                                    }
                                    bw_group1.write(sb.toString());
                                    bw_group1.newLine();
                                }
                                bw_group1.close();
                                br_group1.close();
                            }


                            for(File f:folder_group2.listFiles())
                            {
                                BufferedReader br_group2 = new BufferedReader(new FileReader(f));
                                BufferedWriter bw_group2 = new BufferedWriter(new FileWriter(new File(output_post_group2.getAbsolutePath()+File.separator+f.getName())));

                                String line_group2 = "";
                                while((line_group2=br_group2.readLine())!=null)
                                {
                                    String[] split_line_group2 = line_group2.split("\t");
                                    StringBuilder sb = new StringBuilder();

                                    for(int i = 0; i < split_line_group2.length;i++)
                                    {
                                        if(write_group2.get(i))
                                        {
                                            if(i>0)
                                            {
                                                sb.append("\t");
                                                sb.append(split_line_group2[i]);
                                            }
                                            else
                                            {
                                                sb.append(split_line_group2[i]);
                                            }
                                        }
                                    }
                                    bw_group2.write(sb.toString());
                                    bw_group2.newLine();
                                }
                                bw_group2.close();
                                br_group2.close();
                            }
                            logger.logLine("[TEPIC] TPM filter > 0, end postprocessing of TPM filtered scores");

                            //set to postprocessed TMP filtered data
                            group1_input_dir=output_post_group1.getAbsolutePath();
                            group2_input_dir=output_post_group2.getAbsolutePath();
                        }
                        else
                        {
                            group1_input_dir = folder_pp_input.getAbsolutePath()+File.separator+group1+File.separator+current_HM;
                            group2_input_dir = folder_pp_input.getAbsolutePath()+File.separator+group2+File.separator+current_HM;
                        }


                        File output_mean_affinities = new File(folder_pp_output+File.separator+current_HM+File.separator+current_group_clash+File.separator+options_intern.folder_name_tepic_postprocessing_output_mean_affinities);
                        output_mean_affinities.mkdir();
                        File output_ratios = new File(folder_pp_output+File.separator+current_HM+File.separator+current_group_clash+File.separator+options_intern.folder_name_tepic_postprocessing_output_ratios);
                        output_ratios.mkdir();

                        File output_mean_affinities_group1 = new File(output_mean_affinities.getAbsolutePath()+File.separator+options_intern.file_suffix_tepic_postprocessing_output_mean_affinities+group1+".txt");
                        File output_mean_affinities_group2 = new File(output_mean_affinities.getAbsolutePath()+File.separator+options_intern.file_suffix_tepic_postprocessing_output_mean_affinities+group2+".txt");
                        File output_ratios_group1_group2 = new File(output_ratios.getAbsolutePath()+File.separator+options_intern.file_suffix_tepic_postprocessing_output_ratios+group1+"_"+group2+".txt");

                        HashSet<File> files_to_create = new HashSet<>();
                        files_to_create.add(output_mean_affinities_group1);
                        files_to_create.add(output_mean_affinities_group2);
                        files_to_create.add(output_ratios_group1_group2);

                        for(File f: files_to_create)
                        {
                            BufferedWriter bw = new BufferedWriter(new FileWriter(f));
                            bw.write("");
                            bw.close();
                        }


                        String command_edited = new String(command_base);

                        command_edited+= " " + group1_input_dir+File.separator;
                        command_edited+= " " + group2_input_dir+File.separator;
                        command_edited+= " " + output_mean_affinities_group1.getAbsolutePath();
                        command_edited+= " " + output_mean_affinities_group2.getAbsolutePath();
                        command_edited+= " " + output_ratios_group1_group2.getAbsolutePath();

                        String command_tail = "";

                        if(options_intern.tepic_original_decay)
                        {
                            command_tail += " True";
                        }
                        else
                        {
                            command_tail += " False";
                        }
                        if(options_intern.tepic_not_generated)
                        {
                            command_tail += " True";
                        }
                        else
                        {
                            command_tail+= " False";
                        }
                        if(options_intern.tepic_tpm_cutoff>0)
                        {
                            command_tail+= " True";
                        }
                        else
                        {
                            command_tail+= " False";
                        }

                        command_edited += command_tail;

                        logger.logLine("[TEPIC] execute computeMeanRatioTFAffinities.py with command line: " + command_edited);
                        Process child = Runtime.getRuntime().exec(command_edited);
                        int code = child.waitFor();
                        switch (code){
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
        logger.logLine("Finished postprocessing of TEPIC output");
    }

    /**
     * creates the TEPIC command lines and runs all samples of all histone modification and all timepoints
     * @throws Exception
     */
    public void run_tepic() throws Exception {
        logger.logLine("Start TEPIC.sh");

        String command = "bash";
        String tepic_path = " " + options_intern.path_to_COM2POSE+File.separator+options_intern.directory_for_tepic_scripts_code_tepic_sh;
        command += tepic_path;
        command += " -g "+ options_intern.tepic_input_ref_genome;
        command += " -p "+ options_intern.tepic_path_pwms;


        String command_tail = "";
        if(options_intern.tepic_cores>1)
        {
            command_tail += " -c "+ options_intern.tepic_cores;
        }
        if(!options_intern.tepic_bed_chr_sign.equals(""))
        {
            command_tail += " -d "+ options_intern.tepic_bed_chr_sign;
        }
        if(options_intern.tepic_column_bedfile!=-1)
        {
            command_tail += " -n "+options_intern.tepic_column_bedfile;
        }
        if(!options_intern.tepic_gene_annot.equals(""))
        {
            command_tail += " -a " + options_intern.tepic_gene_annot;
        }
        if(options_intern.tepic_window_size!=50000)
        {
            command_tail += " -w " + options_intern.tepic_window_size;

        }
        if(!options_intern.tepic_onlyDNasePeaks.equals(""))
        {
            command_tail += " -f " + options_intern.tepic_onlyDNasePeaks;
        }
        if(options_intern.tepic_exponential_decay)
        {
            command_tail += " -e TRUE";
        }
        if(options_intern.tepic_not_norm_peak_length)
        {
            command_tail += " -l TRUE";
        }
        if(options_intern.tepic_not_generated)
        {
            command_tail += " -u TRUE";
        }
        if(options_intern.tepic_original_decay)
        {
            command_tail += " -x TRUE";
        }
        if(!options_intern.tepic_psems_length.equals(""))
        {
            command_tail += " -m "+ options_intern.tepic_psems_length;
        }
        if(options_intern.tepic_entire_gene_body)
        {
            command_tail += " -y TRUE";
        }
        if(options_intern.tepic_zipped)
        {
            command_tail += " -z TRUE";
        }
        if(!options_intern.tepic_2bit.equals(""))
        {
            command_tail += " -r " + options_intern.tepic_2bit;
        }
        if(options_intern.tepic_pvalue!=0.05)
        {
            command_tail += " -v " + options_intern.tepic_pvalue;

        }
        if(options_intern.tepic_minutes_per_chr!=3)
        {
            command_tail += " -i " + options_intern.tepic_minutes_per_chr;
        }
        if(options_intern.tepic_chr_prefix)
        {
            command_tail += " -j TRUE";
        }
        if(options_intern.tepic_transcript_based)
        {
            command_tail += " -t TRUE";
        }
        if(!options_intern.tepic_loop_list.equals(""))
        {
            command_tail += " -h " + options_intern.tepic_loop_list;
        }
        if(options_intern.tepic_loop_windows!=5000)
        {
            command_tail += " -s " + options_intern.tepic_loop_windows;
        }
        if(options_intern.tepic_only_peak_features)
        {
            command_tail+= " -q TRUE";
        }
        if(options_intern.tepic_tpm_cutoff>0)
        {
            command_tail += " -T " + options_intern.tepic_tpm_cutoff;
            command_tail += " -E " + options_intern.tepic_ensg_symbol;
            command_tail += " -A " + options_intern.deseq2_input_gene_id;
        }

        File output_TEPIC = new File(options_intern.com2pose_working_directory+File.separator+ options_intern.folder_name_tepic_output_raw);
        output_TEPIC.mkdir();


        File folder = new File(options_intern.tepic_input_directory);
        for(File dirGroup :folder.listFiles())
        {
            if(dirGroup.isDirectory())
            {
                File output_TEPIC_group = new File(output_TEPIC.getAbsolutePath()+File.separator+dirGroup.getName());
                output_TEPIC_group.mkdir();

                logger.logLine("[TEPIC] Start group "+ dirGroup.getName());
                for(File dirHM : dirGroup.listFiles())
                {
                    File output_TEPIC_group_hm = new File(output_TEPIC_group.getAbsolutePath()+File.separator+dirHM.getName());
                    output_TEPIC_group_hm.mkdir();

                    logger.logLine("[TEPIC] Start histone modification: "+ dirHM.getName());
                    for(File sample : dirHM.listFiles())
                    {
                        logger.logLine("[TEPIC] Start sample " + sample.getName());

                        File output_sample = new File(output_TEPIC_group_hm.getAbsolutePath()+File.separator+sample.getName());
                        output_sample.mkdir();

                        String command_sample = new String(command);

                        String output_dir_combined = output_sample.getAbsolutePath()+File.separator+output_sample.getName();

                        command_sample += " -b " + sample.getAbsolutePath();
                        command_sample += " -o " + output_dir_combined;

                        String command_tail_sample = new String(command_tail);
                        if(options_intern.tepic_tpm_cutoff>0)
                        {
                            String n_dir = options_intern.com2pose_working_directory+File.separator+ options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_single+File.separator+dirGroup.getName()+options_intern.file_suffix_deseq2_preprocessing_meanCounts;
                            command_tail_sample += " -G " + n_dir;
                        }

                        String command_execute = command_sample + command_tail_sample;
                        logger.logLine("[TEPIC] execute TEPIC with command line: " + command_execute);
                        Process child = Runtime.getRuntime().exec(command_execute);
                        int code = child.waitFor();
                        switch (code){
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
        logger.logLine("Finished TEPIC.sh");
    }

    /**
     *run created DESeq2 scripts and postprocess for input into DYNAMITE
     */
    public void run_and_postprocess_DESeq2() throws Exception {

        logger.logLine("Start running DESeq2 RScripts");

        File folder = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_R_scripts);

        for (File dir:folder.listFiles())
        {
            if(!dir.isDirectory())
            {
                String command = "Rscript " + dir.getAbsolutePath();
                Process child = Runtime.getRuntime().exec(command);
                logger.logLine("[DESEQ2] Running script " + dir.getName()+": " + command);
                int code = child.waitFor();
                switch (code){
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


        File folder_results = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_output_raw);
        File output_file = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_output);
        output_file.mkdir();

        for(File res:folder_results.listFiles())
        {
            if(!res.isDirectory())
            {
                String name_res = res.getName();
                String[] split_name_res = name_res.split("_");
                String name_out = split_name_res[0]+"_"+split_name_res[1]+options_intern.file_suffix_deseq2_output_DYNAMITE;

                BufferedReader br = new BufferedReader(new FileReader(res));
                BufferedWriter bw = new BufferedWriter(new FileWriter(output_file.getAbsolutePath()+File.separator+name_out));

                bw.write("geneID\tlog2fc");
                bw.newLine();

                String line = br.readLine();

                while((line= br.readLine())!=null)
                {
                    String[] split = line.split("\t");
                    if(!split[2].equals("NA"))
                    {
                        StringBuilder sb = new StringBuilder();
                        sb.append(split[0].substring(1,split[0].length()-1));
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
     *create DESeq2 scripts based on input directory for DESeq2 - each group against each group, save intermediate steps and R Scripts
     */
    public void create_DESeq2_scripts() throws IOException {

        logger.logLine("Start preprocessing nfcore RNA-seq for DESeq2 input");
        HashMap<Integer,String> row_ensg_name=new HashMap<>();

        File ensg_names = new File(options_intern.deseq2_input_gene_id);
        BufferedReader br_ensg_per_line= new BufferedReader(new FileReader(ensg_names));
        String line_ensg_per_line="";
        line_ensg_per_line=br_ensg_per_line.readLine();
        int count_ensg_lines = 0;
        while((line_ensg_per_line=br_ensg_per_line.readLine())!=null)
        {
            row_ensg_name.put(count_ensg_lines,line_ensg_per_line);
            count_ensg_lines++;
        }
        br_ensg_per_line.close();

        File output_intermediate_steps = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_preprocessing);
        if(output_intermediate_steps.exists())
        {
            logger.logLine("Working directory was already used - please us another one or empty this one completely");
            //TODO: after debugging use system exit !!!
            //System.exit(1);
        }
        output_intermediate_steps.mkdir();
        File output_inter_steps_combined=new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_combined);
        output_inter_steps_combined.mkdir();
        File output_inter_steps_single=new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_single);
        output_inter_steps_single.mkdir();
        File output_inter_steps_symbols_ensg_mean_counts = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_preprocessing+File.separator+options_intern.folder_name_deseq2_preprocessing_gene_symbols);
        output_inter_steps_symbols_ensg_mean_counts.mkdir();

        File folder = new File(options_intern.deseq2_input_directory);

        HashMap<String,String> ensg_symbol = new HashMap<>();
        BufferedReader br_ensg_symbol = new BufferedReader(new FileReader(new File(options_intern.tepic_ensg_symbol)));
        String line_ensg_symbol = br_ensg_symbol.readLine();
        while((line_ensg_symbol=br_ensg_symbol.readLine())!=null)
        {
            String[] split = line_ensg_symbol.split("\t");
            if(split.length>1)
            {
                ensg_symbol.put(split[0],split[1]);
            }
        }

        br_ensg_symbol.close();

        HashMap<String, HashSet<String>> timepoints_samples= new HashMap<>();

        //CREATE SINGLES
        for(File fileDir:folder.listFiles())
        {
            if(fileDir.isDirectory())
            {
                HashMap<Integer,Integer> mean_line_counts = new HashMap<>();
                int count_samples = 0;
                File group = new File(output_inter_steps_single.getAbsolutePath()+File.separator+fileDir.getName());
                group.mkdir();

                HashSet<String> x = new HashSet<>();
                timepoints_samples.put(group.getName(),x);
                for(File sample : fileDir.listFiles())
                {
                    String name= sample.getName();
                    timepoints_samples.get(group.getName()).add(name);

                    BufferedReader br = new BufferedReader(new FileReader(sample));
                    BufferedWriter bw = new BufferedWriter(new FileWriter(group.getAbsolutePath()+File.separator+name));
                    String line = br.readLine();
                    bw.write(line);
                    bw.newLine();
                    int count=0;
                    while((line=br.readLine())!=null)
                    {
                        bw.write(row_ensg_name.get(count)+"\t"+line);
                        bw.newLine();
                        int count_line = Integer.parseInt(line);
                        if(mean_line_counts.containsKey(count))
                        {
                            int z = mean_line_counts.get(count);
                            z+= count_line;
                            mean_line_counts.put(count,z);
                        }
                        else
                        {
                            mean_line_counts.put(count,count_line);
                        }
                        count++;
                    }
                    if(count!=count_ensg_lines)
                    {
                        logger.logLine("Error in nfcore RNA-seq data: File "+ sample.getName()+ " has not the same number of rows as in File "+ ensg_names.getName());
                        System.exit(1);
                    }
                    br.close();
                    bw.close();

                    count_samples++;
                }

                BufferedWriter bw_means = new BufferedWriter(new FileWriter(new File(output_inter_steps_single.getAbsolutePath()+File.separator+group.getName()+options_intern.file_suffix_deseq2_preprocessing_meanCounts)));
                bw_means.write(group.getName()+"_MEANS");
                bw_means.newLine();
                for(int i = 0; i < count_ensg_lines; i++)
                {
                    int mean_count = mean_line_counts.get(i)/count_samples;
                    bw_means.write(""+mean_count);
                    bw_means.newLine();

                }
                bw_means.close();

                BufferedReader br_ensg = new BufferedReader(new FileReader(new File(options_intern.deseq2_input_gene_id)));
                BufferedWriter bw_symbol = new BufferedWriter(new FileWriter(output_inter_steps_symbols_ensg_mean_counts.getAbsolutePath()+File.separator+group.getName()+".csv"));
                bw_symbol.write("SYMBOL\tENSG\tMEAN_COUNT");
                bw_symbol.newLine();

                String line = br_ensg.readLine();

                int i = 0;
                while((line=br_ensg.readLine())!=null)
                {
                    StringBuilder sb = new StringBuilder();

                    int mean_count = mean_line_counts.get(i)/count_samples;
                    String gene_symbol_name = "NO_SYMBOL";

                    if(ensg_symbol.containsKey(line))
                    {
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
        HashSet<String> already_combined=new HashSet<>();
        for(String k: timepoints_samples.keySet())
        {
            for(String kk: timepoints_samples.keySet())
            {
                String key1=k+"_"+kk;
                String key2=kk+"_"+k;

                if(already_combined.contains(key1)||already_combined.contains(key2)||k.equals(kk))
                {
                    continue;
                }

                already_combined.add(key1);

                HashSet<String> samples_group1 = timepoints_samples.get(k);
                HashSet<String> samples_group2 = timepoints_samples.get(kk);

                File combined_file = new File(output_inter_steps_combined.getAbsolutePath()+File.separator+key1);
                combined_file.mkdir();

                StringBuilder sb_header = new StringBuilder();
                sb_header.append("geneID");

                ArrayList<BufferedReader> bufferedReaders = new ArrayList<>();
                for(String s: samples_group1)
                {
                    sb_header.append("\t");
                    sb_header.append(s);
                    BufferedReader br = new BufferedReader(new FileReader(new File(output_inter_steps_single+File.separator+k+File.separator+s)));
                    br.readLine();
                    bufferedReaders.add(br);
                }

                for(String s: samples_group2)
                {
                    sb_header.append("\t");
                    sb_header.append(s);
                    BufferedReader br = new BufferedReader(new FileReader(new File(output_inter_steps_single+File.separator+kk+File.separator+s)));
                    br.readLine();
                    bufferedReaders.add(br);
                }

                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(output_inter_steps_combined+File.separator+key1+File.separator+key1+".csv")));
                bw.write(sb_header.toString());
                bw.newLine();

                BufferedReader br_ensg = new BufferedReader(new FileReader(new File(options_intern.deseq2_input_gene_id)));
                String line = "";
                br_ensg.readLine();
                while((line=br_ensg.readLine())!=null)
                {
                    StringBuilder sb = new StringBuilder();
                    sb.append(line);
                    for(BufferedReader br : bufferedReaders)
                    {
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
        File r_scripts = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_R_scripts);
        r_scripts.mkdir();
        for(File dir : output_inter_steps_combined.listFiles())
        {
            if(dir.isDirectory())
            {
                String name_rscript = dir.getName()+".R";

                String group1 = dir.getName().split("_")[0];
                String group2 = dir.getName().split("_")[1];


                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(r_scripts.getAbsolutePath()+File.separator+name_rscript)));

                StringBuilder sb = new StringBuilder();
                sb.append("if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n");
                sb.append("  install.packages(\"BiocManager\")\n");
                sb.append("if (!requireNamespace(\"DESeq2\", quietly = TRUE))\n");
                sb.append("  install.packages(\"DESeq2\")\n");
                sb.append("library(DESeq2)\n");
                sb.append("#"+group1+" VS "+group2+"\n");
                logger.logLine("[DESEQ2] Group " + group1 + " VS " + group2);

                BufferedReader br_header = new BufferedReader(new FileReader(dir.getAbsolutePath()+File.separator+dir.getName()+".csv"));
                String line_header = br_header.readLine();
                br_header.close();

                ArrayList<String> groups = new ArrayList<>();
                ArrayList<String> samples = new ArrayList<>();

                String[] split_header= line_header.split("\t");
                for(int i = 1; i < split_header.length;i++)
                {
                    String[] split_samples = split_header[i].split("-");
                    String[] split_groups = split_samples[0].split("_");

                    samples.add(split_samples[0]);
                    groups.add(split_groups[0]);
                }

                sb.append("metadata_df<-data.frame(sample_id = c(");
                int count_s=0;
                for(String s: samples)
                {
                    if(count_s>0)
                    {
                        sb.append(" ,\"");
                    }
                    else
                    {
                        sb.append("\"");
                    }
                    sb.append(s);
                    sb.append("\"");
                    count_s++;
                }
                sb.append("), group = c(");
                count_s=0;
                for(String s: groups)
                {
                    if(count_s>0)
                    {
                        sb.append(" ,\"");
                    }
                    else
                    {
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
                sb.append(File.separator+dir.getName());
                sb.append(".csv\"\n");
                sb.append("count_df = read.csv(count_path, sep = \"\\t\", header = T, row.names = 1)\n");

                sb.append("dds <- DESeqDataSetFromMatrix(countData=count_df, \n");
                sb.append("                              colData=metadata_df, \n");
                sb.append("                              design=~group)\n");

                if(options_intern.deseq2_count_threshold>0)
                {
                    sb.append("threshold = ");
                    sb.append(options_intern.deseq2_count_threshold);
                    sb.append("\n");
                    sb.append("keep <- rowSums(counts(dds)) >= threshold\n");
                    sb.append("dds <- dds[keep,]\n");
                }

                File output_deseq2 = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_deseq2_output_raw);
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

    /**
     * Filter TEPIC input files for blacklisted regions
     */
    public void filter_blacklist() throws IOException {
        File folder_input = new File(options_intern.tepic_input_directory);
        File output_folder = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_blacklisted_regions);
        File output_folder_new_input = new File(output_folder.getAbsolutePath()+File.separator+options_intern.folder_name_blacklisted_regions_new_input);
        output_folder_new_input.mkdir();

        //set new folder directory for tepic input and save old one
        options_intern.tepic_input_prev=options_intern.tepic_input_directory;
        options_intern.tepic_input_directory = output_folder_new_input.getAbsolutePath();

        logger.logLine("[BLACKLIST] Create chromosome binary trees.");
        //CREATE BINARY TREES
        HashMap<String,BL_binary_tree> chr_tree = new HashMap<>();

        File input_folder_chr = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_blacklisted_regions+File.separator+options_intern.folder_name_blacklisted_regions_preprocessing+File.separator+options_intern.folder_name_blacklisted_regions_preprocessing_sorted);
        for(File fileChr : input_folder_chr.listFiles())
        {
            if(fileChr.isFile())
            {
                String name= fileChr.getName().split("\\.")[0];

                ArrayList<BL_ranges_binary_tree> region = new ArrayList<>();

                BufferedReader br_chr = new BufferedReader(new FileReader(fileChr));
                String line_chr = br_chr.readLine();
                while((line_chr=br_chr.readLine())!=null)
                {
                    String[] split = line_chr.split("\t");

                    BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
                    iu.number=Integer.parseInt(split[0]);
                    iu.chr = split[1];
                    iu.left_border = Integer.parseInt(split[2]);
                    iu.right_border = Integer.parseInt(split[3]);
                    iu.signal = split[4];

                    region.add(iu);
                }
                br_chr.close();

                BL_binary_tree_node root = new BL_binary_tree_node(region.get(0),region.get(0).number);
                BL_binary_tree tree = new BL_binary_tree(root);

                for(int i = 1; i < region.size(); i++)
                {
                    tree.add(region.get(i).number,region.get(i));
                }

                chr_tree.put(name,tree);

            }
        }

        logger.logLine("[BLACKLIST] Filter input files for blacklisted regions.");
        //now filter all files
        int count_matched = 0;

        for(File fileDirTP : folder_input.listFiles())
        {
            if(fileDirTP.isDirectory())
            {
                File output_folder_new_input_TP = new File(output_folder_new_input.getAbsolutePath()+File.separator+fileDirTP.getName());
                output_folder_new_input_TP.mkdir();

                for(File fileDirTP_HM : fileDirTP.listFiles())
                {
                    if(fileDirTP_HM.isDirectory())
                    {
                        File output_folder_new_input_TP_HM = new File(output_folder_new_input_TP.getAbsolutePath()+File.separator+fileDirTP_HM.getName());
                        output_folder_new_input_TP_HM.mkdir();

                        for(File filrDirTP_HM_sample : fileDirTP_HM.listFiles())
                        {
                            if(filrDirTP_HM_sample.isFile())
                            {
                                logger.logLine("[BLACKLIST] Filter: " + fileDirTP.getName() + ": " + fileDirTP_HM.getName() + " - " + filrDirTP_HM_sample.getName());
                                //FILTER HERE WITH BINARY TREE

                                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(output_folder_new_input_TP_HM.getAbsolutePath()+File.separator+filrDirTP_HM_sample.getName())));
                                BufferedReader br = new BufferedReader(new FileReader(filrDirTP_HM_sample));
                                String line ="";
                                int count_line = 0;
                                while((line=br.readLine())!=null)
                                {
                                    String[] split = line.split("\t");

                                    String chr = split[0];
                                    if(!chr.matches(".*chr.*"))
                                    {
                                        chr="chr"+chr;
                                    }

                                    BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
                                    iu.chr = chr;
                                    iu.left_border = Integer.parseInt(split[1]);
                                    iu.right_border = Integer.parseInt(split[2]);

                                    if(!chr_tree.containsKey(chr))
                                    {
                                        count_line++;
                                        continue;
                                    }

                                    BL_binary_tree tree = chr_tree.get(chr);
                                    if(tree.containsNode(iu) == null)
                                    {
                                        bw.write(line);
                                        bw.newLine();
                                    }
                                    else
                                    {
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

        File f_blacklist = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_blacklisted_regions);
        f_blacklist.mkdir();
        File f_blacklist_pre = new File(f_blacklist.getAbsolutePath()+File.separator+options_intern.folder_name_blacklisted_regions_preprocessing);
        f_blacklist_pre.mkdir();
        File f_blacklist_pre_chr = new File(f_blacklist_pre.getAbsolutePath()+File.separator+options_intern.folder_name_blacklisted_regions_preprocessing_perChr);
        f_blacklist_pre_chr.mkdir();

        BufferedReader br_chr = new BufferedReader(new FileReader(new File(options_intern.black_list_dir)));
        ArrayList<BL_ranges_binary_tree> chr_ius= new ArrayList<>();

        String line_chr = "";
        String current_chr = "";
        BufferedWriter bw_chr = new BufferedWriter(new FileWriter(new File(f_blacklist_pre_chr.getAbsolutePath()+File.separator+"test.txt")));
        while((line_chr=br_chr.readLine())!=null)
        {
            String[] split = line_chr.split("\t");

            String chr = split[0];
            if(!chr.equals(current_chr))
            {
                Collections.sort(chr_ius);

                int i = 0;

                for(BL_ranges_binary_tree iu : chr_ius)
                {
                    iu.number=i;
                    bw_chr.write(iu.toString());
                    bw_chr.newLine();
                    i++;
                }

                chr_ius.clear();

                if(!chr.matches(".*chr.*"))
                {
                    chr="chr"+chr;
                }

                bw_chr.close();
                bw_chr = new BufferedWriter(new FileWriter(new File(f_blacklist_pre_chr.getAbsolutePath()+File.separator+chr+".txt")));
                bw_chr.write("#\tCHR\tLEFT_BORDER\tRIGHT_BORDER\tSIGNAL");
                bw_chr.newLine();
                current_chr=split[0];
            }

            BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
            iu.chr=chr;
            iu.left_border = Integer.parseInt(split[1]);
            iu.right_border = Integer.parseInt(split[2]);
            iu.signal = split[3].replace(' ', '_').toUpperCase();

            if(options_intern.black_list_signals.contains(iu.signal))
            {
                chr_ius.add(iu);
            }
        }

        Collections.sort(chr_ius);

        int i = 0;

        for(BL_ranges_binary_tree iu : chr_ius)
        {
            iu.number=i;
            bw_chr.write(iu.toString());
            bw_chr.newLine();
            i++;
        }

        bw_chr.close();
        br_chr.close();

        logger.logLine("[BLACKLIST] prepare chromosomes for binary tree");

        File f_blacklist_pre_sorted = new File(f_blacklist_pre.getAbsolutePath()+File.separator+options_intern.folder_name_blacklisted_regions_preprocessing_sorted);
        f_blacklist_pre_sorted.mkdir();

        for(File fileDir : f_blacklist_pre_chr.listFiles())
        {
            if(fileDir.isFile() && !fileDir.getName().equals("test.txt"))
            {
                ArrayList<BL_ranges_binary_tree> current_ius = new ArrayList<>();
                String header;

                BufferedReader br_sort = new BufferedReader(new FileReader(fileDir));
                header= br_sort.readLine();
                String line_sort ="";
                while((line_sort=br_sort.readLine())!=null)
                {
                    String[] split = line_sort.split("\t");

                    BL_ranges_binary_tree iu = new BL_ranges_binary_tree();
                    iu.number = Integer.parseInt(split[0]);
                    iu.chr = split[1];
                    iu.left_border = Integer.parseInt(split[2]);
                    iu.right_border = Integer.parseInt(split[3]);
                    iu.signal=split[4];

                    current_ius.add(iu);
                }
                br_sort.close();

                ArrayList<BL_ranges_binary_tree> newly_ordered = new ArrayList<>();

                recursive_split_BL(current_ius,newly_ordered);

                BufferedWriter bw_sort = new BufferedWriter(new FileWriter(f_blacklist_pre_sorted.getAbsolutePath()+File.separator+fileDir.getName()));
                bw_sort.write(header);
                bw_sort.newLine();

                for(int j = 0; j < newly_ordered.size(); j++)
                {
                    bw_sort.write(newly_ordered.get(j).toString());
                    bw_sort.newLine();
                }
                bw_sort.close();
            }

        }


        logger.logLine("[BLACKLIST] finished preprocessing blacklist");


    }

    /**
     * preprocess mix histones, search for same peaks and use either the union or the intersection of all
     */
    public void mix_option() throws IOException {
        File file_root_input = new File(options_intern.tepic_input_directory);
        File root_mix_working_dir = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_mix_option);
        root_mix_working_dir.mkdir();

        File f_sample_mix_preprocess = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_sample_mix_preprocessing);
        f_sample_mix_preprocess.mkdir();

        File f_sample_mix_output = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_sample_mix);
        f_sample_mix_output.mkdir();

        if(options_intern.mix_level.equals("SAMPLE_LEVEL"))
        {
            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory=f_sample_mix_output.getAbsolutePath();
        }

        logger.logLine("[MIX] Preprocess input data for sample mix - split chromosomes.");

        for(File fileDir:file_root_input.listFiles())
        {
            if(fileDir.isDirectory())
            {
                String timepoint = fileDir.getName();
                File file_output_tp = new File(f_sample_mix_preprocess + File.separator + timepoint);
                file_output_tp.mkdir();

                for(File fileDirHM : fileDir.listFiles())
                {
                    if (fileDirHM.isDirectory())
                    {
                        File file_output_tp_hm = new File(file_output_tp.getAbsolutePath() + File.separator + fileDirHM.getName());
                        file_output_tp_hm.mkdir();

                        for(File fileDirHM_sample : fileDirHM.listFiles())
                        {
                            if(fileDirHM_sample.isFile())
                            {
                                BufferedWriter bw = new BufferedWriter(new FileWriter(new File(file_output_tp_hm.getAbsolutePath()+File.separator+"test.txt")));
                                BufferedReader br = new BufferedReader(new FileReader(fileDirHM_sample));
                                String line = "";
                                String currentChr ="";
                                while((line=br.readLine())!=null)
                                {
                                    String[] split = line.split("\t");
                                    if(!split[0].equals(currentChr))
                                    {
                                        bw.close();
                                        currentChr=split[0];
                                        File f_output_chr = new File(file_output_tp_hm.getAbsolutePath()+File.separator+currentChr);
                                        f_output_chr.mkdir();

                                        bw = new BufferedWriter(new FileWriter(new File(f_output_chr.getAbsolutePath()+File.separator+fileDirHM_sample.getName())));
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

        logger.logLine("[MIX] Create "+options_intern.mix_option+" of samples");

        for(File fileDir: f_sample_mix_preprocess.listFiles())
        {
            if(fileDir.isDirectory())
            {
                String timepoint = fileDir.getName();
                File f_output_union_samples_tp = new File(f_sample_mix_output.getAbsolutePath()+File.separator+timepoint);
                f_output_union_samples_tp.mkdir();

                for(File fileDirHM : fileDir.listFiles())
                {
                    if(fileDirHM.isDirectory())
                    {
                        String hm = fileDirHM.getName();
                        File f_output_union_samples_tp_hm = new File(f_output_union_samples_tp.getAbsolutePath()+File.separator+hm);
                        f_output_union_samples_tp_hm.mkdir();

                        String file_ending = "";

                        HashMap<String,ArrayList<MIX_Interval>> all_chromosomes = new HashMap();
                        ArrayList<Integer> chr_alpha = new ArrayList<>();
                        ArrayList<String> chr_str = new ArrayList<>();

                        for(File fileDirHM_Chr : fileDirHM.listFiles())
                        {
                            if(fileDirHM_Chr.isDirectory())
                            {
                                String chr = fileDirHM_Chr.getName();

                                ArrayList<MIX_Interval> all_intervals = new ArrayList<>();

                                int sample_number = 0;

                                for(File fileDirHM_Chr_sample : fileDirHM_Chr.listFiles())
                                {
                                    if(fileDirHM_Chr_sample.isFile())
                                    {
                                        String[] f_ending = fileDirHM_Chr_sample.getName().split("\\.");
                                        file_ending=f_ending[f_ending.length-1];
                                        //build array of a union of all peaks
                                        BufferedReader br = new BufferedReader(new FileReader(fileDirHM_Chr_sample));
                                        String line = "";
                                        while((line=br.readLine())!=null)
                                        {
                                            String[] split = line.split("\t");

                                            MIX_Interval_Object mio = new MIX_Interval_Object(split[0],Integer.parseInt(split[1]),Integer.parseInt(split[2]),split[3],Integer.parseInt(split[4]),split[5],Double.parseDouble(split[6]),Double.parseDouble(split[7]),Double.parseDouble(split[8]));
                                            MIX_Interval mi = new MIX_Interval(Integer.parseInt(split[1]),Integer.parseInt(split[2]));
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
                                if(options_intern.mix_occurence_intersection==-1)
                                {
                                    min_occurence=sample_number;
                                }
                                else
                                {
                                    if(sample_number<options_intern.mix_occurence_intersection)
                                    {
                                        min_occurence=sample_number;

                                    }
                                    else
                                    {
                                        min_occurence=options_intern.mix_occurence_intersection;
                                    }
                                }

                                while(!stack_union.isEmpty())
                                {
                                    MIX_Interval t = stack_union.pop();
                                    t.calculate_mean("SAMPLE_LEVEL");

                                    if(options_intern.mix_option.equals("INTERSECTION"))
                                    {
                                        if(t.merged_intervals.size()>=min_occurence)
                                        {
                                            chr_unions.add(t);
                                        }
                                    }
                                    else
                                    {
                                        chr_unions.add(t);
                                    }
                                }

                                Collections.sort(chr_unions);

                                all_chromosomes.put(chr,chr_unions);

                                try {
                                    chr_alpha.add(Integer.parseInt(chr));
                                }
                                catch (Exception e)
                                {
                                    chr_str.add(chr);
                                }
                            }
                        }

                        Collections.sort(chr_alpha);
                        Collections.sort(chr_str);


                        //print Sample_unions => can be used if not HM_Level is used
                        BufferedWriter bw = new BufferedWriter(new FileWriter(f_output_union_samples_tp_hm.getAbsolutePath()+File.separator+timepoint+"_"+hm+"."+file_ending));

                        int peak_counter = 1;

                        for(int chr : chr_alpha)
                        {
                            ArrayList<MIX_Interval> x = all_chromosomes.get(""+chr);
                            for(int i = 0; i < x.size();i++)
                            {
                                bw.write(x.get(i).meanToString(peak_counter));
                                bw.newLine();
                                peak_counter++;
                            }
                        }

                        for(String chr: chr_str)
                        {
                            ArrayList<MIX_Interval> x = all_chromosomes.get(chr);
                            for(int i = 0; i < x.size();i++)
                            {
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

        if(options_intern.mix_level.equals("HM_LEVEL"))
        {

            logger.logLine("[MIX] Preprocess sample unions for HM mix");

            File f_output_preprocessing_hm = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_preprocess_hm_mix);
            f_output_preprocessing_hm.mkdir();

            File f_output_hm = new File(root_mix_working_dir.getAbsolutePath()+File.separator+options_intern.folder_name_mix_option_hm_mix);
            f_output_hm.mkdir();

            options_intern.tepic_input_prev = options_intern.tepic_input_directory;
            options_intern.tepic_input_directory=f_output_hm.getAbsolutePath();

            logger.logLine("[MIX] Identify possible Timepoints with same Histone Modifications");

            HashMap<String,ArrayList<String>> timepoints_histone_modifications = new HashMap<>();
            HashMap<String,ArrayList<String>> deleted_tps = new HashMap<>();
            HashSet<String> available_hms = new HashSet<>();

            //identify possible timepoints
            for(File fileDir:f_sample_mix_output.listFiles())
            {
                if(fileDir.isDirectory())
                {
                    String timepoint = fileDir.getName();
                    ArrayList hm = new ArrayList();

                    for(File fileDirHM : fileDir.listFiles())
                    {
                        if(fileDirHM.isDirectory())
                        {
                            hm.add(fileDirHM.getName());
                            available_hms.add(fileDirHM.getName());
                        }
                    }
                    timepoints_histone_modifications.put(timepoint,hm);
                }
            }


            for(String tp : timepoints_histone_modifications.keySet())
            {
                if(timepoints_histone_modifications.get(tp).size()< available_hms.size())
                {
                    deleted_tps.put(tp,timepoints_histone_modifications.get(tp));
                    timepoints_histone_modifications.remove(tp);
                    continue;
                }

                boolean all_in = true;
                for(String hm : available_hms)
                {
                    if(!timepoints_histone_modifications.get(tp).contains(hm))
                    {
                        all_in=false;
                    }
                }

                if(!all_in)
                {
                    deleted_tps.put(tp,timepoints_histone_modifications.get(tp));
                    timepoints_histone_modifications.remove(tp);
                }
            }

            StringBuilder sb_found_mixing_tps = new StringBuilder();
            sb_found_mixing_tps.append("[MIX] Can perform complete mix for HMs (");
            for(String hm : available_hms)
            {
                sb_found_mixing_tps.append(hm);
                sb_found_mixing_tps.append(" ");
            }
            sb_found_mixing_tps.append(") in timepoints (");
            for(String tp:timepoints_histone_modifications.keySet())
            {
                sb_found_mixing_tps.append(tp);
                sb_found_mixing_tps.append(" ");
            }
            sb_found_mixing_tps.append("). Can perform part-mix or no-mix for timepoints (");
            for(String tp : deleted_tps.keySet())
            {
                sb_found_mixing_tps.append(tp);
                sb_found_mixing_tps.append(" ");
            }
            sb_found_mixing_tps.append(").");
            logger.logLine(sb_found_mixing_tps.toString());


            for(File fileDir : f_sample_mix_output.listFiles())
            {
                if(fileDir.isDirectory())
                {
                    String timepoint = fileDir.getName();
                    File f_output_hm_prepro = new File(f_output_preprocessing_hm.getAbsolutePath()+File.separator+timepoint);
                    f_output_hm_prepro.mkdir();

                    for(File fileDirHM : fileDir.listFiles())
                    {
                        if(fileDirHM.isDirectory())
                        {
                            String hm = fileDirHM.getName();
                            File f_output_hm_prepro_hm = new File(f_output_hm_prepro.getAbsolutePath()+File.separator+"MIX");
                            f_output_hm_prepro_hm.mkdir();

                            for(File fileDirHM_samples : fileDirHM.listFiles())
                            {
                                if(fileDirHM_samples.isFile())
                                {
                                    BufferedWriter bw = new BufferedWriter(new FileWriter(new File(f_output_hm_prepro_hm.getAbsolutePath()+File.separator+"test.txt")));
                                    BufferedReader br = new BufferedReader(new FileReader(fileDirHM_samples));
                                    String line = "";
                                    String currentChr ="";
                                    while((line=br.readLine())!=null)
                                    {
                                        String[] split = line.split("\t");
                                        if(!split[0].equals(currentChr))
                                        {
                                            bw.close();
                                            currentChr=split[0];
                                            File f_output_chr = new File(f_output_hm_prepro_hm.getAbsolutePath()+File.separator+currentChr);
                                            f_output_chr.mkdir();

                                            bw = new BufferedWriter(new FileWriter(new File(f_output_chr.getAbsolutePath()+File.separator+fileDirHM_samples.getName())));
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

            logger.logLine("[MIX] Create "+options_intern.mix_option+" of HMs");

            for(File fileDir: f_output_preprocessing_hm.listFiles())
            {
                if(fileDir.isDirectory())
                {
                    String timepoint = fileDir.getName();
                    File f_output_union_samples_tp = new File(f_output_hm.getAbsolutePath()+File.separator+timepoint);
                    f_output_union_samples_tp.mkdir();

                    for(File fileDirHM : fileDir.listFiles())
                    {
                        if(fileDirHM.isDirectory())
                        {
                            String hm = fileDirHM.getName();
                            File f_output_union_samples_tp_hm = new File(f_output_union_samples_tp.getAbsolutePath()+File.separator+hm);
                            f_output_union_samples_tp_hm.mkdir();

                            String file_ending = "";

                            HashMap<String,ArrayList<MIX_Interval>> all_chromosomes = new HashMap();
                            ArrayList<Integer> chr_alpha = new ArrayList<>();
                            ArrayList<String> chr_str = new ArrayList<>();

                            for(File fileDirHM_Chr : fileDirHM.listFiles())
                            {
                                if(fileDirHM_Chr.isDirectory())
                                {
                                    String chr = fileDirHM_Chr.getName();

                                    ArrayList<MIX_Interval> all_intervals = new ArrayList<>();

                                    int sample_number = 0;

                                    for(File fileDirHM_Chr_sample : fileDirHM_Chr.listFiles())
                                    {
                                        if(fileDirHM_Chr_sample.isFile())
                                        {
                                            String[] f_ending = fileDirHM_Chr_sample.getName().split("\\.");
                                            file_ending=f_ending[f_ending.length-1];
                                            //build array of a union of all peaks
                                            BufferedReader br = new BufferedReader(new FileReader(fileDirHM_Chr_sample));
                                            String line = "";
                                            while((line=br.readLine())!=null)
                                            {
                                                String[] split = line.split("\t");

                                                MIX_Interval_Object mio = new MIX_Interval_Object(split[0],Integer.parseInt(split[1]),Integer.parseInt(split[2]),split[3],Integer.parseInt(split[4]),split[5],Double.parseDouble(split[6]),Double.parseDouble(split[7]),Double.parseDouble(split[8]));
                                                MIX_Interval mi = new MIX_Interval(Integer.parseInt(split[1]),Integer.parseInt(split[2]));
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
                                    if(options_intern.mix_occurence_intersection==-1)
                                    {
                                        min_occurence=sample_number;
                                    }
                                    else
                                    {
                                        if(sample_number<options_intern.mix_occurence_intersection)
                                        {
                                            min_occurence=sample_number;

                                        }
                                        else
                                        {
                                            min_occurence=options_intern.mix_occurence_intersection;
                                        }
                                    }

                                    while(!stack_union.isEmpty())
                                    {
                                        MIX_Interval t = stack_union.pop();
                                        t.calculate_mean("HM_LEVEL");

                                        if(options_intern.mix_option.equals("INTERSECTION"))
                                        {
                                            if(t.merged_intervals.size()>=min_occurence)
                                            {
                                                chr_unions.add(t);
                                            }
                                        }
                                        else
                                        {
                                            chr_unions.add(t);
                                        }
                                    }

                                    Collections.sort(chr_unions);

                                    all_chromosomes.put(chr,chr_unions);

                                    try {
                                        chr_alpha.add(Integer.parseInt(chr));
                                    }
                                    catch (Exception e)
                                    {
                                        chr_str.add(chr);
                                    }
                                }
                            }

                            Collections.sort(chr_alpha);
                            Collections.sort(chr_str);


                            //print Sample_unions => can be used if not HM_Level is used
                            BufferedWriter bw = new BufferedWriter(new FileWriter(f_output_union_samples_tp_hm.getAbsolutePath()+File.separator+timepoint+"_"+hm+"."+file_ending));

                            int peak_counter = 1;

                            for(int chr : chr_alpha)
                            {
                                ArrayList<MIX_Interval> x = all_chromosomes.get(""+chr);
                                for(int i = 0; i < x.size();i++)
                                {
                                    bw.write(x.get(i).meanToString(peak_counter));
                                    bw.newLine();
                                    peak_counter++;
                                }
                            }

                            for(String chr: chr_str)
                            {
                                ArrayList<MIX_Interval> x = all_chromosomes.get(chr);
                                for(int i = 0; i < x.size();i++)
                                {
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
     * read config file
     * @param check_options should options be checked for validity? Should be true if pipeline is run, should be false if analyse programms are run
     */
    public void read_config_file(boolean check_options) throws IOException {

        logger.logLine("Start reading config file at "+ options_intern.config_data_path);

        BufferedReader br = new BufferedReader(new FileReader(new File(options_intern.config_data_path)));
        String line= "";
        while((line=br.readLine())!=null)
        {
            if(line.startsWith("#"))
            {
                continue;
            }

            String[] split=line.split("=");

            switch (split[0])
            {
                case "mix_level":
                    options_intern.mix_level=split[1].substring(1,split[1].length()-1);
                    break;
                case "mix_option":
                    options_intern.mix_option=split[1].substring(1,split[1].length()-1);
                    break;
                case "mix_occurence_intersection":
                    options_intern.mix_occurence_intersection=Integer.parseInt(split[1]);
                    break;
                case "black_list_dir":
                    options_intern.black_list_dir=split[1].substring(1,split[1].length()-1);
                    break;
                case "black_list_signals":
                    options_intern.black_list_signals=new HashSet<>(Arrays.asList(split[1].substring(1,split[1].length()-1).toUpperCase().split(";")));
                    break;
                case "deseq2_input_directory":
                    options_intern.deseq2_input_directory=split[1].substring(1,split[1].length()-1);
                    break;
                case "deseq2_input_gene_id":
                    options_intern.deseq2_input_gene_id=split[1].substring(1,split[1].length()-1);
                    break;
                case "deseq2_count_threshold":
                    options_intern.deseq2_count_threshold=Integer.parseInt(split[1]);
                    break;
                case "tepic_input_directory":
                    options_intern.tepic_input_directory=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_input_ref_genome":
                    options_intern.tepic_input_ref_genome=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_path_pwms":
                    options_intern.tepic_path_pwms=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_cores":
                    options_intern.tepic_cores=Integer.parseInt(split[1]);
                    break;
                case "tepic_bed_chr_sign":
                    options_intern.tepic_bed_chr_sign=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_column_bedfile":
                    options_intern.tepic_column_bedfile=Integer.parseInt(split[1]);
                    break;
                case "tepic_gene_annot":
                    options_intern.tepic_gene_annot=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_window_size":
                    options_intern.tepic_window_size=Integer.parseInt(split[1]);
                    break;
                case "tepic_onlyDNasePeaks":
                    options_intern.tepic_onlyDNasePeaks=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_exponential_decay":
                    options_intern.tepic_exponential_decay=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_not_norm_peak_length":
                    options_intern.tepic_not_norm_peak_length=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_not_generated":
                    options_intern.tepic_not_generated=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_original_decay":
                    options_intern.tepic_original_decay=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_psems_length":
                    options_intern.tepic_psems_length=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_entire_gene_body":
                    options_intern.tepic_entire_gene_body=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_zipped":
                    options_intern.tepic_zipped=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_background_seq":
                    options_intern.tepic_background_seq=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_2bit":
                    options_intern.tepic_2bit=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_pvalue":
                    options_intern.tepic_pvalue=Double.parseDouble(split[1]);
                    break;
                case "tepic_minutes_per_chr":
                    options_intern.tepic_minutes_per_chr=Integer.parseInt(split[1]);
                    break;
                case "tepic_chr_prefix":
                    options_intern.tepic_chr_prefix=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_transcript_based":
                    options_intern.tepic_transcript_based=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_loop_list":
                    options_intern.tepic_loop_list=split[1].substring(1,split[1].length()-1);
                    break;
                case "tepic_loop_windows":
                    options_intern.tepic_loop_windows=Integer.parseInt(split[1]);
                    break;
                case "tepic_only_peak_features":
                    options_intern.tepic_only_peak_features=Boolean.parseBoolean(split[1]);
                    break;
                case "tepic_tpm_cutoff":
                    options_intern.tepic_tpm_cutoff=Double.parseDouble(split[1]);
                    break;
                case "tepic_ensg_symbol":
                    options_intern.tepic_ensg_symbol=split[1].substring(1,split[1].length()-1);
                    break;
                case "tgen_consensus":
                    options_intern.tgen_consensus=Double.parseDouble(split[1]);
                    break;
                case "tgen_consensus_calc":
                    options_intern.tgen_consensus_calc=split[1].substring(1,split[1].length()-1);
                    break;
                case "tgen_no_closest_locus":
                    options_intern.tgen_no_closest_locus=Boolean.parseBoolean(split[1]);
                    break;
                case "tgen_no_closest_tss":
                    options_intern.tgen_no_closest_tss=Boolean.parseBoolean(split[1]);
                    break;
                case "tgen_max_link_distances":
                    options_intern.tgen_max_link_distances=Integer.parseInt(split[1]);
                    break;
                case "tgen_pvalue":
                    options_intern.tgen_pvalue=Double.parseDouble(split[1]);
                    break;
                case "tgen_mt_writing":
                    options_intern.tgen_mt_writing=split[1].substring(1,split[1].length()-1);
                    break;
                case "dynamite_preprocessing_integrate_data_geneIDs":
                    options_intern.dynamite_preprocessing_integrate_data_geneIDs=Integer.parseInt(split[1]);
                    break;
                case "dynamite_preprocessing_integrate_data_log2fc":
                    options_intern.dynamite_preprocessing_integrate_data_log2fc=Integer.parseInt(split[1]);
                    break;
                case "dynamite_preprocessing_integrate_data_consider_geneFile":
                    options_intern.dynamite_preprocessing_integrate_data_consider_geneFile=split[1].substring(1,split[1].length()-1);
                    break;
                case "dynamite_out_var":
                    options_intern.dynamite_out_var=split[1].substring(1,split[1].length()-1);
                    break;
                case "dynamite_cores":
                    options_intern.dynamite_cores=Integer.parseInt(split[1]);
                    break;
                case "dynamite_alpha":
                    options_intern.dynamite_alpha=Double.parseDouble(split[1]);
                    break;
                case "dynamite_testsize":
                    options_intern.dynamite_testsize=Double.parseDouble(split[1]);
                    break;
                case "dynamite_Ofolds":
                    options_intern.dynamite_Ofolds=Integer.parseInt(split[1]);
                    break;
                case "dynamite_Ifolds":
                    options_intern.dynamite_Ifolds=Integer.parseInt(split[1]);
                    break;
                case "dynamite_balanced":
                    options_intern.dynamite_balanced=Boolean.parseBoolean(split[1]);
                    break;
                case "dynamite_performance":
                    options_intern.dynamite_performance=Boolean.parseBoolean(split[1]);
                    break;
                case "dynamite_randomise":
                    options_intern.dynamite_randomise=Boolean.parseBoolean(split[1]);
                    break;
                case "plot_th_coefficient":
                    String[] split_coefficient_ths = split[1].split(";");
                    options_intern.plot_th_coefficient.clear();
                    for(String s: split_coefficient_ths)
                    {
                        options_intern.plot_th_coefficient.add(Double.parseDouble(s));
                    }
                    break;
                case "plot_cutoff_tps":
                    options_intern.plot_cutoff_tps=Integer.parseInt(split[1]);
                    break;
                case "plot_cutoff_hms":
                    options_intern.plot_cutoff_hms=Integer.parseInt(split[1]);
                    break;
                case "plot_cutoff_gcs":
                    options_intern.plot_cutoff_gcs=Integer.parseInt(split[1]);
                    break;
                default:
                    logger.logLine("Misformed cfg file - please use template of: /COM2POSE/config_templates/com2pose_template.cfg");
                    logger.logLine("Do not delete unused parameters in config data!");
                    System.exit(1);
            }
        }
        br.close();

        if(check_options)
        {
            boolean all_set = checkOptions();
            logger.logLine("Check config file parameters for validity");
            if(!all_set)
            {
                logger.logLine("Not all [REQ]uired options set. Please set them in config file");
                logger.logLine("Aborting COM2POSE");
                System.exit(1);
            }
            else
            {
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

        if(!options_intern.mix_level.equals(""))
        {
            if(!options_intern.mix_level.equals("HM_LEVEL")&&!options_intern.mix_level.equals("SAMPLE_LEVEL"))
            {
                logger.logLine("[MIX]: mix_level parameter must be either HM_LEVEL or SAMPLE_LEVEL");
                all_set=false;
            }
            if(!options_intern.mix_option.equals("UNION")&&!options_intern.mix_option.equals("INTERSECTION"))
            {
                logger.logLine("[MIX]: mix_option parameter must be either UNION or INTERSECTION");
                all_set=false;
            }
        }

        /**
         * black list options
         */

        if(!options_intern.black_list_dir.equals(""))
        {
            File f = new File(options_intern.black_list_dir);
            if(!f.exists())
            {
                logger.logLine("[BLACKLIST] blacklist path does not exist!");
                all_set=false;
            }
            if(f.isDirectory())
            {
                logger.logLine("[BLACKLIST] blacklist path is a directory but must be a file");
                all_set=false;
            }

            if(options_intern.black_list_signals.isEmpty())
            {
                logger.logLine("[BLACKLIST] signals must not be empty");
                all_set=false;
            }
        }

        /**
         * DESEQ2 options
         */

        if(options_intern.deseq2_input_directory.equals(""))
        {
            logger.logLine("[DESEQ2] input directory is not given");
            all_set=false;
        }
        else
        {
            File f = new File(options_intern.deseq2_input_directory);
            if(!f.exists())
            {
                logger.logLine("[DESEQ2] input path does not exist!");
                all_set=false;
            }
        }

        if(options_intern.deseq2_input_gene_id.equals(""))
        {
            logger.logLine("[DESEQ2] gene ID file from nfcore RNA-seq is not given");
            all_set=false;
        }
        else
        {
            File f = new File(options_intern.deseq2_input_gene_id);
            if(!f.exists())
            {
                logger.logLine("[DESEQ2] gene ID file from nfcore RNA-seq path does not exist!");
                all_set=false;
            }
        }


        /**
         * TEPIC options
         */
        if(options_intern.tepic_input_ref_genome.equals(""))
        {
            logger.logLine("[TEPIC] Reference genome path is not given");
            all_set=false;
        }
        else
        {
            File f = new File(options_intern.tepic_input_ref_genome);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] Reference genome path does not exist!");
                all_set=false;
            }
        }

        if(options_intern.tepic_gene_annot.equals(""))
        {
            logger.logLine("[TEPIC] Gene annotation path is not given");
            all_set=false;
        }
        else
        {
            File f = new File(options_intern.tepic_gene_annot);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] Gene annotation path does not exist!");
                all_set=false;
            }
        }

        if(options_intern.tepic_input_directory.equals(""))
        {
            logger.logLine("[TEPIC] nfcore ChIP-seq data directory is not given");
            all_set=false;
        }
        else
        {
            File f = new File(options_intern.tepic_input_directory);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] nfcore ChIP-seq data path does not exist!");
                all_set=false;
            }
        }

        if(options_intern.tepic_path_pwms.equals(""))
        {
            logger.logLine("[TEPIC] position specific energy matrix directory is not given");
            all_set=false;
        }
        else
        {
            File f = new File(options_intern.tepic_path_pwms);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] position specific energy matrix path does not exist!");
                all_set=false;
            }
        }

        if(options_intern.tepic_tpm_cutoff>0)
        {
            //check for gene annotation
            if(options_intern.tepic_gene_annot.equals(""))
            {
                logger.logLine("[TEPIC] TPM cutoff set, but no annotation file is given");
                all_set=false;
            }
            else
            {
                File f = new File(options_intern.tepic_gene_annot);
                if(!f.exists())
                {
                    logger.logLine("[TEPIC] TPM cutoff set and gene annotation file path does not exist!");
                    all_set=false;
                }
            }
        }

        //check for map ensg symbol
        if(options_intern.tepic_ensg_symbol.equals(""))
        {
            logger.logLine("[TEPIC] No map of ENSG to Gene Symbol is given!");
        }
        else
        {
            File f = new File(options_intern.tepic_ensg_symbol);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] Map of ENSG to Gene Symbol file path does not exist!");
                all_set=false;
            }
        }

        if(!options_intern.tepic_background_seq.equals("") && !options_intern.tepic_2bit.equals(""))
        {
            logger.logLine("[TEPIC] parameters: tepic_background_seq and tepic_2bit are mutually exclusive");
            all_set=false;
        }

        if(options_intern.tepic_column_bedfile!=-1 && !options_intern.tepic_bed_chr_sign.equals(""))
        {
            logger.logLine("[TEPIC] parameters: tepic_column_bedfile and tepic_bed_chr_sign are mutually exclusive");
            all_set=false;
        }

        if(!options_intern.tepic_bed_chr_sign.equals(""))
        {
            File f = new File(options_intern.tepic_bed_chr_sign);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] tepic_bed_chr_sign file path does not exists!");
                all_set=false;
            }
        }

        if(!options_intern.tepic_psems_length.equals(""))
        {
            File f = new File(options_intern.tepic_psems_length);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] tepic_psems_length file path does not exists!");
                all_set=false;
            }
        }

        if(!options_intern.tepic_loop_list.equals(""))
        {
            File f = new File(options_intern.tepic_loop_list);
            if(!f.exists())
            {
                logger.logLine("[TEPIC] tepic_loop_list file path does not exists!");
                all_set=false;
            }
        }

        /**
         * TGEN options
         */
        if(!options_intern.path_tgen.equals(""))
        {
            File file_tgen = new File(options_intern.path_tgen);
            if(!file_tgen.exists() || !file_tgen.isDirectory())
            {
                logger.logLine("[TGENE] TGene file directory does not exist or is not a directory!");
                all_set=false;
            }

            File tgene_dir = new File(options_intern.path_tgen+File.separator+"bin");

            if(!tgene_dir.exists())
            {
                logger.logLine("[TGENE] TGene binary directory cannot be found: " + tgene_dir.getAbsolutePath());
                all_set=false;
            }

            if(options_intern.tgen_mt_writing.equals(""))
            {
                logger.logLine("[TGENE] Please specify spelling of Mitochondrial DNA, e.g. M or MT (default: MT)");
                all_set=false;
            }

            if(options_intern.tgen_consensus==0.0)
            {
                logger.logLine("[TGENE] tgen_consensus must be in range ]0.0,1.0], it cannot be 0.0, if you do not want to use consensus set path_tgen=\"\"");
                all_set=false;
            }

            if(options_intern.tgen_consensus_calc.equals(""))
            {
                logger.logLine("[TGENE] tgen_consensus_calc must be set.");
                all_set=false;
            }
        }

        /**
         * DYNAMITE OPTIONS
         */

        if(!options_intern.dynamite_preprocessing_integrate_data_consider_geneFile.equals(""))
        {
            File f = new File(options_intern.dynamite_preprocessing_integrate_data_consider_geneFile);
            if(!f.exists())
            {
                logger.logLine("[DYNAMITE] preprocessing consider genes file for integrateData.py does not exist!");
                all_set=false;
            }
        }

        /**
         * PLOT OPTIONS
         */

        if(options_intern.plot_th_coefficient.isEmpty())
        {
            logger.logLine("[PLOTS] plot th coefficients is empty, please use at least one coefficient.");
            all_set=false;
        }


        return all_set;
    }

    private HashMap<String, HashMap<String, HashSet<String>>> checkGroupsTEPIC()
    {
        HashMap<String, HashMap<String, HashSet<String>>> groups = new HashMap<>();

        File folder = new File(options_intern.com2pose_working_directory+File.separator+options_intern.folder_name_tepic_output_raw);

        for(File fileDir : folder.listFiles())
        {
            if(fileDir.isDirectory())
            {
                HashMap<String,HashSet<String>> n_timepoint = new HashMap<>();
                String n_tp_name = fileDir.getName();
                for(File fileDir2 : fileDir.listFiles())
                {
                    String n_hm_name="";
                    HashSet<String> n_hm = new HashSet<>();
                    n_hm_name= fileDir2.getName();

                    for(File fileDir3 : fileDir2.listFiles())
                    {
                        //samples
                        if(fileDir3.isDirectory())
                        {
                            n_hm.add(fileDir3.getName());
                        }
                    }
                    n_timepoint.put(n_hm_name,n_hm);

                }
                groups.put(n_tp_name,n_timepoint);
            }
        }
        return groups;
    }

    /**
     * BLACKLIST preprocessing - creates binary tree friendly files of chromosomes
     */
    private ArrayList<BL_ranges_binary_tree> recursive_split_BL(ArrayList<BL_ranges_binary_tree> region, ArrayList<BL_ranges_binary_tree> newly_ordered)
    {
        if(region.size()>0)
        {
            BL_ranges_binary_tree median = region.get(region.size()/2);
            newly_ordered.add(median);

            ArrayList<BL_ranges_binary_tree> region_left=new ArrayList<>();

            for(int i = 0; i < region.size()/2; i++)
            {
                region_left.add(region.get(i));
            }
            ArrayList<BL_ranges_binary_tree> region_right=new ArrayList<>();
            for(int i = region.size()/2+1; i < region.size(); i++)
            {
                region_right.add(region.get(i));
            }

            recursive_split_BL(region_left,newly_ordered);
            recursive_split_BL(region_right,newly_ordered);

        }
        return newly_ordered;
    }

    /**
     * TGENE preprocessing - creates binary tree friendly files of chromosomes
     */
    private  ArrayList<ENSG_ranges_binary_trees> recursive_split(ArrayList<ENSG_ranges_binary_trees> region, ArrayList<ENSG_ranges_binary_trees> newly_ordered)
    {
        if(region.size()>0)
        {
            ENSG_ranges_binary_trees median = region.get(region.size()/2);
            newly_ordered.add(median);

            ArrayList<ENSG_ranges_binary_trees> region_left=new ArrayList<>();
            for(int i = 0; i < region.size()/2; i++)
            {
                region_left.add(region.get(i));
            }
            ArrayList<ENSG_ranges_binary_trees> region_right=new ArrayList<>();
            for(int i = region.size()/2+1; i < region.size(); i++)
            {
                region_right.add(region.get(i));
            }

            recursive_split(region_left,newly_ordered);
            recursive_split(region_right,newly_ordered);

        }
        return newly_ordered;
    }

    /**
     * mixes the samples of one folder into one file, based on mix_option (UNION or INTERSECTION)
     */
    private static Stack<MIX_Interval> mergeIntervals(ArrayList<MIX_Interval> interval)
    {
        Stack<MIX_Interval> stack = new Stack<>();

        if(stack.empty())
        {
            stack.push(interval.get(0));
        }

        for(int i = 1; i < interval.size(); i++)
        {
            MIX_Interval top = stack.peek();

            if(top.end < interval.get(i).start)
            {
                stack.push(interval.get(i));
            } else if(top.end < interval.get(i).end)
            {
                top.end = interval.get(i).end;
                top.merged_intervals.addAll(interval.get(i).merged_intervals);

                stack.pop();
                stack.push(top);
            }
            else
            {
                top.merged_intervals.addAll(interval.get(i).merged_intervals);

                stack.pop();
                stack.push(top);
            }
        }
        return stack;
    }


}
