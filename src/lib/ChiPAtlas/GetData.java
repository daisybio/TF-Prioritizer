package lib.ChiPAtlas;

import lib.ExecutableStep;
import org.apache.commons.compress.utils.IOUtils;
import tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.TrustAllManager;

import javax.net.ssl.HttpsURLConnection;
import javax.net.ssl.SSLContext;
import javax.net.ssl.TrustManager;
import java.io.*;
import java.net.URL;
import java.nio.file.Files;
import java.util.*;

import static util.FileManagement.extend;
import static util.FileManagement.makeSureFileExists;
import static util.ScriptExecution.executeAndWait;

public class GetData extends ExecutableStep
{
    private final AbstractConfig<File> f_input = TFPRIO.configs.chipAtlas.fileStructure.f_list_csv;
    private final AbstractConfig<File> f_input_dcg = TFPRIO.configs.distributionAnalysis.fileStructure.f_dcg_stats;

    private final AbstractConfig<File> d_output = TFPRIO.configs.chipAtlas.fileStructure.d_peakFiles;

    private final AbstractConfig<String> tissueTypes = TFPRIO.configs.chipAtlas.tissueType;
    private final AbstractConfig<String> genomeVersion = TFPRIO.configs.chipAtlas.genomeVersion;
    private final AbstractConfig<File> igvPath = TFPRIO.configs.igv.pathToIGV;
    // r for regex
    private final AbstractConfig<String> r_geneVersion = TFPRIO.configs.chipAtlas.column_GeneVersion;
    private final AbstractConfig<String> r_antigenClass = TFPRIO.configs.chipAtlas.column_AntigenClass;
    private final AbstractConfig<String> r_antigen = TFPRIO.configs.chipAtlas.column_Antigen;
    private final AbstractConfig<String> r_cellTypeClass = TFPRIO.configs.chipAtlas.column_CellTypeClass;
    private final AbstractConfig<String> r_url = TFPRIO.configs.chipAtlas.column_url;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_input, f_input_dcg));
    }

    @Override protected Set<AbstractConfig<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(
                Arrays.asList(tissueTypes, genomeVersion, r_geneVersion, r_antigenClass, r_antigen, r_cellTypeClass,
                        r_url, igvPath));
    }

    @Override protected void execute()
    {
        HashMap<String, String> tf_url = new HashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(f_input.get())))
        {
            final String header = reader.readLine();
            final String[] splitHeader = header.split(",");

            Integer c_geneVersion = null;
            Integer c_antigenClass = null;
            Integer c_antigen = null;
            Integer c_cellTypeClass = null;
            Integer c_url = null;

            String[] types = tissueTypes.get().split(";");

            ArrayList<String> regexes = new ArrayList<>();
            for (String type : types)
            {
                StringBuilder regex = new StringBuilder();

                for (String subType : type.split(" "))
                {
                    regex.append(".*").append(subType);
                }
                regex.append(".*");

                regexes.add(regex.toString().toUpperCase());
            }

            StringBuilder r_tissueType = new StringBuilder();

            if (regexes.size() > 1)
            {
                r_tissueType.append("(");
                r_tissueType.append(String.join("|", regexes));
                r_tissueType.append(")");
            } else if (regexes.size() == 1)
            {
                r_tissueType.append(regexes.get(0));
            } else
            {
                logger.error("There were no valid tissue types found inside " + tissueTypes.getName() + ".");
            }


            for (int i = 0; i < splitHeader.length; i++)
            {
                if (splitHeader[i].matches(r_geneVersion.get()))
                {
                    c_geneVersion = i;
                } else if (splitHeader[i].matches(r_antigenClass.get()))
                {
                    c_antigenClass = i;
                } else if (splitHeader[i].matches(r_antigen.get()))
                {
                    c_antigen = i;
                } else if (splitHeader[i].matches(r_cellTypeClass.get()))
                {
                    c_cellTypeClass = i;
                } else if (splitHeader[i].matches(r_url.get()))
                {
                    c_url = i;
                }
            }

            String inputLine;
            while ((inputLine = reader.readLine()) != null)
            {
                String[] split = inputLine.split(",");
                if (!split[c_geneVersion].equals(genomeVersion.get()))
                {
                    continue;
                }
                if (!split[c_antigenClass].equals("TFs and others"))
                {
                    continue;
                }
                if (!split[c_cellTypeClass].toUpperCase().matches(r_tissueType.toString()))
                {
                    continue;
                }

                String url = split[split.length - 1];
                String tf = split[c_antigen];
                if (!tf.equals(""))
                {
                    tf_url.put(tf.toUpperCase(), url);
                }

            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        logger.info("Downloading available TF data and creating IGV index files");


        ArrayList<String> tfSymbols = new ArrayList<>();
        HashMap<String, HashSet<String>> tf_splits = new HashMap<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(f_input_dcg.get())))
        {
            String inputLine;
            reader.readLine();
            while ((inputLine = reader.readLine()) != null)
            {
                String tfSymbol = inputLine.split("\t")[1];

                tfSymbols.add(tfSymbol);
                tf_splits.put(tfSymbol.toUpperCase(), new HashSet<>(List.of(tfSymbol.split("\\."))));
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        for (int i = 0; i < tfSymbols.size(); i++)
        {
            String tfSymbol = tfSymbols.get(i);
            boolean found_tf = false;

            ArrayList<String> downloadTfs = new ArrayList<>();

            for (String tfSplit : tf_splits.get(tfSymbol.toUpperCase()))
            {
                if (tf_url.containsKey(tfSplit))
                {
                    found_tf = true;
                    downloadTfs.add(tfSplit);
                }
            }

            if (!found_tf)
            {
                continue;
            }

            int rank = i + 1;

            File d_output_tf = extend(d_output.get(), rank + "_" + tfSymbol);

            for (String tfSplit : downloadTfs)
            {
                executorService.submit(() ->
                {
                    File f_output_tfSplit = extend(d_output_tf, tfSplit + ".bed");

                    TrustManager[] trustAllCerts = new TrustManager[]{new TrustAllManager()};
                    try
                    {
                        SSLContext sc = SSLContext.getInstance("SSL");
                        sc.init(null, trustAllCerts, new java.security.SecureRandom());
                        HttpsURLConnection.setDefaultSSLSocketFactory(sc.getSocketFactory());
                    } catch (Exception e)
                    {
                        logger.error(e.getMessage());
                    }

                    makeSureFileExists(f_output_tfSplit, logger);

                    try (InputStream inputStream = new URL(tf_url.get(tfSplit)).openConnection().getInputStream();
                         OutputStream outputStream = new FileOutputStream(f_output_tfSplit))
                    {
                        IOUtils.copy(inputStream, outputStream);
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }

                    //check if valid bed file, if not, delete
                    long size = 0;
                    try
                    {
                        size = Files.size(f_output_tfSplit.toPath());
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }

                    if (size < 300)
                    {
                        f_output_tfSplit.delete();
                        d_output_tf.delete();
                    } else
                    {
                        //create index file for IGV
                        String command = extend(igvPath.get(), "igvtools").getAbsolutePath() + " index " +
                                f_output_tfSplit.getAbsolutePath();

                        executeAndWait(command, logger);
                    }
                });
            }
        }
    }
}
