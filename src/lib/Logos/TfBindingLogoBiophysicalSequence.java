package lib.Logos;

import lib.ExecutableStep;
import org.apache.commons.compress.utils.IOUtils;
import tfprio.tfprio.TFPRIO;
import util.Configs.ConfigTypes.AbstractConfig;
import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.TrustAllManager;

import javax.net.ssl.*;
import java.io.*;
import java.net.URL;
import java.util.*;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class TfBindingLogoBiophysicalSequence extends ExecutableStep
{
    private final AbstractConfig<File> f_input_dcgResult =
            TFPRIO.configs.distributionAnalysis.fileStructure.f_dcg_stats;
    private final AbstractConfig<File> f_input_pwmFile = TFPRIO.configs.tepic.pathPwms;
    private final AbstractConfig<File> f_scriptTemplate_biophysicalModel =
            TFPRIO.configs.scriptTemplates.f_logos_biophysicalModel;

    private final GeneratedFileStructure d_output_biophysicalModel =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_logos_biophysicalModel;
    private final GeneratedFileStructure d_output_tfSequence =
            TFPRIO.configs.distributionAnalysis.fileStructure.d_logos_tfSequence;
    private final GeneratedFileStructure f_output_jaspar =
            TFPRIO.configs.distributionAnalysis.fileStructure.f_logos_tfSequence_jaspar_pfms;

    private final AbstractConfig<String> s_biophysicalModel_data =
            TFPRIO.configs.distributionAnalysis.fileStructure.s_logos_biophysicalModel_data;
    private final AbstractConfig<String> s_biophysicalModel_script =
            TFPRIO.configs.distributionAnalysis.fileStructure.s_logos_biophysicalModel_script;
    private final AbstractConfig<String> s_biophysicalModel_plot =
            TFPRIO.configs.distributionAnalysis.fileStructure.s_logos_biophysicalModel_image;
    private final AbstractConfig<String> jasparUrl = TFPRIO.configs.jaspar.downloadUrl;
    private final AbstractConfig<String> jasparApi = TFPRIO.configs.jaspar.apiCall;
    private final AbstractConfig<String> jasparLogoDownloadUrl = TFPRIO.configs.jaspar.logoDownloadUrl;
    private final AbstractConfig<String> s_jasparJson =
            TFPRIO.configs.distributionAnalysis.fileStructure.s_logos_tfSequence_jaspar_json;
    private final AbstractConfig<String> s_jasparImage =
            TFPRIO.configs.distributionAnalysis.fileStructure.s_logos_tfSequence_jaspar_image;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        return new HashSet<>(Arrays.asList(f_input_dcgResult, f_input_pwmFile, f_scriptTemplate_biophysicalModel));
    }

    @Override public Set<GeneratedFileStructure> getCreatedFileStructure()
    {
        return new HashSet<>(Arrays.asList(d_output_biophysicalModel, d_output_tfSequence, f_output_jaspar));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(
                Arrays.asList(s_biophysicalModel_data, s_biophysicalModel_script, s_biophysicalModel_plot, jasparUrl,
                        jasparApi, jasparLogoDownloadUrl, s_jasparJson, s_jasparImage));
    }

    @Override protected void execute()
    {
        String scriptTemplate = null;
        try
        {
            scriptTemplate = readFile(f_scriptTemplate_biophysicalModel.get());
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
        assert scriptTemplate != null;

        ArrayList<String> tfSymbols = new ArrayList<>();

        try (BufferedReader reader = new BufferedReader(new FileReader(f_input_dcgResult.get())))
        {
            String inputLine;
            reader.readLine();
            while ((inputLine = reader.readLine()) != null)
            {
                String[] split = inputLine.split("\t");
                tfSymbols.add(split[1]);
            }
        } catch (IOException e)
        {
            e.printStackTrace();
        }

        logger.info("Preparing energies of biophysical model for logo creation");

        //retrieve binding energy
        HashMap<String, ArrayList<String>> tfSymbol_bindingProfile = new HashMap<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(f_input_pwmFile.get())))
        {
            String inputLine;
            reader.readLine();
            String tfName = null;

            while ((inputLine = reader.readLine()) != null)
            {
                if (inputLine.startsWith(">"))
                {
                    tfName = inputLine.split("\t")[1].split("_")[0].toUpperCase().replace(":", ".");

                    tfSymbol_bindingProfile.put(tfName, new ArrayList<>());
                } else
                {
                    if (tfName != null)
                    {
                        tfSymbol_bindingProfile.get(tfName).add(inputLine);
                    }
                }
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        for (int i = 0; i < tfSymbols.size(); i++)
        {
            int rank = i + 1;
            String tfName = tfSymbols.get(i);

            File d_biophysicalModel_tf = extend(d_output_biophysicalModel.get(), rank + "_" + tfName);

            tfName = tfName.replaceAll("\\.", ":");
            ArrayList<String> tfBindingEnergies;
            if (tfName.contains("::"))
            {
                tfBindingEnergies = new ArrayList<>();
                String[] split = tfName.split("::");
                for (String singleTf : split)
                {
                    if (tfSymbol_bindingProfile.containsKey(singleTf))
                    {
                        tfBindingEnergies.addAll(tfSymbol_bindingProfile.get(singleTf));
                    }
                }
                if (tfSymbol_bindingProfile.containsKey(tfSymbols.get(i)))
                {
                    tfBindingEnergies.addAll(tfSymbol_bindingProfile.get(tfSymbols.get(i)));
                }
            } else
            {
                tfBindingEnergies = tfSymbol_bindingProfile.get(tfName);
            }

            File f_output_data = extend(d_biophysicalModel_tf, s_biophysicalModel_data.get());
            makeSureFileExists(f_output_data, logger);
            try (BufferedWriter bw_binding_energies = new BufferedWriter(new FileWriter(f_output_data)))
            {
                for (String line_energy : tfBindingEnergies)
                {
                    bw_binding_energies.write(line_energy);
                    bw_binding_energies.newLine();
                }
            } catch (IOException e)
            {
                logger.error(e.getMessage());
            }

            File f_output_script = extend(d_biophysicalModel_tf, s_biophysicalModel_script.get());
            File f_output_plot = extend(d_biophysicalModel_tf, s_biophysicalModel_plot.get());
            makeSureFileExists(f_output_plot, logger);

            String script = scriptTemplate;
            script = script.replace("{INPUTFILE}", f_output_data.getAbsolutePath());
            script = script.replace("{OUTPUTFILE}", f_output_plot.getAbsolutePath());

            try
            {
                writeFile(f_output_script, script);
            } catch (IOException e)
            {
                logger.error(e.getMessage());
            }
            executorService.submit(() -> executeAndWait(f_output_script, logger));
        }
        finishAllQueuedThreads();
        logger.info("Download JASPAR matrix ids from: " + jasparUrl.get());

        // Create a new trust manager that trust all certificates
        TrustManager[] trustAllCerts = new TrustManager[]{new TrustAllManager()};
        // Activate the new trust manager
        try
        {
            SSLContext sc = SSLContext.getInstance("SSL");
            sc.init(null, trustAllCerts, new java.security.SecureRandom());
            HttpsURLConnection.setDefaultSSLSocketFactory(sc.getSocketFactory());
        } catch (Exception e)
        {
            logger.error(e.getMessage());
        }

        makeSureFileExists(f_output_jaspar.get(), logger);

        try (InputStream is = new URL(jasparUrl.get()).openConnection().getInputStream();
             OutputStream outputStream = new FileOutputStream(f_output_jaspar.get()))
        {
            IOUtils.copy(is, outputStream);
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        logger.info("JASPAR matrix id download complete.");

        HashMap<String, HashSet<String>> tf_to_matrixID = new HashMap<>();
        try (BufferedReader reader = new BufferedReader(new FileReader(f_output_jaspar.get())))
        {
            String line_jaspar;
            while ((line_jaspar = reader.readLine()) != null)
            {
                if (line_jaspar.startsWith(">"))
                {
                    String[] split = line_jaspar.split("\t");
                    String matrixID = split[0].substring(1);
                    String tf = split[1].toUpperCase();

                    if (tf_to_matrixID.containsKey(tf))
                    {
                        tf_to_matrixID.get(tf).add(matrixID);
                    } else
                    {
                        HashSet<String> x = new HashSet<>();
                        x.add(matrixID);
                        tf_to_matrixID.put(tf, x);
                    }
                }
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        logger.info("Retrieve additional information from JASPAR API: " + jasparApi.get());
        logger.info("Retrieve sequence logos from JASPAR: " + jasparLogoDownloadUrl.get());

        for (int i = 0; i < tfSymbols.size(); i++)
        {
            int rank = i + 1;
            String tfSymbol = tfSymbols.get(i);

            File d_tf = extend(d_output_tfSequence.get(), rank + "_" + tfSymbol);

            if (tf_to_matrixID.containsKey(tfSymbol))
            {
                HashSet<String> matrix_ids = tf_to_matrixID.get(tfSymbol);

                for (String matrixId : matrix_ids)
                {
                    executorService.submit(() ->
                    {
                        //get json via curl
                        String command = jasparApi.get();
                        command = command.replace("{MATRIXID}", matrixId);
                        File f_output = extend(d_tf, matrixId + s_jasparJson.get());
                        makeSureFileExists(f_output, logger);
                        command = command.replace("{OUTPUTFILE}", f_output.getAbsolutePath());

                        try
                        {
                            logger.debug("Excuting command: " + command);
                            Process curl  = Runtime.getRuntime().exec(command);
                            int returnCode = curl.waitFor();
                            if(returnCode != 0)
                            {
                                System.out.println(new String(curl.getErrorStream().readAllBytes()));
                                throw new IOException("Received return code: " + returnCode);
                            }
                        } catch (IOException | InterruptedException e)
                        {
                            logger.error(e.getMessage());
                        }

                        //get logos via static file download
                        // Create a new trust manager that trust all certificates
                        File f_matrix_svg = extend(d_tf, matrixId + s_jasparImage.get());
                        String matrix_svg_download_url = jasparLogoDownloadUrl.get() + matrixId + s_jasparImage.get();

                        try
                        {
                            SSLContext sc = SSLContext.getInstance("SSL");
                            sc.init(null, trustAllCerts, new java.security.SecureRandom());
                            HttpsURLConnection.setDefaultSSLSocketFactory(sc.getSocketFactory());
                        } catch (Exception e)
                        {
                            logger.error(e.getMessage());
                        }

                        try (InputStream is = new URL(matrix_svg_download_url).openConnection().getInputStream();
                             OutputStream outputStream = new FileOutputStream(f_matrix_svg))
                        {
                            IOUtils.copy(is, outputStream);
                        } catch (IOException e)
                        {
                            logger.error(e.getMessage());
                        }
                    });
                }
            }
        }
    }
}
