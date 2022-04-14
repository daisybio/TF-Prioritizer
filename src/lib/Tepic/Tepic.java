package lib.Tepic;

import lib.ExecutableStep;
import tfprio.TFPRIO;
import tfprio.Workflow;
import util.Configs.ConfigTypes.AbstractConfig;
import util.FileFilters.Filters;

import java.io.*;
import java.util.*;

import static util.FileManagement.*;
import static util.ScriptExecution.executeAndWait;

public class Tepic extends ExecutableStep
{
    private AbstractConfig<File> d_input;
    private final AbstractConfig<File> d_output = TFPRIO.configs.tepic.fileStructure.d_outputRaw;

    private final AbstractConfig<File> f_referenceGenome = TFPRIO.configs.tepic.inputReferenceGenome;
    private final AbstractConfig<File> f_pwms = TFPRIO.configs.tepic.pathPwms;
    private final AbstractConfig<Integer> tpmCutoff = TFPRIO.configs.tepic.tpmCutoff;
    private final AbstractConfig<File> tepicExecutable = TFPRIO.configs.tepic.executable;
    private final AbstractConfig<String> tfBindingSiteSearch = TFPRIO.configs.tepic.tfBindingSiteSearch;
    private final AbstractConfig<Boolean> mixMutuallyExclusive = TFPRIO.configs.mixOptions.mutuallyExclusive;
    private final AbstractConfig<Boolean> tgeneTargetGenes = TFPRIO.configs.tepic.tgeneTargetGenes;
    private final AbstractConfig<String> s_tgene_links = TFPRIO.configs.tgene.fileStructure.s_output_links;
    private final AbstractConfig<File> d_deseq2_meanCounts =
            TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_meanCounts;
    private final AbstractConfig<File> d_tgene_output = TFPRIO.configs.tgene.fileStructure.d_output;
    private final AbstractConfig<String> s_outputRaw_trapSequences =
            TFPRIO.configs.tepic.fileStructure.s_outputRaw_trapSequences;

    // if tpmCut > 0
    private final AbstractConfig<File> ensgSymbolFile = TFPRIO.configs.tepic.ensgSymbolFile;
    private final AbstractConfig<File> inputGeneID = TFPRIO.configs.deSeq2.inputGeneID;

    // Optional configs
    private final AbstractConfig<Integer> threadLimit = TFPRIO.configs.general.threadLimit;
    private final AbstractConfig<File> bedChromatinSignal = TFPRIO.configs.tepic.bedChromatinSignal;
    private final AbstractConfig<Integer> columnBedfile = TFPRIO.configs.tepic.columnBedfile;
    private final AbstractConfig<File> geneAnnotationFile = TFPRIO.configs.tepic.geneAnnotationFile;
    private final AbstractConfig<Integer> windowSize = TFPRIO.configs.tepic.windowSize;
    private final AbstractConfig<File> onlyDNasePeaks = TFPRIO.configs.tepic.onlyDNasePeaks;
    private final AbstractConfig<Boolean> exponentialDecay = TFPRIO.configs.tepic.exponentialDecay;
    private final AbstractConfig<Boolean> doNotNormalizePeakLength = TFPRIO.configs.tepic.doNotNormalizePeakLength;
    private final AbstractConfig<Boolean> doNotGenerate = TFPRIO.configs.tepic.doNotGenerate;
    private final AbstractConfig<Boolean> originalDecay = TFPRIO.configs.tepic.originalDecay;
    private final AbstractConfig<File> psemsLengthFile = TFPRIO.configs.tepic.psemsLengthFile;
    private final AbstractConfig<Boolean> entireGeneBody = TFPRIO.configs.tepic.entireGeneBody;
    private final AbstractConfig<Boolean> doZip = TFPRIO.configs.tepic.doZip;
    private final AbstractConfig<File> twoBitFile = TFPRIO.configs.tepic.twoBitFile;
    private final AbstractConfig<Double> pValue = TFPRIO.configs.tepic.pValue;
    private final AbstractConfig<Integer> maxMinutesPerChromosome = TFPRIO.configs.tepic.maxMinutesPerChromosome;
    private final AbstractConfig<Boolean> chromosomePrefix = TFPRIO.configs.tepic.chromosomePrefix;
    private final AbstractConfig<Boolean> transcriptBased = TFPRIO.configs.tepic.transcriptBased;
    private final AbstractConfig<File> loopListFile = TFPRIO.configs.tepic.loopListFile;
    private final AbstractConfig<Integer> loopWindows = TFPRIO.configs.tepic.loopWindows;
    private final AbstractConfig<Boolean> onlyPeakFeatures = TFPRIO.configs.tepic.onlyPeakFeatures;

    @Override protected Set<AbstractConfig<File>> getRequiredFileStructure()
    {
        Set<AbstractConfig<File>> requirements = new HashSet<>(
                Arrays.asList(d_input, f_referenceGenome, f_pwms, tepicExecutable, d_deseq2_meanCounts,
                        d_tgene_output));

        if (tpmCutoff.get() > 0)
        {
            requirements.add(ensgSymbolFile);
            requirements.add(inputGeneID);
        }
        return requirements;
    }

    @Override protected Set<AbstractConfig<File>> getCreatedFileStructure()
    {
        return new HashSet<>(List.of(d_output));
    }

    @Override protected Set<AbstractConfig<?>> getRequiredConfigs()
    {
        return new HashSet<>(
                Arrays.asList(tpmCutoff, tfBindingSiteSearch, mixMutuallyExclusive, tgeneTargetGenes, s_tgene_links,
                        s_outputRaw_trapSequences));
    }

    @Override protected Set<AbstractConfig<?>> getOptionalConfigs()
    {
        return new HashSet<>(
                Arrays.asList(threadLimit, bedChromatinSignal, columnBedfile, geneAnnotationFile, windowSize,
                        onlyDNasePeaks, exponentialDecay, doNotNormalizePeakLength, doNotGenerate, originalDecay,
                        psemsLengthFile, entireGeneBody, doZip, twoBitFile, pValue, maxMinutesPerChromosome,
                        chromosomePrefix, transcriptBased, loopListFile, loopWindows, onlyPeakFeatures));
    }

    @Override protected void updateInputDirectory()
    {
        d_input = Workflow.getLatestInputDirectory();
    }

    @Override protected void execute()
    {
        //check_tepic_input_with_options();

        File output_TEPIC = d_output.get();

        for (File d_group : Objects.requireNonNull(d_input.get().listFiles(Filters.directoryFilter)))
        {
            for (File d_hm : Objects.requireNonNull(d_group.listFiles(Filters.directoryFilter)))
            {
                for (File f_sample : Objects.requireNonNull(d_hm.listFiles(Filters.fileFilter)))
                {
                    Map<String, String> sampleConfigs = new HashMap<>();
                    String sampleName;
                    if (f_sample.getName().contains("."))
                    {
                        sampleName = f_sample.getName().substring(0, f_sample.getName().lastIndexOf("."));
                    } else
                    {
                        sampleName = f_sample.getName();
                    }

                    File d_output = extend(output_TEPIC, d_group.getName(), d_hm.getName(), sampleName);
                    File d_output_combined = extend(d_output, sampleName);

                    try
                    {
                        makeSureDirectoryExists(d_output);
                        makeSureDirectoryExists(d_output_combined);
                    } catch (IOException e)
                    {
                        logger.error(e.getMessage());
                    }

                    sampleConfigs.put("b", f_sample.getAbsolutePath());
                    sampleConfigs.put("o", d_output_combined.getAbsolutePath());

                    if (tpmCutoff.get() > 0)
                    {
                        File n_dir;

                        if (!mixMutuallyExclusive.get())
                        {
                            n_dir = extend(d_deseq2_meanCounts.get(), d_group.getName() + ".tsv");

                        } else
                        {
                            String[] names_sample = f_sample.getName().split("_");
                            String name_sample = names_sample[0];
                            n_dir = extend(d_deseq2_meanCounts.get(), name_sample + ".tsv");
                        }

                        sampleConfigs.put("G", n_dir.getAbsolutePath());
                    }

                    if (TFPRIO.configs.tgene.pathToExecutable.isSet() && tgeneTargetGenes.get())
                    {
                        File f_tgene_input;

                        if (!mixMutuallyExclusive.get())
                        {
                            String sample_name = f_sample.getName();

                            if (sample_name.contains("."))
                            {
                                sample_name = sample_name.substring(0, sample_name.lastIndexOf("."));
                            }
                            f_tgene_input = extend(d_tgene_output.get(), d_group.getName(), d_hm.getName(), sample_name,
                                    s_tgene_links.get());
                        } else
                        {
                            String[] names_sample = f_sample.getName().split("_");
                            String name_sample = names_sample[0];

                            f_tgene_input = extend(d_tgene_output.get(), d_group.getName(), d_hm.getName(),
                                    name_sample + "_" + d_hm.getName(), s_tgene_links.get());
                        }
                        sampleConfigs.put("L", f_tgene_input.getAbsolutePath());
                    }

                    sampleConfigs.put("S", extend(d_output, s_outputRaw_trapSequences.get()).getAbsolutePath());

                    executorService.submit(() -> executeAndWait(getCommand(sampleConfigs), logger));
                }
            }
        }
    }

    private void updateTPM(File d_group, File d_output_combined)
    {
        //correct TPMs
        HashMap<String, String> ensgTF_tpm = new HashMap<>();

        File d_output_tpms = d_output_combined.getParentFile();

        for (File f_dir : Objects.requireNonNull(d_output_tpms.listFiles(Filters.fileFilter)))
        {
            String name = f_dir.getName();
            if (name.endsWith("_TPM_values.txt"))
            {
                try (BufferedReader reader = new BufferedReader(new FileReader(f_dir)))
                {
                    String line_tpms;
                    reader.readLine();
                    while ((line_tpms = reader.readLine()) != null)
                    {
                        String[] split = line_tpms.split("\t");
                        ensgTF_tpm.put(split[0], split[2]);
                    }
                } catch (IOException e)
                {
                    logger.error(e.getMessage());
                }
            }
        }

        //get correct TPM file
        File f_tpm_toBeChanged = extend(TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_tpm_tpmResults.get(),
                d_group.getName() + ".tsv");
        File f_tpmChanged = extend(TFPRIO.configs.deSeq2.fileStructure.d_preprocessing_tpm_updated.get(),
                d_group.getName() + ".tsv");
        try
        {
            makeSureFileExists(f_tpmChanged);
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }

        try (BufferedReader reader = new BufferedReader(new FileReader(f_tpm_toBeChanged));
             BufferedWriter writer = new BufferedWriter(new FileWriter(f_tpmChanged)))
        {
            String inputLine;
            writer.write(reader.readLine());
            writer.newLine();

            while ((inputLine = reader.readLine()) != null)
            {
                String[] split = inputLine.split("\t");
                String ensg = split[0];

                if (ensgTF_tpm.containsKey(ensg))
                {
                    String tpm_exchanged = ensgTF_tpm.get(ensg);
                    split[3] = tpm_exchanged;

                    StringBuilder sb_tf_line = new StringBuilder();
                    sb_tf_line.append(ensg);

                    for (int i = 1; i < split.length; i++)
                    {
                        sb_tf_line.append("\t");
                        sb_tf_line.append(split[i]);
                    }

                    writer.write(sb_tf_line.toString());
                    writer.newLine();
                } else
                {
                    writer.write(inputLine);
                    writer.newLine();
                }
            }
        } catch (IOException e)
        {
            logger.error(e.getMessage());
        }
    }

    private String getCommand(Map<String, String> sampleConfigs)
    {
        StringBuilder sb_command = new StringBuilder();
        sb_command.append(tepicExecutable.get().getAbsolutePath());

        Map<String, AbstractConfig<?>> configs = new HashMap<>()
        {{
            put("g", f_referenceGenome);
            put("p", f_pwms);
            put("d", bedChromatinSignal);
            put("n", columnBedfile);
            put("a", geneAnnotationFile);
            put("w", windowSize);
            put("f", onlyDNasePeaks);
            put("e", exponentialDecay);
            put("l", doNotNormalizePeakLength);
            put("u", doNotGenerate);
            put("x", originalDecay);
            put("m", psemsLengthFile);
            put("y", entireGeneBody);
            put("z", doZip);
            put("r", twoBitFile);
            put("v", pValue);
            put("i", maxMinutesPerChromosome);
            put("j", chromosomePrefix);
            put("t", transcriptBased);
            put("h", loopListFile);
            put("s", loopWindows);
            put("q", onlyPeakFeatures);

            if (tpmCutoff.get() > 0)
            {
                put("T", tpmCutoff);
                put("E", ensgSymbolFile);
                put("A", inputGeneID);
            }

            put("B", tfBindingSiteSearch);
        }};

        Map<String, String> stringConfigs = new HashMap<>();

        for (Map.Entry<String, AbstractConfig<?>> entry : configs.entrySet())
        {
            if (entry.getValue().isSet())
            {
                if (entry.getValue().get().getClass().equals(Boolean.class))
                {
                    if ("eq".contains(entry.getKey()))
                    {
                        // Boolean parameters
                        stringConfigs.put(entry.getKey(), entry.getValue().toString().toUpperCase());
                    } else
                    {
                        // Flags
                        if ((boolean) entry.getValue().get())
                        {
                            stringConfigs.put(entry.getKey(), "");
                        }
                    }
                } else
                {
                    stringConfigs.put(entry.getKey(), entry.getValue().toString());
                }
            }
        }

        stringConfigs.putAll(sampleConfigs);
        stringConfigs.put("c", "1");

        for (Map.Entry<String, String> entry : stringConfigs.entrySet())
        {
            sb_command.append(" -").append(entry.getKey()).append(" ").append(entry.getValue());
        }

        return sb_command.toString();
    }

    @Override protected int getShutDownTimeOutMinutes()
    {
        return Integer.MAX_VALUE;
    }
}
