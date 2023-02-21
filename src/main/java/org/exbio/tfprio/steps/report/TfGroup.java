package org.exbio.tfprio.steps.report;

import org.json.JSONObject;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Path;
import java.util.*;

import static org.exbio.pipejar.util.FileManagement.extend;
import static org.exbio.pipejar.util.FileManagement.softLink;

public class TfGroup {
    private final String symbol;
    private final Collection<TranscriptionFactor> transcriptionFactors;
    private final Map<String, Double> hmDcg = new HashMap<>();
    private final File outDir;
    private final Path base;
    private final Map<String, Map<String, String>> hmPairingHeatmap = new HashMap<>();
    private final Map<String, Map<String, Map<String, String>>> hmPairingTgIgv = new HashMap<>();
    private final Map<String, String> tfSequence = new HashMap<>();
    private final Map<String, String> distributionPlots = new HashMap<>();
    private final Map<String, Map<String, Double>> hmPairingRegressionCoefficients = new HashMap<>();

    private String biophysicalLogo;
    private JSONObject confusionMatrix;

    public TfGroup(String symbol, Collection<TranscriptionFactor> transcriptionFactors, File srcDir) {
        this.symbol = symbol;
        this.transcriptionFactors = transcriptionFactors;
        this.outDir = extend(srcDir, "assets", "tfs", symbol);
        this.base = srcDir.toPath();
    }

    public void setBiophysicalLogo(File biophysicalLogo) throws IOException {
        File out = extend(outDir, "biophysicalLogo.png");
        softLink(out, biophysicalLogo);
        this.biophysicalLogo = base.relativize(out.toPath()).toString();
    }

    public void setTfSequence(Collection<File> tfSequence) {
        File outDir = extend(this.outDir, "tfSequence");
        tfSequence.stream().forEach(f -> {
            String name = f.getName().substring(0, f.getName().lastIndexOf('.'));
            File out = extend(outDir, f.getName());
            try {
                softLink(out, f);
            } catch (IOException e) {
                throw new UncheckedIOException(e);
            }
            this.tfSequence.put(name, base.relativize(out.toPath()).toString());
        });
    }

    public void setHeatmap(String pairing, String hm, File file) {
        File out = extend(this.outDir, "heatmaps", hm, pairing, "heatmap.png");
        try {
            softLink(out, file);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        hmPairingHeatmap.computeIfAbsent(hm, k -> new HashMap<>()).put(pairing,
                base.relativize(out.toPath()).toString());
    }

    public void setIgv(String pairing, String hm, File tgDirectory) {
        File outCurrent = extend(this.outDir, "igv", hm, pairing);
        Arrays.stream(Objects.requireNonNull(
                tgDirectory.listFiles(f -> f.getName().endsWith(".svg") || f.getName().endsWith(".png")))).forEach(
                tgFile -> {
                    String tg = tgFile.getName().substring(0, tgFile.getName().lastIndexOf('.'));
                    File out = new File(outCurrent, tgFile.getName());
                    try {
                        softLink(out, tgFile);
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                    hmPairingTgIgv.computeIfAbsent(hm, k -> new HashMap<>()).computeIfAbsent(pairing,
                            k -> new HashMap<>()).put(tg, base.relativize(out.toPath()).toString());
                });
    }

    public void setDistributionPlot(String hm, File file) {
        File out = extend(this.outDir, "distributionPlots", hm + ".png");
        try {
            softLink(out, file);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        distributionPlots.put(hm, base.relativize(out.toPath()).toString());
    }

    public void setConfusionMatrix(JSONObject confusionMatrix) {
        this.confusionMatrix = confusionMatrix;
    }


    public void setRank(String hm, Double rank) {
        hmDcg.put(hm, rank);
    }

    public Collection<TranscriptionFactor> getTranscriptionFactors() {
        return transcriptionFactors;
    }

    public String getSymbol() {
        return symbol;
    }

    public JSONObject toJSON() {
        return new JSONObject() {{
            put("symbol", symbol);
            put("hmDcg", hmDcg);
            put("transcriptionFactors", transcriptionFactors.stream().map(TranscriptionFactor::toJSON).toList());
            put("biophysicalLogo", biophysicalLogo);
            put("tfSequence", tfSequence);
            put("heatmaps", hmPairingHeatmap);
            put("igv", hmPairingTgIgv);
            put("distributionPlots", distributionPlots);
            put("regressionCoefficients", hmPairingRegressionCoefficients);
            put("confusionMatrix", confusionMatrix);
        }};
    }

    public void setRegressionCoefficient(String hm, String pairing, Double value) {
        hmPairingRegressionCoefficients.computeIfAbsent(hm, k -> new HashMap<>()).put(pairing, value);
    }
}
