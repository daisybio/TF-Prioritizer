package lib.Report;

import java.io.File;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.net.URLEncoder;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.util.*;

import tfprio.TFPRIO;
import org.json.JSONObject;
import util.FileManagement;
import util.Comparators.IntegerStringComparator;

public class StructureElements
{
    private static final DecimalFormat formatter = new DecimalFormat("0.###");

    static String getBasicData(TranscriptionFactor transcriptionFactor) throws IOException
    {
        String template = FileManagement.readFile(TFPRIO.configs.report.inputStructure.f_basicData.get());


        String log2fcTemplate = FileManagement.readFile(TFPRIO.configs.report.inputStructure.f_basicDataEntry.get());

        log2fcTemplate = log2fcTemplate.replace("{NAME}", "LOG2FC");

        log2fcTemplate = log2fcTemplate.replace("{DATA}",
                getTabularData(transcriptionFactor.getName(), transcriptionFactor.getLog2fc()));

        log2fcTemplate = log2fcTemplate.replace("{INFO-ID}", transcriptionFactor.getName() + "-" + "log2fc");

        template = template.replace("{GENEID}",
                FileManagement.readFile(TFPRIO.configs.report.inputStructure.f_basicDataGeneID.get()));
        template = template.replace("{GENE-ID}", transcriptionFactor.getGeneID());
        template = template.replace("{GENECARDS-LINK}",
                TFPRIO.configs.report.genecardsUrl.get().replace("{GENE}", transcriptionFactor.getName()));
        template = template.replace("{TF-NAME}", transcriptionFactor.getName());

        template = template.replace("{LOG2FC}",
                getBasicDataEntry(transcriptionFactor.getName(), "LOG2FC", transcriptionFactor.getLog2fc()));

        template = template.replace("{TPM}", getBasicDataEntry(transcriptionFactor.getName(), "TPM", new HashMap<>()
        {{
            put("", transcriptionFactor.getTpm());
        }}));

        template = template.replace("{NORMEX}",
                getBasicDataEntry(transcriptionFactor.getName(), "Normalized expression", new HashMap<>()
                {{
                    put("", transcriptionFactor.getNormex());
                }}));

        return template;
    }

    static String getBasicData(TranscriptionFactorGroup tfGroup) throws IOException
    {
        if (!tfGroup.realGroup)
        {
            return getBasicData(tfGroup.getTranscriptionFactors().get(0));
        }

        StringBuilder basicData = new StringBuilder();
        String tfTemplate = FileManagement.readFile(TFPRIO.configs.report.inputStructure.f_home_tf.get());

        for (TranscriptionFactor tf : tfGroup.getTranscriptionFactors())
        {
            String tfString = tfTemplate.replace("{BASICDATA}", getBasicData(tf));
            tfString = tfString.replace("{ID}", String.valueOf(tf.getName().hashCode()));
            tfString = tfString.replace("{BUTTONBAR}", "");
            tfString = tfString.replace("{TF_NAME}", tf.getName());
            basicData.append(tfString);
        }

        return basicData.toString();
    }

    private static String getBasicDataEntry(String tfName, String name, Map<String, Map<String, Number>> data)
            throws IOException
    {
        if (data.size() == 0)
        {
            return "";
        }

        String template = FileManagement.readFile(TFPRIO.configs.report.inputStructure.f_basicDataEntry.get());

        template = template.replace("{ANALYSIS-NAME}", name);

        template = template.replace("{FILENAME}", tfName + "_" + name.replace(" ", "-"));

        template = template.replace("{DATA}", getTabularData(String.valueOf((tfName + " " + name).hashCode()), data));

        template = template.replace("{TF-NAME}", tfName);

        template = template.replace("{CSV-DATA}", getCsvData(data));

        return template;
    }

    static String getCsvData(Map<String, Map<String, Number>> data) throws UnsupportedEncodingException
    {
        if (data.size() == 0)
        {
            return "";
        }

        StringBuilder sb_data = new StringBuilder();

        List<String> rows = new ArrayList<>(data.keySet());
        Set<String> columnsSet = new HashSet<>();
        Collections.sort(rows);

        for (String row : rows)
        {
            columnsSet.addAll(data.get(row).keySet());
        }

        List<String> columns = new ArrayList<>(columnsSet);
        Collections.sort(columns);

        if (!(rows.size() == 1 && rows.get(0).isEmpty()))
        {
            sb_data.append(",");
        }
        for (String column : columns)
        {
            sb_data.append(column + ",");
        }
        sb_data.deleteCharAt(sb_data.toString().length() - 1);
        sb_data.append("\n");

        for (String row : rows)
        {
            if (!(rows.size() == 1 && rows.get(0).isEmpty()))
            {
                sb_data.append(row + ", ");
            }

            for (String column : columns)
            {
                if (!column.equals(row) && data.get(row).get(column) != null)
                {
                    sb_data.append(formatter.format(data.get(row).get(column)));
                } else
                {
                    sb_data.append("-");
                }
                sb_data.append(",");
            }
            sb_data.append("\n");
        }

        return URLEncoder.encode(sb_data.toString(), StandardCharsets.UTF_8.toString());
    }

    static String getTabularData(String id, Map<String, Map<String, Number>> data)
    {
        if (data.size() == 0)
        {
            return "";
        }

        StringBuilder sb_data = new StringBuilder();

        List<String> rows = new ArrayList<>(data.keySet());
        Set<String> columnsSet = new HashSet<>();
        Collections.sort(rows);

        for (String row : rows)
        {
            columnsSet.addAll(data.get(row).keySet());
        }

        List<String> columns = new ArrayList<>(columnsSet);
        Collections.sort(columns);

        sb_data.append("<table>");
        sb_data.append("<tr>");
        if (!(rows.size() == 1 && rows.get(0).isEmpty()))
        {
            sb_data.append("<th></th>");
        }
        int i = 0;
        for (String column : columns)
        {
            sb_data.append("<th id='").append(id).append("-col-").append(i).append("'>");
            sb_data.append(column);
            sb_data.append("</th>");
            i++;
        }
        sb_data.append("</tr>");


        i = 0;
        for (String row : rows)
        {
            sb_data.append("<tr>");
            if (!(rows.size() == 1 && rows.get(0).isEmpty()))
            {
                sb_data.append("<th id='").append(id).append("-row-").append(i).append("'>");
                sb_data.append(row);
                sb_data.append("</th>");
            }

            int j = 0;
            for (String column : columns)
            {
                String parameters = "\"" + id + "\", " + j + ", " + i;
                sb_data.append("<td onmouseover='tableMouseOver(").append(parameters)
                        .append(")' onmouseout" + "='tableMouseOut(").append(parameters).append(")'>");
                if (!column.equals(row) && data.get(row).get(column) != null)
                {
                    sb_data.append(formatter.format(data.get(row).get(column)));
                } else
                {
                    sb_data.append("-");
                }
                sb_data.append("</td>");

                j++;
            }

            sb_data.append("</tr>");

            i++;
        }

        sb_data.append("</table>");

        return sb_data.toString();
    }

    static String getButtonBar(TranscriptionFactorGroup tfGroup) throws IOException
    {
        String buttonbar = FileManagement.readFile(TFPRIO.configs.report.inputStructure.f_home_buttonbar.get());

        buttonbar = buttonbar.replace("{VALIDATION}",
                "VALIDATION" + File.separator + tfGroup.getName() + File.separator + tfGroup.getName() + ".html");
        buttonbar = buttonbar.replace("{DISTRIBUTION}",
                "DISTRIBUTION" + File.separator + tfGroup.getName() + File.separator + tfGroup.getName() + ".html");
        buttonbar = buttonbar.replace("{REGRESSION}",
                "REGRESSION" + File.separator + tfGroup.getName() + File.separator + tfGroup.getName() + ".html");

        buttonbar = buttonbar.replace("{HASVALIDATION}", tfGroup.hasValidation() ? "" : "disabled");
        buttonbar = buttonbar.replace("{HASREGRESSION}", tfGroup.hasRegression() ? "" : "disabled");
        buttonbar = buttonbar.replace("{HASDISTRIBUTION}", tfGroup.hasDistribution() ? "" : "disabled");

        return buttonbar;
    }

    static String getFrame(String title, File bodyFile) throws IOException
    {
        String frame = FileManagement.readFile(TFPRIO.configs.report.inputStructure.f_frame.get());

        frame = frame.replace("{TITLE}", title);

        String body = FileManagement.readFile(bodyFile);

        frame = frame.replace("{BODY}", body);

        return frame;
    }

    static String setGeneCardLinks(String text, TranscriptionFactorGroup tfGroup)
    {
        StringBuilder sb_links = new StringBuilder();

        if (tfGroup.getTranscriptionFactors().size() > 1)
        {
            sb_links.append("<div class='dropdown container'>");
            sb_links.append("<button onclick='toggleDropdown(\"geneCardsDropdown\")' id='geneCardsDropdown-dropdown'>");
            sb_links.append("GeneCards");
            sb_links.append("</button>");
            sb_links.append("<div class='geneCards dropdown content' id='geneCardsDropdown-dropdown-content'>");
            for (TranscriptionFactor tf : tfGroup.getTranscriptionFactors())
            {
                sb_links.append("<a class='dropdown geneCards' href='" +
                        TFPRIO.configs.report.genecardsUrl.get().replace("{GENE}", tf.getName()) + "', target" +
                        "='_blank'>");
                sb_links.append(tf.getName());
                sb_links.append("</a>");
            }
            sb_links.append("</div>");
            sb_links.append("</div>");
            sb_links.append("<script>add_dropdown_closing('geneCardsDropdown');</script>");
        } else
        {
            sb_links.append(
                    "<a href='" + TFPRIO.configs.report.genecardsUrl.get().replace("{GENE}", tfGroup.getName()) +
                            "' target" + "='_blank'>");
            sb_links.append("GeneCards");
            sb_links.append("</a>");
        }

        return text.replace("{GENECARD_LINKS}", sb_links.toString());
    }

    static String setBasicData(String text, TranscriptionFactorGroup tfGroup) throws IOException
    {
        text = text.replace("{BASICDATA}", getBasicData(tfGroup));
        return text;
    }

    static String setBasicData(String text, TranscriptionFactor tf) throws IOException
    {
        text = text.replace("{BASICDATA}", getBasicData(tf));
        return text;
    }

    static String generateImageSelector(String id, File sourceDir, List<SelectorTypes> types)
    {
        return generateImageSelector(id, sourceDir, types, false, new JSONObject());
    }

    static String generateImageSelector(String id, File sourceDir, List<SelectorTypes> types, JSONObject data)
    {
        return generateImageSelector(id, sourceDir, types, false, data);
    }

    static String generateImageSelector(String id, File sourceDir, List<SelectorTypes> types,
                                        boolean enableFilterOptions, JSONObject data)
    {
        List<List<String>> options = new ArrayList<>();

        for (SelectorTypes type : types)
        {
            options.add(Report.existingValues.get(type));
        }

        return generateImageSelector(sourceDir, id, options, enableFilterOptions, data);
    }

    static String generateImageSelector(File sourceDir, String id, List<List<String>> options, JSONObject data)
    {
        return generateImageSelector(sourceDir, id, options, false, data);
    }

    static String generateImageSelector(File sourceDir, String id, List<List<String>> options,
                                        boolean enableFilterOptions, JSONObject data)
    {
        StringBuilder sb_imageSelector = new StringBuilder();

        sb_imageSelector.append("<div class='panel' id='{ID}'>");

        if (enableFilterOptions)
        {
            sb_imageSelector.append("<div class='buttonbar'>");
            sb_imageSelector.append(
                    "<button value='MIR' id='" + id + "-filterMiRNA' class='filterOption " + id + "'>miRNA</button>");
            sb_imageSelector.append(
                    "<button value='ENSG' id='" + id + "-filterENSG' class='filterOption " + id + "'>ENSG</button>");
            sb_imageSelector.append(
                    "<button value='RNA' id='" + id + "-filterRNA' class='filterOption " + id + "'>RNA</button>");
            sb_imageSelector.append(
                    "<button value='PSEUDOGENE' id='" + id + "-filterPseudogenes' class='filterOption " + id +
                            "'>Pseudogenes</button>");
            sb_imageSelector.append(
                    "<button value='' id='" + id + "-filterSymbols' class='filterOption " + id + "'>Symbols</button>");
            sb_imageSelector.append("</div>");
        }

        int i = 0;

        for (List<String> levelOptions : options)
        {
            if (levelOptions.size() > 0)
            {
                sb_imageSelector.append("<div class='buttonbar'>");

                for (String entry : levelOptions)
                {
                    String value;
                    if (i == options.size() - 1)
                    {
                        if (entry.matches(".+\\.(gif|jpg|jpeg|tiff|png|svg)$"))
                        {
                            value = entry;
                            entry = entry.substring(0, entry.lastIndexOf("."));
                        } else
                        {
                            value = entry + ".png";
                        }

                    } else
                    {
                        value = entry;
                    }
                    sb_imageSelector.append("<button class='{ID} selector' value='").append(value).append("' id='{ID}-")
                            .append(i).append("-").append(value).append("' " +
                                    ((options.size() == 1 && levelOptions.size() == 1) ? "style='display: none'" : "") +
                                    " onclick='update_selection(this, \"{ID}\")'>");
                    sb_imageSelector.append(entry);
                    sb_imageSelector.append("</button>");
                }
                sb_imageSelector.append("</div>");
            } else
            {
                sb_imageSelector.append("<div class='buttonbar centered'>");

                sb_imageSelector.append(
                        "<button style='width: 15%' class='selector' onclick='downloadActive(\"{ID}\")'>");
                sb_imageSelector.append("Download data");
                sb_imageSelector.append("</button>");

                sb_imageSelector.append(
                        "<button class='selector narrow' id='{ID}-next-option' onclick='move_lowest_level(\"{ID}\", " +
                                "-1, {ID" + "}Combinations)'>");
                sb_imageSelector.append("<");
                sb_imageSelector.append("</button>");

                sb_imageSelector.append("<div class='dropdown container'>");
                sb_imageSelector.append(
                        "<button class='selector dropdown imageSelector' id='{ID}-dropdown' onclick='toggleDropdown" +
                                "(\"{ID}\")'></button>");
                sb_imageSelector.append("<div class=\"dropdown content\" id=\"{ID}-dropdown-content\"></div>");
                sb_imageSelector.append("<script>add_dropdown_closing('{ID}')</script>");

                sb_imageSelector.append("</div>");

                sb_imageSelector.append(
                        "<button class='selector narrow' id='{ID}-previous-option' onclick='move_lowest_level" +
                                "(\"{ID}\", 1)'>");
                sb_imageSelector.append(">");
                sb_imageSelector.append("</button>");

                sb_imageSelector.append(
                        "<button style='width: 15%' class='selector' onclick='openImageInTab(\"{ID}-image\")'>");
                sb_imageSelector.append("Open in new tab");
                sb_imageSelector.append("</button>");

                sb_imageSelector.append("</div>");
            }
            i++;
        }

        if (options.get(options.size() - 1).size() > 0)
        {
            {
                sb_imageSelector.append("<div class='buttonbar centered'>");

                sb_imageSelector.append("<button class='selector' onclick='openModal(\"{ID}-modal\")'>");
                sb_imageSelector.append("Zoom in");
                sb_imageSelector.append("</button>");

                sb_imageSelector.append("<button class='selector' onclick='downloadActive(\"{ID}\")'>");
                sb_imageSelector.append("Download data");
                sb_imageSelector.append("</button>");

                sb_imageSelector.append("<button class='selector' onclick='openImageInTab(\"{ID}-image\")'>");
                sb_imageSelector.append("Open in new tab");
                sb_imageSelector.append("</button>");

                sb_imageSelector.append("</div>");
            }
        }

        {
            sb_imageSelector.append("<img class='modal source' onclick='openModal(\"{ID}-modal\")'  id='{ID}-image'>");
        }

        {
            sb_imageSelector.append("<div id='{ID}-modal' class='modal background'>");

            sb_imageSelector.append("<div class='buttonbar centered'>");
            sb_imageSelector.append("<div id='{ID}-modal-caption' class='modal caption'></div>");
            sb_imageSelector.append(
                    "<button class=\"modal button close\" onclick=\"closeModal('{ID}-modal')" + "\">&times;</button>");
            sb_imageSelector.append("</div>");

            sb_imageSelector.append("<div class='buttonbar'>");
            sb_imageSelector.append(
                    "<button class=\"modal button\" id=\"{ID}-modal-leftarrow\"onclick=\"move_lowest_level('{ID}', -1)\">\n" +
                            "            &lt;</button>");
            sb_imageSelector.append("<img class=\"modal content\" id=\"{ID}-modal-image\">");
            sb_imageSelector.append(
                    "<button class=\"modal button\" id=\"{ID}-modal-rightarrow\" onclick=\"move_lowest_level('{ID}', 1)\">\n" +
                            "            &gt;</button>");
            sb_imageSelector.append("</div>");

            sb_imageSelector.append("</div>");
        }

        sb_imageSelector.append("</div>");

        sb_imageSelector.append("<script>var {ID}CombinationsUnfiltered = {COMBINATIONS};</script>");
        sb_imageSelector.append("<script>var {ID}Combinations = {COMBINATIONS};</script>");
        sb_imageSelector.append("<script>var {ID}DataCombinations = " + data.toString() + ";</SCRIPT>");
        sb_imageSelector.append("<script>init_filterOptions(\"{ID}\")</script>");
        sb_imageSelector.append("<script>init_selection(\"{ID}\")</script>");

        String imageSelector = sb_imageSelector.toString().replace("{ID}", id);

        return imageSelector.replace("{COMBINATIONS}", getFileStructureAsString(sourceDir));
    }

    private static String getFileStructureAsString(File sourceDir)
    {
        boolean onlyFiles = true;
        boolean onlyDirs = true;
        StringBuilder sb_result = new StringBuilder();

        if (sourceDir.listFiles() == null)
        {
            return "";
        }

        for (File entry : Objects.requireNonNull(sourceDir.listFiles()))
        {
            onlyFiles = onlyFiles && entry.isFile();
            onlyDirs = onlyDirs && entry.isDirectory();
        }

        if (onlyFiles)
        {
            ArrayList<String> fileNames = new ArrayList<>();

            for (File entry : Objects.requireNonNull(sourceDir.listFiles()))
            {
                String name = entry.getName();
                if (!name.equals("ALL.png"))
                {
                    fileNames.add(name);
                }
            }

            Collections.sort(fileNames);
            fileNames.sort(new IntegerStringComparator());

            sb_result.append("[");

            for (String fileName : fileNames)
            {
                sb_result.append("\"" + fileName + "\",");
            }

            sb_result.deleteCharAt(sb_result.length() - 1);
            sb_result.append("]");
        }

        if (onlyDirs)
        {
            sb_result.append("{");

            for (File entry : Objects.requireNonNull(sourceDir.listFiles()))
            {
                sb_result.append("\"" + entry.getName() + "\":");
                sb_result.append(getFileStructureAsString(entry));
                sb_result.append(",");
            }

            sb_result.deleteCharAt(sb_result.length() - 1);
            sb_result.append("}");
        }

        return sb_result.toString();
    }
}
