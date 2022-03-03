package util.Report;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

import util.FileManagement;

public class StructureElements
{
    private static final DecimalFormat formatter = new DecimalFormat("0.###");

    static String getBasicData(TranscriptionFactor transcriptionFactor) throws IOException
    {
        String template = FileManagement.loadFile(Report.options_intern.path_to_COM2POSE + File.separator +
                Report.options_intern.f_report_resources_basicdata_html);


        String log2fcTemplate = FileManagement.loadFile(Report.options_intern.path_to_COM2POSE + File.separator +
                Report.options_intern.f_report_resources_basicdata_entry_html);

        log2fcTemplate = log2fcTemplate.replace("{NAME}", "LOG2FC");

        log2fcTemplate = log2fcTemplate.replace("{DATA}",
                getTabularData(transcriptionFactor.getName(), transcriptionFactor.getLog2fc()));

        log2fcTemplate = log2fcTemplate.replace("{INFO-ID}", transcriptionFactor.getName() + "-" + "log2fc");

        template = template.replace("{GENEID}", FileManagement.loadFile(
                Report.options_intern.path_to_COM2POSE + File.separator + "ext" + File.separator + "REPORT" +
                        File.separator + "BASICDATA" + File.separator + "GENEID.html"));
        template = template.replace("{GENE-ID}", transcriptionFactor.getGeneID());
        template = template.replace("{GENECARDS-LINK}",
                Report.options_intern.link_report_genecards.replace("{GENE}", transcriptionFactor.getName()));

        template = template.replace("{LOG2FC}",
                getBasicDataEntry(transcriptionFactor.getName(), "LOG2FC", transcriptionFactor.getLog2fc()));

        template = template.replace("{TPM}", getBasicDataEntry(transcriptionFactor.getName(), "TPM", new HashMap<>()
        {{
            put("", transcriptionFactor.getTpm());
        }}));

        template = template.replace("{NORMEX}",
                getBasicDataEntry(transcriptionFactor.getName(), "Norm. expression", new HashMap<>()
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
        String tfTemplate = FileManagement.loadFile(Report.options_intern.path_to_COM2POSE + File.separator +
                Report.options_intern.f_report_resources_home_tf_html);

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

        String template = FileManagement.loadFile(Report.options_intern.path_to_COM2POSE + File.separator +
                Report.options_intern.f_report_resources_basicdata_entry_html);

        template = template.replace("{NAME}", name);

        template = template.replace("{DATA}", getTabularData(String.valueOf((tfName + " " + name).hashCode()), data));

        template = template.replace("{INFO-ID}", tfName + "-" + name);

        return template;
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
        String buttonbar = FileManagement.loadFile(Report.options_intern.path_to_COM2POSE + File.separator +
                Report.options_intern.f_report_resources_home_buttonbar_html);

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

    static String getFrame(String title, String bodyPath) throws IOException
    {
        String frame = FileManagement.loadFile(Report.options_intern.path_to_COM2POSE + File.separator +
                Report.options_intern.f_report_resources_frame_html);

        frame = frame.replace("{TITLE}", title);

        String body = FileManagement.loadFile(bodyPath);

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
                        Report.options_intern.link_report_genecards.replace("{GENE}", tf.getName()) + "', target" +
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
                    "<a href='" + Report.options_intern.link_report_genecards.replace("{GENE}", tfGroup.getName()) +
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
        List<List<String>> options = new ArrayList<>();

        for (SelectorTypes type : types)
        {
            options.add(Report.existingValues.get(type));
        }

        return generateImageSelector(sourceDir, id, options);
    }

    static String generateImageSelector(File sourceDir, String id, List<List<String>> options)
    {
        StringBuilder sb_imageSelector = new StringBuilder();

        sb_imageSelector.append("<div class='panel' id='{ID}'>");

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
                                    " onclick='update_selection(this, \"{ID}\", {ID}Combinations)'>");
                    sb_imageSelector.append(entry);
                    sb_imageSelector.append("</button>");
                }
                sb_imageSelector.append("</div>");
            } else
            {
                sb_imageSelector.append("<div class='buttonbar centered'>");

                sb_imageSelector.append("<button class='selector' onclick='openModal(\"{ID}-modal\")'>");
                sb_imageSelector.append("Zoom in");
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
                                "(\"{ID}\", 1, " + "{ID}Combinations)'>");
                sb_imageSelector.append(">");
                sb_imageSelector.append("</button>");

                sb_imageSelector.append("<button class='selector' onclick='openImageInTab(\"{ID}-image\")'>");
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
                    "<button class=\"modal button\" id=\"{ID}-modal-leftarrow\"onclick=\"move_lowest_level('{ID}', -1, {ID}Combinations)\">\n" +
                            "            &lt;</button>");
            sb_imageSelector.append("<img class=\"modal content\" id=\"{ID}-modal-image\">");
            sb_imageSelector.append(
                    "<button class=\"modal button\" id=\"{ID}-modal-rightarrow\" onclick=\"move_lowest_level('{ID}', 1, {ID}Combinations)\">\n" +
                            "            &gt;</button>");
            sb_imageSelector.append("</div>");

            sb_imageSelector.append("</div>");
        }

        sb_imageSelector.append("</div>");

        sb_imageSelector.append("<script>let {ID}Combinations = {COMBINATIONS};</script>");
        sb_imageSelector.append("<script>init_selection(\"{ID}\", {ID}Combinations)</script>");

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
            fileNames.sort(new StringComparator());

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

    static class StringComparator implements Comparator<String>
    {
        @Override public int compare(String a, String b)
        {
            return prefixNum(a) - prefixNum(b);
        }

        private int prefixNum(String a)
        {
            if (a.matches("[0-9]+_.*"))
            {
                return Integer.parseInt(a.split("_")[0]);
            } else
            {
                return 0;
            }
        }
    }
}
