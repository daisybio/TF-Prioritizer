package util.Report;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

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

        template = template.replace("{LOG2FC}", log2fcTemplate);

        template = template.replace("{TPM}", getBasicDataEntry("TPM", transcriptionFactor.getTpm()));

        template = template.replace("{NORMEX}", getBasicDataEntry("Norm. expression", transcriptionFactor.getNormex()));

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
            tfString = tfString.replace("{GENEID}", tf.getGeneID());
            tfString = tfString.replace("{BUTTONBAR}", "");
            tfString = tfString.replace("{TF_NAME}", tf.getName());
            basicData.append(tfString);
        }

        return basicData.toString();
    }

    private static String getBasicDataEntry(String name, Map<String, Number> data) throws IOException
    {
        if (data.size() == 0)
        {
            return "";
        }

        String template = FileManagement.loadFile(Report.options_intern.path_to_COM2POSE + File.separator +
                Report.options_intern.f_report_resources_basicdata_entry_html);

        template = template.replace("{NAME}", name);

        StringBuilder sb_data = new StringBuilder();

        sb_data.append("<div class='keyvaluepaircontainer'>");

        for (Map.Entry<String, Number> kvPair : data.entrySet())
        {
            sb_data.append("<div class=\"keyvaluepair\"><h4>").append(kvPair.getKey()).append("</h4><p>")
                    .append(formatter.format(kvPair.getValue())).append("</p></div>");
        }

        sb_data.append("</div>");

        template = template.replace("{DATA}", sb_data.toString());

        return template;
    }

    static String getTabularData(String id, Map<String, Map<String, Number>> data)
    {
        if (data.size() == 0)
        {
            return "";
        }

        StringBuilder sb_data = new StringBuilder();

        List<String> columns = new ArrayList<>(data.keySet());
        Set<String> rowsSet = new HashSet<>();
        Collections.sort(columns);

        sb_data.append("<table>");
        sb_data.append("<tr>");
        sb_data.append("<th></th>");
        int i = 0;
        for (String column : columns)
        {
            sb_data.append("<th id='").append(id).append("-col-").append(i).append("'>");
            sb_data.append(column);
            sb_data.append("</th>");
            i++;

            rowsSet.addAll(data.get(column).keySet());
        }
        sb_data.append("</tr>");

        List<String> rows = new ArrayList<>(rowsSet);
        Collections.sort(rows);

        i = 0;
        for (String row : rows)
        {
            sb_data.append("<tr>");

            sb_data.append("<th id='").append(id).append("-row-").append(i).append("'>");
            sb_data.append(row);
            sb_data.append("</th>");

            int j = 0;
            for (String column : columns)
            {
                String parameters = "\"" + id + "\", " + j + ", " + i;
                sb_data.append("<td onmouseover='tableMouseOver(").append(parameters)
                        .append(")' onmouseout" + "='tableMouseOut(").append(parameters).append(")'>");
                if (!column.equals(row) && data.get(column).get(row) != null)
                {
                    sb_data.append(formatter.format(data.get(column).get(row)));
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

        for (TranscriptionFactor tf : tfGroup.getTranscriptionFactors())
        {
            String command =
                    "window.open('" + Report.options_intern.link_report_genecards.replace("{GENE}", tf.getName()) +
                            "');\n";
            sb_links.append(command);
        }

        return text.replace("{GENECARD_BUTTON_ACTION}", sb_links.toString());
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
        StringBuilder sb_imageSelector = new StringBuilder();

        sb_imageSelector.append("<div class='panel' id='{ID}'>");

        int i = 0;

        for (SelectorTypes type : types)
        {
            if (type != SelectorTypes.EMPTY_DROPDOWN)
            {
                sb_imageSelector.append("<div class='buttonbar'>");

                for (String entry : Report.existingValues.get(type))
                {
                    String value;
                    if (i == types.size() - 1 && !entry.matches(".+\\.(gif|jpg|jpeg|tiff|png)$"))
                    {
                        value = entry + ".png";
                    } else
                    {
                        value = entry;
                    }
                    sb_imageSelector.append("<button class='{ID} selector' value='").append(value).append("' id='{ID}-")
                            .append(i).append("-").append(value)
                            .append("' onclick='update_selection(this, \"{ID}\", {ID}Combinations)'>");
                    sb_imageSelector.append(entry);
                    sb_imageSelector.append("</button>");
                }
                sb_imageSelector.append("</div>");
            } else
            {
                sb_imageSelector.append("<div class='buttonbar centered'>");

                sb_imageSelector.append("<button onclick='openModal(\"{ID}-modal\")'>");
                sb_imageSelector.append("Zoom in");
                sb_imageSelector.append("</button>");

                sb_imageSelector.append(
                        "<button class='narrow' onclick='move_lowest_level(\"{ID}\", -1, " + "{ID}Combinations)'>");
                sb_imageSelector.append("<");
                sb_imageSelector.append("</button>");

                sb_imageSelector.append("<div class='dropdown container'>");
                sb_imageSelector.append(
                        "<button class='dropdown button' id='{ID}-dropdown' onclick='toggleDropdown(\"{ID}\")" +
                                "'></button>");
                sb_imageSelector.append("<div class=\"dropdown content\" id=\"{ID}-dropdown-content\"></div>");
                sb_imageSelector.append("<script>add_dropdown_closing('{ID}')</script>");

                sb_imageSelector.append("</div>");

                sb_imageSelector.append(
                        "<button class='narrow' onclick='move_lowest_level(\"{ID}\", 1, " + "{ID}Combinations)'>");
                sb_imageSelector.append(">");
                sb_imageSelector.append("</button>");

                sb_imageSelector.append("<button onclick='openImageInTab(\"{ID}-image\")'>");
                sb_imageSelector.append("Open in new tab");
                sb_imageSelector.append("</button>");

                sb_imageSelector.append("</div>");
            }
            i++;
        }

        if (types.get(types.size() - 1) != SelectorTypes.EMPTY_DROPDOWN)
        {
            {
                sb_imageSelector.append("<div class='buttonbar centered'>");

                sb_imageSelector.append("<button onclick='openModal(\"{ID}-modal\")'>");
                sb_imageSelector.append("Zoom in");
                sb_imageSelector.append("</button>");

                sb_imageSelector.append("<button onclick='openImageInTab(\"{ID}-image\")'>");
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

        for (File entry : Objects.requireNonNull(sourceDir.listFiles()))
        {
            onlyFiles = onlyFiles && entry.isFile();
            onlyDirs = onlyDirs && entry.isDirectory();
        }

        if (onlyFiles)
        {
            sb_result.append("[");
            for (File entry : Objects.requireNonNull(sourceDir.listFiles()))
            {
                sb_result.append("\"" + entry.getName() + "\",");
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
