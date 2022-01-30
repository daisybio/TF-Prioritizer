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
        return getButtonBar(tfGroup.getName(), tfGroup.hasValidation(), tfGroup.hasDistribution(),
                tfGroup.hasRegression());
    }

    private static String getButtonBar(String name, boolean hasValidation, boolean hasDistribution,
                                       boolean hasRegression) throws IOException
    {
        String buttonbar = FileManagement.loadFile(Report.options_intern.path_to_COM2POSE + File.separator +
                Report.options_intern.f_report_resources_home_buttonbar_html);

        buttonbar = buttonbar.replace("{VALIDATION}",
                "VALIDATION" + File.separator + name + File.separator + name + ".html");
        buttonbar = buttonbar.replace("{DISTRIBUTION}",
                "DISTRIBUTION" + File.separator + name + File.separator + name + ".html");
        buttonbar = buttonbar.replace("{REGRESSION}",
                "REGRESSION" + File.separator + name + File.separator + name + ".html");

        buttonbar = buttonbar.replace("{HASVALIDATION}", hasValidation ? "" : "disabled");
        buttonbar = buttonbar.replace("{HASREGRESSION}", hasRegression ? "" : "disabled");
        buttonbar = buttonbar.replace("{HASDISTRIBUTION}", hasDistribution ? "" : "disabled");

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

    static String generateThreeLevelImageSelector(String name, File sourceDir, File targetDir, boolean keepFileNameAsIs)
            throws IOException
    {
        return generateThreeLevelImageSelector(name, sourceDir, targetDir, new ArrayList<>(), keepFileNameAsIs, true);
    }

    static String generateThreeLevelImageSelector(String name, File sourceDir, File targetDir,
                                                  ArrayList<String> specialRemovables, boolean keepFileNameAsIs,
                                                  boolean compression) throws IOException
    {
        String suffix = ".png";
        HashMap<String, HashMap<String, ArrayList<String>>> combinations = new HashMap<>();

        for (File group_dir : Objects.requireNonNull(sourceDir.listFiles()))
        {
            if (group_dir.isFile())
            {
                continue;
            }

            String group = group_dir.getName();

            if (group.equals("A_SESSIONS"))
            {
                continue;
            }

            combinations.put(group, new HashMap<>());

            for (File subgroup_dir : Objects.requireNonNull(group_dir.listFiles()))
            {
                if (subgroup_dir.isFile())
                {
                    continue;
                }

                String subgroup = subgroup_dir.getName();

                if (subgroup.split("_").length > 1)
                {
                    subgroup = subgroup.split("_")[1];
                }

                combinations.get(group).put(subgroup, new ArrayList<>());

                for (File image_file : Objects.requireNonNull(subgroup_dir.listFiles()))
                {
                    ArrayList<String> removables = new ArrayList<>(List.of(group, subgroup, "threshold"));
                    removables.addAll(specialRemovables);

                    if (image_file.isDirectory() || !image_file.getName().endsWith(suffix))
                    {
                        continue;
                    }

                    String relevantFileName = image_file.getName();
                    relevantFileName = relevantFileName.replace(suffix, "");

                    if (!keepFileNameAsIs)
                    {
                        for (String entry : removables)
                        {
                            relevantFileName = relevantFileName.replace(entry, "");
                        }


                        relevantFileName = relevantFileName.replaceAll("_+", "_");

                        while (relevantFileName.endsWith("_"))
                        {
                            relevantFileName = relevantFileName.substring(0, relevantFileName.length() - 1);
                        }

                        while (relevantFileName.startsWith("_"))
                        {
                            relevantFileName = relevantFileName.substring(1);
                        }
                    }

                    if (targetDir != null)
                    {
                        FileManagement.copyFile(image_file, new File(
                                targetDir.getAbsolutePath() + File.separator + group + File.separator + subgroup +
                                        File.separator + relevantFileName + suffix), compression);
                    }

                    combinations.get(group).get(subgroup).add(relevantFileName);
                }
            }
        }

        String three_level_image_selector =
                FileManagement.loadFile(Report.options_intern.f_report_resources_three_level_image_selector_html);

        Set<String> groups = combinations.keySet();
        Set<String> subgroups = new HashSet<>();

        StringBuilder sb_groups = new StringBuilder();

        for (String group : groups)
        {
            subgroups.addAll(combinations.get(group).keySet());

            sb_groups.append("<button class=\"{ID} group-selector\" onclick=\"select_group" + "('{ID}', this, " +
                    "{ID}Combinations)\" " + "value=\"" + group + "\">" + group + "</button>");
        }
        three_level_image_selector = three_level_image_selector.replace("{GROUPS}", sb_groups.toString());

        StringBuilder sb_subgroups = new StringBuilder();
        for (String subgroup : subgroups)
        {
            sb_subgroups.append(
                    "<button class=\"{ID} subgroup-selector\" onclick=\"select_subgroup" + "('{ID}', " + "this, " +
                            "{ID}Combinations)\" " + "value=\"" + subgroup + "\">" + subgroup + "</button>");
        }

        three_level_image_selector = three_level_image_selector.replace("{SUBGROUPS}", sb_subgroups.toString());

        String json;
        {
            HashMap<String, HashMap<String, String>> lv1 = new HashMap<>();
            HashMap<String, String> lv2 = new HashMap<>();

            for (String hm : combinations.keySet())
            {
                lv1.put(hm, new HashMap<>());
                for (String group : combinations.get(hm).keySet())
                {
                    StringBuilder sb_genes = new StringBuilder("[");
                    ArrayList<String> genes = combinations.get(hm).get(group);
                    genes.sort(new StringComparator());

                    for (String gene : genes)
                    {
                        sb_genes.append("\"");
                        sb_genes.append(gene);
                        sb_genes.append("\",");
                    }
                    sb_genes.setLength(sb_genes.length() - 1);
                    sb_genes.append("]");
                    lv1.get(hm).put(group, sb_genes.toString());
                }
                lv2.put(hm, mapToJson(lv1.get(hm)));
            }

            json = mapToJson(lv2);
        }

        three_level_image_selector = three_level_image_selector.replace("{ID}", name);

        three_level_image_selector = three_level_image_selector.replace("{COMBINATIONS}", json);

        return three_level_image_selector;
    }

    static String generateTwoLevelImageSelector(String name, File sourceDir, File targetDir, boolean compression)
            throws IOException
    {
        String suffix = ".png";
        HashMap<String, ArrayList<String>> combinations = new HashMap<>();

        for (File group_dir : Objects.requireNonNull(sourceDir.listFiles()))
        {
            if (group_dir.isFile())
            {
                continue;
            }

            String group = group_dir.getName();

            combinations.put(group, new ArrayList<>());


            for (File image_file : Objects.requireNonNull(group_dir.listFiles()))
            {
                if (image_file.isDirectory() || !image_file.getName().endsWith(suffix))
                {
                    continue;
                }

                String relevantFileName = image_file.getName();
                relevantFileName = relevantFileName.replace(suffix, "");

                if (targetDir != null)
                {
                    FileManagement.copyFile(image_file, new File(
                            targetDir.getAbsolutePath() + File.separator + group + File.separator + relevantFileName +
                                    suffix), compression);
                }

                combinations.get(group).add(relevantFileName);
            }

        }

        String two_level_image_selector =
                FileManagement.loadFile(Report.options_intern.f_report_resources_two_level_image_selector_html);

        Set<String> groups = combinations.keySet();
        Set<String> subgroups = new HashSet<>();

        {
            StringBuilder sb_groups = new StringBuilder();

            for (String group : groups)
            {
                subgroups.addAll(combinations.get(group));

                sb_groups.append("<button class=\"{ID} group-selector\" onclick=\"select_group" + "('{ID}', this, " +
                        "{ID}Combinations)\" " + "value=\"" + group + "\">" + group + "</button>");
            }
            two_level_image_selector = two_level_image_selector.replace("{GROUPS}", sb_groups.toString());
        }

        {
            StringBuilder sb_subgroups = new StringBuilder();
            for (String subgroup : subgroups)
            {
                sb_subgroups.append(
                        "<button class=\"{ID} subgroup-selector\" onclick=\"select_subgroup" + "('{ID}', " + "this, " +
                                "{ID}Combinations)\" " + "value=\"" + subgroup + "\">" + subgroup + "</button>");
            }

            two_level_image_selector = two_level_image_selector.replace("{SUBGROUPS}", sb_subgroups.toString());
        }

        String json;
        {
            HashMap<String, String> lv1 = new HashMap<>();

            for (String group : combinations.keySet())
            {
                StringBuilder sb_subgroups = new StringBuilder("[");
                for (String subgroup : combinations.get(group))
                {
                    sb_subgroups.append("\"");
                    sb_subgroups.append(subgroup);
                    sb_subgroups.append("\",");
                }
                sb_subgroups.setLength(sb_subgroups.length() - 1);
                sb_subgroups.append("]");
                lv1.put(group, sb_subgroups.toString());
            }
            json = mapToJson(lv1);
        }

        two_level_image_selector = two_level_image_selector.replace("{ID}", name);

        two_level_image_selector = two_level_image_selector.replace("{COMBINATIONS}", json);

        return two_level_image_selector;
    }

    private static String mapToJson(Map<String, String> map)
    {
        StringBuilder sb_output = new StringBuilder("{");

        for (Map.Entry<String, String> entry : map.entrySet())
        {
            boolean valueIsJson = (entry.getValue().startsWith("{") && entry.getValue().endsWith("}")) ||
                    (entry.getValue().startsWith("[") && entry.getValue().endsWith("]"));
            sb_output.append("\"");
            sb_output.append(entry.getKey());
            sb_output.append("\":");
            if (!valueIsJson)
            {
                sb_output.append("\"");
            }
            sb_output.append(entry.getValue());
            if (!valueIsJson)
            {
                sb_output.append("\"");
            }
            sb_output.append(",");
        }
        sb_output.setLength(sb_output.length() - 1);
        sb_output.append("}");

        return sb_output.toString();
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
}
