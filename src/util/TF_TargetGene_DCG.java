package util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class TF_TargetGene_DCG {
    public String TF = "";
    public HashMap<String, TargetGene_DCG> target_gene_affinity_values = new HashMap<>();

    public ArrayList<TargetGene_DCG> get_ordered_target_gene_list() {
        ArrayList<TargetGene_DCG> result = new ArrayList<>();

        for (String key_target_gene : target_gene_affinity_values.keySet()) {
            TargetGene_DCG dcg = target_gene_affinity_values.get(key_target_gene);
            dcg.set_mean_affnity_value();
            result.add(dcg);
        }

        Collections.sort(result);


        return result;
    }
}
