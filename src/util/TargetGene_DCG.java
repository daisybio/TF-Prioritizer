package util;

import java.lang.annotation.Target;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class TargetGene_DCG implements Comparable {
    public String target_gene = "";
    public ArrayList<Double> affinity_values = new ArrayList<>();
    public double mean_affinity_value = 0.0;

    public void set_mean_affnity_value() {
        mean_affinity_value = calculateAverage(affinity_values);
    }

    private double calculateAverage(List<Double> marks) {
        Double sum = 0.0;
        if (!marks.isEmpty()) {
            for (Double mark : marks) {
                sum += mark;
            }
            return sum.doubleValue() / marks.size();
        }
        return sum;
    }

    @Override public int compareTo(Object o) {
        TargetGene_DCG other = (TargetGene_DCG) o;

        if (other.mean_affinity_value > this.mean_affinity_value) {
            return 1;
        }
        if (other.mean_affinity_value < this.mean_affinity_value) {
            return -1;
        }

        return 0;
    }
}
