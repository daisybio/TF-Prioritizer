package org.exbio.tfprio.util;

import java.util.Comparator;

public class NumericStringComparator implements Comparator<String> {
    public int compare(String a, String b) {
        if (isAnyInteger(a) && isAnyInteger(b)) {
            return stringToInteger(a) - stringToInteger(b);
        }
        return a.compareTo(b);
    }

    private int stringToInteger(String str) {
        if (isInteger(str)) {
            return Integer.parseInt(str);
        }
        if (isIntegerPrefix(str)) {
            return Integer.parseInt(str.split("_")[0]);
        }
        return 0;
    }

    private boolean isAnyInteger(String str) {
        return isInteger(str) || isIntegerPrefix(str);
    }

    private boolean isInteger(String str) {
        return str.matches("[0-9]+");
    }

    private boolean isIntegerPrefix(String str) {
        return str.matches("[0-9]+_.*");
    }
}
