package org.exbio.tfprio.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;

public class Helpers {
    public static String getDatasetVersion(String species) throws IOException {
        URL url = new URL("http://www.ensembl.org/biomart/martservice?type=datasets&mart=ENSEMBL_MART_ENSEMBL");
        HttpURLConnection conn = (HttpURLConnection) url.openConnection();
        conn.setRequestMethod("GET");
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(conn.getInputStream()))) {
            for (String line; (line = reader.readLine()) != null; ) {
                if (line.contains(species)) {
                    String[] split = line.split("\t");
                    return split[4];
                }
            }
        }
        return null;
    }

}
