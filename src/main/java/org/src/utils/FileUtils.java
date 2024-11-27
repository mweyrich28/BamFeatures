package org.src.utils;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

public class FileUtils {
    public static String parseGTFAttributes(String[] attributeEntries, String attributeName) {
        for (int i = 0; i < attributeEntries.length; i++) {
            String trimmedEntry = attributeEntries[i].trim();

            int posSpace = trimmedEntry.indexOf(' ');
            String attributeKey = trimmedEntry.substring(0, posSpace);

            String attributeVal;

            if (trimmedEntry.endsWith("\"")) {
                attributeVal = trimmedEntry.substring(posSpace + 2, trimmedEntry.length() - 1);
            }
            // there are entries like exon_number, which are not encapsulated in quotes
            else {
                attributeVal = trimmedEntry.substring(posSpace + 1, trimmedEntry.length() - 1);
            }

            if (attributeKey.equals(attributeName)) {
                return attributeVal;
            }
        }
        return null;
    }

    public static void writeFile(String fileName, ArrayList<String> content) throws IOException {
        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(fileName)));

        for (int i = 0; i < content.size(); i++) {
            bw.write(content.get(i));
            bw.newLine();
        }
        bw.close();
    }


    public static ArrayList<String> readLines(File r) throws IOException {
        ArrayList<String> fileLines = new ArrayList<>();
        BufferedReader buff = new BufferedReader(new FileReader(r));
        String line;
        while((line = buff.readLine()) != null){
            fileLines.add(line);
        }
        return fileLines;
    }
}



