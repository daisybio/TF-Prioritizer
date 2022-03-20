package util;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

import static util.FileManagement.readFile;
import static util.FileManagement.writeFile;

public class Hashing
{
    public static String hash(String text)
    {
        MessageDigest digest = null;
        try
        {
            digest = MessageDigest.getInstance("SHA-256");
        } catch (NoSuchAlgorithmException ignore)
        {
        }
        assert digest != null;
        byte[] hashBytes = digest.digest(text.getBytes(StandardCharsets.UTF_8));
        return new String(hashBytes, StandardCharsets.UTF_8);
    }

    public static String hash(File file) throws IOException
    {
        String content = readFile(file);
        return hash(content);
    }

    public static void hashFileToFile(File source, File target) throws IOException
    {
        String hash = hash(source);
        writeFile(target, hash);
    }

    public static boolean checkHashInFileForFile(File dataFile, File hashFile) throws IOException
    {
        String necessaryHash = hash(dataFile);
        String presentHash = readFile(hashFile);
        return presentHash.equals(necessaryHash);
    }
}
