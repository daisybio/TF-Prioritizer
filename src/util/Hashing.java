package util;

import util.Comparators.FileComparator;
import util.Configs.ConfigTypes.AbstractConfig;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.security.DigestInputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.*;

public class Hashing
{
    public static String hash(String text)
    {
        MessageDigest digest = getMessageDigest();
        byte[] hashBytes = digest.digest(text.getBytes(StandardCharsets.UTF_8));
        return bytesToString(hashBytes);
    }

    private static String hashHashes(List<String> hashes)
    {
        return hash(String.join("", hashes));
    }

    public static String hashConfigs(Iterable<AbstractConfig<?>> configs)
    {
        ArrayList<String> hashes = new ArrayList<>();

        for (AbstractConfig<?> config : configs)
        {
            String hash = hash(config.toString());
            hashes.add(hash);
        }
        Collections.sort(hashes);
        return hashHashes(hashes);
    }

    public static String hashFiles(List<File> files) throws IOException
    {
        return hashFiles(true, files);
    }

    public static String hashFiles(boolean includeHiddenFiles, List<File> files) throws IOException
    {
        files.sort(new FileComparator());
        Vector<FileInputStream> fileStreams = new Vector<>();

        for (File file : files)
        {
            collectInputStreams(file, fileStreams, includeHiddenFiles);
        }

        MessageDigest md = getMessageDigest();
        try (SequenceInputStream seqStream = new SequenceInputStream(fileStreams.elements());
             DigestInputStream stream = new DigestInputStream(seqStream, md))
        {
            byte[] buffer = new byte[1024 * 8];
            while (stream.read(buffer) != -1)
            {
            }
        }

        return bytesToString(md.digest());
    }

    private static void collectInputStreams(File file, List<FileInputStream> foundStreams, boolean includeHiddenFiles)
            throws FileNotFoundException
    {


        if (file.isFile() && (includeHiddenFiles || !file.isHidden()))
        {
            foundStreams.add(new FileInputStream(file));
        } else if (file.isDirectory())
        {
            File[] fileList = file.listFiles();
            assert fileList != null;
            Arrays.sort(fileList, new FileComparator());

            for (File f : fileList)
            {
                if (!includeHiddenFiles && f.isHidden())
                {
                    continue;
                }
                collectInputStreams(f, foundStreams, includeHiddenFiles);
            }
        }
    }

    private static MessageDigest getMessageDigest()
    {
        MessageDigest digest = null;
        try
        {
            digest = MessageDigest.getInstance("MD5");
        } catch (NoSuchAlgorithmException ignore)
        {
        }
        assert digest != null;

        return digest;
    }

    private static String bytesToString(byte[] bytes)
    {
        byte[] HEX_ARRAY = "0123456789ABCDEF".getBytes(StandardCharsets.UTF_8);
        byte[] hexChars = new byte[bytes.length * 2];

        for (int j = 0; j < bytes.length; j++)
        {
            int v = bytes[j] & 0xFF;
            hexChars[j * 2] = HEX_ARRAY[v >>> 4];
            hexChars[j * 2 + 1] = HEX_ARRAY[v & 0x0F];
        }
        return new String(hexChars, StandardCharsets.UTF_8);
    }
}
