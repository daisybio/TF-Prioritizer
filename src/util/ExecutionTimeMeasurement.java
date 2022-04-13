package util;

import java.util.concurrent.TimeUnit;

public class ExecutionTimeMeasurement
{
    private final long startTimeMillis;
    private long stopTimeMillis;
    private boolean running;

    public ExecutionTimeMeasurement()
    {
        startTimeMillis = System.currentTimeMillis();
        running = true;
    }

    public void stop()
    {
        stopTimeMillis = System.currentTimeMillis();
        running = false;
    }

    public double stopAndGetDeltaSeconds()
    {
        stop();
        return getDeltaSeconds();
    }

    public double getDeltaSeconds()
    {
        return (double) getDeltaMillis() / 1e3;
    }

    public long getDeltaMillis()
    {
        long compareTime = running ? System.currentTimeMillis() : stopTimeMillis;
        return compareTime - startTimeMillis;
    }

    public String stopAndGetDeltaFormatted()
    {
        stop();
        return getDeltaFormatted();
    }

    public String getDeltaFormatted()
    {
        long delta = getDeltaMillis();
        return formatMillis(delta);
    }

    public static String formatMillis(long millis)
    {
        String format = "%02d %s ";

        long hours = TimeUnit.MILLISECONDS.toHours(millis);

        long minutes = TimeUnit.MILLISECONDS.toMinutes(millis) - TimeUnit.HOURS.toMinutes(hours);

        long seconds = TimeUnit.MILLISECONDS.toSeconds(millis) - TimeUnit.MINUTES.toSeconds(minutes) -
                TimeUnit.HOURS.toSeconds(hours);

        StringBuilder sb_output = new StringBuilder();
        if (hours > 0)
        {
            sb_output.append(String.format(format, hours, "h"));
        }
        if (minutes > 0 || hours > 0)
        {
            sb_output.append(String.format(format, minutes, "min"));
        }

        sb_output.append(String.format(format, seconds, "sec"));

        return sb_output.toString();
    }
}
