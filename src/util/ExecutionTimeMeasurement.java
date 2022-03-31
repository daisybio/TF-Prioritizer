package util;

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
}
