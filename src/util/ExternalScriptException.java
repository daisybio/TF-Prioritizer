package util;

import java.io.IOException;

public class ExternalScriptException extends IOException
{
    private final int returnCode;

    public ExternalScriptException(int returnCode, String message)
    {
        super(message);
        this.returnCode = returnCode;
    }

    public int getReturnCode()
    {
        return returnCode;
    }
}
