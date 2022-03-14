package util.Configs;

public class RequiredConfig<T> extends Config<T>
{
    public RequiredConfig(Class<T> configClass)
    {
        super(configClass);
    }

    public boolean isWriteable()
    {
        return true;
    }

    public boolean isValid()
    {
        return isSet();
    }
}
