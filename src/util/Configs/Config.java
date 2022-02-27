package util.Configs;

public class Config<Type>
{
    private final Type defaultValue;
    private Type actualValue;

    public Config(Type defaultValue)
    {
        this.actualValue = this.defaultValue = defaultValue;
    }

    public void setValue(Object value)
    {
        try
        {
            this.actualValue = (Type) value;
        } catch (Exception e)
        {
            e.printStackTrace();
        }
    }

    public Type get()
    {
        return actualValue;
    }

    public String toString()
    {
        if (actualValue != null)
        {
            return actualValue.toString();
        }
        return "{NULL}";
    }
}
