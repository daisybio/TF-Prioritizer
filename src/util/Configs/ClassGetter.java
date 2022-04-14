package util.Configs;

import java.util.List;
import java.util.Map;

public class ClassGetter
{
    public static Class<List<String>> getStringList()
    {
        return (Class<List<String>>) ((Class) List.class);
    }

    public static Class<List<Double>> getDoubleList()
    {
        return (Class<List<Double>>) ((Class) List.class);
    }

    public static Class<Map<String, String>> getStringStringMap()
    {
        return (Class<Map<String, String>>) ((Class) Map.class);
    }
}
