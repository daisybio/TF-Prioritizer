package util.Configs.ConfigValidators;

public class PositiveIntegerValidator extends IntegerRangeValidator
{
    public PositiveIntegerValidator()
    {
        super(1, Integer.MAX_VALUE);
    }
}
