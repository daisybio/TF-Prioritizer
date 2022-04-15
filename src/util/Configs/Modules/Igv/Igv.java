package util.Configs.Modules.Igv;

import util.Configs.ConfigTypes.GeneratedFileStructure;
import util.Configs.ConfigTypes.InputConfig;
import util.Configs.ConfigTypes.InputFileStructure;
import util.Configs.ClassGetter;
import util.Configs.ConfigTypes.SourceDirectoryFileStructure;
import util.Configs.ConfigValidators.IntegerRangeValidator;
import util.Configs.ConfigValidators.ListNotEmptyValidator;
import util.Configs.ConfigValidators.PositiveIntegerValidator;
import util.Configs.Modules.AbstractModule;
import util.Logger;

import java.lang.reflect.InvocationTargetException;
import java.util.List;
import java.util.Map;

import static util.FileManagement.extend;

public class Igv extends AbstractModule
{
    public FileStructure fileStructure;

    public final InputFileStructure pathToIGV = extend(sourceDirectory, "IGV_2.11.2");
    public final InputConfig<Integer> topLog2fc = new InputConfig<>(Integer.class, new PositiveIntegerValidator());
    public final InputConfig<Boolean> topLog2fcIncludeLncRnaPseudogenes = new InputConfig<>(Boolean.class);
    public final InputConfig<List<String>> includePredictionData = new InputConfig<>(ClassGetter.getStringList());
    public final InputConfig<List<String>> importantLociAllPrioTf =
            new InputConfig<>(ClassGetter.getStringList(), new ListNotEmptyValidator<>());
    public final InputFileStructure pathToTfChipSeq = new InputFileStructure();
    public final InputFileStructure pathToTdf = new InputFileStructure();
    public final InputConfig<List<String>> enhancerDatabases = new InputConfig<>(ClassGetter.getStringList());
    public final InputConfig<Map<String, String>> grcSynonymDict = new InputConfig<>(ClassGetter.getStringStringMap());
    public final InputConfig<String> speciesReferenceGenome = new InputConfig<>(String.class);

    public Igv(GeneratedFileStructure workingDirectory, SourceDirectoryFileStructure sourceDirectory, Logger logger)
            throws InvocationTargetException, IllegalAccessException, NoSuchMethodException, InstantiationException
    {
        super(workingDirectory, sourceDirectory, logger);
        init();
    }
}
