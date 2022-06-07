from time import sleep
import pandas as pd
from pyliftover import LiftOver

MAX_LIFTOVER_ATTEMPTS = 10

def convert(c, x, y, converter):
    first = converter.convert_coordinate(c, int(x))
    second = converter.convert_coordinate(c, int(y))

    if(first is None or second is None):
        return None, None

    if len(first) == 0 or len(second) == 0:
        return None, None

    return str(first[0][1]), str(second[0][1])


def main():
    path_to_X = "{INPUTFILE}"
    path_to_newSave = "{OUTPUTFILE}"

    version_of_our_dat = "{REFERENCE_GENOME}"
    version_of_their_dat = "{FOUND_GENOME}"

    please_convert = True
    if(version_of_our_dat == version_of_their_dat):
        please_convert = False

    df = pd.read_csv(path_to_X, sep="\t")

    data_output = []

    column_names = []

    for col in df.columns:
        column_names.append(col)

    attempt = 0
    worked = False

    while attempt < MAX_LIFTOVER_ATTEMPTS and not worked:
        try:
            converter = LiftOver(version_of_their_dat, version_of_our_dat, search_dir=None, cache_dir=None,
                                 use_web=True)
            worked = True
        except AttributeError:
            attempt += 1
            sleep(5)

    if not worked:
        raise IOError("Could not get liftover file")

    for _, row in df.iterrows():
        mgi_symbol = row[column_names.__getitem__(3)]
        chromosome = str(row[column_names.__getitem__(0)])
        start_position = row[column_names.__getitem__(1)]
        end_position = row[column_names.__getitem__(2)]

        x_converted = start_position
        y_converted = end_position

        if please_convert:
            x_converted, y_converted = convert(
                chromosome, start_position, end_position, converter)
            if(x_converted is None or y_converted is None):
                continue

        start_position = x_converted
        end_position = y_converted

        data_output.append([chromosome, start_position,
                           end_position, mgi_symbol, version_of_our_dat])

    df_converted = pd.DataFrame(data_output, columns=column_names)
    df_converted.to_csv(path_to_newSave, sep="\t", index=False)


main()
