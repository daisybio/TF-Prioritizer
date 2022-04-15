from pyliftover import LiftOver
import pandas as pd

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

    converter = LiftOver(version_of_their_dat, version_of_our_dat)

    for index, row in df.iterrows():
        mgi_symbol = row[column_names.__getitem__(3)]
        chromosome = str(row[column_names.__getitem__(0)])
        start_position = row[column_names.__getitem__(1)]
        end_position = row[column_names.__getitem__(2)]

        ref_genome = row[column_names.__getitem__(4)]

        search_string = chromosome+":" + \
            str(start_position)+"-"+str(end_position)
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
