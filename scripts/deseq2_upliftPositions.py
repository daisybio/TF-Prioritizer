import pandas as pd
from pyliftover import LiftOver


def convert(c, x, y, s, converter):
    first = converter.convert_coordinate(c, int(x), s)
    second = converter.convert_coordinate(c, int(y), s)

    if (first is None or second is None):
        return None, None

    if len(first) == 0 or len(second) == 0:
        return None, None

    return str(first[0][1]), str(second[0][1])


def main():
    path_to_X = "{INPUTFILE}"
    version = "{VERSION}"
    path_to_newSave = "{OUTPUTFILE}"

    version_of_our_dat = "{REFERENCE_GENOME}"

    grc_dict = {{REFERENCE_GENOME_DICT}}

    version_convert_from = grc_dict.get(version)
    please_convert = True
    if (version_of_our_dat == version_convert_from):
        please_convert = False

    df = pd.read_csv(path_to_X, sep="\t")

    data_output = []

    column_names = []

    for col in df.columns:
        column_names.append(col)
    if please_convert:
        converter = LiftOver(version_convert_from, version_of_our_dat)
    for _, row in df.iterrows():
        mgi_symbol = row[column_names.__getitem__(0)]
        chromosome_not_edited = str(row[column_names.__getitem__(1)])
        chromosome = "chr" + str(row[column_names.__getitem__(1)])
        start_position = row[column_names.__getitem__(2)]
        end_position = row[column_names.__getitem__(3)]
        strand = row[column_names.__getitem__(4)]
        band = row[column_names.__getitem__(5)]

        ensg = row[column_names.__getitem__(6)]

        x_converted = start_position
        y_converted = end_position

        if please_convert:
            x_converted, y_converted = convert(chromosome, start_position, end_position, strand, converter)
            if (x_converted is None or y_converted is None):
                continue

        start_position = x_converted
        end_position = y_converted

        data_output.append([mgi_symbol, chromosome_not_edited, start_position, end_position, strand, band, ensg])

    df_converted = pd.DataFrame(data_output, columns=column_names)
    df_converted.to_csv(path_to_newSave, sep="\t", index=False)


main()
