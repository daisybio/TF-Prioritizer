import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_context("notebook")
color = "#A6CEE3"
sns.set_context("talk")
sns.set(font_scale=2)
sns.set_style("whitegrid")
plt.figure(figsize=(26, 20))


def generate(sourcePath: str, threshold: float, targetPath: str, hm: str, group1: str, group2: str):
    df = pd.read_table(sourcePath).sort_values(['value'], ascending=False)
    name = hm + ": " + group1 + " VS " + group2

    # Remove suffix
    df['TF'] = df['TF'].str.split('_', expand=True)[0]
    # Sort and filter
    df_temp1 = df[df['value'] > threshold].set_index('TF')
    df_temp2 = df[df['value'] < -threshold].set_index('TF')
    df = pd.concat([df_temp1, df_temp2], axis=0)
    df.columns = [name]

    if ('Peak' in df.index):
        df = pd.DataFrame.drop(df, index='Peak')

    if not df.empty:  # Bar Plot
        sns.set(font_scale=3)
        if (len(df) > 61):
            sns.set(font_scale=2)
        if (len(df) > 81):
            sns.set(font_scale=1)
        sns.set_style("whitegrid")
        ax = sns.barplot(x=df.index, y=name, data=df, color=color)
        ax.set_title(hm + "\n" + group1 + ":" + group2)
        ax.set_ylabel('Normalized feature value\n(regression coefficient)')
        ax.set_xlabel('Transcription Factor')
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(targetPath)
        plt.clf()
    plt.figure(figsize=(26, 20))
    return df

def stages(dataPath: str, plotPath: str, df: pd.DataFrame):
    df.to_csv(dataPath, index=True, header=True)

    if not df.empty:
        change_colnames = True
        sns.set(font_scale=2)
        if (len(df) > 41):
            sns.set(font_scale=1)
            change_colnames = False
        if (len(df) > 61):
            sns.set(font_scale=0.8)
            change_colnames = False
        sns.set_style("whitegrid")
        column_names = list(df.columns)
        n_column_names = []
        for column_name in column_names:
            split1 = column_name.split(":")
            hm = split1[0]
            split2 = split1[1].split("VS")
            group1 = split2[0]
            group2 = split2[1]
            if change_colnames:
                n_column_names.append(hm + "\n" + group1 + ":" + group2)
            else:
                n_column_names.append(hm + " -" + group1 + ":" + group2)
        df.columns = n_column_names
        plot = sns.heatmap(df.transpose(), cmap="Paired", square=True, vmin=1, vmax=1, cbar=False, linewidths=0.5,
                        linecolor='black', xticklabels=True)
        plot.set_xlabel('Transcription Factor')
        plt.savefig(plotPath)

{CALLS}

{DIFFERENT_STAGES}

{SAME_STAGES}
