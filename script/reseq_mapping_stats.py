import sys
import glob
import fire
import pandas as pd
from io import StringIO
from pathlib import Path


STATS_LABEL = {
    'mapping': 'SN',
    'coverage': 'COV',
}

USE_COLS = {
    'mapping': [1, 2],
    'coverage': [2, 3],
}

SKIP_ROWS = {
    'mapping': [0, 1, 3, 4, 5, 11, 13, 16, 17, 20, 21, 25, 26, 27, 28, 29, 32, 33, 34, 35, 36],
    'coverage': None,
}

STATS_FILE = {
    'mapping': ['genome'],
    'coverage': ['genome', 'exon', 'cds'],
}


def extract_label_df(samtools_stats, label,
                     usecols, skiprows=None):
    # stats name pattern sample_id.[genome|cds|exon].stat
    samtools_stats = Path(samtools_stats)
    sample_id = '.'.join(samtools_stats.name.split('.')[:-2])
    content_list = []
    with open(samtools_stats) as stats_inf:
        for eachline in stats_inf:
            eachline = eachline.strip()
            if eachline.startswith(label):
                if '#' in eachline:
                    eachline_inf = eachline.strip().split('\t')
                    content_list.append('\t'.join(eachline_inf[:-1]))
                else:
                    content_list.append(eachline)
    stats_stringio = StringIO('\n'.join(content_list))
    sn_df = pd.read_table(stats_stringio, header=None,
                          index_col=0, names=[
                              'Item', sample_id],
                          usecols=usecols, skiprows=skiprows)
    # coverage > 1000 row name is 1000, duplicated index cause concat trouble
    if label == 'COV':
        sn_df = sn_df[sn_df.index < 1000]
    return sn_df


def check_stats(stats):
    if stats not in STATS_LABEL:
        print('unsupported stats!')
        supported_stats = ', '.join(STATS_LABEL.keys())
        print(f'supported stats [{supported_stats}].')
        sys.exit(1)


def mapping_stats_summary(mapping_stats_dir, stats, out_dir):
    check_stats(stats)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stats_file = STATS_FILE[stats]
    stats_label = STATS_LABEL[stats]
    stats_usecols = USE_COLS[stats]
    stats_skiprows = SKIP_ROWS[stats]
    for stats_i in stats_file:
        stats_i_files = glob.glob(f'{mapping_stats_dir}/*.{stats_i}.stat')
        if not stats_i_files:
            sys.exit(f'Wrong stats dir: {mapping_stats_dir}')
        stats_df_list = [extract_label_df(
            file_i, stats_label, stats_usecols, stats_skiprows)
            for file_i in stats_i_files]
        stats_df = pd.concat(stats_df_list, axis=1)
        stats_df.fillna(0, inplace=True)
        if len(stats_file) == 1:
            out_file = out_dir / f'{stats}.summary.txt'
        else:
            out_file = out_dir / f'{stats_i}.{stats}.summary.txt'
        out_file2 = out_file.with_suffix('.csv')
        stats_df.to_csv(out_file, sep='\t')
        stats_df.to_csv(out_file2)


if __name__ == '__main__':
    fire.Fire(mapping_stats_summary)
