import glob
import fire
import pandas as pd
from io import StringIO
from pathlib import Path


SNP_STATS = {
    'varType': '# Variantss by type',
    'varImpact': '# Effects by impact',
    'varClass': '# Effects by functional class',
    'varEffects': '# Count by effects',
    'varRegion': '# Count by genomic region',
    'varHoHe': '# Hom/Het table',
}


def extract_snpeff_stats(snpeff_stats_file, stats_name, suffix_pattern):
    # stats name pattern sample_id.hq.vcf.stat.csv
    snpeff_stats_file = Path(snpeff_stats_file)
    sample_id = snpeff_stats_file.name.replace(suffix_pattern, '')

    stats_flag_content = SNP_STATS[stats_name]
    stats_flag = 0
    content_list = []
    with open(snpeff_stats_file) as stats_inf:
        for eachline in stats_inf:
            if stats_flag_content in eachline:
                stats_flag = 1
                continue
            if stats_flag:
                if eachline.strip():
                    if eachline[0] == '#':
                        stats_flag = 0
                    elif 'Missense_Silent_ratio' in eachline:
                        pass
                    else:
                        content_list.append(eachline)
    stats_stringio = StringIO(''.join(content_list))
    sn_df = pd.read_csv(stats_stringio, sep=' , ',
                        engine='python')
    if stats_name == 'varHoHe':
        sn_df = sn_df.loc[[1, 2]]
        sn_df.columns = ['Type', 'Count']
        sn_df.loc[:, 'Percent'] = sn_df.Count / sn_df.Count.sum()
        sn_df.loc[:, 'Percent'] = [f'{count_i * 100}%' for count_i in sn_df.Percent]
    sn_df.loc[:, 'sample_id'] = sample_id
    return sn_df


def snp_stats_plot_and_show(snp_stats_dir, out_dir, suffix_pattern=".hq.vcf.stat.csv"):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    snpeff_files = glob.glob(f'{snp_stats_dir}/*{suffix_pattern}')
    for snp_stats_i in SNP_STATS:
        snp_stats_list = [extract_snpeff_stats(snpeff_file_i, snp_stats_i, suffix_pattern)
                          for snpeff_file_i in snpeff_files]
        snp_stats_df = pd.concat(snp_stats_list)
        snp_stats_df = snp_stats_df.set_index(['sample_id', 'Type'])
        show_df = snp_stats_df.unstack()
        show_df.columns.names = [None, None]
        show_df.index.name = ''
        # TODO 0% for percent value
        show_df.fillna(0, inplace=True)
        show_stats_file = out_dir / f'{snp_stats_i}.summary.csv'
        show_df.to_csv(show_stats_file)

        plot_df = snp_stats_df.unstack(level=0).loc[:, 'Count']
        plot_df.columns.name = ''
        plot_df.fillna(0, inplace=True)
        plot_file = out_dir / f'{snp_stats_i}.count.summary.txt'
        plot_df.to_csv(plot_file, sep='\t')


if __name__ == '__main__':
    fire.Fire(snp_stats_plot_and_show)
