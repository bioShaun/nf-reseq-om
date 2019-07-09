import fire
import gzip
import pandas as pd
from pathlib import PurePath


def merged_chr_size_inf(split_bed_df):
    merge_chr_size = split_bed_df.groupby(['chrom'])['end'].max()
    merge_chr_size_list = []
    for n, size_i in enumerate(merge_chr_size):
        eachline = '##contig=<ID={0},length={1}>'.format(
            merge_chr_size.index[n], size_i)
        merge_chr_size_list.append(eachline)
    return ('\n'.join(merge_chr_size_list))


def merge_vcf_split_chr(vcf_file, split_chr_inf, outfile):
    split_bed_df = pd.read_csv(split_chr_inf,
                               index_col=3,
                               header=None,
                               names=['chrom', 'start', 'end'],
                               sep='\t')
    merge_chr_size_str = merged_chr_size_inf(split_bed_df)
    vcf_file = PurePath(vcf_file)
    out_inf = open(outfile, 'a')
    is_gz_file = vcf_file.suffix == '.gz'
    if is_gz_file:
        split_vcf_inf = gzip.open(vcf_file, 'rt')
    else:
        split_vcf_inf = open(vcf_file)

    # extract header and merge chr info
    chr_header_flag = 1
    for eachline in split_vcf_inf:
        eachline = eachline.strip()
        # split chrom size -> merged chrom size
        if eachline.startswith('##contig='):
            if chr_header_flag:
                print(merge_chr_size_str, file=out_inf)
                chr_header_flag = 0
            continue
        elif not eachline.startswith('#'):
            break
        else:
            print(eachline, file=out_inf)

    # change vcf body chr-pos info
    vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
    out_cols = vcf_df.columns
    vcf_df = vcf_df.merge(split_bed_df, left_on=[0], right_index=True)
    vcf_df.loc[:, 0] = vcf_df.chrom
    vcf_df.loc[:, 1] = vcf_df.loc[:, 1] + vcf_df.start
    vcf_df = vcf_df.loc[:, out_cols]
    vcf_df.to_csv(out_inf, sep='\t', index=False, header=False)

    out_inf.close()


if __name__ == '__main__':
    fire.Fire(merge_vcf_split_chr)
