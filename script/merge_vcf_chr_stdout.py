import fire
import pandas as pd
from pathlib import PurePath
from Bio import bgzf


def merged_chr_size_inf(split_bed_df):
    merge_chr_size = split_bed_df.groupby(['chrom'])['end'].max()
    merge_chr_size_list = []
    for n, size_i in enumerate(merge_chr_size):
        eachline = '##contig=<ID={0},length={1}>'.format(
            merge_chr_size.index[n], size_i)
        merge_chr_size_list.append(eachline)
    return ('\n'.join(merge_chr_size_list))


# def merged_chr_vcf_file(vcf_file):
#     # a.xxx.[vcf|table].[gz|txt] -> a.xxx.catChr.[vcf|table].[gz|txt]
#     merge_chr_vcf_name = vcf_file.name
#     name_prefix = '.'.join(merge_chr_vcf_name.split('.')[:-2])
#     name_suffix = '.'.join(merge_chr_vcf_name.split('.')[-2:])
#     merge_chr_vcf_name = '{0}.catChr.{1}'.format(name_prefix, name_suffix)
#     merge_chr_vcf_file = str(vcf_file.with_name(merge_chr_vcf_name))
#     return merge_chr_vcf_file


def merge_vcf_split_chr(vcf_file, split_chr_inf):
    split_bed_df = pd.read_csv(split_chr_inf,
                               index_col=3,
                               header=None,
                               names=['chrom', 'start', 'end'],
                               sep='\t')
    merge_chr_size_str = merged_chr_size_inf(split_bed_df)
    vcf_file = PurePath(vcf_file)
    is_gz_file = vcf_file.suffix == '.gz'

    if is_gz_file:
        split_vcf_inf = bgzf.BgzfReader(vcf_file)
    else:
        split_vcf_inf = open(vcf_file)

    chr_header_flag = 1
    # TODO add merge chr command information in vcf header
    for eachline in split_vcf_inf:
        eachline = eachline.strip()
        eachline_inf = eachline.split('\t')
        chrom = eachline_inf[0]
        # split chrom size -> merged chrom size
        if eachline.startswith('##contig='):
            if chr_header_flag:
                print(merge_chr_size_str)
                chr_header_flag = 0
            continue
        elif chrom in split_bed_df.index:
            new_chrom, offset = split_bed_df.loc[chrom, ['chrom', 'start']]
            eachline_inf[0] = new_chrom
            eachline_inf[1] = str(offset + int(eachline_inf[1]))
            eachline = '\t'.join(eachline_inf)
        print(eachline)
    split_vcf_inf.close()


if __name__ == '__main__':
    fire.Fire(merge_vcf_split_chr)
