import fire
import gtfparse
import pandas as pd
from pathlib import PurePath


def file_exist(file_path):
    file_path = Path(file_path)
    err_msg = f'{file_path} not exist!'
    assert file_path.exists(), err_msg


def check_reseq_refs(genome, gtf, bwa_idx, split_bed_dir, bedfiles):
    # bedfiles sep with ,
    bedfiles = bedfiles.split(',')
    map(file_exist, bedfiles)
    # check ref file existance
    genome = Path(genome)
    genome_dict = genome.with_suffix('.dict')
    genome_fai = genome.with_suffix(f'{genome.suffix}.fai')
    map(file_exist, [genome, genome_dict, genome_fai,
                     bwa_idx, split_bed_dir])


if __name__ == '__main__':
    fire.Fire(check_reseq_refs)
