#!/bin/env python
# -*- encoding=utf8 -*-
# import ConfigParser
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vcf_file", type=str,
                    help="vcf form snpEff",
                    required=True)
parser.add_argument("-o", "--out_name", type=str,
                    help="extract table name",
                    required=True)
args = parser.parse_args()

# column 8: sub-column:2,6,10
# cut -f8 | cut -f2,6,10 -d"|"


def extractAnnotaion(strs):
    AlleType = []
    AlleGene = []
    AlleLocu = []
    for annStr in strs.split(","):
        annDetailArray = strs.split("|")
        AlleType.append(annDetailArray[1])
        AlleGene.append(annDetailArray[3])
        AlleLocu.append(annDetailArray[9])

    AlleTypeRE = '|'.join(list(set(AlleType)))
    AlleGeneRE = '|'.join(list(set(AlleGene)))
    AlleLocuRE = '|'.join(list(set(AlleLocu)))
    return([AlleTypeRE, AlleGeneRE, AlleLocuRE])

# column 1,2,4,5,9-


def extractAlleFreq(lists):
    alle_freq = []
    for each_record in lists:
        if each_record.count('./.') > 0:
            alle_freq.append("0")
        else:
            each_sample = each_record.split(":")[1]
            alle_freq.append(each_sample)
    return(alle_freq)


vcf_file = gzip.open(args.vcf_file)
out_table = open(args.out_name, "w")

# find header line
cnt = 0
while cnt < 100000:
    cnt = cnt + 1
    line = vcf_file.readline()
    if line[0:4] == "#CHR":
        break
header = line.split("\t")[9:]
header_out = ["Chrom", "Pos", "Refer", "Alter",
              "AlleNum", "Feature", "Gene", "Alle"] + header
out_table.write('\t'.join(header_out))

while 1:
    lines = vcf_file.readlines(10000)
    if not lines:
        break
    for line in lines:
        line = line.decode()
        line_column = line.strip().split("\t")
        line_header = [line_column[i] for i in [0, 1, 3, 4]]
        Alle_number = [str(line_column[4].count(",") + 1)]
        line_ann = extractAnnotaion(line_column[7])
        line_all = extractAlleFreq(line_column[9:])
        out_inf = line_header + Alle_number + line_ann + line_all
        out_table.write('\t'.join(out_inf) + "\n")

vcf_file.close()
out_table.close()
