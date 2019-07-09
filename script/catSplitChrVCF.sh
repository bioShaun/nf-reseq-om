#!/bin/bash

htslib=/public/software/samtools/samtools-1.9/htslib-1.9/

    inputFile=$1
    outFile=$2

    grep "^#" $inputFile > $outFile
	
	sed -i '/##contig=<ID=chr/d' $outFile
	sed -i '/#CHROM/d' $outFile

	echo "##contig=<ID=chr1A,length=594102056>
##contig=<ID=chr1B,length=689851870>
##contig=<ID=chr1D,length=495453186>
##contig=<ID=chr2A,length=780798557>
##contig=<ID=chr2B,length=801256715>
##contig=<ID=chr2D,length=651852609>
##contig=<ID=chr3A,length=750843639>
##contig=<ID=chr3B,length=830829764>
##contig=<ID=chr3D,length=615552423>
##contig=<ID=chr4A,length=744588157>
##contig=<ID=chr4B,length=673617499>
##contig=<ID=chr4D,length=509857067>
##contig=<ID=chr5A,length=709773743>
##contig=<ID=chr5B,length=713149757>
##contig=<ID=chr5D,length=566080677>
##contig=<ID=chr6A,length=618079260>
##contig=<ID=chr6B,length=720988478>
##contig=<ID=chr6D,length=473592718>
##contig=<ID=chr7A,length=736706236>
##contig=<ID=chr7B,length=750620385>
##contig=<ID=chr7D,length=638686055>
##contig=<ID=chrUn,length=480980714>" >> $outFile

	grep  "#CHROM" $inputFile >> $outFile

    while read line
    do
        chrName=$(echo $line|cut -f1 -d" ")
        splitName=$(echo $line|cut -f4 -d" ")
        splitPos=$(echo $line|cut -f2 -d" ")
        grep "^$splitName" $inputFile |awk -v chrName="$chrName" -v splitPos="$splitPos" '{$1=chrName;$2=$2+splitPos;print $0}'|sed 's/ /\t/g' >> $outFile
    done < /public/database/exomeCapture/wheat/tcuni_180718/split.cat.bed

    grep "^chrUn" $inputFile >> $outFile

    $htslib/bgzip $outFile
    $htslib/tabix --csi $outFile.gz

