#!/usr/bin/env nextflow

def helpMessage() {
    log.info """

    Usage:

    References If not specified in the configuration file or you wish to overwrite any of the references.
      --genome                      Name of Genomes reference
      --fasta                       Path to reference fasta
      --gtf                         Path to reference gtf
      --bwa_index                   Path to reference bwa index
      --star_index                  Path to reference star index
      --hisat_index                 Path to reference hisat index
      --exon_bed                    Path to reference exon bed file
      --cds_bed                     Path to reference cds bed file
      --padded_bed                  Path to reference pedded bed file
      --split_bed                   Path to split bed directory
      --known_vcf                   Path to known vcf file

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --outdir                      The output directory where the results will be saved

    Other options:
      --snpEff_db                   snpEff database name
      --quality                     GATK reads quality threshold
      --depth                       GATK reads depth threshold
      --data_type                   reseq, exome or rnaseq data
      --merge_chr_bed               bedfile to guide merge split chr information
      --aligner                     rnaseq mapping software star or hisat
      --skip_merge                  skip merge vcf for each sample
      --run_fastp                   run pipeline to fastp end
      --run_align                   run pipeline to alignment end
      --run_snp                     run whole snp pipeline

    """.stripIndent()
}

// workflow internal path&files
script_dir = file("$baseDir/script/")

/*
 * SET UP CONFIGURATION VARIABLES
 */

 // Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


def check_ref_exist = {file_path, file_type ->
    if (file_path) {
        file_path = file(file_path)
        if( !file_path.exists() ) exit 1, "${file_type} file not found: ${file_path}"
        return file_path
    } else {
        exit 1, "No reference genome ${file_type} specified!"
    }
}

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the Genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}


// default parameters
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bwa_index = params.genome ? params.genomes[ params.genome ].bwa_index ?: false : false
params.star_index = params.genome ? params.genomes[ params.genome ].star_index ?: false : false
params.hisat_index = false
params.cds_bed = params.genome ? params.genomes[ params.genome ].cds_bed ?: false : false
params.exon_bed = params.genome ? params.genomes[ params.genome ].exon_bed ?: false : false
params.padded_bed = params.genome ? params.genomes[ params.genome ].padded_bed ?: false : false
params.split_bed = params.genome ? params.genomes[ params.genome ].split_bed ?: false : false
params.snpEff_db = params.genome ? params.genomes[ params.genome ].snpEff_db ?: false : false
params.merge_chr_bed = params.genome ? params.genomes[ params.genome ].merge_chr_bed ?: false : false

// reference files
cds_bed_file = file(params.cds_bed)
exon_bed_file = file(params.exon_bed)
genome_fa = file(params.fasta)
genome_path = genome_fa.getParent()
genome_fai = file("${params.fasta}.fai")
genome_dict = file("${genome_path}/${genome_fa.baseName}.dict")


if (params.known_vcf) {
    known_vcf = file(params.known_vcf)
    known_vcf_index = file("${params.known_vcf}.tbi")
} else {
    known_vcf = false
    known_vcf_index = false
}

if (params.data_type == 'exome') {
    if (params.padded_bed) {
        split_bed_dir  = "${params.split_bed}/padded"
    } else {
        split_bed_dir  = "${params.split_bed}/exon"
    }
    
} else {
    split_bed_dir  = "${params.split_bed}/genome"
}

if (params.padded_bed) {
    padded_bed_file = file(params.padded_bed)
} else {
    padded_bed_file = exon_bed_file
}

// hisat index
if (params.hisat_index) {
    gtf = check_ref_exist(params.gtf, 'gtf file')
    hisat_index = Channel
                        .fromPath("${params.hisat_index}/*.ht2*")
                        .ifEmpty { exit 1, "HISAT2 index not found: ${params.hisat_index}" }

    alignment_splicesites = file("${params.hisat_index}/${gtf.baseName}.hisat2_splice_sites.txt")
}



// prepare split bed files
Channel
    .fromPath("${split_bed_dir}/*bed")
    .ifEmpty { exit 1, "Cannot find any bed file in directory: ${split_bed_dir}\n!" }
    .into { haplotype_beds; combine_gvcf_beds }


// Prepare analysis fastq files
Channel
    .fromFilePairs("${params.reads}")
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n!" }
    .set { raw_fq_files }


// pipeline control params

if (params.run_snp) {
    params.run_align = true
}

if (params.run_align) {
    params.run_fastp = true
}

if (! params.run_fastp) {
    exit 1, "Nothing to run! One of --run_fastp, --run_align, --run_snp is needed! "
}

/*
 * Fastp
 */

process fastp {

    tag "${name}"

    when:
    params.run_fastp

    publishDir "${params.outdir}/fastp_trimmed_reads/${name}", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".fq.gz") > 0 ? null: "$filename"}

    input:
    set name, file(reads) from raw_fq_files

    output:
    file "*trimmed.R*.fq.gz" into trimmed_reads
    file "${name}.json" into fastp_json
    file "${name}.html" into fastp_html

    cpus = 4

    script:
    """
    fastp \\
        --in1 ${reads[0]} \\
        --in2 ${reads[1]} \\
        --out1 ${name}.trimmed.R1.fq.gz \\
        --out2 ${name}.trimmed.R2.fq.gz \\
        --json ${name}.json \\
        --html ${name}.html    
    """
}  

process fastp_summary {

    publishDir "${params.outdir}/summary/", mode: 'copy'

    when:
    params.run_fastp

    input:
    file 'fastp/*' from fastp_json.collect()
    
    output:
    file 'qc' into qc_dir
    
    script:
    """
    python ${script_dir}/extract_fastp_info.py \\
        --fastp-dir fastp \\
        --outdir qc
    
    """
}


/*
* Alignment
*/
if (params.data_type == 'rnaseq') {

    if (params.aligner == 'star') {
        star_index_file = check_ref_exist(params.star_index, 'STAR index')
        process star_mapping {
            tag "${sample_name}"

            when:
            params.run_align

            publishDir "${params.outdir}/alignment/${sample_name}", mode: 'copy'

            queue 'bm'

            input:
            file reads from trimmed_reads
            file star_index_file from star_index_file

            output:
            file "${sample_name}.bam" into unsort_bam

            cpus = 40

            script:
            sample_name = reads[0].toString() - '.trimmed.R1.fq.gz'

            """
            STAR \\
                --genomeDir ${star_index_file} \\
                --readFilesIn ${reads}  \\
                --runThreadN 40 \\
                --twopassMode Basic \\
                --outSAMtype BAM Unsorted  \\
                --readFilesCommand zcat \\
                --outSAMattrRGline  ID:${sample_name} CN:TCuni SM:${sample_name} LB:${sample_name} PI:350 PL:Illumina \\
                --outFileNamePrefix ${sample_name}.

            mv ${sample_name}.Aligned.out.bam ${sample_name}.bam
            """

        }
    } else {
        process hisat_mapping {
            tag "${sample_name}"

            when:
            params.run_align

            publishDir "${params.outdir}/alignment/${sample_name}", mode: 'copy'

            input:
            file reads from trimmed_reads
            file index from hisat_index.collect()
            
            output:
            file "${sample_name}.bam" into unsort_bam

            cpus = 20
            
            script:
            sample_name = reads[0].toString() - '.trimmed.R1.fq.gz'
            index_base = index[0].toString() - ~/.\d.ht2l?/
            """
            hisat2 -x ${index_base} \\
                    -1 ${reads[0]} \\
                    -2 ${reads[1]} \\
                    --known-splicesite-infile ${alignment_splicesites} \\
                    --no-mixed \\
                    --no-discordant \\
                    -p ${task.cpus} \\
                    --met-stderr \\
                    --new-summary \\
                    --summary-file ${sample_name}.hisat2_summary.txt \\
                    --rg-id ${sample_name} \\
                    --rg SM:${sample_name}\\
                    --rg LB:${sample_name} \\
                    --rg PI:350 \\
                    --rg PL:Illumina \\
                    --rg CN:TCuni \\
                    | samtools view -bS -F 4 -F 8 -F 256 - > ${sample_name}.bam           
            """
        }
    }


} else {
    bwa_index_file = check_ref_exist(params.bwa_index, 'bwa index')
    process bwa_mapping {
        tag "${sample_name}"

        when:
        params.run_align

        input:
        file reads from trimmed_reads
        file bwa_index_file from bwa_index_file
        file genome_fa from genome_fa

        output: 
        file "${sample_name}.bam" into unsort_bam

        cpus = 40

        script:
        sample_name = reads[0].toString() - '.trimmed.R1.fq.gz'
        """
        bwa mem -M -a \\
            -R \"@RG\\tID:${sample_name}\\tSM:${sample_name}\\tLB:${sample_name}\\tPI:350\\tPL:Illumina\\tCN:TCuni\" \\
            -t ${task.cpus} \\
            -K 10000000 \\
            ${bwa_index_file}/${genome_fa.getName()} \\
            ${reads[0]} \\
            ${reads[1]} \\
            | \\
        samtools view -O bam \\
            --threads ${task.cpus} \\
            -o ${sample_name}.bam
        """
    }
}


/*
* Sort bam
*/
process sort_bam {
    tag "${sample_name}"

    when:
    params.run_align

    input:
    file bam from unsort_bam

    output: 
    file "${sample_name}.sort.bam" into samtools_stats_bam, to_rmdup_bam

    cpus = 8

    script:
    sample_name = bam.baseName
    """	
    samtools fixmate --threads ${task.cpus} \\
        -m ${bam} ${sample_name}.fixmate.bam

    samtools sort -m 2400M --threads ${task.cpus} \\
	    -o ${sample_name}.sort.bam \\
	    ${sample_name}.fixmate.bam
    """
}

/*
* Reads coverage stats
*/
process reads_cov_stats {

    tag "${sample_name}"

    when:
    params.run_align

    publishDir "${params.outdir}/alignment/${sample_name}", mode: 'copy'

    input:
    file samtools_stats_bam from samtools_stats_bam
    file cds_bed_file from cds_bed_file
    file exon_bed_file from exon_bed_file

    output:
    file "${sample_name}.*.stat" into reads_cov_results

    cpus = 8

    script:
    sample_name = samtools_stats_bam.baseName - '.sort'
    """
    samtools stats \\
	    --threads ${task.cpus} \\
	    --target-regions  ${cds_bed_file} \\
	    ${samtools_stats_bam} \\
	    > ${sample_name}.cds.stat

    samtools stats \\
	    --threads ${task.cpus} \\
	    --target-regions  ${exon_bed_file} \\
	    ${samtools_stats_bam} \\
	    > ${sample_name}.exon.stat

    samtools stats \\
	    --threads ${task.cpus} \\
	    ${samtools_stats_bam} \\
	    > ${sample_name}.genome.stat 
    """
}


process align_summary {

    when:
    params.run_align

    publishDir "${params.outdir}/summary/", mode: 'copy'

    input:
    file 'alignment_stats/*' from reads_cov_results.collect()
    
    output:
    file "alignment" into mapping_dir
    
    script:
    """
    python ${script_dir}/reseq_mapping_stats.py \\
        --mapping-stats-dir alignment_stats \\
        --stats mapping \\
        --out-dir alignment 

    python ${script_dir}/reseq_mapping_stats.py \\
        --mapping-stats-dir alignment_stats \\
        --stats coverage \\
        --out-dir alignment
    """
}


/*
* Reads mark duplication and recalibration
*/
process bam_remove_duplicate {
    tag "${sample_name}"

    when:
    params.run_snp

    input:
    file bam from to_rmdup_bam
  
    output:
    file "${sample_name}.rmdup.bam" into br_rmdup_bam
    file "${sample_name}.rmdup.bai" into br_rmdup_bam_idx
  
    cpus = 8

    script:
    sample_name = bam.baseName - '.sort'
    """
    samtools markdup -r \\
	    --threads ${task.cpus} \\
	    ${bam} \\
	    ${sample_name}.rmdup.bam

    samtools index ${sample_name}.rmdup.bam ${sample_name}.rmdup.bai -@ 8
    """
}



if (params.data_type == 'rnaseq') {

    /*
    * Split N
    */

    process SplitNCigarReads {

        tag "${sample_name}"

        when:
        params.run_snp
            
        input:
        file bam from br_rmdup_bam
        file fasta from genome_fa
        file fa_dict from genome_dict
        file fai_idx from genome_fai
        
        output:
        file "${sample_name}.qc.bam" into qc_bam, sample_bam
        file "${sample_name}.qc.bai" into qc_bam_idx

        cpus = 8
        
        script:
        sample_name = bam.baseName - '.rmdup'

        """
        gatk SplitNCigarReads \\
            --reference ${fasta} \\
            --input ${bam} \\
            --output ${sample_name}.qc.bam \\
            --max-reads-in-memory 1000000  \\
        """
    } 

} else {
    /*
    * bam BaseRecalibrator
    */
    process bam_BaseRecalibrator {
        tag "${sample_name}"

        // publishDir "${params.outdir}/alignment/${sample_name}", mode: 'copy'

        when:
        params.known_vcf && params.run_snp

        input:
        file bam from br_rmdup_bam
        file refer from genome_fa
        file refer_fai from genome_fai
        file refer_dict from genome_dict
        file known_vcf from known_vcf
        file known_vcf_index from known_vcf_index
    
        output:
        file "${sample_name}.recal.table" into recal_table
        file bam into bqsr_rmdup_bam
    
        cpus = 8

        script:
        sample_name = bam.baseName - '.rmdup'
        """
        gatk BaseRecalibrator \\
            --reference ${refer} \\
            --input ${bam} \\
            --output ${sample_name}.recal.table \\
            --known-sites ${known_vcf}
        """
    }

    /*
    * bam ApplyBQSR
    */
    process bam_ApplyBQSR {
        tag "${sample_name}"

        publishDir "${params.outdir}/alignment/${sample_name}", mode: 'copy'

        when:
        params.known_vcf && params.run_snp

        input:
        file bam from bqsr_rmdup_bam
        file recal_table_file from recal_table
        file refer from genome_fa
        file refer_fai from genome_fai
        file refer_dict from genome_dict
    
        output:
        file "${sample_name}.qc.bam" into qc_bam, sample_bam
        file "${sample_name}.qc.bai" into qc_bam_idx
    
        cpus = 8

        script:
        sample_name = bam.baseName - '.rmdup'
        """       
        gatk ApplyBQSR \\
            --bqsr-recal-file ${recal_table_file} \\
            --input ${bam} \\
            --output ${sample_name}.qc.bam \\
            --read-filter AmbiguousBaseReadFilter \\
            --read-filter MappingQualityReadFilter \\
            --read-filter NonZeroReferenceLengthAlignmentReadFilter \\
            --read-filter ProperlyPairedReadFilter \\
            --minimum-mapping-quality 30 

        """
    }

}



/*
*  GATK HaplotypeCaller
*/
process gatk_HaplotypeCaller {
    tag "${sample_name}|${chr_name}"

    when:
    params.known_vcf && params.run_snp

    input:
    file bam from qc_bam
    each file(bed) from haplotype_beds
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict    

    output:
    file "${sample_name}.${chr_name}.hc.g.vcf.gz" into sample_gvcf, chr_gvcf
    file "${sample_name}.${chr_name}.hc.g.vcf.gz.tbi" into sample_gvcf_index, chr_gvcf_index
    
    cpus = 8

    script:
    sample_name = bam.baseName - '.qc'

    chr_name = bed.baseName
    """
    gatk HaplotypeCaller  \\
        --input ${bam} \\
        --output ${sample_name}.${chr_name}.hc.g.vcf.gz \\
        --reference ${refer} \\
        --intervals ${bed} \\
        --emit-ref-confidence GVCF \\
        --read-filter AmbiguousBaseReadFilter \\
        --read-filter MappingQualityReadFilter \\
        --read-filter NonZeroReferenceLengthAlignmentReadFilter \\
        --read-filter ProperlyPairedReadFilter \\
        --minimum-mapping-quality 30         
    """
}

/*
* GATK CombineGVCFs
*/
process gatk_CombineGVCFs {
    tag "Chrom: ${chr_name}"

    publishDir "${params.outdir}/gvcf/by_chr/${chr_name}", mode: 'copy'

    when:
    params.known_vcf && params.run_snp

    input:
    file ('gvcf/*') from chr_gvcf.collect()
    file ('gvcf/*') from chr_gvcf_index.collect()
    each file(bed) from combine_gvcf_beds
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict    
    
    output:
    file "${chr_name}.g.vcf.gz" into merged_sample_gvcf
    file "${chr_name}.g.vcf.gz.tbi" into merged_sample_gvcf_index
    
    script:
    chr_name = bed.baseName
    """
    ls gvcf/*.${chr_name}.hc.g.vcf.gz > ${chr_name}.gvcf.list

    /public/software/GATK/gatk-4.0.10.1/gatk CombineGVCFs \\
    	--output ${chr_name}.g.vcf.gz \\
	    --reference ${refer} \\
	    --variant ${chr_name}.gvcf.list
    """
}

/*
* GATK GenotypeGVCFs
*/
process gatk_GenotypeGVCFs {
    tag "Chrom: ${chr_name}"

    publishDir "${params.outdir}/vcf/by_chr/${chr_name}", mode: 'copy'

    when:
    params.known_vcf && params.run_snp && !params.skip_merge

    input:
    file gvcf from merged_sample_gvcf
    file gvcf_index from merged_sample_gvcf_index
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict    
    
    output:
    file "${chr_name}.vcf.gz" into merged_sample_vcf
    file "${chr_name}.vcf.gz.tbi" into merged_sample_vcf_index

    cpus = 8
    
    script:
    chr_name = gvcf.baseName - '.g.vcf'

    """
    gatk GenotypeGVCFs \\
        --reference ${refer} \\
        --variant ${gvcf} \\
        --intervals ${split_bed_dir}/${chr_name}.bed \\
        --output ${chr_name}.vcf.gz
    """
}

/*
* Concat vcf
*/
process concat_vcf {

    when:
    params.known_vcf && params.run_snp && !params.skip_merge

    input:
    file ('vcf/*') from merged_sample_vcf.collect()
    file ('vcf/*') from merged_sample_vcf_index.collect()

    output:
    file "raw.vcf.gz" into all_sample_raw_vcf, m_all_sample_raw_vcf
    file "raw.vcf.gz.tbi" into all_sample_raw_vcf_idx, m_all_sample_raw_vcf_idx

    cpus = 8
    
    script:
    """
    bcftools concat \\
	    vcf/*.vcf.gz | \\
        bgzip > raw.vcf.gz

    tabix -p vcf raw.vcf.gz
    """
}

/*
* Basic quality filter
*/
process vcf_base_qual_filter {

    when:
    params.known_vcf && params.run_snp && !params.skip_merge

    input:
    file raw_vcf from all_sample_raw_vcf
    file raw_vcf_idx from all_sample_raw_vcf_idx
    
    output:
    file "hq.vcf.gz" into all_hq_vcf, all_hq_vcf_table
    file "hq.vcf.gz.tbi" into all_hq_vcf_idx, all_hq_vcf_idx_table

    cpus = 8
    
    script:
    """
    bcftools filter -e '%QUAL<${params.quality} || INFO/DP<${params.depth}' \\
	    raw.vcf.gz | \\
        bgzip > hq.vcf.gz

    tabix -p vcf hq.vcf.gz
    """
}


/*
* SNP Table 1
*/
process snp_gatk_table {

    publishDir "${params.outdir}/vcf/all", mode: 'copy'

    when:
    params.known_vcf && params.snpEff_db && params.run_snp && !params.skip_merge

    input:
    file vcf from all_hq_vcf_table
    file vcf_idx from all_hq_vcf_idx_table
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict 

    output:
    file "vcf.gatk.table.txt" into gatk_vcf_table

    script:
    """
    java -jar /public/software/GATK/GATK.3.8/GenomeAnalysisTK.jar \\
        -T VariantsToTable \\
        -R ${refer} \\
        -V ${vcf} \\
        -F CHROM -F POS -F REF -F ALT \\
        -GF AD -GF DP -GF GQ -GF PL \\
        -o vcf.gatk.table.txt
    """
}

/*
* snpeff for combined vcf
*/
process snpEff_for_all {

    publishDir "${params.outdir}/vcf/all", mode: 'copy'

    when:
    params.known_vcf && params.snpEff_db && params.run_snp && !params.skip_merge

    input:
    file raw_vcf from m_all_sample_raw_vcf
    file raw_vcf_idx from m_all_sample_raw_vcf_idx
    file vcf from all_hq_vcf
    file vcf_idx from all_hq_vcf_idx
    
    output:
    file "raw.vcf.gz" into final_all_raw_vcf
    file "raw.vcf.gz.*" into final_all_raw_vcf_idx    
    file "hq.vcf.gz" into final_all_hq_vcf
    file "hq.vcf.gz.*" into final_all_hq_vcf_idx
    file "hq.vcf.stat.csv"
    file "hq.vcf.stat.html"
    file "hq.vcf.stat.genes.txt"
    file "hq.ann.vcf.gz" into all_sample_anno_vcf, all_sample_anno_split_vcf
    file "hq.ann.vcf.gz.*" into all_sample_anno_vcf_idx, all_sample_anno_split_vcf_idx
    
    cpus = 8  

    script:
    if (params.merge_chr_bed)
        """
        mv raw.vcf.gz raw.split.vcf.gz
        mv raw.vcf.gz.tbi raw.split.vcf.gz.tbi

        python ${script_dir}/merge_vcf_chr_pd.py \\
            --vcf-file raw.split.vcf.gz \\
            --split-chr-inf ${params.merge_chr_bed} \\
            --outfile raw.vcf            
        bgzip raw.vcf
        tabix --csi raw.vcf.gz

        mv hq.vcf.gz hq.split.vcf.gz
        mv hq.vcf.gz.tbi hq.split.vcf.gz.tbi
        python ${script_dir}/merge_vcf_chr_pd.py \\
            --vcf-file hq.split.vcf.gz \\
            --split-chr-inf ${params.merge_chr_bed} \\
            --outfile hq.vcf 
        bgzip hq.vcf
        tabix --csi hq.vcf.gz

        snpEff -Xmx10g \\
            -c ${params.snpEff}/snpEff.config \\
            -csvStats hq.vcf.stat.csv \\
            -htmlStats hq.vcf.stat.html \\
            -v ${params.snpEff_db}  \\
            -quiet \\
            hq.vcf.gz \\
            | bgzip > hq.ann.vcf.gz
        
        tabix --csi hq.ann.vcf.gz
        """
    else
        """
        snpEff -Xmx10g \\
            -c ${params.snpEff}/snpEff.config \\
            -csvStats hq.vcf.stat.csv \\
            -htmlStats hq.vcf.stat.html \\
            -v ${params.snpEff_db}  \\
            -quiet \\
            ${vcf} \\
            | bgzip > hq.ann.vcf.gz
        
        tabix -p vcf hq.ann.vcf.gz        
        """
}


/*
* extract vcf file for each sample
*/

process gatk_CombineGVCFs_by_sample {
    tag "${sample_name}"

    publishDir "${params.outdir}/gvcf/by_sample/${sample_name}", mode: 'copy'

    when:
    params.known_vcf  && params.run_snp

    input:
    file ('gvcf/*') from sample_gvcf.collect()
    file ('gvcf/*') from sample_gvcf_index.collect()
    file padded_bed_file from padded_bed_file
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict    
    file bam from sample_bam

    output:
    file "${sample_name}.g.vcf.gz" into merged_sample_chr_gvcf
    file "${sample_name}.g.vcf.gz.tbi" into merged_sample_chr_gvcf_idx

    cpus = 8
    
    script:
    sample_name = bam.baseName - '.qc'
    """
    ls gvcf/${sample_name}.*.hc.g.vcf.gz > ${sample_name}.gvcf.list

    gatk CombineGVCFs \\
    	--output ${sample_name}.g.vcf.gz \\
	    --reference ${refer} \\
	    --variant ${sample_name}.gvcf.list
    """
}

process gatk_GenotypeGVCFs_by_sample {
    tag "${sample_name}"

    when:
    params.known_vcf && params.run_snp

    input:
    file gvcf from merged_sample_chr_gvcf
    file gvcf_index from merged_sample_chr_gvcf_idx
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict    
    
    output:
    file "${sample_name}.raw.vcf.gz" into merged_sample_chr_vcf, m_merged_sample_chr_vcf
    file "${sample_name}.raw.vcf.gz.tbi" into merged_sample_chr_vcf_index, m_merged_sample_chr_vcf_index

    cpus = 8
    
    script:
    sample_name = gvcf.baseName - '.g.vcf'
    if (params.data_type == 'exome')
        """
        gatk GenotypeGVCFs \\
            --reference ${refer} \\
            --variant ${gvcf} \\
            --intervals ${padded_bed_file} \\
            --output ${sample_name}.raw.vcf.gz
        """
    else 
        """
        gatk GenotypeGVCFs \\
            --reference ${refer} \\
            --variant ${gvcf} \\
            --output ${sample_name}.raw.vcf.gz        
        """
}

process vcf_base_qual_filter_by_sample {
    tag "${sample_name}"

    when:
    params.known_vcf && params.run_snp

    input:
    file raw_vcf from merged_sample_chr_vcf
    file raw_vcf_idx from merged_sample_chr_vcf_index
    
    output:
    file "${sample_name}.hq.vcf.gz" into sample_hq_vcf
    file "${sample_name}.hq.vcf.gz.tbi" into sample_hq_vcf_idx

    cpus = 8
    
    script:
    sample_name = raw_vcf.baseName - '.raw.vcf'
    """
    bcftools filter -e '%QUAL<${params.quality} || INFO/DP<${params.depth}' \\
	    ${raw_vcf} | \\
        bgzip > ${sample_name}.hq.vcf.gz

    tabix -p vcf ${sample_name}.hq.vcf.gz
    """
}

/*
* snpEff for each sample
*/
process snpEff_for_sample {

    tag "${sample_name}"

    publishDir "${params.outdir}/vcf/by_sample/${sample_name}", mode: 'copy'

    when:
    params.known_vcf && params.snpEff_db && params.run_snp

    input:
    file "raw_vcf/*" from m_merged_sample_chr_vcf.collect()
    file "raw_vcf/*" from m_merged_sample_chr_vcf_index.collect()
    file vcf from sample_hq_vcf
    file vcf_idx from sample_hq_vcf_idx
    
    output:
    file "${sample_name}.raw.vcf.gz"
    file "${sample_name}.raw.vcf.gz.*"
    file "${sample_name}.hq.vcf.gz"
    file "${sample_name}.hq.vcf.gz.*"
    file "${sample_name}.hq.vcf.stat.csv" into sample_vcf_stats
    file "${sample_name}.hq.vcf.stat.html"
    file "${sample_name}.hq.vcf.stat.genes.txt"
    file "${sample_name}.ann.vcf.gz" into single_sample_anno_vcf
    file "${sample_name}.ann.vcf.gz.*" into single_sample_anno_vcf_idx

    cpus = 10
    
    script:
    sample_name = vcf.baseName - '.hq.vcf'
    if (params.merge_chr_bed)
        """
        mv raw_vcf/${sample_name}.raw.vcf.gz raw_vcf/${sample_name}.raw.split.vcf.gz

        python ${script_dir}/merge_vcf_chr_pd.py \\
            --vcf-file raw_vcf/${sample_name}.raw.split.vcf.gz \\
            --split-chr-inf ${params.merge_chr_bed} \\
            --outfile ${sample_name}.raw.vcf 
        bgzip ${sample_name}.raw.vcf 
        tabix --csi ${sample_name}.raw.vcf.gz

        mv ${sample_name}.hq.vcf.gz ${sample_name}.hq.split.vcf.gz
        mv ${sample_name}.hq.vcf.gz.tbi ${sample_name}.hq.split.vcf.gz.tbi

        python ${script_dir}/merge_vcf_chr_pd.py \\
            --vcf-file ${sample_name}.hq.split.vcf.gz \\
            --split-chr-inf ${params.merge_chr_bed} \\
            --outfile ${sample_name}.hq.vcf 
        bgzip ${sample_name}.hq.vcf 
        tabix --csi ${sample_name}.hq.vcf.gz

        snpEff -Xmx10g \\
            -c ${params.snpEff}/snpEff.config \\
            -csvStats ${sample_name}.hq.vcf.stat.csv \\
            -htmlStats ${sample_name}.hq.vcf.stat.html \\
            -v ${params.snpEff_db}  \\
            -quiet \\
            ${sample_name}.hq.vcf.gz \\
            | bgzip > ${sample_name}.ann.vcf.gz

        tabix --csi ${sample_name}.ann.vcf.gz        
        """
    else
        """
        snpEff -Xmx10g \\
            -c ${params.snpEff}/snpEff.config \\
            -csvStats ${sample_name}.hq.vcf.stat.csv \\
            -htmlStats ${sample_name}.hq.vcf.stat.html \\
            -v ${params.snpEff_db}  \\
            -quiet \\
            ${vcf} \\
            | bgzip > ${sample_name}.ann.vcf.gz

        tabix -p vcf ${sample_name}.ann.vcf.gz
        """
}

/*
* SNP Table 2
*/
process snp_inhouse_table {

    publishDir "${params.outdir}/vcf/all", mode: 'copy'

    when:
    params.known_vcf && params.snpEff_db && params.run_snp && !params.skip_merge

    input:
    file vcf from all_sample_anno_vcf
    file vcf_idx from all_sample_anno_vcf_idx
    file refer from genome_fa
    file refer_fai from genome_fai
    file refer_dict from genome_dict 

    output:
    file "vcf.table.txt" into om_vcf_table

    script:
    """   
    python ${script_dir}/extractTableFromsnpEff.py \\
        -v ${vcf}  \\
        -o vcf.table.txt
    """
}


/*
* snp summary and plot
*/
process snp_summary {
    publishDir "${params.outdir}/summary/", mode: 'copy'

    when:
    params.known_vcf && params.snpEff_db && params.run_snp

    input:
    file qc_dir from qc_dir
    file mapping_dir from mapping_dir
    file 'snp_stats/*' from sample_vcf_stats.collect()

    output:
    file "snp/*csv" into snp_summary_table
    file "plot" into plot_dir

    script:
    """
    python ${script_dir}/reseq_snpeff_summary.py \\
        --snp-stats-dir snp_stats \\
        --out-dir snp

    /public/software/R/R-3.5.1/executable/bin/Rscript \\
        ${script_dir}/reseq.R \\
        --stats_dir . \\
        --out_dir ./plot \\
        --mapping \\
        --genome_cov \\
        --variant
    """
}
