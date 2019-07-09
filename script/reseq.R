suppressMessages(library(ggplot2))
suppressMessages(library(omplotr))
suppressMessages(library(argparser))
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(scales))
suppressMessages(library(plyr))

options(stringsAsFactors = F)

p <- arg_parser('Reseq stats plot collection.')
p <- add_argument(p, '--stats_dir',
                  help = 'Directory store plot stats files.')
p <- add_argument(p, '--out_dir',
                  help = 'Plot output directory.')
p <- add_argument(p, '--big_sample_number',
                  help = 'Sample number limit to show each sample detail in plot.',
                  default = 15,
                  type = 'numeric')
p <- add_argument(p, '--depth_limit',
                  help = 'Sequencing depth limit to plot.',
                  default = 100,
                  type = 'numeric')
p <- add_argument(p, '--exome',
                  help = 'Exome sequencing results.',
                  flag = TRUE)
p <- add_argument(p, '--mapping',
                  help = 'Generate mapping plot.',
                  flag = TRUE)
p <- add_argument(p, '--genome_cov',
                  help = 'Generate genome coverage plot.',
                  flag = TRUE)
p <- add_argument(p, '--variant',
                  help = 'Generate variant stats plot.',
                  flag = TRUE)
argv <- parse_args(p)

# # for test
# stats_dir <- './reseq_stats//'
# out_dir <- './reseq_plots/'
# BIG_SAMPLE_NUM <- 15
# MAX_DEPTH <- 100
# mapping_flag <- TRUE
# genome_cov_flag <- TRUE
# exome_flag <- FALSE
# variant_flag <- TRUE

# read parameters
stats_dir <- argv$stats_dir
out_dir <- argv$out_dir
BIG_SAMPLE_NUM <- argv$big_sample_number
mapping_flag <- argv$mapping
genome_cov_flag <- argv$genome_cov
exome_flag <- argv$exome
MAX_DEPTH <- argv$depth_limit
variant_flag <- argv$variant

# mapping part
if (mapping_flag) {
  bwa_mapping_stats <- file.path(stats_dir, 'alignment',
                                 'mapping.summary.txt')
  if (! file.exists(bwa_mapping_stats)) {
    stop('Mapping stats file not exists!')
  }
  
  mapping_plot_prefix <- file.path(out_dir, 'alignment', 'Mapping_stats')
  om_bwa_mapping_plot(bwa_mapping_stats, BIG_SAMPLE_NUM, mapping_plot_prefix)  
}
  
# genome cov part
if (genome_cov_flag) {
  cov_types <- c('genome', 'cds', 'exon')
  if (exome_flag) {
    cov_types <- cov_types[2:3]
  }
  coverage_tables <- file.path(stats_dir, 'alignment',
                               paste(cov_types, 'coverage.summary.txt', sep='.'))
  if (! all(file.exists(coverage_tables))) {
    stop('Coverage files missing!')
  }
  for (n in seq(length(cov_types))) {
    each_reg <- cov_types[n]
    each_table <- coverage_tables[n]
    cov_plot_prefix <- file.path(out_dir, 'alignment', 
                                 paste('Reads_coverage', each_reg, sep = '_'))
    om_reads_cov_plot(each_table, BIG_SAMPLE_NUM, MAX_DEPTH, cov_plot_prefix)
  }
}

# variant stats
if (variant_flag) {
  impact_map_file <- system.file("extdata", "variant_stats", "snpeff_varEffects.csv", package = "omplotr")
  variant_dir <- file.path(stats_dir, 'snp')
  var_stats <- lapply(om_const_reseq_variant$var_file_labs,
                      om_var_pie_stats,
                      stats_dir=variant_dir,
                      impact_map_file=impact_map_file,
                      varRegion_order=om_const_reseq_variant$varRegion_order)

  
  for (m_p_var_df in var_stats) {
    samples <- unique(m_p_var_df$variable)
    stats_name <- unique(m_p_var_df$fig)
    for (each_sample in samples) {
      sample_var_stats <- dplyr::filter(m_p_var_df, variable == each_sample)
      if (! grepl('varEffects', stats_name)[1]) {
        color_pal <- om_const_reseq_variant[['color_pal_map']][stats_name]
        var_plot_prefix <- file.path(out_dir, 'variants', each_sample, 
                                     paste(each_sample, stats_name, sep = '_'))      
        om_var_pie_plot(sample_var_stats, color_pal=color_pal, out_prefix=var_plot_prefix)          
      } else {
        for (each_impact in stats_name) {
          sample_ip_var_stats <- dplyr::filter(sample_var_stats, fig == each_impact)
          sample_ip_var_stats$value <- sample_ip_var_stats$value / sum(sample_ip_var_stats$value)
          sample_ip_var_stats$Percentage <- str_percent(sample_ip_var_stats$value)
          var_plot_prefix <- file.path(out_dir, 'variants', each_sample,
                                       paste(each_sample, each_impact, sep = '_'))      
          om_var_pie_plot(sample_ip_var_stats, out_prefix=var_plot_prefix)             
        }
      }
    }
  }
  
  var_stats_df <- plyr::ldply(var_stats, data.frame)
  var_stats_df$fig <- factor(var_stats_df$fig,
                             levels = om_const_reseq_variant$var_plot_order)
  var_summary_plot_prefix <- file.path(out_dir, 'variants', 
                                       'Variant_stats_summary')
  om_var_summary_plot(var_stats_df, var_summary_plot_prefix)
}

