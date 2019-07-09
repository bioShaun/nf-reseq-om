import jinja2
import fire
from pathlib import Path, PurePath
import os
import glob
import sys

script_dir, _ = os.path.split(os.path.abspath(__file__))

env = jinja2.Environment(loader=jinja2.FileSystemLoader(
    searchpath='{}/template'.format(script_dir)
)
)
template = env.get_template('index.html')
# code that fills in display_dictionary with the values to send to the template


def table2dict(table_file, name, sep='\t'):
    table_dict = dict()
    if table_file.exists():
        with open(table_file) as table_inf:
            for n, eachline in enumerate(table_inf):
                eachline_inf = eachline.split(sep)
                if n == 0:
                    label = '{}_header'.format(name)
                    table_dict[label] = eachline_inf
                else:
                    label = '{}_body'.format(name)
                    table_dict.setdefault(label, []).append(eachline_inf)
    if 'snp' in name:
        table_dict['snp'] = True
    return table_dict


def plot2report(plot_path, outpath, plot_flag, plot_name=None):
    plot_dict = dict()
    plots = glob.glob(str(plot_path))
    outpath = PurePath(outpath)
    if plots:
        plot = plots[0]
        plot_path = PurePath(plot)
        if plot_name is None:
            plot_name = plot_path.stem
        outfile_path = outpath / f'{plot_name}{plot_path.suffix}'
        os.system(f'cp {plot_path} {outfile_path}')
        plot_dict[plot_flag] = True
    return plot_dict


def exom_report(result_dir, proj_name, report_dir=None):
    result_dir = Path(result_dir)
    if report_dir is None:
        report_dir = result_dir / 'report'
    else:
        report_dir = Path(report_dir)
    if report_dir.is_dir():
        os.system(f'rm -r {report_dir}')
    display_dictionary = {}
    display_dictionary['project_name'] = proj_name

    # add fastqc table
    qc_table = result_dir / 'reads_qc/data.summary.csv'
    display_dictionary.update(
        table2dict(qc_table, 'seq', sep=','))

    # add aligment table
    align_table = result_dir / 'alignment/mapping.summary.csv'
    display_dictionary.update(
        table2dict(align_table, 'align', sep=','))

    # snp stats
    # summary
    snp_summary_table = result_dir / 'snp/overall.varSummary.txt'
    display_dictionary.update(
        table2dict(snp_summary_table, 'snp_summary'))

    snp_number_table = result_dir / 'snp/overall.varNum.txt'
    display_dictionary.update(
        table2dict(snp_number_table, 'snp_number'))

    snp_impact_table = result_dir / 'snp/overall.varImpact.txt'
    display_dictionary.update(
        table2dict(snp_impact_table, 'snp_impact'))

    snp_effect_table = result_dir / 'snp/overall.varEffects.txt'
    display_dictionary.update(
        table2dict(snp_effect_table, 'snp_effect'))

    snp_region_table = result_dir / 'snp/overall.varRegion.txt'
    display_dictionary.update(
        table2dict(snp_region_table, 'snp_region'))

    report_dir.mkdir(parents=True, exist_ok=True)
    os.system('cp -r {script_dir}/template/* {report_dir}'.format(
        script_dir=script_dir,
        report_dir=report_dir
    ))

    # plots
    report_plot_path = report_dir / 'imgs'
    mapping_plot = result_dir / 'plot/alignment/Mapping_stats.png'
    display_dictionary.update(
        plot2report(mapping_plot, report_plot_path, 'mapping_plot'))

    # genome_cov_plot = result_dir / 'plot/alignment/Reads_coverage_genome.png'
    # display_dictionary.update(
    #     plot2report(genome_cov_plot, report_plot_path, 'genome_cov_plot')
    # )

    # cds_cov_plot = result_dir / 'plot/alignment/Reads_coverage_cds.png'
    # display_dictionary.update(
    #     plot2report(cds_cov_plot, report_plot_path, 'cds_cov_plot')
    # )

    variant_summary_plot = result_dir / \
        'plot/variants/Variant_stats_summary.png'
    if variant_summary_plot.exists():
        display_dictionary.update(
            plot2report(variant_summary_plot,
                        report_plot_path, 'variant_summary')
        )
        variant_summary_plot_dir = result_dir / 'plot/variants/'
        for dir_i in variant_summary_plot_dir.iterdir():
            if dir_i.is_dir():
                example_sample = dir_i.name
        varType_plot = result_dir / \
            f'plot/variants/{example_sample}/{example_sample}_varType.png'
        display_dictionary.update(
            plot2report(varType_plot, report_plot_path,
                        'variant_type', 'varType'))

        varRegion_plot = result_dir / \
            f'plot/variants/{example_sample}/{example_sample}_varRegion.png'
        display_dictionary.update(
            plot2report(varRegion_plot, report_plot_path,
                        'variant_region', 'varRegion'))

        varEffects_high_plot = result_dir / \
            f'plot/variants/{example_sample}/{example_sample}_varEffects-HIGH.png'
        display_dictionary.update(
            plot2report(varEffects_high_plot, report_plot_path,
                        'variant_effect_high', 'varEffects-HIGH'))

        varEffects_moderate_plot = result_dir / \
            f'plot/variants/{example_sample}/{example_sample}_varEffects-MODERATE.png'
        display_dictionary.update(
            plot2report(varEffects_moderate_plot,
                        report_plot_path,
                        'variant_effect_moderate', 'varEffects-MODERATE'))

        varEffects_low_plot = result_dir / \
            f'plot/variants/{example_sample}/{example_sample}_varEffects-LOW.png'
        display_dictionary.update(
            plot2report(varEffects_low_plot, report_plot_path,
                        'variant_effect_low', 'varEffects-LOW'))

        varEffects_modifier_plot = result_dir / \
            f'plot/variants/{example_sample}/{example_sample}_varEffects-MODIFIER.png'
        display_dictionary.update(
            plot2report(varEffects_modifier_plot,
                        report_plot_path,
                        'variant_effect_modifier', 'varEffects-MODIFIER'))
        varImpact_plot = result_dir / \
            f'plot/variants/{example_sample}/{example_sample}_varImpact.png'
        display_dictionary.update(plot2report(varImpact_plot,
                                              report_plot_path, 'variant_impact', 'varImpact'))

    # deltaSNP_plot = result_dir / 'mapping/*deltaSNP.png'
    # Gprime_plot = result_dir / 'mapping/*Gprime.png'
    # negLog10Pval_plot = result_dir / 'mapping/*negLog10Pval.png'
    # plot2report(deltaSNP_plot, report_plot_path, 'deltaSNP')
    # plot2report(Gprime_plot, report_plot_path, 'Gprime')
    # plot2report(negLog10Pval_plot, report_plot_path, 'negLog10Pval')

    # display_dictionary.update({'pca': True, 'snp_index': True})
    display_html = template.render(display_dictionary)
    report_html = report_dir / 'index.html'
    with open(report_html, 'w') as out_inf:
        out_inf.write(display_html)
    os.system(f'tar -zcf {report_dir}.tar.gz {report_dir}')


if __name__ == '__main__':
    fire.Fire(exom_report)
