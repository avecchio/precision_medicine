from scipy.stats import chisquare, chi2_contingency

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

import argparse
import sys

def read_phenotype_file(filename):
    df = pd.read_csv(filename, sep='\t', header=0)
    records = df.to_dict('records')
    phenotype_map = {}
    for record in records:
        phenotype_map[record['Sample']] = record['Phenotype']
    return phenotype_map

def read_vcf_file(filename):
    with open(filename, 'r') as fp:
        lines = fp.readlines()
    start_index = 0
    for line in lines:
        if line[0:2] != '##':
            break
        start_index += 1
    columns = lines[start_index].replace("\n", "").split('\t')
    entries = lines[start_index + 1:]
    split_entries = [entry.replace("\n", "").split('\t') for entry in entries]
    df = pd.DataFrame(split_entries, columns=columns)
    return df.to_dict('records')

def genotype_to_allele(genotype):
    if genotype == '0|0':
        return 'aa'
    elif genotype in ['1|0', '0|1']:
        return 'Aa'
    elif genotype == '1|1':
        return 'AA'
    else:
        print('error')

phenotype_map = read_phenotype_file("PRS_phen.txt")
snps = read_vcf_file('GWAS_data.vcf')

genotypes = {}
snp_pvals = []
for snp in snps:
    snp_id = snp['ID']

    genotype_freq_table = {
        'cases' : {'aa': 0, 'Aa': 0, 'AA': 0},
        'controls': {'aa': 0, 'Aa': 0, 'AA': 0}
    }

    for sample in phenotype_map:
        phenotype = phenotype_map[sample]
        genotype = snp[sample].split(":")[0]
        allele = genotype_to_allele(genotype)
        if phenotype == 1:
            genotype_freq_table['cases'][allele] += 1
        else:
            genotype_freq_table['controls'][allele] += 1

    contingency_df = pd.DataFrame(genotype_freq_table).transpose()
    #print(contingency_df)
    #values = contingency_df.to_numpy().flatten()
    try:
        res = chi2_contingency(contingency_df, correction=False)
        snp_pvals.append({'id': snp_id, 'pval': res.pvalue})
    except:
        snp_pvals.append({'id': snp_id, 'pval': 1})
        #print(contingency_df)


    #if genotype not in genotypes:
    #    genotypes[genotype] = 0
    #genotypes[genotype] += 1


'''

#############################################
# Core processing or statistical calculations
#############################################


        # input: frequeny of alleles: aa (int), Aa (int) AA (int)
#        where the column matches the column in the VCF file
# output: the hwe chi square (int) and p_value (float) statistics
#         calculate hwe chi square and p value statistic for each variant
#         If any column (aa, Aa, AA) sums to 0 in the contingency table, the resulting values
#         will be 1 (chi square) and 0 (p value) respectively
def calculate_hwe(aa, Aa, AA):

    # calculate the expected values
    total = aa + Aa + AA
    p = (Aa + (AA*2)) / (total * 2)
    q = 1 - p

    expected_aa = (q * q) * total
    expected_Aa = 2 * (p * q) * total
    expected_AA = (p * p) * total

    observed_values = [aa, Aa, AA]
    expected_values = [
        expected_aa,
        expected_Aa,
        expected_AA
    ]

    # if any of the columns (aa, Aa, or AA) sum to zero,
    # chi square is set to 0 and the p value is set to 1 as the actual
    # chi square equation cannot properly calcuate the real value
    for column_index in range(0, len(observed_values)):
        if observed_values[column_index] + expected_values[column_index] == 0:
            return 0, 1

    # using scipy chisquare to determine if there is a statistical difference
    # between all observed vs expected values
    res = chisquare(observed_values, f_exp=expected_values)
    chi_square, p_value = res

    return chi_square, p_value


###################################
#Output generators
###################################

def record_significant_p_values(variants, id_key, p_value_key, hwe_p_value_key, significance, output_filepath):
    record_entries = []
    record_entries.append(f'VAR_ID P_VALUE HWE_P_VALUE')
    for variant in variants:
        if variant[p_value_key] < significance:
            p_val = variant[p_value_key]
            hwe_p_value = variant[hwe_p_value_key]
            id = variant[id_key]
            record_entries.append(f'{id} {p_val} {hwe_p_value}')

    with open(output_filepath, 'w') as f:
        f.write('\n'.join(record_entries))

def generate_manhattan_plot(
        data, chromosome_key, p_value_key, export_path = None, significance = 5e-8,
        colors = ['#E24E42', '#008F95'], refSNP = False
    ):

    data['-log10(p_value)'] = -np.log10(data[p_value_key])
    data[chromosome_key] = data[chromosome_key].astype('category')
    data['ind'] = range(len(data))
    data_grouped = data.groupby((chromosome_key))

    fig = plt.figure()
    ax = fig.add_subplot(111)

    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(data_grouped):
        group.plot(kind='scatter', x='ind', y='-log10(p_value)', color=colors[num % len(colors)], ax=ax, s= 10000/len(data))
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))

    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(data)])
    ax.set_ylim([0, data['-log10(p_value)'].max() + 1])
    ax.set_xlabel('Chromosome')
    plt.axhline(y=significance, color='black', linestyle='-', linewidth = 1)
    plt.xticks(fontsize=8, rotation=60)
    plt.yticks(fontsize=8)

    if refSNP:
        for index, row in data.iterrows():
            if row['-log10(p_value)'] >= significance:
                ax.annotate(str(row[refSNP]), xy = (index, row['-log10(p_value)']))

    if export_path:
        plt.savefig(export_path)

    plt.show()

# Input:
# variants:
# phenotype_responses:
# This is considered to be the main method
def evaluate_variant_phenotype_significance(variants, phenotype_responses):

    # calculate statistics for Question 1
    alt_alleles = count_multiple_alternate_alleles(variants)
    print(f'The number of variants with >1 alt alleles: {alt_alleles}')
    passing_alleles = count_passing_variants(variants)
    print(f'The number of variants that passed all filters: {passing_alleles}')
    snps_count = count_snp_variants(variants)
    print(f'The number of variants that are SNPs: {snps_count}')
    indels_count = count_indel_variants(variants)
    print(f'The number of variants that are INDELs: {indels_count}')
    deletions_count = count_deletion_indels(variants)
    print(f'The number of variants that are deletions (subset of INDELS): {deletions_count}')
    
    snp_phenotype_p_values = []
    
    # calculate the 
    for variant in variants:
        snp_phenotype_frequencies = count_phenotype_responses_for_variant(variant, phenotype_responses)
        chi_square, p_value = variant_chi_square_analysis(snp_phenotype_frequencies['genotype_freq_table'])

        gft = snp_phenotype_frequencies['genotype_freq_table']['controls']

        hwe_chi_square, hwe_p_value = calculate_hwe(gft['aa'], gft['Aa'], gft['AA'])

        snp_phenotype_p_value = snp_phenotype_frequencies.copy()
        del snp_phenotype_p_value['genotype_freq_table']
        snp_phenotype_p_value['p_value'] = p_value
        snp_phenotype_p_value['hwe_p_value'] = hwe_p_value
        snp_phenotype_p_values.append(snp_phenotype_p_value)

    chromosome_key = 'CHROM'
    p_value_key = 'p_value'
    hwe_p_value_key = "hwe_p_value"

    statistical_significance_threshold = 5e-8

    hwe_significant_variants_count = count_hwe_significant_variants(snp_phenotype_p_value, statistical_significance_threshold)
    print(f'The number of variants that are HWE significant are: {hwe_significant_variants_count}')

    record_significant_p_values(
        snp_phenotype_p_values,
        chromosome_key,
        p_value_key,
        hwe_p_value_key,
        statistical_significance_threshold,
        "significant_variants.txt"
    )
    
    generate_manhattan_plot(
        pd.DataFrame(snp_phenotype_p_values),
        chromosome_key=chromosome_key,
        p_value_key=p_value_key,
        significance= -math.log10(statistical_significance_threshold),
        export_path="genotype_significance.png"
    )


###################################
# MAIN
###################################

if __name__ == '__main__':
    # setup the arg parser
    parser = argparse.ArgumentParser(description='GWAS Analysis')

    # add two required arguments. The file paths for the phenotype responses and VCF, respectively
    parser.add_argument('--phen_responses', dest='phen', required=True,
                    help='a two column (no header) text file containing phenotype responses per subject')
    parser.add_argument('--vcf', dest='vcf', required=True,
                    help='The VCF file that contains genotype data for variants')
    args = parser.parse_args()

    # read the phenotype response file and extract all entries
    phenotypes = read_phenotype_responses_from_file(args.phen)
    # read the vcf file and extract all entries
    vcf_data = extract_vcf_entries(args.vcf)
    # calculate all of the statistics and determine statistically significant variants
    evaluate_variant_phenotype_significance(vcf_data, phenotypes)
'''