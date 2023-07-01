from scipy.stats import chi2_contingency
from mne.stats import fdr_correction, bonferroni_correction
import pandas as pd

# parse phenotype file and return phenotype dictionary
# where key:value is sample_id:phenotype
def read_phenotype_file(filename):
    df = pd.read_csv(filename, sep='\t', header=0)
    records = df.to_dict('records')
    phenotype_map = {}
    for record in records:
        phenotype_map[record['Sample']] = record['Phenotype']
    return phenotype_map

# read vcf file into a list of dictionaries using pandas
def read_vcf_file(filename):
    with open(filename, 'r') as fp:
        lines = fp.readlines()
    # traverse file ignoring comments and find start of data table
    start_index = 0
    for line in lines:
        if line[0:2] != '##':
            break
        start_index += 1
    # get columns
    columns = lines[start_index].replace("\n", "").split('\t')
    # get entries which is the next line after the column headers
    entries = lines[start_index + 1:]
    # split each line into array by tab delimiter
    split_entries = [entry.replace("\n", "").split('\t') for entry in entries]
    # construct dataframe and return array of dictionaries
    df = pd.DataFrame(split_entries, columns=columns)
    return df.to_dict('records')

# method split out to convert diploid genotype into representative allele
def genotype_to_allele(genotype):
    if genotype == '0|0':
        return 'AA'
    elif genotype in ['1|0', '0|1']:
        return 'Aa'
    elif genotype == '1|1':
        return 'aa'
    else:
        print('error')

# start of main code, read in all of the files
phenotype_map = read_phenotype_file("PRS_phen.txt")
snps = read_vcf_file('GWAS_data.vcf')

genotypes = {}
snp_pvals = []

# calculate pvalue for each snp
for snp in snps:
    snp_id = snp['ID']

    # construct base contingency table
    genotype_freq_table = {
        'cases' : {'aa': 0, 'Aa': 0, 'AA': 0},
        'controls': {'aa': 0, 'Aa': 0, 'AA': 0}
    }

    # get allele for each sample and store in contingency table
    for sample in phenotype_map:
        # get genotype and phenotype of sample
        phenotype = phenotype_map[sample]
        genotype = snp[sample].split(":")[0]

        # convert genotype to allele
        allele = genotype_to_allele(genotype)

        # store allele in contingency table
        cases = 1
        if phenotype == cases:
            genotype_freq_table['cases'][allele] += 1
        # else is control
        else:
            genotype_freq_table['controls'][allele] += 1

    # create dataframe from dictionary
    contingency_df = pd.DataFrame(genotype_freq_table).transpose()

    # attempt to perform chi square
    try:
        res = chi2_contingency(contingency_df, correction=False)
        # if successful, get pvalue and record
        snp_pvals.append({'id': snp_id, 'pval': res.pvalue})
    # else, columns or rows sum to zero. Set pvalue for snp to 1
    except:
        snp_pvals.append({'id': snp_id, 'pval': 1})

# extract pvalues for each snp
pvals = [snp_pval['pval'] for snp_pval in snp_pvals]

# perform bonferroni correction
bon_p_adjusted = bonferroni_correction(pvals, alpha=0.05)
# construct dictionary to record count of pvalues that are significant after correction
bon_tests = {}
for p in bon_p_adjusted[0]:
    pass_test = str(p)
    if pass_test not in bon_tests:
        bon_tests[pass_test] = 0
    bon_tests[pass_test] += 1

# print results
print('bonferonni')
print(bon_tests)

# perform false discovery rate correction
fdr_p_adjusted = fdr_correction(pvals, alpha=0.05, method='indep')
# construct dictionary to record count of pvalues that are significant after correction
fdr_tests = {}
for p in fdr_p_adjusted[0]:
    pass_test = str(p)
    if pass_test not in fdr_tests:
        fdr_tests[pass_test] = 0
    fdr_tests[pass_test] += 1

# print results
print('fdr')
print(fdr_tests)

