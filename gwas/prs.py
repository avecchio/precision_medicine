import pandas as pd
import numpy as np
import math
from sklearn import metrics

# helper method to read tsv file into a list of dictionaries
def read_tsv(filename):
    df = pd.read_csv(filename, sep='\t', header=0)
    return df.to_dict('records')

# parse phenotype file and return phenotype dictionary
# where key:value is sample_id:phenotype
def read_phenotype_file(filename):
    records = read_tsv(filename)
    phenotype_map = {}
    for record in records:
        phenotype_map[record['Sample']] = record['Phenotype']
    return phenotype_map

# read stats file and return dictionary
# where the key is the sample_id and the value is the stats as an embedded dictionary
def read_stats_file(filename):
    records = read_tsv(filename)
    stats_map = {}
    for record in records:
        stats_map[record['ID']] = record

    return stats_map

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

# method split out to convert diploid genotype into dosing value
def genotype_to_dosage(genotype):
    if genotype == '0|0':
        return 0
    elif genotype in ['1|0', '0|1']:
        return 1
    elif genotype == '1|1':
        return 2
    else:
        print('error')

# start of main code, read in all of the files
phenotype_map = read_phenotype_file("PRS_phen.txt")
snps = read_vcf_file('GWAS_data.vcf')
snp_stats = read_stats_file("PRS_GWAS_summ_stats.txt")
total_snps = len(snp_stats.keys())


# determine the PRS for each significance cutoff
all_prs_data = []
for cutoff in [0.01, 0.05, 0.1, 0.5]:
    prs_stats = {}

    filtered_snp_stats = [snp_stat for snp_stat in snp_stats.values() if snp_stat['p-values'] <= cutoff]
    
    filtered_snp_ids = [snp['ID'] for snp in filtered_snp_stats]
    filtered_snps = [snp for snp in snps if snp['ID'] in filtered_snp_ids]

    # for each sample, calculate the PRS
    for sample in phenotype_map:
        # get the phenotype for the sample
        phenotype = phenotype_map[sample]
        effect_sums = []
        # for each snp, multiply dosage by natural log of effect size
        for snp in filtered_snps:
            genotype = snp[sample].split(":")[0]
            dosage = genotype_to_dosage(genotype)
            snp_id = snp['ID']
            effect_size = snp_stats[snp_id]['OR']
            effect_sums.append(math.log(effect_size) * dosage)
        # calculate prs
        ploidy = 2
        prs = sum(effect_sums) / (total_snps * ploidy)
        # store prs and associated phenotype
        prs_stats[sample] = {
            'phenotype': phenotype,
            'prs': prs
        }

    # save phenotype and prs scores for each sample to later calculate auc
    phenotypes = []
    prs_scores = []
    for sample in prs_stats:
        phenotypes.append(prs_stats[sample]['phenotype'])
        prs_scores.append(prs_stats[sample]['prs'])
        # also save all prs scores and sample as per instructions
        all_prs_data.append({
            'sample': sample,
            'prs_score': prs_stats[sample]['prs']
        })

    # calculate AUC for each set of prs scores and associated phenotypes
    y = np.array(phenotypes)
    pred = np.array(prs_scores)
    fpr, tpr, thresholds = metrics.roc_curve(y, pred)
    auc = metrics.auc(fpr, tpr)
    print(f'{str(cutoff)}:{str(auc)}')

    # if the cutoff is 0.05, find the PRS value for the avatar and convert back to OR
    if cutoff == 0.05:
        avatar = 'genome_Isaac_Winters_v5_Full_20180326104258_GENOTYPE'
        avatar_or_risk = math.exp(prs_stats[avatar]['prs'])
        print('Avatar OR Risk:' + str(avatar_or_risk))

# write all prs data to a csv per the instructions
prs_scores_df = pd.DataFrame(all_prs_data)
prs_scores_df.to_csv('prs_scores.csv', index=False)