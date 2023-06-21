from scipy.stats import chisquare, chi2_contingency

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

import argparse
import sys

def read_tsv(filename):
    df = pd.read_csv(filename, sep='\t', header=0)
    return df.to_dict('records')

def read_phenotype_file(filename):
    records = read_tsv(filename)
    phenotype_map = {}
    for record in records:
        phenotype_map[record['Sample']] = record['Phenotype']
    return phenotype_map

def read_stats_file(filename):
    records = read_tsv(filename)
    stats_map = {}
    for record in records:
        stats_map[record['ID']] = record

    return stats_map

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

def genotype_to_dosage(genotype):
    if genotype == '0|0':
        return 0
    elif genotype in ['1|0', '0|1']:
        return 1
    elif genotype == '1|1':
        return 2
    else:
        print('error')

phenotype_map = read_phenotype_file("PRS_phen.txt")
snps = read_vcf_file('GWAS_data.vcf')
snp_stats = read_stats_file("PRS_GWAS_summ_stats.txt")
total_snp_stats = len(snp_stats.keys())

genotypes = {}

prs_stats = {}
for cutoff in [0.01, 0.05, 0.1, 0.5]:
    prs_scores = []
    filtered_snp_stats = [snp_stat for snp_stat in snp_stats.values() if snp_stat['p-values'] <= cutoff]
    filtered_snp_ids = [snp['ID'] for snp in filtered_snp_stats]
    filtered_snps = [snp for snp in snps if snp['ID'] in filtered_snp_ids]
    for sample in phenotype_map:
        phenotype = phenotype_map[sample]
        snp_sums = []
        for snp in filtered_snps:
            genotype = snp[sample].split(":")[0]
            dosage = genotype_to_dosage(genotype)
            snp_id = snp['ID']
            effect_size = snp_stats[snp_id]['OR']
            snp_sums.append(math.log(effect_size) * dosage)
        prs = sum(snp_sums) / (total_snp_stats * 2)
        prs_stats[sample] = {
            'phenotype': phenotype,
            'prs': prs
        }
        prs_scores.append(prs)    
    print(sorted(prs_scores))
    #print(prs_stats)
    #genotype = snp[sample].split(":")[0]
#for snp in snps:
