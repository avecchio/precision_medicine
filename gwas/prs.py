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

phenotype_map = read_phenotype_file("PRS_phen.txt")
snps = read_vcf_file('GWAS_data.vcf')

genotypes = {}
for snp in snps:

    for sample in phenotype_map:
        phenotype = phenotype_map[sample]
        genotype = snp[sample].split(":")[0]
