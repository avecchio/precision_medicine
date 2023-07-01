import pandas as pd

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


interested_snps = [
    '16:48258198',
    #'rs17822931',
    '16:31107689',
    #'rs9923231',
    '19:45411941',
    #'rs429358',
    '19:45412079'
    #'rs7412'
]

avatar = 'genome_Isaac_Winters_v5_Full_20180326104258_GENOTYPE'
snps = read_vcf_file('GWAS_data.vcf')
print('ID REF ALT GENOTYPE')
for snp in snps:
    if snp['ID'] in interested_snps:
        id = snp['ID']
        ref = snp['REF']
        alt = snp['ALT']
        genotype = snp[avatar].split(":")[0]
        print(f'{id} {ref} {alt} {genotype}')
