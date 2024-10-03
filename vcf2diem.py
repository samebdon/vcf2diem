#!/usr/bin/env python3

"""vcf2diem.py

Usage: 
 vcf2diem.py -v <FILE> [-h -n -l <INT>]

Options:
 -v, --vcf <FILE>                       VCF file
 -n, --no_location                      Turn off SNP location printing
 -l, --non_callable_limit <INT>         Maximum number of noncallable genotypes per site (default = no limit)
 -h, --help                             Print this message
"""

# Example Command
# mamba install -c conda-forge -c bioconda numpy pandas scikit-allel tqdm docopt
# python vcf2diem.py -v vcf_file.vcf.gz

import numpy as np
import pandas as pd
import sys, os
import allel
import random
from tqdm import tqdm
from timeit import default_timer as timer
from docopt import docopt

class GenotypeData:
    def __init__(
        self,
        chromosome,
        vcf_dict,
    ):
        self.chromosome = chromosome
        self.vcf_dict = vcf_dict
        self.mask_array = None
        self.genotype_array = None
        self.most_common_alleles = None
        self.second_most_common_alleles = None
        self.positions = None

    def get_mask_array(self, biallelic=False):
        is_chrom_array = self.vcf_dict["variants/CHROM"] == self.chromosome
        is_SNP_array = self.vcf_dict["variants/is_snp"]

        mask_array = (is_SNP_array == True) & (is_chrom_array == True)
        if isinstance(self.vcf_dict["variants/NUMALT"][0], int) & biallelic:
            numalt_array = self.vcf_dict["variants/NUMALT"]
            mask_array = mask_array & (numalt_array == 1) # here i am taking sites with only 1 alt allele
            
        self.mask_array = mask_array

    def get_genotype_array(self):
        snp_gts = self.vcf_dict["calldata/GT"][self.mask_array]
        self.genotype_array = allel.GenotypeArray(snp_gts)

    def get_allele_order(self):
        """
        Currently doesnt account for ties just uses what numpy picks
        Tossing a coin would be preferable
        Maybe possible with a sorted key?
        """
        acs = self.genotype_array.count_alleles()

        #can i use the sorted function to get an output like argsort?
        self.acs = np.flip(acs.argsort(), axis=1) 
        #print(acs)
        #print(acs.argsort())
        #print(np.lexsort((acs, lambda v: (v, random.random())), axis=1))
        #sorted(acs, key=lambda v: (v, random.random()))


    def map_alleles(self):
        mapping = self.acs
        if mapping.shape[1] > 2:
            mapping[:,2:] = -2
        self.genotype_array = self.genotype_array.map_alleles(mapping)

    def get_pos(self):
        self.positions =  self.vcf_dict["variants/POS"][self.mask_array]
    
    def convert_to_diem_df(self, limit=None):
        """
        Assumes the genotype array is coded so that 0 is the most frequent allele
        and 1 is the second most frequent allele. 

        All other alleles are excluded!
        """
        snp_ga = self.genotype_array
        snp_ga[snp_ga > 1] = -2
        snp_ga[snp_ga < 0] = -2
        summed_snp_ga = np.sum(snp_ga, axis=2)

        df = pd.DataFrame(summed_snp_ga)
        df.replace([i for i in range(-1,-10,-1)], "_", inplace=True)
        S = pd.DataFrame(["S"] * df.shape[0], columns=['S'])
        pos = pd.DataFrame(self.positions, columns=['pos'])
        df = pd.concat([pos, S, df], axis=1)
        exclusions = get_exclusions(df, limit)

        return df.loc[~exclusions, :]

def load_vcf(vcf_f):
    query_fields = [
        "samples",
        "calldata/GT",
        "variants/CHROM",
        "variants/POS",
        "variants/NUMALT",
        "variants/is_snp",
    ]

    vcf_dict = allel.read_vcf(vcf_f, fields=query_fields)

    return vcf_dict

def write_diem(df, chromosome_name, print_pos=True):
    if print_pos:
        np.savetxt(
            f"./diem_files/snp_pos/{str(chromosome_name)}.snp_pos.diem.txt",
            df['pos'].values,
            fmt="%s",
            delimiter="",
        )
    
    df.drop(columns=['pos'], inplace=True)
    np.savetxt(
        f"./diem_files/diem_input/{str(chromosome_name)}.diem.txt",
        df.values,
        fmt="%s",
        delimiter="",
    )
    print(f"Chromosome: {str(chromosome_name)} written")

def get_exclusions(df, limit=None):
    new_df = df.drop(columns=['pos'])
    new_df.drop(columns=['S'], inplace=True)
    unc_count_arr = (new_df.to_numpy() == '_').sum(axis=1)
    new_df = new_df.replace("_", 0)

    singletons = get_singletons_list(new_df, unc_count_arr, limit)
    no_homozygotes = np.sum(new_df==0, axis=1) + np.sum(new_df==2, axis=1) == 0
    invariants = new_df.sum(axis=1) == new_df.shape[1]*2

    exclusions = np.logical_and(singletons, no_homozygotes)
    exclusions = np.logical_and(exclusions, invariants)

    if limit:
        limit = int(limit)
        unc_bool_arr = unc_count_arr <= limit
        exclusions = np.logical_and(exclusions, unc_bool_arr)

    return exclusions

def get_singletons_list(df, unc_count_arr, limit=None):
    row_sum_arr = df.sum(axis=1)
    max_sum_arr = [2*(df.shape[1]-unc_count) for unc_count in unc_count_arr]
    return [True if (x<=1) or (x>=max_sum-1) else False for x, max_sum in zip(row_sum_arr, max_sum_arr)]

def main():
    args = docopt(__doc__)

    vcf_f = args["--vcf"]

    try:
        start_time = timer()

        path = os.path.join("./diem_files/diem_input/")
        os.makedirs(path, exist_ok=True)

        if args["--no_location"]:
            print_pos = False
        else:
            print_pos = True
            path = os.path.join("./diem_files/snp_pos/")
            os.makedirs(path, exist_ok=True)

        vcf_dict = load_vcf(vcf_f=vcf_f)
        print(f'Loaded VCF file: {vcf_f}')

        all_contigs = np.unique(vcf_dict["variants/CHROM"])
        chromosome_names = [x for x in all_contigs if 'CAM' not in x] #need better contig filter

        for chromosome in chromosome_names:

            chromosome_genotype_data = GenotypeData(chromosome, vcf_dict)
            chromosome_genotype_data.get_mask_array()
            chromosome_genotype_data.get_genotype_array()
            chromosome_genotype_data.get_allele_order()
            chromosome_genotype_data.map_alleles() #sets most common and second most common alleles to 0 and 1, and everything else to -2
            chromosome_genotype_data.get_pos()
            diem_df = chromosome_genotype_data.convert_to_diem_df(args['--non_callable_limit'])
            write_diem(diem_df, chromosome, print_pos)

    except KeyboardInterrupt:
        sys.stderr.write(
            "\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time)
        )
        sys.exit(-1)

if __name__ == "__main__":
    main()
