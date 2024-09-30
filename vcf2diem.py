"""vcf2diem.py

Usage: 
 vcf2diem.py -v <FILE> [-h -l -n <INT>]

Options:
 -v, --vcf <FILE>                       VCF file
 -l, --location                         Print SNP location
 -n, --non_callable_limit <INT>         Maximum number of noncallable genotypes per site (default = no limit)
 -h, --help                             Print this message
"""

# Example Command
# mamba install -c conda-forge -c bioconda numpy pandas scikit-allel tqdm docopt
# python vcf2diem.py -v vcf_file.vcf.gz

import numpy as np
import pandas as pd
import sys, os
import allel
from tqdm import tqdm
from timeit import default_timer as timer
from docopt import docopt

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


def get_chromosome_data(vcf_dict, chromosome_name):

    is_chrom_array = vcf_dict["variants/CHROM"] == chromosome_name

    is_SNP_array = vcf_dict["variants/is_snp"]

    if isinstance(vcf_dict["variants/NUMALT"][0], int):
        numalt_array = vcf_dict["variants/NUMALT"]
        mask_array = (
            (numalt_array == 1) & (is_SNP_array == True) & (is_chrom_array == True)
        )
    else:
        mask_array = (is_SNP_array == True) & (is_chrom_array == True)

    snp_pos = vcf_dict["variants/POS"][mask_array]
    snp_gts = vcf_dict["calldata/GT"][mask_array]
    snp_ga = allel.GenotypeArray(snp_gts)
    snp_ga[snp_ga > 1] = -2
    snp_ga[snp_ga < 0] = -2
    summed_snp_ga = np.sum(snp_ga, axis=2)
    df = pd.DataFrame(summed_snp_ga)
    df = df.replace([i for i in range(-1,-10,-1)], "_")
    S = pd.DataFrame(["S"] * df.shape[0])
    pos = pd.DataFrame(snp_pos)
    df = pd.concat([pos, S, df], axis=1)
    colnames = list(df.columns)
    colnames[0]='pos'
    colnames[1]='S'
    df.columns=colnames
    non_singletons_list = non_singletons(df)
    df['non_singletons'] = non_singletons_list
    df = df.loc[df.non_singletons, :]
    return df


def write_diem(df, chromosome_name):
    df.drop(columns=['non_singletons'], inplace=True)
    if print_pos:
        np.savetxt(
            "./diem_files/snp_pos/" + str(chromosome_name) + ".snp_pos.diem.txt",
            df['pos'].values,
            fmt="%s",
            delimiter="",
        )

    
    df.drop(columns=['pos'], inplace=True)
    np.savetxt(
        "./diem_files/diem_input/" + str(chromosome_name) + ".diem.txt",
        df.values,
        fmt="%s",
        delimiter="",
    )
    print("Chromosome: " + str(chromosome_name) + "...")

def non_singletons(df):
    new_df = df.drop(columns=['pos'])
    new_df.drop(columns=['S'], inplace=True)
    unc_count_arr = (new_df.to_numpy() == '_').sum(axis=1)
    new_df = new_df.replace("_", 0)
    row_sum_arr = new_df.sum(axis=1)
    max_sum_arr = [2*(new_df.shape[1]-unc_count) for unc_count in unc_count_arr]
    non_singletons = [False if (x<=1) or (x>=max_sum-1) else True for x, max_sum in zip(row_sum_arr, max_sum_arr)]

    if args["--non_callable_limit"]:
        limit = int(args["--non_callable_limit"])
        unc_bool_arr = unc_count_arr <= limit
        non_singletons_pass = np.logical_and(non_singletons, unc_bool_arr)
        return non_singletons_pass

    else:
        return non_singletons

if __name__ == "__main__":

    args = docopt(__doc__)

    vcf_f = args["--vcf"]

    if args["--location"]:
        print_pos = True
    else:
        print_pos = False

    try:
        start_time = timer()

        path = os.path.join("./diem_files/diem_input/")
        os.makedirs(path, exist_ok=True)

        if print_pos:
            path = os.path.join("./diem_files/snp_pos/")
            os.makedirs(path, exist_ok=True)

        vcf_dict = load_vcf(vcf_f=vcf_f)

        all_contigs = np.unique(vcf_dict["variants/CHROM"])
        chromosome_names = np.array([x for x in all_contigs if 'contig' not in x])

        for chromosome in chromosome_names:
            chromosome_data = get_chromosome_data(
                vcf_dict=vcf_dict, chromosome_name=chromosome
            )
            write_diem(df=chromosome_data, chromosome_name=chromosome)

        print("[+] Total runtime: %.3fs" % (timer() - start_time))

    except KeyboardInterrupt:
        sys.stderr.write(
            "\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time)
        )
        sys.exit(-1)
