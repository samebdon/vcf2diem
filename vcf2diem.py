#!/usr/bin/env python3

"""vcf2diem.py

Usage: 
 vcf2diem.py -v <FILE> [-h -n -l <INT> -f <STR> -m -c <INT>]

Options:
 -v, --vcf <FILE>                       VCF file
 -n, --no_annotations                   Suppress writing SNP annotations
 -l, --non_callable_limit <INT>         Maximum number of noncallable genotypes allowed per site (default = no limit)
 -f, --contig-filter-string <STR>       String identifying contigs to ignore (default = None)
 -m, --missing-homs                     Include sites missing a ref and/or alt homozygote (default = False)
 -c, --chunks <INT>                      Split diem_files/diem_input/per_chromosome into chunked files (used /bin/bash)
 -h, --help                             Print this message
"""

# Example Command
# mamba install -c conda-forge -c bioconda numpy pandas scikit-allel tqdm docopt
# ./vcf2diem.py -v vcf_file.vcf.gz

import numpy as np
import pandas as pd
import sys, os
import allel
import random
import subprocess
from tqdm import tqdm
from timeit import default_timer as timer
from docopt import docopt

"""
Feature creeps:

- Ref and alt allele in annotation file
- Reason for exclusion column
- User-defined chunk number
"""

class GenotypeData:
    def __init__(
        self,
        chromosome,
        vcf_dict,
    ):
        self.chromosome = chromosome
        self.vcf_dict = vcf_dict
        self.mask_array = None
        self.nucleotide_array = None
        self.most_frequent_nucleotides = None
        self.genotype_array = None
        self.allele_order = None
        self.positions = None
        self.qual = None

    def get_mask_array(self, biallelic=False):
        """
        Only included SNPs
        Non SNPs are not written to excluded
        """

        is_chrom_array = self.vcf_dict["variants/CHROM"] == self.chromosome
        is_SNP_array = self.vcf_dict["variants/is_snp"]

        mask_array = (is_SNP_array == True) & (is_chrom_array == True)
        if isinstance(self.vcf_dict["variants/NUMALT"][0], int) & biallelic:
            numalt_array = self.vcf_dict["variants/NUMALT"]
            mask_array = mask_array & (
                numalt_array == 1
            )  # here i am taking sites with only 1 alt allele

        self.mask_array = mask_array

    def get_nucleotide_array(self):
        ref = self.vcf_dict['variants/REF'][:, None]
        alt = self.vcf_dict['variants/ALT']
        self.nucleotide_array = np.append(ref, alt, axis=1)[self.mask_array]

    def get_genotype_array(self):
        snp_gts = self.vcf_dict["calldata/GT"][self.mask_array]
        self.genotype_array = allel.GenotypeArray(snp_gts)

    def get_allele_order(self):
        """
        For each pos gets the indices of the sorted allele counts in descending order.
        Ties are broken by values drawn from a random normal distribution.
        """
        acs = self.genotype_array.count_alleles()
        rand_acs = np.random.randn(*acs.shape)
        self.allele_order = np.flip(np.lexsort((rand_acs, acs), axis=1))

    def get_most_frequent_nucleotides(self):
        self.most_frequent_nucleotides = np.take_along_axis(self.nucleotide_array,self.allele_order[:,:2], axis=1)

    def map_alleles(self):
        mapping = self.allele_order
        if mapping.shape[1] > 2:
            mapping[:, 2:] = -2
        self.genotype_array = self.genotype_array.map_alleles(mapping)

    def get_pos(self):
        self.positions = self.vcf_dict["variants/POS"][self.mask_array]

    def get_qual(self):
        self.qual = self.vcf_dict["variants/QUAL"][self.mask_array]

    def convert_to_diem_df(self, limit, exclude_missing_homs):
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
        df.replace([i for i in range(-1, -10, -1)], "_", inplace=True)
        exclusions = get_exclusions(df, limit, exclude_missing_homs)
        # Maybe exclusions could be a column with these codes

        new_cols = pd.DataFrame(
            {
                "pos": self.positions,
                "qual": self.qual,
                "S": ["S"] * df.shape[0],
                "exclusions": exclusions,
                "ref_allele": self.most_frequent_nucleotides[:,0],
                "alt_allele": self.most_frequent_nucleotides[:,1]
            }
        )

        df = pd.concat([new_cols, df], axis=1)

        return df


def load_vcf(vcf_f):
    query_fields = [
        "samples",
        "calldata/GT",
        "variants/CHROM",
        "variants/POS",
        "variants/NUMALT",
        "variants/is_snp",
        "variants/QUAL",
        "variants/REF",
        "variants/ALT"
    ]

    vcf_dict = allel.read_vcf(vcf_f, fields=query_fields)
    return vcf_dict


def write_samples(samples):
    np.savetxt(
        "./diem_files/sampleNames.diem.txt",
        samples,
        fmt="%s",
        delimiter="",
    )


def write_diem(df, chromosome_name, write_annotations=True):
    np.savetxt(
        f"./diem_files/diem_input/per_chromosome/{str(chromosome_name)}.diem_input.txt",
        df.loc[df["exclusions"] == False]
        .drop(columns=["pos", "qual", "exclusions", "ref_allele","alt_allele"])
        .values,
        fmt="%s",
        delimiter="",
    )
    print(f"Chromosome: {str(chromosome_name)} written")

    if write_annotations:
        df["chrom"] = chromosome_name
        df["start"] = df["pos"] - 1  ## VCF 1 based, BED 0 based

        df.loc[df["exclusions"] == False][["chrom", "start", "pos", "qual", "ref_allele","alt_allele"]].to_csv(
            f"./diem_files/annotations/included/per_chromosome/{str(chromosome_name)}.included.annotations.bed",
            sep="\t",
            header=None,
            index=None,
        )

        df.loc[df["exclusions"] == True][["chrom", "start", "pos", "qual"]].to_csv(
            f"./diem_files/annotations/excluded/per_chromosome/{str(chromosome_name)}.excluded.annotations.bed",
            sep="\t",
            header=None,
            index=None,
        )


def write_empty_files(chromosome_name, write_annotations=True):
    open(
        f"./diem_files/diem_input/per_chromosome/{str(chromosome_name)}.diem_input.txt",
        "a",
    ).close()
    print(f"Chromosome: {str(chromosome_name)} written (EMPTY)")

    if write_annotations:
        open(
            f"./diem_files/annotations/included/per_chromosome/{str(chromosome_name)}.included.annotations.bed",
            "a",
        ).close()
        open(
            f"./diem_files/annotations/excluded/per_chromosome/{str(chromosome_name)}.excluded.annotations.bed",
            "a",
        ).close()


def get_exclusions(df, limit=None, exclude_missing_homs=None):
    """
    Excludes:
        Singleton sites
        Sites where there is not one homozygote of each variant
        Invariant sites (all alt hom slips through in very few cases)
        Sites where the number of missing genotypes is above the user-defined threshold

    Invalid lines:
        All 0
        All 2
        All 0 and a 1
        All 2 and a 1

    If limit:
        Count of _ above limit

    If exclude missing homs:
        No 0
        No 2
    """

    df = df.copy(deep=True)
    unc_count_arr = (df.to_numpy() == "_").sum(axis=1)
    df = df.replace("_", np.nan)
    row_sum_arr = df.sum(axis=1)
    max_sum_arr = [2 * (df.shape[1] - unc_count) for unc_count in unc_count_arr]

    singletons = get_singletons(max_sum_arr, row_sum_arr)  # 1
    invariants = get_invariants(df, max_sum_arr)  # 1
    exclusions = np.logical_or(singletons, invariants)

    if exclude_missing_homs:
        no_ref_hom = np.sum(df == 0, axis=1) == 0
        no_alt_hom = np.sum(df == 2, axis=1) == 0
        missing_homs = np.logical_or(no_ref_hom, no_alt_hom)  # 2
        exclusions = np.logical_or(exclusions, missing_homs)

    if limit:
        limit = int(limit)
        unc_bool_arr = unc_count_arr > limit
        exclusions = np.logical_or(exclusions, unc_bool_arr)

    return exclusions


def get_singletons(max_sum_arr, row_sum_arr):
    return pd.Series(
        [
            True if (x == 1) or (x == max_sum - 1) else False
            for x, max_sum in zip(row_sum_arr, max_sum_arr)
        ]
    )


def get_invariants(df, max_sum_arr):
    return np.logical_or(df.sum(axis=1) == max_sum_arr, df.sum(axis=1) == 0)


def chunk(chr_path, chunk_path, inc_path, inc_chunk_path, num_chunks):
    num_total_sites = int(
        subprocess.run(
            f"wc -l {chr_path}/* | tail -n1 |  cut -f2 -d' '",
            shell=True,
            executable="/bin/bash",
            capture_output=True,
            text=True,
        ).stdout
    )
    chunksize = get_chunksize(num_chunks, num_total_sites)
    for in_path, out_path in [(chr_path, chunk_path), (inc_path, inc_chunk_path)]:
        subprocess.run(
            f"split -l {chunksize} -d <(cat {in_path}/*) {out_path}/chunk_",
            shell=True,
            executable="/bin/bash",
        )

def get_chunksize(c, n):
    return int(np.divide(n, c) + np.sum(np.remainder(n, c) > 0))

def main():
    args = docopt(__doc__)

    vcf_f = args["--vcf"]
    start_time = timer()

    try:

        chr_path = os.path.join("./diem_files/diem_input/per_chromosome")
        os.makedirs(chr_path, exist_ok=True)

        if args["--no_annotations"]:
            write_annotations = False
        else:
            write_annotations = True
            inc_path = os.path.join("./diem_files/annotations/included/per_chromosome")
            excl_path = os.path.join("./diem_files/annotations/excluded/per_chromosome")
            os.makedirs(inc_path, exist_ok=True), os.makedirs(excl_path, exist_ok=True)

        if args["--missing-homs"]:
            exclude_missing_homs = False
        else:
            exclude_missing_homs = True

        vcf_dict = load_vcf(vcf_f=vcf_f)
        print(f"Loaded VCF file: {vcf_f}")

        all_samples = np.array(vcf_dict["samples"])
        write_samples(all_samples)

        all_contigs = np.unique(vcf_dict["variants/CHROM"])
        if args["--contig-filter-string"]:
            chromosome_names = [
                x for x in all_contigs if str(args["--contig-filter-string"]) not in x
            ]
        else:
            chromosome_names = all_contigs.tolist()

        for chromosome in chromosome_names:
            chromosome_genotype_data = GenotypeData(chromosome, vcf_dict)
            chromosome_genotype_data.get_mask_array()
            chromosome_genotype_data.get_nucleotide_array()
            chromosome_genotype_data.get_genotype_array()
            try:
                chromosome_genotype_data.get_allele_order()
            except ValueError:
                write_empty_files(chromosome, write_annotations)
                continue
            chromosome_genotype_data.get_most_frequent_nucleotides()
            chromosome_genotype_data.map_alleles()  # sets most common and second most common alleles to 0 and 1, and everything else to -2
            chromosome_genotype_data.get_pos()
            chromosome_genotype_data.get_qual()
            diem_df = chromosome_genotype_data.convert_to_diem_df(
                args["--non_callable_limit"], exclude_missing_homs
            )
            write_diem(diem_df, chromosome, write_annotations)

        if args["--chunks"]:
            print("Chunking...")

            chunk_path = os.path.join("./diem_files/diem_input/per_chunk")
            inc_chunk_path = os.path.join("./diem_files/annotations/included/per_chunk")
            os.makedirs(chunk_path, exist_ok=True), os.makedirs(
                inc_chunk_path, exist_ok=True
            )

            chunk(chr_path, chunk_path, inc_path, inc_chunk_path, num_chunks=args["--chunks"])

    except KeyboardInterrupt:
        sys.stderr.write(
            "\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time)
        )
        sys.exit(-1)


if __name__ == "__main__":
    main()

