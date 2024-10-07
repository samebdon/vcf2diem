Usage
```
  vcf2diem.py -v <FILE> [-h -n -l <INT> -f <STR> -m -c]

Options:
 -v, --vcf <FILE>                       VCF file
 -n, --no_annotations                   Suppress writing SNP annotations
 -l, --non_callable_limit <INT>         Maximum number of noncallable genotypes allowed per site (default = no limit)
 -f, --contig-filter-string <STR>       String identifying contigs to ignore (default = None)
 -m, --missing-homs                     Include sites missing a ref and/or alt homozygote (default = False)
 -c, --chunks                           EXPERIMENTAL: Split diem_files/diem_input/per_chromosome into chunked files (requires split_into_chunks.sh in the same directory as vcf2diem.py)
 -h, --help                             Print this message
"""

mamba install -c conda-forge -c bioconda numpy pandas scikit-allel tqdm docopt
./vcf2diem.py -v your_vcf.gz
```
