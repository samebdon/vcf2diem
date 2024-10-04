Usage
```
  vcf2diem.py -v <FILE> [-h -n -l <INT> -f <STR> -r]

Options:
 -v, --vcf <FILE>                       VCF file
 -n, --no_location                      Suppress SNP location
 -l, --non_callable_limit <INT>         Maximum number of noncallable genotypes allowed per site (default = no limit)
 -f, --contig-filter-string <STR>       String identifying contigs to ignore (default = None)
 -r, --require-homs                     Ignore sites without a ref and alt homozygote (default = False)
 -h, --help                             Print this message

mamba install -c conda-forge -c bioconda numpy pandas scikit-allel tqdm docopt
./vcf2diem.py -v your_vcf.gz
```
