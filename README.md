# vcf-preprocessing-rvas
## Step 1: Filter variants and annotate
run `./variant_filter_annotate.sh input.vcf.gz filter1_output.vcf.gz snpeff_output.vcf.gz snpsift_output.vcf.gz`
## Step 2: Use GATK best practices filter to further filter
run `./hard_filter.py input.vcf output.vcf`
## Step 3: Determine sample ancestry and remove samples of non European ancestry
check instructions in ancestry_analysis_w_admixture.txt file
## Step 4: Generate intermediate mutations.tsv file
run `./gzvcf_to_mutations.py input.vcf.gz famfile outfile`
## Step 5: Generate QC stats and remove outlier samples
to generate QC stats run `./generate_qc_stats.py input_mutations_file cpgfile famfile outfile`
remove samples which are outside +/- three standard deviations for ts/tv, het/hom, perhaploidSNV, etc.
## Step 6: Perform SDT test
run `./sdt_test.py input_mutations_file famfile outfile`
