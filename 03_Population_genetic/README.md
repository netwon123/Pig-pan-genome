# Population Genetic Analyses

## PCA
```bash
plink --allow-extra-chr --threads 20   -bfile testacc_qc   --pca 2   --out testacc_qc
```

## Phylogeny Analysis
```bash
# Generate distance matrix
plink2 --bfile ont203.snp_QC   --out Distance_snp.out   --recode   --distance-matrix

# Convert IDs for PHYLIP format
awk '{print $2"           "}' 58inds_maf005_geno01_phased.mdist.id | cut -c1-10 > 58inds_maf005_geno01_phased_phylip.id 

# Merge IDs and distance matrix
paste -d '' 58inds_maf005_geno01_phased_phylip.id 58inds_maf005_geno01_phased.mdist | sed '1i 58' > 58inds_maf005_geno01_phased_phylip.mdist

# Tree construction
neighbor   # from phylip-3.697
```

## Admixture Analysis
```bash
for K in 2 3 4; do
  admixture -j16 --cv admix.bed $K | tee log${K}.out
done
```

## Fst Calculation
```bash
vcftools --vcf pop_miss03_maf005.recode.vcf   --weir-fst-pop ead.txt   --weir-fst-pop eaw.txt   --out fst_ead_eaw   --fst-window-size 1000000   --fst-window-step 500000
```
