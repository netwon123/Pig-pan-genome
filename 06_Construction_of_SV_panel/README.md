#  SNP & SV Processing and Imputation Pipeline

This pipeline describes the workflow for **SNP & SV data processing, filtering, merging, phasing, and imputation**.  
It is based on `bcftools`, `vcftools`, `plink`, and `beagle` with Slurm job submission.  

---

## 1. SNP Processing

### 1.1 Simplify VCF
```bash
# Remove QUAL, FILTER, INFO, keep only FORMAT/GT
bcftools annotate --threads 2 --remove QUAL,FILTER,INFO,^FORMAT/GT \
    bq_298.genome.snp.bqsr_snps.vcf.gz \
    -Oz -o bq_298.genome.snp.simplfied.vcf.gz
```

### 1.2 Convert single haploid genotypes `.` → diploid `./.`
```bash
zcat bq_298.genome.snp.simplfied.vcf.gz | \
    awk '{for(i = 1; i <= NF; i++) if(i >= 10 && i <= 2755 && $i == ".") {$i = "./."} {print $0}}' \
    | tr -s ' ' '\t' | bgzip -c > bq_298.genome.snp.bi-allele.vcf.gz
tabix -p vcf bq_298.genome.snp.bi-allele.vcf.gz
```

### 1.3 Merge Multiple SNP VCFs
```
bcftools merge --threads 10 bq_298.genome.snp.bi-allele.vcf.gz \
    bq_42.genome.bi-allele.vcf.gz \
    ygq_953.genome.snps.vcf.gz \
    yl_gte_1602.genome.bi-allele.vcf.gz \
    zls_207.genome.bi-allele.vcf.gz \
    jx273_snpindel_v2.vcf.gz \
    -Oz -o ./40.original_vcf_data/panel_3376.genome.snps.vcf.gz
tabix -p vcf panel_3376.genome.snps.vcf.gz
```

### 1.4 Update Header and Filter Samples
```bash
# Update VCF header
bcftools reheader -h header.txt panel_3376.genome.snps.vcf.gz -o panel_3376.genome.snps.reheadered.vcf.gz
tabix -p vcf panel_3376.genome.snps.reheadered.vcf.gz

# Exclude jx_273 samples
bcftools view -S ^jx_273.sample.list --threads 10 panel_3376.genome.snps.reheadered.vcf.gz -Oz -o panel_3103.genome.snps.vcf.gz

# Sort samples
bcftools query -l panel_3103.genome.snps.vcf.gz | sort > panel_3103.samples.sorted.txt
bcftools view -S panel_3103.samples.sorted.txt --threads 10 panel_3103.genome.snps.vcf.gz -Oz -o panel_3103.genome.snps.sorted.vcf.gz
```

### 1.5 PLINK Conversion and QC
```bash
# Convert to PLINK
plink --vcf panel_3103.genome.snps.sorted.vcf.gz --vcf-half-call m \
    --allow-extra-chr --double-id --make-bed --out panel_3103.genome.snps --threads 10

# Heterozygosity
plink --bfile panel_3103.genome.snps --allow-extra-chr --het --out panel_3103.genome.snps.het --threads 10 --memory 600000

# IBD
plink --bfile panel_3103.genome.snps --genome --out panel_3103.genome.snps.ibd_output --threads 10 --memory 600000
```

---

## 2. SV Processing

### 2.1 Extract INS and DEL Variants
```bash
bcftools view -i 'SVTYPE="INS" || SVTYPE="DEL"' panel_all3061_sv.rname2.1-18.anno.vcf.gz -Oz -o panel_3061.chr1-18.sv.ins.del.vc.gz
tabix -p vcf panel_3061.chr1-18.sv.ins.del.vc.gz
```

### 2.2 Filter SVs
```bash
bcftools view -i 'SVLEN>=50' panel_3061.chr1-18.sv.ins.del.vc.gz | \
bcftools view -i 'INFO/SVLEN <= 30000000' | \
bcftools view -e 'AC==1 || AC==2' -S 2337_panel.sample.txt -o panel_2337.chr1-18.sv.filtered.vc.gz

vcftools --gzvcf panel_2337.chr1-18.sv.filtered.vc.gz --max-missing 0.8 --recode --recode-INFO-all --out panel_2337.chr1-18.sv.filtered.miss
bgzip panel_2337.chr1-18.sv.filtered.miss.recode.vcf
tabix -p vcf panel_2337.chr1-18.sv.filtered.miss.recode.vcf.gz

bcftools view -m2 -M2 -i 'INFO/SVLEN >50' panel_2337.chr1-18.sv.filtered.miss.recode.vcf.gz -Oz -o panel_2337.chr1-18.sv.only_bi.vcf.gz
tabix -p vcf panel_2337.chr1-18.sv.only_bi.vcf.gz
```

---

## 3. Merge SNP and SV Per Chromosome

```bash
#!/bin/bash
#SBATCH --job-name=panel_2337
#SBATCH --partition=low,big,amd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=200G
#SBATCH --error=chr%a_%j.err
#SBATCH --output=chr%a_%j.out
#SBATCH --array=14

CHROM=$SLURM_ARRAY_TASK_ID

# SNP
bcftools view --threads 10 -v snps -m2 -M2 -S 2337_panel.sample.sorted.txt \
    -r ${CHROM} panel_3376.genome.snps.reheadered.vcf.gz \
    -Oz -o panel_2337.chr${CHROM}.snp.sample_sorted.vcf.gz

# SV
bcftools view --threads 10 -r $CHROM ../39.new_sv/panel_2337.chr1-18.sv.only_bi.vcf.gz -Oz -o panel_2337.chr${CHROM}.sv.vcf.gz
tabix -p vcf panel_2337.chr${CHROM}.sv.vcf.gz

# Merge
bcftools concat --threads 10 -a -d snps panel_2337.chr${CHROM}.snp.sample_sorted.vcf.gz \
    panel_2337.chr${CHROM}.sv.vcf.gz -Oz -o panel_2337.chr${CHROM}.snp.sv.vcf.gz
bcftools sort panel_2337.chr${CHROM}.snp.sv.vcf.gz -Oz -o panel_2337.chr${CHROM}.snp.sv.sorted.vcf.gz
tabix -p vcf panel_2337.chr${CHROM}.snp.sv.sorted.vcf.gz
```

---

## 4. Phasing & Self-Imputation with Beagle

```bash
java -Xmx190g -jar beagle.jar \
    gt=panel_2337.chr${CHROM}.snp.sv.sorted.vcf.gz out=panel_2337.chr${CHROM}.snp.sv.phased impute=true window=4.0 nthreads=20 ne=100
```

Reheader:
```bash
bcftools reheader -h header_panel_2337.txt -o panel_2337.chr${CHROM}.snp.sv.phased.reheadered.vcf.gz panel_2337.chr${CHROM}.snp.sv.phased.vcf.gz
tabix -p vcf panel_2337.chr${CHROM}.snp.sv.phased.reheadered.vcf.gz
```

SNP filtering:
```bash
vcftools --gzvcf panel_2337.chr${CHROM}.snp.phased.vcf.gz \
    --min-alleles 2 --max-alleles 2 \
    --max-missing 0.9 --maf 0.01 \
    --recode --recode-INFO-all --out panel_2337.chr${CHROM}.snp.phased
bgzip panel_2337.chr${CHROM}.snp.phased.recode.vcf
tabix -p vcf panel_2337.chr${CHROM}.snp.phased.recode.vcf.gz
```

Merge SNP+SV again:
```bash
bcftools concat -a -d snps --threads 20 panel_2337.chr${CHROM}.snp.phased.recode.vcf.gz \
    panel_2337.chr${CHROM}.sv.phased.vcf.gz \
    -Oz -o panel_2337.chr${CHROM}.snp.sv.phased.merged.vcf.gz
bcftools sort panel_2337.chr${CHROM}.snp.sv.phased.merged.vcf.gz -Oz -o panel_2337.chr${CHROM}.snp.sv.phased.final_sorted.vcf.gz
tabix -p vcf panel_2337.chr${CHROM}.snp.sv.phased.final_sorted.vcf.gz
```

---

## 5. HiFi SNP Imputation

```bash
#!/bin/bash
#SBATCH --job-name=hifi_imputed
#SBATCH --partition=low,big,amd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --error=chr%a_%j.err
#SBATCH --output=chr%a_%j.out
#SBATCH --array=16,17,18

CHROM=$SLURM_ARRAY_TASK_ID

java -jar beagle.jar \
    gt=$hifi_dir/chip.hifi_18.snp.chr${CHROM}.vcf.gz \
    ref=$dir1/panel_2337.chr${CHROM}.snp.sv.phased.final_sorted.vcf.gz \
    out=$dir3/hifi_18.panel_2337.chr${CHROM}.snp.sv.imputed \
    impute=true ap=true gp=true seed=1234 ne=1000 nthreads=10
```

Post-processing:
```bash
# Reheader
bcftools reheader -h header.txt --threads 10 hifi_18.panel_2337.chr${CHROM}.snp.sv.imputed.vcf.gz -o hifi_18.panel_2337.chr${CHROM}.snp.sv.imputed.reheadered.vcf.gz

# Add MAF
bcftools +fill-tags hifi_18.panel_2337.chr${CHROM}.snp.sv.imputed.reheadered.vcf.gz --threads 10 -Oz -o hifi_18.panel_2337.chr${CHROM}.snp.sv.imputed.maf.vcf.gz -- -t MAF

# Extract DR2
bcftools query -f '%INFO/DR2\n' hifi_18.panel_2337.chr${CHROM}.snp.sv.imputed.maf.vcf.gz > DR2.hifi_18.panel_2337.chr${CHROM}.snp.sv.imputed.txt
```
