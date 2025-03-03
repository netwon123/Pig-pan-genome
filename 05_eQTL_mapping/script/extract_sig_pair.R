library(data.table)
library(dplyr)
ARGS <- commandArgs(trailingOnly = TRUE)
tissue = ARGS[1]
cat("Processing tissue: ", tissue, "...\n")
signif_pair <- function(tissue){
  fdrdf=fread(paste0("joint.cis_qtl_mt.txt.gz"))
  fdrdf=fdrdf[,c("pheno_id", "pval_g1_threshold")]
  for (chr in c(1:18)){
    df=fread(paste0(tissue,".cis_qtl_pairs.", chr,".txt.gz"))
    DF = merge(df, fdrdf, by="pheno_id")
    DF=filter(DF, pval_g1 < pval_g1_threshold)
    if (nrow(DF) > 0){
      fwrite(DF, paste0( tissue,".signifpairs",chr,".txt.gz"), sep = "\t", row.names=F, compress="gzip", quote=F)
    }else{
      cat("Chromosome ", chr, "does not have significant QTLs...\n")
    }
    rm(df)
    rm(DF)
  }
}
signif_pair(tissue)
DF=data.frame()
for (chr in c(1:18)){
  signif_file=paste0( tissue,".signifpairs",chr,".txt.gz")
  if (file.exists(signif_file)){
    df = fread(signif_file, sep = "\t", header=T)
    DF=rbind(DF, df)
    system(paste0("rm -f ", signif_file))
  }
}
fwrite(DF, paste0( tissue,".signifpairs.txt.gz"), sep = "\t", row.names=F, compress="gzip", quote=F)

cat("Done!\n")
