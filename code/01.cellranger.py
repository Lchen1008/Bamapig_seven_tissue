## Enter the working directory
cd /WORK/Project/zengbo/Sc_RNA_pig

## Set the software path
export PATH=/WORK/Project/zengbo/software/cellranger-5.0.1:$PATH

spaceranger count --help

# Create a reference database file
nohup cellranger mkref \
--genome=ref_SC_pig_ensembl111_add_MYH24 \
--fasta=/WORK/Project/zengbo/spatial/genome-database/pig-ref/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
--genes=/WORK/Project/zengbo/spatial/genome-database/pig-ref/ensembl-pig-add-MYH24.gtf

# This step is very time-consuming, taking approximately 8 hours or more with a single thread!
cellranger count --help

nohup cellranger count --id=counts_SC_liver \
 --fastqs=/WORK/Project/zengbo/Sc_RNA_pig/fastq/liver \
  --transcriptome=/WORK/Project/zengbo/Sc_RNA_pig/ref/ref_SC_pig_ensembl111_add_MYH24 \
 --localcores=40 --localmem=400 --no-bam \
 --sample=P_Liver \
 --expect-cells=10000

