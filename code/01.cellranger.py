##进入工作目录
cd /WORK/Project/zengbo/Sc_RNA_pig
##加载软件路径
export PATH=/WORK/Project/zengbo/software/cellranger-5.0.1:$PATH
spaceranger count --help

#制作参考数据库文件
nohup cellranger mkref \
--genome=ref_SC_pig_ensembl111_add_MYH24 \
--fasta=/WORK/Project/zengbo/spatial/genome-database/pig-ref/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
--genes=/WORK/Project/zengbo/spatial/genome-database/pig-ref/ensembl-pig-add-MYH24.gtf
#这一步非常耗时，单线程大概8小时以上！
cellranger count --help
nohup cellranger count --id=counts_SC_liver \
 --fastqs=/WORK/Project/zengbo/Sc_RNA_pig/fastq/liver \
  --transcriptome=/WORK/Project/zengbo/Sc_RNA_pig/ref/ref_SC_pig_ensembl111_add_MYH24 \
 --localcores=40 --localmem=400 --no-bam \
 --sample=P_Liver \
 --expect-cells=10000
