#! /usr/bin/bash

###eRNA Quantification for rareRCC patient samples###
#This script looks for enhancer regions and then aligns them against mRNA-seq reads to quantify transcribed eRNAs.

### Procedure ###

#1) Retrieve enhancer region in bed file from FANTOM5. 
curl -o F5.hg38.enhancers.bed.gz https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/F5.hg38.enhancers.bed.gz 
gzip -dc F5.hg38.enhancers.bed.gz > F5.hg38.enhancers.bed
#-----

#2) Retrieve GTF file from ESEMBL.
curl -o Homo_sapiens.GRCh38.104.gtf.gz http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz   
gzip -dc Homo_sapiens.GRCh38.104.gtf.gz > Homo_sapiens.GRCh38.104.gtf
#-----

#5) Add 'chr' in gtf file to make it consistent with bed.
grep -v "^#" Homo_sapiens.GRCh38.104.gtf | awk '{print "chr"$0}' > Homo_sapiens.GRCh38.104.chr.gtf
#-----

#6) Retain the coding regions in gtf file.
grep -E '#|transcript_biotype "protein_coding"' Homo_sapiens.GRCh38.104.chr.gtf > Homo_sapiens.GRCh38.104.tb.chr.gtf
grep -E '#|gene_biotype "protein_coding"' Homo_sapiens.GRCh38.104.chr.gtf > Homo_sapiens.GRCh38.104.gb.chr.gtf

grep -E '#|gene_biotype "lncRNA"' Homo_sapiens.GRCh38.104.chr.gtf > Homo_sapiens.GRCh38.104.LNC.gb.chr.gtf
grep -E '#|transcript_biotype "lncRNA"' Homo_sapiens.GRCh38.104.chr.gtf > Homo_sapiens.GRCh38.104.LNC.tb.chr.gtf

awk '$3=="transcript"' Homo_sapiens.GRCh38.104.tb.chr.gtf > Homo_sapiens.GRCh38.104.T.tb.chr.gtf
awk '$3=="transcript"' Homo_sapiens.GRCh38.104.gb.chr.gtf > Homo_sapiens.GRCh38.104.T.gb.chr.gtf



curl -o gencode.v38.long_noncoding_RNAs.gtf.gz http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.long_noncoding_RNAs.gtf.gz   
gzip -dc gencode.v38.long_noncoding_RNAs.gtf.gz > gencode.v38.long_noncoding_RNAs.gtf

curl -o ENCFF356LFX.bed.gz https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz   
gzip -dc ENCFF356LFX.bed.gz > ENCFF356LFX.bed
#GRCh38_unified_blacklist.bed

curl -o hg38.knownGene.gtf.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz
gzip -dc hg38.knownGene.gtf.gz > hg38.knownGene.gtf


# Refgene files
curl -o hg38.refGene.gtf.gz  http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz  
gzip -dc hg38.refGene.gtf.gz > hg38.refGene.gtf
# 1893860

curl -o MANE.GRCh38.v0.95.select_ensembl_genomic.gtf.gz https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v0.95.select_ensembl_genomic.gtf.gz
gzip -dc MANE.GRCh38.v0.95.select_ensembl_genomic.gtf.gz > MANE.GRCh38.v0.95.select_ensembl_genomic.gtf
# 507114

cat Homo_sapiens.GRCh38.104.chr.gtf | cut -f 1-5 | grep -w "gene" | cut -f 1,4,5 > genes.bed

cat Homo_sapiens.GRCh38.104.chr.gtf | cut -f 1-5 | grep -w "exon" | cut -f 1,4,5 > exons.bed

bedtools subtract -a genes.bed -b exons.bed > introns.bed

# for further filters
#https://www.gencodegenes.org/human/

# for RefSeq Reference Genome Annotation
# https://www.ncbi.nlm.nih.gov/genome/guide/human/
# 208040
#-----

#7) Remove regions from bed file 
bedtools intersect -a F5.hg38.enhancers.bed \
      -b Homo_sapiens.GRCh38.104.T.gb.chr.gtf \
      Homo_sapiens.GRCh38.104.T.tb.chr.gtf \
      Homo_sapiens.GRCh38.104.LNC.gb.chr.gtf \
      Homo_sapiens.GRCh38.104.LNC.tb.chr.gtf \
      gencode.v38.long_noncoding_RNAs.gtf \
      ENCFF356LFX.bed \
      hg38.knownGene.gtf \
      introns.bed \
      hg38.refGene.gtf \
      MANE.GRCh38.v0.95.select_ensembl_genomic.gtf \
      GRCh38_latest_genomic.gff \
      gencode.v38.chr_patch_hapl_scaff.basic.annotation.gtf \
      gencode.v38.chr_patch_hapl_scaff.annotation.gtf \
      chess2.2.gtf \
      chess2.2_assembly.gtf chess2.2_and_refseq.gtf -wa -v > Homo_sapiens.GRCh38.104.T-LNC.gb-tb_ALL_BL_knowngenes_Introns1_RefGene3_add_chess.chr.bed

#9) Convert bed to saf file.
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}'  Homo_sapiens.GRCh38.104.T-LNC.gb-tb_ALL_BL_knowngenes_Introns1_RefGene3_add_chess.chr.bed > Homo_sapiens.GRCh38.104.T-LNC.gb-tb_ALL_BL_knowngenes_Introns1_RefGene3_add_chess.chr.bed.saf
#-----

#10) Using FeatureCounts
featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/Homo_sapiens.GRCh38.104.T-LNC.gb-tb_ALL_BL_knowngenes_Introns1_RefGene3_add_chess.chr.bed.saf \
-F SAF -o counts_filterd_enhancer_15159_TCGA_chRCCpatients_v3.txt 098a03fe-6a1b-4248-834a-9fa70ea02000_gdc_realn_rehead.bam
#-----
# -p specifies that fragment will be counted instead of reads.
#"A read is said to overlap a feature if at least one read base is found to overlap the feature. For paired-end data, a fragment (or template) is said to overlap a feature if any of the two reads from that fragment is found to overlap the feature" 

#11) Take unnecessary columns out
#awk '{$2=$3=$4=$5=$6=""; print $0}' counts_filterd_only_TPC_rep1_and_2.txt
#-----

#12) Remove the first line which had the path to featureCounts command
sed -i '1d' counts_filterd_only_TPC_rep1_and_2.txt
#-----
#Use vim to remove path: 


#13) Run R script for using EdgeR
#-----

########################
coverageBed -abam C3L-01469-CPT0389960009-alig.bam -hist -b /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/Homo_sapiens.GRCh38.104.T-LNC.gb-tb_ALL_BL_knowngenes_Introns1_RefGene3_add_chess.chr.bed > eRNAs.bed.coverage
bedtools intersect -abam C3L-01469-CPT0389960009-alig.bam -b /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/Homo_sapiens.GRCh38.104.T-LNC.gb-tb_ALL_BL_knowngenes_Introns1_RefGene3_add_chess.chr.bed > eRNAs.bed.bam

awk '{print $3}' Homo_sapiens.GRCh38.104.chr.gtf | sort | uniq
CDS
exon
five_prime_utr
gene
Selenocysteine
start_codon
stop_codon
three_prime_utr
transcript

cat Homo_sapiens.GRCh38.104.chr.gtf | cut -f 1-5 | grep -w "CDS" | cut -f 1,4,5 > CDS.bed

cat Homo_sapiens.GRCh38.104.chr.gtf | cut -f 1-5 | grep -w "exon" | cut -f 1,4,5 > exon.bed

cat Homo_sapiens.GRCh38.104.chr.gtf | cut -f 1-5 | grep -w "five_prime_utr" | cut -f 1,4,5 > five_prime_utr.bed

cat Homo_sapiens.GRCh38.104.chr.gtf | cut -f 1-5 | grep -w "gene" | cut -f 1,4,5 > gene.bed

cat Homo_sapiens.GRCh38.104.chr.gtf | cut -f 1-5 | grep -w "Selenocysteine" | cut -f 1,4,5 > Selenocysteine.bed

cat Homo_sapiens.GRCh38.104.chr.gtf | cut -f 1-5 | grep -w "start_codon" | cut -f 1,4,5 > start_codon.bed

cat Homo_sapiens.GRCh38.104.chr.gtf | cut -f 1-5 | grep -w "stop_codon" | cut -f 1,4,5 > stop_codon.bed

cat Homo_sapiens.GRCh38.104.chr.gtf | cut -f 1-5 | grep -w "three_prime_utr" | cut -f 1,4,5 > three_prime_utr.bed

cat Homo_sapiens.GRCh38.104.chr.gtf | cut -f 1-5 | grep -w "transcript" | cut -f 1,4,5 > transcript.bed

awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}'  CDS.bed > CDS.bed.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}'  exon.bed > exon.bed.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}'  gene.bed > gene.bed.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}'  five_prime_utr.bed > five_prime_utr.bed.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}'  Selenocysteine.bed > Selenocysteine.bed.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}'  start_codon.bed > start_codon.bed.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}'  stop_codon.bed > stop_codon.bed.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}'  three_prime_utr.bed > three_prime_utr.bed.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}'  transcript.bed > transcript.bed.saf
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}'  introns.bed > introns.bed.saf



featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/CDS.bed.saf \
-F SAF -o counts_filterd_CDS_TCGA_chRCCpatients_v3.txt 098a03fe-6a1b-4248-834a-9fa70ea02000_gdc_realn_rehead.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/exon.bed.saf \
-F SAF -o counts_filterd_exon_TCGA_chRCCpatients_v3.txt 098a03fe-6a1b-4248-834a-9fa70ea02000_gdc_realn_rehead.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/gene.bed.saf \
-F SAF -o counts_filterd_gene_TCGA_chRCCpatients_v3.txt 098a03fe-6a1b-4248-834a-9fa70ea02000_gdc_realn_rehead.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/five_prime_utr.bed.saf \
-F SAF -o counts_filterd_five_prime_utr_TCGA_chRCCpatients_v3.txt 098a03fe-6a1b-4248-834a-9fa70ea02000_gdc_realn_rehead.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/Selenocysteine.bed.saf \
-F SAF -o counts_filterd_Selenocysteine_TCGA_chRCCpatients_v3.txt 098a03fe-6a1b-4248-834a-9fa70ea02000_gdc_realn_rehead.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/start_codon.bed.saf \
-F SAF -o counts_filterd_start_codon_TCGA_chRCCpatients_v3.txt 098a03fe-6a1b-4248-834a-9fa70ea02000_gdc_realn_rehead.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/stop_codon.bed.saf \
-F SAF -o counts_filterd_stop_codon_TCGA_chRCCpatients_v3.txt 098a03fe-6a1b-4248-834a-9fa70ea02000_gdc_realn_rehead.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/three_prime_utr.bed.saf \
-F SAF -o counts_filterd_three_prime_utr_TCGA_chRCCpatients_v3.txt 098a03fe-6a1b-4248-834a-9fa70ea02000_gdc_realn_rehead.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/transcript.bed.saf \
-F SAF -o counts_filterd_transcript_TCGA_chRCCpatients_v3.txt 098a03fe-6a1b-4248-834a-9fa70ea02000_gdc_realn_rehead.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/introns.bed.saf \
-F SAF -o counts_filterd_transcript_TCGA_chRCCpatients_v3.txt 098a03fe-6a1b-4248-834a-9fa70ea02000_gdc_realn_rehead.bam

##############################################################################################################

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/CDS.bed.saf \
-F SAF -o counts_filterd_CDS_CPTAC_chRCCpatients_v3.txt C3L-01469-CPT0389960009-alig.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/exon.bed.saf \
-F SAF -o counts_filterd_exon_CPTAC_chRCCpatients_v3.txt C3L-01469-CPT0389960009-alig.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/gene.bed.saf \
-F SAF -o counts_filterd_gene_CPTAC_chRCCpatients_v3.txt C3L-01469-CPT0389960009-alig.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/five_prime_utr.bed.saf \
-F SAF -o counts_filterd_five_prime_utr_CPTAC_chRCCpatients_v3.txt C3L-01469-CPT0389960009-alig.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/Selenocysteine.bed.saf \
-F SAF -o counts_filterd_Selenocysteine_CPTAC_chRCCpatients_v3.txt C3L-01469-CPT0389960009-alig.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/start_codon.bed.saf \
-F SAF -o counts_filterd_start_codon_CPTAC_chRCCpatients_v3.txt C3L-01469-CPT0389960009-alig.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/stop_codon.bed.saf \
-F SAF -o counts_filterd_stop_codon_CPTAC_chRCCpatients_v3.txt C3L-01469-CPT0389960009-alig.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/three_prime_utr.bed.saf \
-F SAF -o counts_filterd_three_prime_utr_CPTAC_chRCCpatients_v3.txt C3L-01469-CPT0389960009-alig.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/transcript.bed.saf \
-F SAF -o counts_filterd_transcript_CPTAC_chRCCpatients_v3.txt C3L-01469-CPT0389960009-alig.bam

featureCounts -p -a /mctp/share/users/gondal/cptac3-rarercc/rna-seq-v5-GRCh38/01_Shell_Script/introns.bed.saf \
-F SAF -o counts_filterd_transcript_TCGA_chRCCpatients_v3.txt C3L-01469-CPT0389960009-alig.bam
