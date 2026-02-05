# Genome-Assembly-Using-Canu
 This repository contains code and workflows for performing de novo genome assembly using Canu, a long-read assembler optimized for PacBio and Oxford Nanopore sequencing data.

# Key Features:

Long-Read Error Correction: Corrects sequencing errors inherent to long-read technologies.

Read Trimming: Removes low-quality regions and adapters prior to assembly.

De Novo Genome Assembly: Constructs contiguous genome assemblies using Canu’s adaptive k-mer and overlap-based approach.

Scalable Workflow: Supports large genomes and can be executed on local machines or HPC clusters.

Assembly Evaluation: Compatible with downstream tools for assessing assembly quality (e.g., contiguity and completeness).

# Requirements:

Canu

Java

Long-read FASTQ/FASTA files (PacBio or Oxford Nanopore)

# Code

(base) [suresh@node3 canu_assembly]$ nohup /opt/canu-2.0/Linux-amd64/bin/canu -p asm -d . genomeSize=1.35g -pacbio-hifi cattle.ccs.fastq.gz correctedErrorRate=0.015 corMinCoverage=3 minOverLapLength=1000 useGrid=false maxThreads=192 &


mito assembly:
(base) [suresh@node3 mito_genome_assembly]$ minimap2 -ax map-pb -t 40 ref_mito_seq.fasta ../Chilika_PFY2021N04Denovo1055_cattle.ccs.fastq.gz -o mito_reads.sam
(base) [suresh@node3 mito_genome_assembly]$ samtools view -bS mito_reads.sam > mito_reads.bam
(base) [suresh@node3 mito_genome_assembly]$ samtools view -h -F 4 mito_reads.bam > mito_reads_mapped.bam
(base) [suresh@node3 mito_genome_assembly]$ samtools fastq mito_reads_mapped.bam > mito_reads_mapped.fastq

(base) [suresh@node3 mito_genome_assembly]$ nohup /opt/canu-2.0/Linux-amd64/bin/canu -p asm -d . genomeSize=16.341m -pacbio-hifi mito_reads_mapped.fastq corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 stopOnLowCoverage=0.1 minInputCoverage=0.1 useGrid=false maxThreads=40 &

(base) [suresh@node3 mito_genome_assembly]$ nohup blastn -db /Analysis/Cattle_PFY2021N04Denovo1054/mito_genome_assembly/ref_mito_seq.fasta -query asm.contigs.fasta -outfmt 6 &

(base) [suresh@node3 mito_genome_assembly]$ samtools faidx asm.contigs.fasta tig00000001:16423-32780 >tig00000001:16423-32780_mito.fasta

#Again Mito reads from Buffalo (>NC_006295.1 Bubalus carabanensis mitochondrion, complete genome)
(base) [suresh@node3 mito_assembly_buffalo_ref]$ minimap2 -ax map-pb -t 40 Mito_Buffalo_ref.fasta ../../Chilika_PFY2021N04Denovo1055_cattle.ccs.fastq.gz -o mito_reads.sam

#Busco
(base) [suresh@node4 assembly-results]$ nohup python /apps_940/busco/scripts/run_BUSCO.py -i final.p_ctg.fasta -o busco_analysis_mammalian -l /apps_940/busco/lineage/mammalia_odb9/ -m genome -c 192 -z &

BUSCO version is: 3.1.0
The lineage dataset is: mammalia_odb9 (Creation date: 2016-02-13, number of species: 50, number of BUSCOs: 4104)
To reproduce this run: python /apps_940/busco/scripts/run_BUSCO.py -i final.p_ctg.fasta -o busco_analysis_mammalian -l /apps_940/busco/lineage/mammalia_odb9/ -m genome -c 192 -z -sp human

Summarized benchmarking in BUSCO notation for file final.p_ctg.fasta
BUSCO was run in mode: genome

        C:93.0%[S:91.9%,D:1.1%],F:3.0%,M:4.0%,n:4104

        3820    Complete BUSCOs (C)
        3773    Complete and single-copy BUSCOs (S)
        47      Complete and duplicated BUSCOs (D)
        123     Fragmented BUSCOs (F)
        161     Missing BUSCOs (M)
        4104    Total BUSCO groups searched


Chilika stats
stats for final.p_ctg.fasta
sum = 2647072203, n = 452, ave = 5856354.43, largest = 57802038
N50 = 15889895, n = 46
N60 = 12853537, n = 64
N70 = 10099134, n = 88
N80 = 7006491, n = 119
N90 = 4583560, n = 163
N100 = 27843, n = 452
N_count = 0
Gaps = 0
(base) [suresh@node4 assembly-results]$ python Assembly_Falcon/2-asm-falcon/pb-assembly/scripts/get_asm_stats.py final.p_ctg.fasta 
{
 "asm_contigs": 452,
 "asm_esize": 21085941,
 "asm_max": 57802038,
 "asm_mean": 5856354,
 "asm_median": 1610929,
 "asm_min": 27843,
 "asm_n50": 15889895,
 "asm_n90": 4583560,
 "asm_n95": 2509104,
 "asm_total_bp": 2647072203
}

 
Meryl
(/Analysis/tools/merqury) [suresh@node3 Merqury]$ meryl count k=21 cattle.ccs.fastq.gz output O


mammalia_odb9

Optical Map-hybrid scaffolding
(base) [suresh@node4 hybrid_scaffolds]$ perl5.16.3 /Analysis/vgp-pipeline/tools/pipeline/Solve3.6.1_11162020/HybridScaffold/1.0/hybridScaffold.pl -c /Analysis/vgp-pipeline/tools/pipeline/1.0/HybridScaffold/1.0/hybridScaffold_DLE1_config.xml -b /Analysis/Chilika_PFY2021N04Denovo1055_cattle/OM_data/PFY2021N04Denovo1055_OM/PFY2021N04Denovo1055_OM_-_De_novo_pipeline_results/output/contigs/exp_refineFinal1/EXP_REFINEFINAL1.cmap -n ../IPA/RUN/assembly-results/final.p_ctg.fasta -u CTTAAG -z out.zip -w status_out.txt -B 2 -N 2 -f -g -r /Analysis/vgp-pipeline/tools/pipeline/Solve3.6.1_11162020/RefAligner/1.0/RefAligner -p /Analysis/vgp-pipeline/tools/pipeline/Solve3.6.1_11162020/Pipeline/1.0/ -o output

(base) [suresh@node4 hybrid_scaffolds]$ assembly-stats EXP_REFINEFINAL1_bppAdjust_cmap_final_p_ctg_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta
stats for EXP_REFINEFINAL1_bppAdjust_cmap_final_p_ctg_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta
sum = 2645144328, n = 61, ave = 43363021.77, largest = 156782546
N50 = 79336693, n = 13
N60 = 71730516, n = 16
N70 = 64536426, n = 20
N80 = 52387028, n = 25
N90 = 42534419, n = 30
N100 = 127000, n = 61
N_count = 10647162
Gaps = 1826
"asm_contigs": 61,
 "asm_esize": 83962810,
 "asm_max": 156782546,
 "asm_mean": 43363021,
 "asm_median": 34346819,
 "asm_min": 127000,
 "asm_n50": 79336693,
 "asm_n90": 42534419,
 "asm_n95": 23907480,
 "asm_total_bp": 2645144328


pbmm2 align cns.ctg.fasta /data/analysis/dr.shekhar/analyses/bams/m54333_181117_013313.subreads.bam m54333_181117_013313.subreads.sorted.bam --sort -j 40 -J 40


