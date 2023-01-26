from pathlib import Path
import random

print(""" USAGE"
1. Open screen/tmux so that the job will continue even after closing: e.g. tmux
2. Activate CONDA environment with Snakemake: e.g. conda activate /data/umcg-sandreusanchez/PhyloAdaptationPipeline/CONDA_ENVS/Snake
3. Run Snakemke with a profile for sending HPC jobs / (might require also tracking, check: https://github.com/jdblischak/smk-simple-slurm ). e.g. snakemake --profile Snake_profile/ --snakefile GUTS_MGS_pipeline.smk --cluster-status status_Script.sh
""")



Dic_files = {}
#I = "V2_3G064_FKDN220284340-1A_HJT23DSX3_L1" #V5G071_FKDN220284427-1A_HKGMMDSX3_L2
with open("FILES/Files_list.txt") as I:
	for line in I:
		l = line.rstrip().split()
		#if l[0] != I: continue
		Dic_files[l[0]] = l[1:]

#print(Dic_files)

def search(ID, direction):
	Found = Dic_files[ID]
	if direction == "forward": return(Found[1])
	else: return(Found[2])


rule all:
	input:
		"Output/Merged_taxonomy.tsv",
		"Output/All_kneaddata_reads_stats_F.csv",
		"Output/Pathway/GBM_quant/"),
		"Output/Subspecies/metaSNV/")
	message: "Pipeline complete"




rule KO_to_GBM:
	input: "Output/Pathway/Merged_genefamilies_KO.tsv"
	output: directory("Output/Pathway/GBM_quant")
	params:
		GBM = "GBMs/GBM.inputfile.txt"
	shell:
		"""
		java -jar omixer-rpm-1.1.jar -a 1  -d {parms.GBM} -o {output} -i {input}
		"""
rule MergePathways:
        input: expand("Output/Pathway/{Sample}_kneaddata/{Sample}_kneaddata_cleaned_paired_merged_genefamilies.tsv", zip,Sample=Dic_files.keys() )
        output: "Output/Pathway/Merged_genefamilies.tsv"
        shell:
                """
                ml Miniconda3/4.7.10
                set +u; /data/umcg-tifn/rgacesa/conda_dag3_mp3; set -u
                humann_join_tables --input Output/Pathway/ --output Output/Pathway/Merged_genefamilies.tsv
		"""

rule Pathway_to_KO:
	input: "Output/Pathway/{Sample}_kneaddata/{Sample}_kneaddata_cleaned_paired_merged_genefamilies.tsv"
	output:	temp("Output/Pathway/{Sample}_kneaddata_cleaned_paired_merged_genefamiliesKO.tsv")
	shell:
		"""
		ml Miniconda3/4.7.10
		set +u; /data/umcg-tifn/rgacesa/conda_dag3_mp3; set -u
		humann_regroup_table --input {input} --groups uniref90_ko --output {output}
		"""

rule QuantifyPathways:
	input: "Reads_clean/{Sample}_kneaddata_cleaned_paired_merged.fastq" 
	output: temp("Output/Pathway/{Sample}_kneaddata/{Sample}_kneaddata_cleaned_paired_merged_genefamilies.tsv")
	shell:
		"""
		module load Miniconda3/4.7.10
		set +u; source activate /data/umcg-tifn/rgacesa/conda_dag3_mp3 ; set -u
		humann --input {input} --output Output/Pathway/{wildcards.Sample}
		"""
rule metaSNV:
	input: expand("Output/Subspecies/Alignment/{Sample}.bam",zip,Sample=Dic_files.keys()),
	output: directory("Output/Subspecies/metaSNV")
	shell:
		"""
		ml Anaconda3/2020.11 SAMtools/1.10-GCC-9.3.0 BWA/0.7.17-GCC-10.2.0
		ENV=/scratch/umcg-sandreusanchez/BIOM_project/environment
		set +u; source activate $ENV  ; set -u
		motus snv_call -d Output/Subspecies/Alignments/ -o {output}
		"""

rule Align_for_calling:
	input:
		Forward = "Reads_clean/{Sample}_kneaddata_cleaned_pair_1.fastq",
		Reverse = "Reads_clean/{Sample}_kneaddata_cleaned_pair_2.fastq"
	output:
		temp("Output/Subspecies/Alignment/{Sample}.bam")
	shell:
		"""
		ml Anaconda3/2020.11 SAMtools/1.10-GCC-9.3.0 BWA/0.7.17-GCC-10.2.0
		ENV=/scratch/umcg-sandreusanchez/BIOM_project/environment
		set +u; source activate $ENV  ; set -u
		motus map_snv -f {input.Forward} -r {input.Reverse} -t 2 > {output}
		"""

rule merge_motus:
	input: expand('Output/taxonomy/taxonomy_profile.{Sample}.txt',zip,Sample=Dic_files.keys()),
	output: "Output/Merged_taxonomy.tsv"
	shell:
		"""
		ml Anaconda3/2020.11 SAMtools/1.10-GCC-9.3.0 BWA/0.7.17-GCC-10.2.0
		ENV=/scratch/umcg-sandreusanchez/BIOM_project/environment
		set +u; source activate $ENV  ; set -u
		motus merge -d Output/taxonomy/ -o {output}
"""
rule motus:
	input:
		Forward = "Reads_clean/{Sample}_kneaddata_cleaned_pair_1.fastq",
		Reverse = "Reads_clean/{Sample}_kneaddata_cleaned_pair_2.fastq", 
	output:
		temp('Output/taxonomy/taxonomy_profile.{Sample}.txt')
	params:
		Specificity=4
	resources:
		mem = "40gb",
		time = "23:59:00",
		threads = 3
	shell:
		"""
		ml Anaconda3/2020.11
ml SAMtools/1.10-GCC-9.3.0
ml BWA/0.7.17-GCC-10.2.0

ENV=/scratch/umcg-sandreusanchez/BIOM_project/environment
set +u; source activate $ENV ; set -u
MoTus=/scratch/umcg-sandreusanchez/BIOM_project/environment/bin/motus

# -c : Report counts instead of relative abundance
# -t : thread number
# -l : minimum read match (default 75)
# -g : Sensitivity (0) vs specificity (6)
# -A : print all taxonomic levels
# -q: Give complete taxonomy
$MoTus profile -f {input.Forward} -r {input.Reverse} -t {resources.threads} -g {params.Specificity} -c -n {wildcards.Sample} -q -A  > {output}
"""

rule merge_stats:
	input:
		expand("Output/{Sample}_kneaddata_reads_stats.csv", zip, Sample=Dic_files.keys())
	output:
		Merged = "Output/All_kneaddata_reads_stats_F.csv"
	run:
		Add = ""
		First = True
		for i in Path("Output/").glob("*_kneaddata_reads_stats.csv"):
			with open(i) as F:
				for line in F:
					if First == True:
						First = False
						Add += line
					if "Sample," in line: continue
					Add += line

		with open(output.Merged, "w") as O:
  			O.write(Add)


rule kneaddata:
	input:
		Forward = lambda wildcards: search(wildcards.Sample, "forward"),
		Reverse = lambda wildcards: search(wildcards.Sample, "reverse"),
	output:
		Paired1 = temp("Reads_clean/{Sample}_kneaddata_cleaned_pair_1.fastq"),
		Paired2 = temp("Reads_clean/{Sample}_kneaddata_cleaned_pair_2.fastq"),
		Merged = temp("Reads_clean/{Sample}_kneaddata_cleaned_paired_merged.fastq"),
		Reads_stats = temp("Output/{Sample}_kneaddata_reads_stats.csv")
	resources:
		mem = "24gb",
                time = "48:00:00"
	shell:
		"""
module load Miniconda3/4.7.10
set +u; source activate /data/umcg-tifn/rgacesa/conda_dag3_mp3 ; set -u
kneaddata --input {input.Forward} --input {input.Reverse} --threads 1 --processes 1 --output Reads_clean/{wildcards.Sample}/ --log logs/{wildcards.Sample}_kneaddata.log -db /data/umcg-tifn/rgacesa/dag3_pipeline_v3_dbs/hg37dec_v0.1  --trimmomatic /data/umcg-tifn/rgacesa/conda_dag3_mp3/share/trimmomatic-0.39-2/ --fastqc fastqc --sequencer-source none --trimmomatic-options "LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50" --bypass-trf
mv Reads_clean/{wildcards.Sample}/*_kneaddata_paired_1.fastq {output.Paired1}
mv Reads_clean/{wildcards.Sample}/*_kneaddata_paired_2.fastq {output.Paired2}
cat {output.Paired1} {output.Paired2} > {output.Merged}
rm -r Reads_clean/{wildcards.Sample}/
# -- PARSE KNEADDATA RESULTS --
python /data/umcg-tifn/rgacesa/dag3_pipeline_3c/utils/knead3parser.py --infile logs/{wildcards.Sample}_kneaddata.log --outfile {output.Reads_stats}
"""






 		









