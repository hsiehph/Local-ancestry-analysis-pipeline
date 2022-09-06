import os, re
import pandas as pd

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

shell.executable("/bin/bash")
shell.prefix("source %s/env.cfg; set -eo pipefail; " % SNAKEMAKE_DIR)

if not os.path.exists("log"):
	os.makedirs("log")

if config == {}:
	configfile: "%s/config.yaml" % SNAKEMAKE_DIR

BED_removed = config["BED_removed"]
ListSample = config["sampleFile4RFMix"].keys()
sampleID = config["phasedSNV_BEDs"].keys()
list_chroms = eval(config["list_chroms"])
cutoff_genLength = config["cutoff_genLength"]
print(cutoff_genLength)


def _get_sampleBEDfile(wildcards):
	return config["phasedSNV_BEDs"][wildcards.sampleID]

def _get_geneticMapfile(wildcards):
	return config["geneticMaps"][wildcards.chrID]

def _get_refVCFfile(wildcards):
	return config["refPanel_VCFs"][wildcards.chrID]


rule all:
	input:	expand("plots_PH/{sampleID}.{ListSample}.msp.tsv.LAIplots.PH.pdf", sampleID=sampleID, ListSample=ListSample),
			expand("plots_DP/{sampleID}.{ListSample}.msp.tsv.LAIplots.DP.pdf", sampleID=sampleID, ListSample=ListSample),
			expand("plots_PH_denoise/{sampleID}.{ListSample}.{cutoff_genLength}.msp.tsv.LAIplots.PH.pdf", sampleID=sampleID, ListSample=ListSample, cutoff_genLength=cutoff_genLength),
			expand("plots_DP_denoise/{sampleID}.{ListSample}.{cutoff_genLength}.msp.tsv.LAIplots.DP.pdf", sampleID=sampleID, ListSample=ListSample, cutoff_genLength=cutoff_genLength)
	wildcard_constraints:
			ListSample="|".join(config["sampleFile4RFMix"].keys())
			


rule plotHap_PH_denoise:
	input: expand("{{sampleID}}/{{ListSample}}_denoise/{{sampleID}}.{chrID}.{{ListSample}}.{{cutoff_genLength}}.msp.tsv", chrID=list_chroms)
	output: "plots_PH_denoise/{sampleID}.{ListSample}.{cutoff_genLength}.msp.tsv.LAIplots.PH.pdf"
	wildcard_constraints:
			ListSample="|".join(config["sampleFile4RFMix"].keys())
	params: sge_opts="-l mfree=5G -l h_rt=2:00:00"
	shell:
		" x=`basename {output} .PH.pdf`; Rscript scripts/karyoplotR.PH.r {wildcards.sampleID}/{wildcards.ListSample}_denoise/ plots_PH_denoise/${{x}} {wildcards.cutoff_genLength} ; touch {output} "


rule plotHap_DP_denoise:
	input: expand("{{sampleID}}/{{ListSample}}_denoise/{{sampleID}}.{chrID}.{{ListSample}}.{{cutoff_genLength}}.msp.tsv", chrID=list_chroms)
	output: "plots_DP_denoise/{sampleID}.{ListSample}.{cutoff_genLength}.msp.tsv.LAIplots.DP.pdf"
	wildcard_constraints:
			ListSample="|".join(config["sampleFile4RFMix"].keys())
	params: sge_opts="-l mfree=5G -l h_rt=2:00:00"
	shell:
		" Rscript scripts/karyoplotR.DP.r {wildcards.sampleID}/{wildcards.ListSample}_denoise/ {output} "


rule denoiseMSPfile:
	input: "{sampleID}/{ListSample}/{sampleID}.{chrID}.{ListSample}.msp.tsv"
	output: "{sampleID}/{ListSample}_denoise/{sampleID}.{chrID}.{ListSample}.{cutoff_genLength}.msp.tsv"
	wildcard_constraints:
			ListSample="|".join(config["sampleFile4RFMix"].keys())
	params: sge_opts="-l mfree=4G -l h_rt=12:00:00"
	shell:
		" python scripts/denoise_rfmixMSPoutput.py {input} {wildcards.cutoff_genLength} > {output} "


rule plotHap_PH:
	input: expand("{{sampleID}}/{{ListSample}}/{{sampleID}}.{chrID}.{{ListSample}}.msp.tsv", chrID=list_chroms)
	output: "plots_PH/{sampleID}.{ListSample}.msp.tsv.LAIplots.PH.pdf"
	wildcard_constraints:
			ListSample="|".join(config["sampleFile4RFMix"].keys())
	params: sge_opts="-l mfree=5G -l h_rt=2:00:00"
	shell: 
		" x=`basename {output} .PH.pdf`; Rscript scripts/karyoplotR.PH.r {wildcards.sampleID}/{wildcards.ListSample}/ plots_PH/${{x}} ; touch {output} "


rule plotHap_DP:
	input: expand("{{sampleID}}/{{ListSample}}/{{sampleID}}.{chrID}.{{ListSample}}.msp.tsv", chrID= list_chroms)
	output: "plots_DP/{sampleID}.{ListSample}.msp.tsv.LAIplots.DP.pdf"
	wildcard_constraints:
			ListSample="|".join(config["sampleFile4RFMix"].keys())
	params: sge_opts="-l mfree=5G -l h_rt=2:00:00"
	shell:
		" Rscript scripts/karyoplotR.DP.r {wildcards.sampleID}/{wildcards.ListSample}/ {output} "


rule RFMix:
	input: vcf="{sampleID}/VCFbyChr/{sampleID}.{chrID}.rmdup.masked.vcf.gz",
			geneMap= lambda wildcards: config["geneticMaps"][wildcards.chrID],
			refVCF= lambda wildcards: config["refPanel_VCFs"][wildcards.chrID],
			ListSample= lambda wildcards: config["sampleFile4RFMix"][wildcards.ListSample]
	output: "{sampleID}/{ListSample}/{sampleID}.{chrID}.{ListSample}.msp.tsv"
	wildcard_constraints:
			ListSample="|".join(config["sampleFile4RFMix"].keys())
	params: sge_opts="-l mfree=6G -l h_rt=72:00:00 -pe serial 5 -l testing=true ", 
			prefix=lambda wildcards: wildcards.sampleID + "/" + wildcards.ListSample + "/" + wildcards.sampleID + "." + wildcards.chrID + "." + wildcards.ListSample
	shell:
		" rfmix -f {input.vcf} -r {input.refVCF} -g {input.geneMap} -m {input.ListSample} -o {params.prefix} --chromosome={wildcards.chrID} --n-threads=5 -G 15 -e 5 -w 0.4 -n 5 -c 50 -s 150 --rf-minimum-snps=100 --reanalyze-reference"


rule VCFbyChr:
	input: "{sampleID}/rawVCF/{sampleID}.raw.vcf.gz" 
	output: "{sampleID}/VCFbyChr/{sampleID}.{chrID}.rmdup.masked.vcf.gz" 
	params: sge_opts="-l mfree=4G -l h_rt=12:00:00 -pe serial 4 "
	priority: 10
	shell:
		" bcftools view -r {wildcards.chrID} {input} -Ov | bcftools norm -m + - -Ov | bcftools view -m2 -M2 -v snps | bedtools intersect -a - -b {BED_removed} -v -header | bgzip -c > {output}"


rule BED2VCF:
	input: inBED = _get_sampleBEDfile
	output: "{sampleID}/rawVCF/{sampleID}.raw.vcf.gz" 
	params: sge_opts="-l mfree=10G -l h_rt=12:00:00"
	priority: 10
	shell:
		" cat {input.inBED} | python ~/bin/scripts/scripts_uw/HGSVC_LAI/PeteBED2VCF.py {wildcards.sampleID} | bgzip -c > {output} ; tabix -p vcf {output} "


