import subprocess

###################################################################################################
## Merge tools
###################################################################################################
rule decompress_for_merging:
	input:
		vcf = OUT_FOLDER + "/sv_discovery/{tool}/{sample}/{sample}.{tool}.vcf.gz"
	output:
		vcf = OUT_FOLDER + "/sv_discovery/merged/{sample}/{sample}.{tool}.vcf"
	shell:	
		"zcat {input.vcf} > {output.vcf}; "

rule merge_tools:
	input:
		vcf = expand(OUT_FOLDER + "/sv_discovery/merged/{{sample}}/{{sample}}.{tool}.vcf", tool=TOOLS)
	output:
		vcf = OUT_FOLDER + "/sv_discovery/merged/{sample}/{sample}.vcf"
	conda:
		SNAKEDIR + "envs/survivor.yaml"
	params:
		vcf_list = OUT_FOLDER + "/sv_discovery/merged/{sample}/vcf.list",
		breakpoint_dist = "100",
		min_num_calls = "1",
		use_type = "1",
		use_strand = "1",
		dist_based = "0",
		min_sv_size = "50"
	shell:
		"ls {input.vcf} > {params.vcf_list}; "
		"SURVIVOR merge {params.vcf_list} \
			{params.breakpoint_dist} \
			{params.min_num_calls} \
			{params.use_type} \
			{params.use_strand} \
			{params.dist_based} \
			{params.min_sv_size} \
			{output.vcf}; "

# rule compress_merge_tools:
# 	input:
# 		vcf = OUT_FOLDER + "/sv_discovery/merged/{sample}/{sample}.vcf"
# 	output:
# 		vcfgz = OUT_FOLDER + "/sv_discovery/merged/{sample}/{sample}.vcf.gz",
# 		tbi = OUT_FOLDER + "/sv_discovery/merged/{sample}/{sample}.vcf.gz.tbi"
# 	conda:
# 		SNAKEDIR + "envs/bcftools.yaml"
# 	shell:
# 		"bcftools sort -Oz -o {output.vcfgz} {input.vcf}; "
# 		"tabix -p vcf {output.vcfgz}; "

###################################################################################################
## Merge samples for Joint-genotyping
###################################################################################################
rule merge_samples_01:
	input:
		vcf = expand(OUT_FOLDER + "/sv_discovery/merged/{sample}/{sample}.vcf", sample=participant_id)
	output:
		vcf_list = OUT_FOLDER + "/merged_cohort/vcf.list"
	run:
		shell("ls {input.vcf} > {output.vcf_list}; ")
		if( "SV_PANEL" in config ):
				for sv_ref in config["SV_PANEL"]:
					shell("ls {sv_ref} >> {output.vcf_list}; ")

# Options of using SURVIVOR or Jasmine. Jasmine by default.
useJAMINE=False
if (useJAMINE):
	checkpoint merge_samples_02:
		input:
			vcf_list = OUT_FOLDER + "/merged_cohort/vcf.list"
		output:
			vcf = OUT_FOLDER + "/merged_cohort/raw_merged.vcf"
		conda:
			SNAKEDIR + "envs/jasmine.yaml"
		params:
			vcf_list = OUT_FOLDER + "/merged_cohort/vcf.list",
		shell:
			"jasmine file_list={params.vcf_list} out_file={output.vcf}; " 
			"jasmine --dup_to_ins --postprocess_only out_file={output.vcf}; "
else:
	checkpoint merge_samples_02:
		input:
			vcf_list = OUT_FOLDER + "/merged_cohort/vcf.list"
		output:
			vcf = OUT_FOLDER + "/merged_cohort/raw_merged.vcf"
		conda:
			SNAKEDIR + "envs/survivor.yaml"
		params:
			vcf_list = OUT_FOLDER + "/merged_cohort/vcf.list",
			breakpoint_dist = "100",
			min_num_calls = "1",
			use_type = "1",
			use_strand = "1",
			dist_based = "0",
			min_sv_size = "50"
		shell:
			"SURVIVOR merge {params.vcf_list} \
				{params.breakpoint_dist} \
				{params.min_num_calls} \
				{params.use_type} \
				{params.use_strand} \
				{params.dist_based} \
				{params.min_sv_size} \
				{output.vcf}; "

checkpoint compress_merge_samples_02:
	input:
		vcf = OUT_FOLDER + "/merged_cohort/raw_merged.vcf"
	output:
		vcfgz = OUT_FOLDER + "/merged_cohort/raw_merged.vcf.gz",
		tbi = OUT_FOLDER + "/merged_cohort/raw_merged.vcf.gz.tbi"
	conda:
		SNAKEDIR + "envs/bcftools.yaml"
	shell:
		"bcftools sort -Oz -o {output.vcfgz} {input.vcf}; "
		"tabix -p vcf {output.vcfgz}; "

###################################################################################################
## Joint-genotyping
###################################################################################################
def get_chr(VCF):
	try:
		CHROM = subprocess.run(["tabix","-l", VCF], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout.decode('utf-8').splitlines()
		zipbObj = zip(CHROM, CHROM)
		return(dict(zipbObj))
	except:
		return(range(1,22))

rule avg_cov_by_readlen:
	input:
	  vcf = OUT_FOLDER + "/merged_cohort/raw_merged.vcf.gz",
		bam = OUT_FOLDER + "/input/{sample}.bam"
	output:
		avg_cov_by_readlen = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.avg_cov_by_readlen" 
	conda:
		SNAKEDIR + "envs/graphtyper.yaml"
	shell:
		"samtools idxstats $(readlink {input.bam}) | head -n -1 | awk \'{{sum+=$3+$4; ref+=$2}} END{{print sum/ref}}\' > {output.avg_cov_by_readlen}; "

rule genotype_01:
	input:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		gvcf = OUT_FOLDER + "/merged_cohort/raw_merged.vcf.gz",
		avg_cov_by_readlen = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.avg_cov_by_readlen" 
	output:
		vcf = temp(OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.{chrom}.vcf")
	conda:
		SNAKEDIR + "envs/graphtyper.yaml"
	params:
		outdir = directory(OUT_FOLDER + "/sv_genotyping/{sample}"),
		chrom = lambda wildcards: get_chr(OUT_FOLDER + "/merged_cohort/raw_merged.vcf.gz")[wildcards.chrom]
	priority: 50
	shell:
		"graphtyper genotype_sv {REFERENCE_FASTA} {input.gvcf} \
			--force_no_filter_zero_qual \
			--avg_cov_by_readlen={input.avg_cov_by_readlen} \
			--log={params.outdir}/{params.chrom}/{wildcards.sample}.log \
			--sam={input.bam} \
			--region={params.chrom} \
			--max_files_open=1 \
			--threads=1 \
			--output={params.outdir} && "
		"graphtyper vcf_concatenate {params.outdir}/{params.chrom}/*.vcf.gz > {output.vcf}.tmp; "
		"rm -rf {params.outdir}/{params.chrom}; "
		"mv {output.vcf}.tmp {output.vcf}; "

rule compress_genotype_01:
	input:
		vcf = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.{chrom}.vcf"
	output:
		vcf = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.{chrom}.vcf.gz"
	conda:
		SNAKEDIR + "envs/bcftools.yaml"
	params:
		chrom = lambda wildcards: get_chr(OUT_FOLDER + "/merged_cohort/raw_merged.vcf.gz")[wildcards.chrom]
	shell:
		"bcftools sort -Oz -o {output.vcf} {input.vcf}; "
		"tabix -p vcf {output.vcf}; "

rule genotype_02:
	input:
		tbi = OUT_FOLDER + "/merged_cohort/raw_merged.vcf.gz.tbi",
		vcf = lambda wildcards: expand(OUT_FOLDER + "/sv_genotyping/{{sample}}/{{sample}}.{chrom}.vcf.gz", chrom=get_chr(OUT_FOLDER + "/merged_cohort/raw_merged.vcf.gz"))
	output:
		vcf = temp(OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.vcf")
	conda:
		SNAKEDIR + "envs/graphtyper.yaml"
	shell:
		"graphtyper vcf_concatenate {input.vcf} > {output.vcf}; "

rule compress_and_filter_genotype_02:
	input:
		vcf = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.vcf"
	output:
		vcfgz = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.vcf.gz",
		tbi = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.vcf.gz.tbi",
		filt_vcfgz = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}_filt.vcf.gz",
		filt_tbi = OUT_FOLDER + "/sv_genotyping/{sample}/{sample}_filt.vcf.gz.tbi"
	conda:
		SNAKEDIR + "envs/bcftools.yaml"
	shell:
		"bcftools view --include \"SVMODEL='AGGREGATED'\" {input.vcf} | bcftools sort -Oz -o {output.vcfgz}; "
		"tabix -p vcf {output.vcfgz}; "
		"bcftools view -e QUAL==0 -Oz -o {output.filt_vcfgz} {output.vcfgz}; "
		"tabix -p vcf {output.filt_vcfgz}; "

rule merge_genotyped_samples:
	input:
		vcf = expand(OUT_FOLDER + "/sv_genotyping/{sample}/{sample}.vcf.gz", sample=participant_id)
	output:
		vcf_full = OUT_FOLDER + "/merged_cohort/gt_merged_full.vcf"
	conda:
		SNAKEDIR + "envs/graphtyper.yaml"		
	shell:
		"graphtyper vcf_merge {input.vcf} --sv | grep -E -v \"VarType=XG|VarType=IG\" > {output.vcf_full}; "

rule compress_merge_genotyped_samples:
	input:
		vcf = OUT_FOLDER + "/merged_cohort/gt_merged_full.vcf"
	output:
		vcfgz_full = OUT_FOLDER + "/merged_cohort/gt_merged_full.vcf.gz",
		tbi_full = OUT_FOLDER + "/merged_cohort/gt_merged_full.vcf.gz.tbi",
		vcf = OUT_FOLDER + "/merged_cohort/gt_merged.vcf",
		vcfgz = OUT_FOLDER + "/merged_cohort/gt_merged.vcf.gz",
		tbi = OUT_FOLDER + "/merged_cohort/gt_merged.vcf.gz.tbi"
	conda:
		SNAKEDIR + "envs/bcftools.yaml"		
	shell:
		"bcftools sort -Oz -o {output.vcfgz_full} {input.vcf}; "
		"tabix -p vcf {output.vcfgz_full}; "
		"bcftools view -e QUAL==0 {output.vcfgz_full} > {output.vcf}; "
		"bgzip -c {output.vcf} > {output.vcfgz}; "
		"tabix -p vcf {output.vcfgz}; "

###################################################################################################
## For benchmarking genotyping individual tools
###################################################################################################
# Deactivated in the pipeline.
rule genotype_discovery:
	input:
		vcf = OUT_FOLDER + "/sv_discovery/{tool}/{sample}/{sample}.{tool}.vcf.gz"
	output:
		avg_cov_by_readlen = OUT_FOLDER + "/sv_discovery/{tool}/{sample}/{sample}.{tool}.avg_cov_by_readlen",
		vcf = OUT_FOLDER + "/sv_discovery/{tool}/{sample}/{sample}.{tool}.gt.vcf",
		vcfgz = OUT_FOLDER + "/sv_discovery/{tool}/{sample}/{sample}.{tool}.gt.vcf.gz",
		tbi = OUT_FOLDER + "/sv_discovery/{tool}/{sample}/{sample}.{tool}.gt.vcf.gz.tbi"
	conda:
		SNAKEDIR + "envs/graphtyper.yaml"
	params:
		bam = OUT_FOLDER + "/input/{sample}.bam",
		outdir = directory(OUT_FOLDER + "/sv_discovery/{tool}/{sample}/gt")
	shell:
		"samtools idxstats $(readlink {params.bam}) | head -n -1 | awk \'{{sum+=$3+$4; ref+=$2}} END{{print sum/ref}}\' > {output.avg_cov_by_readlen}; "
		"for chr in $(tabix -l {input.vcf}); do \
			graphtyper genotype_sv {REFERENCE_FASTA} {input.vcf} \
			--force_no_filter_zero_qual \
			--avg_cov_by_readlen={output.avg_cov_by_readlen} \
			--log={params.outdir}/$chr/{wildcards.sample}.log \
			--sam={params.bam} \
			--region=$chr \
			--max_files_open=1 \
			--threads=1 \
			--output={params.outdir}; \
		done; "
		"graphtyper vcf_concatenate --no_sort {params.outdir}/*/*.vcf.gz | \
			bcftools view --include \"SVMODEL='AGGREGATED'\" | bcftools view -e QUAL==0 | bcftools sort > {output.vcf};" 
		"bgzip -c {output.vcf} > {output.vcfgz}; "
		"tabix -p vcf {output.vcfgz}; "
