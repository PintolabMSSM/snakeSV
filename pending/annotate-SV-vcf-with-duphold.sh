
####################
# RUN DUPHOLD      #
####################

# Download the static binary for duphold
wget https://github.com/brentp/duphold/releases/download/v0.2.3/duphold
chmod a+x duphold


# Run duphold. To annotate the SVs with SNVs we use the echtvar VCF files that have been normalized and are prefiltered for SNVs
# The echtvar VCF files are generated in the WGS-Sample-contamination-checks.sh script
DUPHOLD_BIN=/sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/duphold
REFERENCE_FASTA=/sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/ref-bkp/GRCh38.primary_assembly.genome.fa
ECHTVAR_VCF_FOLDER=/sc/arion/projects/pintod02c/Xiao/EPIASD_WGS/vcf-echtvar
BAM_FOLDER=/sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/input

for INPUTVCF in results/sv_discovery/*/*/*.vcf results/sv_discovery/*/*/*.vcf.gz
do
   # Format sample and output prefix
   SAMPLE_PREFIX=`basename ${INPUTVCF} .vcf`
   SAMPLE_PREFIX=`basename ${SAMPLE_PREFIX} .vcf.gz`
   SAMPLE_PREFIX=`basename ${SAMPLE_PREFIX} .delly`
   SAMPLE_PREFIX=`basename ${SAMPLE_PREFIX} .manta`
   SAMPLE_PREFIX=`basename ${SAMPLE_PREFIX} .smoove`
   OUTPUT_PREFIX=${INPUTVCF%.vcf.gz}
   OUTPUT_PREFIX=${OUTPUT_PREFIX%.vcf}
   
   # Submit duphold jobs
   submitjob 12 -c 12 -m 25 -q express -P acc_PVI \
   ${DUPHOLD_BIN} \
      --threads 8 \
      --snp     ${ECHTVAR_VCF_FOLDER}/${SAMPLE_PREFIX}.vcf \
      --vcf     ${INPUTVCF} \
      --bam     ${BAM_FOLDER}/${SAMPLE_PREFIX}.bam \
      --fasta   ${REFERENCE_FASTA} \
      --output  ${OUTPUT_PREFIX}.duphold.bcf
done

