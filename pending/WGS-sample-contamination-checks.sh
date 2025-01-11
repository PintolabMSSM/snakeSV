
################
# RUN SOMALIER #
################

# Input parameters
BAMINPUT=/sc/arion/projects/pintod02c/Xiao/EPIASD_WGS/bam/staging
REFGEN=/sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/ref-bkp/GRCh38.primary_assembly.genome.fa
REFVCF=/sc/arion/work/pintod02/conda/envs/snakesv_env/opt/snakeSV/workflow/../resources//somalier_sites/sites.hg38.vcf.gz
SOMALIERPED=pedfile_for_somalier_ready-to-use.csv
SOMALIERBIN=/sc/arion/work/pintod02/conda/envs/snakesv_env/opt/snakeSV/workflow/../resources//somalier
SOMALIEREXT=/sc/arion/projects/pintod02c/snakesv_test/WGS-QC-checks/results/somalier_extracted
SOMALIERAGG=/sc/arion/projects/pintod02c/snakesv_test/WGS-QC-checks/results/somalier

# Run somalier extract
for INPUT in ${BAMINPUT}/*.bam
do
   ${SOMALIERBIN} extract \
      -d results/somalier_extracted/ \
      -s ${REFVCF} \
      -f ${REFGEN} \
      ${INPUT}
done

# Aggregate ancestry data
${SOMALIERBIN} ancestry \
   --n-pcs=20 \
   -o ${SOMALIERAGG} \
   --labels /sc/arion/work/pintod02/conda/envs/snakesv_env/opt/snakeSV/resources/somalier_sites/ancestry-labels-1kg.tsv \
   /sc/arion/work/pintod02/conda/envs/snakesv_env/opt/snakeSV/resources/somalier_sites/1kg-somalier/*.somalier ++ ${SOMALIEREXT}/*.somalier;

# Aggregate relatedness data
${SOMALIERBIN} relate \
   --infer \
   --ped ${SOMALIERPED} \
   -o ${SOMALIERAGG}/somalier.somalier-relate \
   ${SOMALIEREXT}/*.somalier \
   > ${SOMALIERAGG}/somalier.somalier-relate.log 2>&1



####################
# RUN VERIFYBAMID2 #
####################

# Input parameters
BAMINPUT=/sc/arion/projects/pintod02c/Xiao/EPIASD_WGS/bam/staging
REFGEN=/sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/ref-bkp/GRCh38.primary_assembly.genome.fa
VBID2EXT=/sc/arion/projects/pintod02c/snakesv_test/WGS-QC-checks/results/verifyBamID2_extracted/
VBID2AGG=/sc/arion/projects/pintod02c/snakesv_test/WGS-QC-checks/results/verifyBamID2/

# Make sure the output dirs exist
mkdir -p ${VBID2EXT}
mkdir -p ${VBID2AGG}

# Activate verifybamid2 environment
conda activate verifybamid2

# Run within/between ancestry analysis with verifybamid2
for INPUT in ${BAMINPUT}/*.bam
do
   for SNPREF in hgdp.10k.b38 1000g.phase3.10k.b38 hgdp.100k.b38 1000g.phase3.100k.b38
   do
      OUT=`basename $INPUT .bam`
      # Within-ancestry analysis
      submitjob 12 -c 4 -m 12 -q express -P acc_PVI \
      verifybamid2 \
         --SVDPrefix /sc/arion/work/pintod02/conda/envs/verifybamid2/share/verifybamid2-2.0.1-12/resource/${SNPREF}.vcf.gz.dat \
         --BamFile   ${INPUT}\
         --Output    ${VBID2EXT}${OUT}.${SNPREF}.withinAncestry.check \
         --Reference ${REFGEN} \
         --WithinAncestry
      
      # Between-ancestry analysis
      submitjob 12 -c 4 -m 12 -q express -P acc_PVI \
      verifybamid2 \
         --SVDPrefix /sc/arion/work/pintod02/conda/envs/verifybamid2/share/verifybamid2-2.0.1-12/resource/${SNPREF}.vcf.gz.dat \
         --BamFile   ${INPUT}\
         --Output    ${VBID2EXT}${OUT}.${SNPREF}.betweenAncestry.check \
         --Reference ${REFGEN}
   done
done

# Aggregate results
grep '' ${VBID2EXT}*.hgdp.100k.b38.withinAncestry.check.selfSM          > ${VBID2AGG}/hgdp.100k.b38.withinAncestry.selfSM.txt
grep '' ${VBID2EXT}*.1000g.phase3.100k.b38.withinAncestry.check.selfSM  > ${VBID2AGG}/1000g.phase3.100k.b38.withinAncestry.selfSM.txt
grep '' ${VBID2EXT}*.hgdp.100k.b38.betweenAncestry.check.selfSM         > ${VBID2AGG}/hgdp.100k.b38.betweenAncestry.selfSM.txt
grep '' ${VBID2EXT}*.1000g.phase3.100k.b38.betweenAncestry.check.selfSM > ${VBID2AGG}/1000g.phase3.100k.b38.betweenAncestry.selfSM.txt


############################################
# MAKE VERIFYBAMID2 PCA PLOTS FOR ANCESTRY #
############################################

# Parameters
VBID2EXT=/sc/arion/projects/pintod02c/snakesv_test/WGS-QC-checks/results/verifyBamID2_extracted/
VBID2AGG=/sc/arion/projects/pintod02c/snakesv_test/WGS-QC-checks/results/verifyBamID2/
PLOTDIR=/sc/arion/work/pintod02/opt/VerifyBamID/bin

# Transpose PCA data for plots
cd ${VBID2EXT}
for i in *.Ancestry; do transpose $i > $i.transposed; done

# Aggregate PCA data for plots
cd ${VBID2EXT}
for SET in hgdp.10k.b38.betweenAncestry 1000g.phase3.10k.b38.betweenAncestry hgdp.10k.b38.withinAncestry 1000g.phase3.10k.b38.withinAncestry
do
   grep -iH 'IntendedSample' *.${SET}.check.Ancestry.transposed  | perl -pe 's/:/\t/' | cut -f 1,3,4 | perl -pe 's/\.hgdp.*transposed//; s/\.1000g.*transposed//' \
      > ${VBID2AGG}/${SET}.ancestry-PCs.txt
done

# Create hgdp plots
cd ${PLOTDIR}
module load R
for INPUT in ${VBID2AGG}/*hgdp.10k*.ancestry-PCs.txt
do
   PREFIX=`basename ${INPUT} .txt`
   sh ./run.plot.sh \
      -r hgdp \
      -i ${VBID2AGG}/${PREFIX}.txt \
      -o ${VBID2AGG}/${PREFIX}.plot 
done

# Create 1000g plots
cd ${PLOTDIR}
module load R
for INPUT in ${VBID2AGG}/*1000g.phase3.10k*.ancestry-PCs.txt
do
   PREFIX=`basename ${INPUT} .txt`
   sh ./run.plot.sh \
      -r 1000g \
      -i ${VBID2AGG}/${PREFIX}.txt \
      -o ${VBID2AGG}/${PREFIX}.plot 
done


############################
# RUN ECHTVAR ON VCF files #
############################

IPATH=/sc/arion/projects/pintod02c/Xiao/EPIASD_WGS/vcf/staging
OPATH=/sc/arion/projects/pintod02c/Xiao/EPIASD_WGS/vcf-echtvar/staging
RSPATH=/sc/arion/projects/pintod02c/snakesv_test/WGS-QC-checks/results/sceVCF_extracted/
AGPATH=/sc/arion/projects/pintod02c/snakesv_test/WGS-QC-checks/results/sceVCF
REFGEN=/sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/ref-bkp/GRCh38.primary_assembly.genome.fa
ECHTREF=/sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/gnomad.v3.1.2.echtvar.v2.zip
LIBDIR=/sc/arion/work/pintod02/opt/snakeSV/resources

# Make sure the output folders exist
mkdir -p ${RSPATH}
mkdir -p ${AGPATH}

# Load environment
module load bcftools

# Normalize the vcf files first with bcftools and filter for PASS biallelic snps
for i in  ${IPATH}/*.vcf.gz
do
   NAME=`basename $i .vcf.gz`
   submitjob 12 -c 4 -m 15 -q express -P acc_PVI \
   bcftools view ${IPATH}/${NAME}.vcf.gz -f PASS,. -m2 -M2 -v snps -O b \
      \| bcftools norm -m -both -w 10000 -f ${REFGEN} -O b \
      \| ${LIBDIR}/echtvar anno - ${OPATH}/${NAME}.vcf -e ${ECHTREF}
done

# Run echtvar with all gnomad 'AF' types
for ALLELEFREQ in gnomad_controls_and_biobanks_af gnomad_af gnomad_popmax_af
do
   for i in  ${OPATH}/*.vcf
   do
      NAME=`basename $i .vcf`
      submitjob 12 -c 2 -m 15 -q express -P acc_PVI \
      ${LIBDIR}/sceVCF --af_field=${ALLELEFREQ} -o ${RSPATH}/${NAME}.${ALLELEFREQ}.txt ${OPATH}/${NAME}.vcf
   done
done

# Aggregate the results
for ALLELEFREQ in gnomad_controls_and_biobanks_af gnomad_af gnomad_popmax_af
do
   grep '' ${RSPATH}/*.${ALLELEFREQ}.txt > ${AGPATH}/sceVCF-aggregated-results_${ALLELEFREQ}.txt
done


