
####################
# RUN VERIFYBAMID2 #
####################

for INPUT in /sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/input/*.bam
do
   OUT=`basename $INPUT .bam`
   submitjob 12 -c 12 -m 50 -q express -P acc_PVI \
   verifyBamID \
   --vcf /sc/arion/work/pintod02/conda/envs/snakesv_env/opt/snakeSV/resources/somalier_sites/sites.hg38.vcf.gz \
   --bam  ${INPUT}\
   --out /sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/verifyBamID/${OUT}.verifyBamID.check \
   --verbose \
   --ignoreRG
done


for INPUT in /sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/input/*.bam
do
   for REF in hgdp.10k.b38 1000g.phase3.10k.b38 hgdp.100k.b38 1000g.phase3.100k.b38
   do
      OUT=`basename $INPUT .bam`
      submitjob 12 -c 12 -m 50 -q express -P acc_PVI \
      verifybamid2 \
      --SVDPrefix /sc/arion/work/pintod02/conda/envs/verifybamid2/share/verifybamid2-2.0.1-12/resource/${REF}.vcf.gz.dat \
      --BamFile   ${INPUT}\
      --Output    /sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/verifyBamID2/${OUT}.${REF}.withinAncestry.check \
      --Reference /sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/ref-bkp/GRCh38.primary_assembly.genome.fa \
      --WithinAncestry
   done
done

for INPUT in /sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/input/*.bam
do
   for REF in hgdp.10k.b38 1000g.phase3.10k.b38 hgdp.100k.b38 1000g.phase3.100k.b38
   do
      OUT=`basename $INPUT .bam`
      submitjob 12 -c 12 -m 50 -q express -P acc_PVI \
      verifybamid2 \
      --SVDPrefix /sc/arion/work/pintod02/conda/envs/verifybamid2/share/verifybamid2-2.0.1-12/resource/${REF}.vcf.gz.dat \
      --BamFile   ${INPUT}\
      --Output    /sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/verifyBamID2/${OUT}.${REF}.betweenAncestry.check \
      --Reference /sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/ref-bkp/GRCh38.primary_assembly.genome.fa
   done
done



grep '' /sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/verifyBamID2/*.hgdp.100k.b38.withinAncestry.check.selfSM > ~/hgdp.100k.b38.withinAncestry.selfSM.txt
grep '' /sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/verifyBamID2/*.1000g.phase3.100k.b38.withinAncestry.check.selfSM > ~/1000g.phase3.100k.b38.withinAncestry.selfSM.txt

grep '' /sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/verifyBamID2/*.hgdp.100k.b38.betweenAncestry.check.selfSM > ~/hgdp.100k.b38.betweenAncestry.selfSM.txt
grep '' /sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/verifyBamID2/*.1000g.phase3.100k.b38.betweenAncestry.check.selfSM > ~/1000g.phase3.100k.b38.betweenAncestry.selfSM.txt


############################################
# MAKE VERIFYBAMID2 PCA PLOTS FOR ANCESTRY #
############################################

# Transpose PCA data for plots
cd /sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/verifyBamID2
for i in *.Ancestry; do transpose $i > $i.transposed; done

# Aggregate PCA data for plots
for SET in hgdp.10k.b38.betweenAncestry 1000g.phase3.10k.b38.betweenAncestry hgdp.10k.b38.withinAncestry 1000g.phase3.10k.b38.withinAncestry
do
   grep -iH 'IntendedSample' *.${SET}.check.Ancestry.transposed  | perl -pe 's/:/\t/' | cut -f 1,3,4 | perl -pe 's/\.hgdp.*transposed//; s/\.1000g.*transposed//' \
      > /sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/verifyBamID2_aggregated/${SET}.ancestry-PCs.txt
done

# Create hgdp plots
cd /sc/arion/work/pintod02/opt/VerifyBamID/bin
module load R
FOLDER=/sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/verifyBamID2_aggregated
for INPUT in ${FOLDER}/*hgdp.10k*.ancestry-PCs.txt
do
   PREFIX=`basename ${INPUT} .txt`
   sh ./run.plot.sh \
      -r hgdp \
      -i ${FOLDER}/${PREFIX}.txt \
      -o ${FOLDER}/${PREFIX}.plot 
done

# Create 1000g plots
cd /sc/arion/work/pintod02/opt/VerifyBamID/bin
module load R
FOLDER=/sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/results/verifyBamID2_aggregated
for INPUT in ${FOLDER}/*1000g.phase3.10k*.ancestry-PCs.txt
do
   PREFIX=`basename ${INPUT} .txt`
   sh ./run.plot.sh \
      -r 1000g \
      -i ${FOLDER}/${PREFIX}.txt \
      -o ${FOLDER}/${PREFIX}.plot 
done


############################
# RUN ECHTVAR ON VCF files #
############################

IPATH=/sc/arion/projects/pintod02c/Xiao/EPIASD_WGS/vcf
OPATH=/sc/arion/projects/pintod02c/Xiao/EPIASD_WGS/vcf-echtvar
GENREF=/sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/ref-bkp/GRCh38.primary_assembly.genome.fa
ECHTREF=/sc/arion/projects/pintod02c/snakesv_test/2024-12-22_results-759-samples/gnomad.v3.1.2.echtvar.v2.zip

# Normalize the vcf files first with bcftools and filter for PASS biallelic snps
for i in  ${IPATH}/*.vcf.gz
do
   NAME=`basename $i .vcf.gz`
   submitjob 12 -c 4 -m 15 -q express -P acc_PVI \
   bcftools view ${IPATH}/${NAME}.vcf.gz -f PASS,. -m2 -M2 -v snps -O b \
      \| bcftools norm -m -both -w 10000 -f ${GENREF} -O b \
      \| ./echtvar anno - ${OPATH}/${NAME}.vcf -e ${REF}
done

# Run echtvar with all gnomad 'AF' types
for ALLELEFREQ in gnomad_controls_and_biobanks_af gnomad_af gnomad_popmax_af
do
   for i in  ${OPATH}/*.vcf
   do
      NAME=`basename $i .vcf`
      submitjob 12 -c 2 -m 15 -q express -P acc_PVI \
      ./sceVCF --af_field=${ALLELEFREQ} -o results/sceVCF/${NAME}.${ALLELEFREQ}.txt ${OPATH}/${NAME}.vcf
   done
done

# Aggregate the results
for ALLELEFREQ in gnomad_controls_and_biobanks_af gnomad_af gnomad_popmax_af
do
   grep '' results/sceVCF/*.${ALLELEFREQ}.txt > results/sceVCF_aggregated/sceVCF-aggregated-results_${ALLELEFREQ}.txt
done


