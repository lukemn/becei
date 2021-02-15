
function doAln(){
    # make aln over junction for viewing 
    
    L=$1
    R=$2
    
    LH=$(head -1 $L | cut -c2-)
    RH=$(head -1 $R | cut -c2-)
    H=${LH}_${RH}
    HREF=$H.j.fa
    FLANKING=${H}_flank.fa
    SPR=$H.spanning.fa
    PAF=$H.paf	
    BAM=$H.bam
    
    echo ">${H}" > $HREF    
    grep -v ">" $FLANKING | tr -d "\n" >> $HREF
    samtools faidx $HREF

    # for samtools tview
    minimap2 -t $NP --cs -c -ax asm5 -f 0.00001 $HREF $SPR |\
        samtools sort - | samtools view -F4 -q20 -b - > $BAM
    samtools index $BAM
    samtools view $BAM | cut -f1 | sort | uniq > tmp
    grep -f tmp $PAF | cut -f1-12 | sort -k1,1 -k6,6 |\
        cut -f2,3 -d"/" | tr ":" "\t" | cut -f1-6,9,10-13 > $H.alns
    samtools tview --reference $HREF $BAM
    cat $H.alns	| column -t
}

function LGconsensus() {

    # map reads, call variants with deppvariant
    # apply homozygous PASS variants
    local H=$1
    local REF=$2
    
    BAM=$H.bam
    VCF=$H.dpv.vcf
    
    minimap2 -t $NP --cs -c -ax asm5 -f 0.00001 $REF $FQ |\
        samtools sort -m 4G -u -@ $((NP/2)) - |\
        samtools view -@ $((NP/2)) -F4 -q10 -b - > $BAM
    samtools index $BAM
    samtools faidx $REF
    samtools stats $BAM | grep "^SN" | cut -f2- > $BAM.stats    
    samtools view $BAM | cut -f5 | sump.py > $H.cons.log

    SIF=/scratch/work/public/singularity/deepvariant-1.1.0.sif
    singularity exec --bind `pwd` $SIF \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=PACBIO \
        --ref=$REF \
        --reads=$BAM \
        --output_vcf=$VCF  \
        --output_gvcf=$H.g.vcf.gz \
        --logging_dir=$H \
        --num_shards=$NP
    
    echo "Deepvariant" >> $H.cons.log    
    bgzip -f $VCF && tabix -p vcf $VCF.gz
    bcftools filter -i 'FILTER="PASS" && GT="AA"' -Ou $VCF.gz |\
        bcftools norm -f $REF - > $H.norm.vcf
    bgzip -f $H.norm.vcf && tabix -p vcf $H.norm.vcf.gz
    bcftools consensus -f $REF -o $H.cons.fa $H.norm.vcf.gz &>> $H.cons.log    

#     echo "Bcftools" >> $H.cons.log    
#     bcftools mpileup --threads $((NP/2)) -Ou -d 60 -f $REF $BAM |\
#         bcftools call --threads $((NP/2)) -mv -Ov > $H.bcf.vcf
#     bgzip -f $H.bcf.vcf && tabix -p vcf $H.bcf.vcf.gz
#     bcftools norm --threads $((NP/2)) -f $REF $H.bcf.vcf.gz > $H.bcf.norm.vcf
#     bgzip -f $H.bcf.norm.vcf && tabix -p vcf $H.bcf.norm.vcf.gz
#     bcftools consensus -f $REF $H.bcf.norm.vcf.gz -o $H.bcf.cons.fa \
#         -i 'GT=AA' &>> $H.cons.log    
        
}


function mapJunc(){
    # check for reads spanning two (+) single seq
    L=$1
    R=$2
    READS=$3
    
    if [ ! -s $L ]; then echo "no L: $L"; exit; fi
    if [ ! -s $R ]; then echo "no R: $R"; exit; fi
    if [ ! -s $READS ]; then echo "no READS: $READS"; exit; fi        

    FLANK=20000
    MINAL=1000
    MQ=0
    
    samtools faidx $L
    samtools faidx $R	
    
    LH=$(head -1 $L | cut -c2-)
    RH=$(head -1 $R | cut -c2-)
    LL=$(awk '{print $2}' $L.fai)
    LS=$((LL-FLANK))
    
    H=${LH}_${RH}
    FLANKING=${H}_flank.fa
    PAF=$H.paf
    BAM=$H.bam
    SPR=$H.spanning
    
    echo -e "$LH\t$LS\t$LL" | bedtools getfasta -fi $L -bed - > $FLANKING
    echo -e "$RH\t1\t$FLANK" | bedtools getfasta -fi $R -bed - >> $FLANKING
    minimap2 -t $NP --cs -cx asm5 $FLANKING $READS > $PAF
    
    awk -v M=$MINAL -v Q=$MQ '( $11 >= M && $12 >= Q )' $PAF |\
        cut -f1,6 | sort | uniq | cut -f1 | sort | uniq -c |\
        grep -w 2 | trst | cut -f3 > $SPR
    NSPAN=$(cat $SPR | wc -l)
    echo -e "\t$H $NSPAN candidate spanning reads"    
    if [[ $NSPAN -gt 0 ]]; then
        # get candidate reads
        fastaSubset.py $SPR $READS $SPR.fa &
    else
        rm -f $SPR ${H}*
    fi
}

# but only ~6x in 15kb reads? (LG6)
# READS=/scratch/lmn3/becei/hifi/Q20_15kb.fa
READS=/scratch/lmn3/becei/hifi/Q20_12kb.fa
# for consensus
FQ=/scratch/lmn3/becei/hifi/Q20_10kb.fq
NP=12
source ~/bin/quickmerge/.quickmergerc 

LG=1
# Hifiasm: 17(-) <1cM> 30(-) <0> 54(?) <0> 6(+)
# Flye m10 scaffold_8 spans 54-6
J=1
# flip 17 & 30
for i in 17 30; do ~/bin/seqtk/seqtk seq -r ptg0000${i}l.fa > ptg0000${i}l_rc.fa; done
# skip 17-30 junction
# map 30-54 junction. check both ori
for i in 54; do ~/bin/seqtk/seqtk seq -r ptg0000${i}l.fa > ptg0000${i}l_rc.fa; done
# none
mapJunc ptg000030l_rc.fa ptg000054l.fa $READS
# 3 candidate
mapJunc ptg000030l_rc.fa ptg000054l_rc.fa $READS
# no reads directly span junction, but 3 show gapped alignment
doAln ptg000030l_rc.fa ptg000054l_rc.fa
# samtools tview --reference ptg000030l_ptg000054l.j.fa ptg000030l_ptg000054l.bam
# cat ptg000030l_ptg000054l.alns
# 3 good reads strand consistent
125043886/ccs  17349  2863   14770  +  ptg000030l  8096   19997  11836  11907  60
125043886/ccs  17349  20     2845   +  ptg000054l  14814  17640  2794   2827   60
134547784/ccs  15260  11025  15239  +  ptg000030l  8096   12308  4165   4215   60
134547784/ccs  15260  5230   11007  +  ptg000054l  11868  17640  5716   5777   60
134547784/ccs  15260  919    6683   +  ptg000054l  11877  17553  435    5765   60
150276163/ccs  15342  23     6945   -  ptg000030l  8096   14999  6699   6923   60
150276163/ccs  15342  6963   12742  -  ptg000054l  11868  17640  5676   5779   60
# 23267483/ccs   16648  20     13541  -  ptg000030l  7      13514  13351  13521  46
# 23267483/ccs   16648  3395   15251  -  ptg000054l  7      11851  11722  11856  0
# 65863714/ccs   16923  19     7789   -  ptg000030l  7      7780   7740   7773   0
# 65863714/ccs   16923  19     9497   -  ptg000054l  7      9488   9448   9481   60
# merge? nope
# hifi 30=8592421
merge_wrapper.py -l 5000 -pre LG${LG}_J${J} ptg000030l_rc.fa ptg000054l_rc.fa
# chew back first 11.86 kb of ptg000054l_rc, then stitch and take consensus  
grep -v ">" ptg000030l_rc.fa | tr -d "\n"  > tmp
echo -e "ptg000054l\t11860\t489597" |\
    bedtools getfasta -fi ptg000054l_rc.fa -bed - | grep -v ">" >> tmp
echo ">ptg000030l_ptg000054l_stitch" > LG${LG}_J${J}.stitch.fa    
fold tmp >> LG${LG}_J${J}.stitch.fa
# iterate consensus
# variants: 43, 21, 9, 1
# reads mapped: 58646, 58890, 58835, 58845
# bases mapped (cigar): 532883171, 533190109, 533210405, 533266866
cat *log
LGconsensus LG${LG}_J${J}.stitch.1 LG${LG}_J${J}.stitch.fa
for i in {1..5}; do
    LGconsensus LG${LG}_J${J}.stitch.$((i+1)) LG${LG}_J${J}.stitch.${i}.cons.fa
done
    
J=2
# 0 candidates, as expected
mapJunc ptg000054l.fa ptg000006l.fa $READS
# 214 candidates
mapJunc ptg000054l_rc.fa ptg000006l.fa $READS
doAln ptg000054l_rc.fa ptg000006l.fa
# stitch with m10 scaffold_8
# hifi 54=489597, 6=4773269
# merge=5117367
merge_wrapper.py -l 5000 -pre LG${LG}_J${J} scaffold_8.fa hifiasm.fa
ln -s merged_LG${LG}_J${J}.fasta LG${LG}_J${J}.stitch.fa
LGconsensus LG${LG}_J${J}.stitch.1 LG${LG}_J${J}.stitch.fa
for i in {1..5}; do
    LGconsensus LG${LG}_J${J}.stitch.$((i+1)) LG${LG}_J${J}.stitch.${i}.cons.fa
done
# variants: 49, 10, 4, 1, 0, 0

# pseudochromosome
# ignore J1 for now
grep -v ">" ptg000017l_rc.fa > tmp
echo "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" >> tmp
grep -v ">" ptg000030l_rc.fa >> tmp
echo "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" >> tmp
grep -v ">" LG1_J2.stitch.6.cons.fa >> tmp
echo ">LG$LG" > lg$LG.pchrom.fa
cat tmp | tr -d "\n" | fold >> lg$LG.pchrom.fa



LG=2
# Hifiasm: 37(+) <0> 9(+)
# Flye m10 contig_68 spans these two
J=1
cat ptg0000* > hifiasm.fa
merge_wrapper.py -l 5000 -pre LG${LG}_J${J} hifiasm.fa contig_68.fa
# hifiasm sum 13988657, contig_68 = 13584573
# merged 13981387
H=LG${LG}_J${J}
ln -s merged_${H}.fasta ${H}.stitch.fa
LGconsensus ${H}.stitch.1 ${H}.stitch.fa
for i in {1..5}; do LGconsensus ${H}.stitch.$((i+1)) ${H}.stitch.${i}.cons.fa; done
# pseudochromosome
echo ">LG$LG" > lg$LG.pchrom.fa
grep -v ">" LG2_J1.stitch.6.cons.fa | tr -d "\n" | fold >> lg$LG.pchrom.fa


LG=3
# m10: 2(-) <0.53> 72(-) <0> scaffold_77(+) <0.53> 36(+)
J=1
# thread m10 with hifiasm HiC 3,4
merge_wrapper.py -l 5000 -pre LG${LG}_J${J} hic.fa m10.fa
# hifiasm hic sum=12471910, m10 sum=13550355
# merged=13447135
H=LG${LG}_J${J}
ln -fs merged_${H}.fasta ${H}.stitch.fa
LGconsensus ${H}.stitch.1 ${H}.stitch.fa
for i in {1..8}; do LGconsensus ${H}.stitch.$((i+1)) ${H}.stitch.${i}.cons.fa; done
# doesn't converge to 0 varianst added, and q stats rise slightly
# pseudochromosome
echo ">LG$LG" > lg$LG.pchrom.fa
grep -v ">" LG3_J1.stitch.5.cons.fa | tr -d "\n" | fold >> lg$LG.pchrom.fa


LG=4
# hifiasm 28(-) <0.53> 68(+) <0> 41(-) <0> 107(?) <0> 35(+)
# but use hifiasm hic 11 ~= m10 hic 1
# no joins to do, but take consensus for both and cf qual scores
# hifiasm=24292982, m10 sum=23482171
J=0
ln -s HiC_scaffold_1_m10HiC.fa m10.fa
ln -s HiC_scaffold_11_hifiasmHiC.fa hifiasm.fa
for A in hifiasm m10 ; do
    H=$A.LG${LG}_J${J}
    ln -fs $A.fa $H.stitch.fa
    LGconsensus $H.stitch.1 $H.stitch.fa        
    for i in {1..2}; do LGconsensus ${H}.stitch.$((i+1)) ${H}.stitch.${i}.cons.fa; done
done
# hifiasm on top
hifiasm.LG4_J0.stitch.1.cons.log:Applied 140 variants
hifiasm.LG4_J0.stitch.2.cons.log:Applied 52 variants
m10.LG4_J0.stitch.1.cons.log:Applied 766 variants
m10.LG4_J0.stitch.2.cons.log:Applied 96 variants
hifiasm.LG4_J0.stitch.1.bam.stats:reads mapped: 118301
hifiasm.LG4_J0.stitch.2.bam.stats:reads mapped: 119243
m10.LG4_J0.stitch.1.bam.stats:reads mapped:     116646
m10.LG4_J0.stitch.2.bam.stats:reads mapped:     117624
hifiasm.LG4_J0.stitch.1.bam.stats:error rate:   5.371619e-03    hifiasm.LG4_J0.stitch.2.bam.stats:error rate:   4.027733e-03
m10.LG4_J0.stitch.1.bam.stats:error rate:       7.546640e-03
m10.LG4_J0.stitch.2.bam.stats:error rate:       6.187933e-03

A=hifiasm; H=$A.LG${LG}_J${J}
for i in {3..6}; do LGconsensus ${H}.stitch.$((i+1)) ${H}.stitch.${i}.cons.fa; done

# pseudochromosome
echo ">LG$LG" > lg$LG.pchrom.fa
grep -v ">" hifiasm.LG4_J0.stitch.5.cons.fa | tr -d "\n" | fold >> lg$LG.pchrom.fa


LG=5
# m10: 190(?) <0> 3(-)
# or hifiasm 87 in place of m10 190
J=1
# nope
merge_wrapper.py -l 5000 -pre LG${LG}_J${J} contig_190.fa contig_166.fa
# yes
merge_wrapper.py -l 5000 -pre LG${LG}_J${J} ptg000087l.fa contig_166.fa
# contig_166==hifiasm 24,38: 12549184, merged 12682650
H=LG${LG}_J${J}
ln -s merged_${H}.fasta $H.stitch.fa
LGconsensus ${H}.stitch.1 ${H}.stitch.fa
for i in {1..5}; do LGconsensus ${H}.stitch.$((i+1)) ${H}.stitch.${i}.cons.fa; done
# variants: 99, 31, 7, 1, 1, 0
# mean Q: 58.408 > 58.458
# pseudochromosome
echo ">LG$LG" > lg$LG.pchrom.fa
grep -v ">" LG5_J1.stitch.6.cons.fa | tr -d "\n" | fold >> lg$LG.pchrom.fa


LG=6
# hifiasm 43(+) <0> 3(+) <0> 62(-)
# thread with hifiasm hic 9 and 10
J=1
cat ptg0000* > hifiasm.fa
cat HiC_scaffold_* > hic.fa
merge_wrapper.py -l 5000 -pre LG${LG}_J${J} hic.fa hifiasm.fa
# hifiasm sum: 13746859, hic sum 12331126
# merged 13737257
H=LG${LG}_J${J}
ln -s merged_${H}.fasta $H.stitch.fa
LGconsensus ${H}.stitch.1 ${H}.stitch.fa
for i in {1..5}; do LGconsensus ${H}.stitch.$((i+1)) ${H}.stitch.${i}.cons.fa; done
# variants: 89, 28, 21, 13, 2, 1
cat *log
grep "reads mapped:" *stats
grep "bases mapped (cigar)" *stats
# pseudochromosome
echo ">LG$LG" > lg$LG.pchrom.fa
grep -v ">" LG6_J1.stitch.6.cons.fa | tr -d "\n" | fold >> lg$LG.pchrom.fa

# merge
cat $(find . -name *pchrom* | sort) | sed 's/>/\n>/' | tail -n+2 > stitch1.fa

#########################################################################################
# stitch2, first 5.4Mb disco (the rest is hifiasm 28)
# instead try threading hifiasm with HiC 11
# nope, this gives colinearity regardless of ref/query 
LG=4; J=1
cd stitch2/lg4
for i in ptg000028l ptg000068l ptg000041l ptg000107l ptg000035l; do 
    fastaSubsetSingle.py $i ../../hifiasm/hap/QG2082.p_ctg.fa > $i.fa
done
cat ptg* > hifiasm.fa
merge_wrapper.py -l 5000 -pre LG${LG}_J${J} ../../stitch1/lg4//HiC_scaffold_11_hifiasmHiC.fa hifiasm.fa
# try manual stitch. all oriented seqs are concordant
# hifiasm 28(-) <0.53> 68(+) <0> 41(-) <0> 107(?) <0> 35(+)
# orientations differ with hic 11: ignore
# 107 is - strand on hic 11 (28, 41, 35 are all +)
# 107 is SUPER repetitive, with a few islands of complexity
# dupe

# 68 is ordered by 4 markers spanning 2cM
# 41 ordered by 3 markers spanning 2cM

# 35 is ordered by just two markers, ~0.25cM apart. 5' (+) is full of telomeric repeats
# 3' is also repetitive, but not telomeric
# so will assume -

~/bin/seqtk/seqtk seq -r ptg000041l.fa > ptg000041l_rc.fa
mapJunc ptg000068l.fa ptg000041l_rc.fa $READS
# 3 candidates. 2 aln Q20 (samtools). and none for 41 +
doAln ptg000068l.fa ptg000041l_rc.fa
12583429/ccs  12445  11494  12425  -  ptg000068l  15730  16810  686  1103  0
12583429/ccs  12445  11494  12425  -  ptg000068l  16163  17175  666  1024  0
12583429/ccs  12445  11494  12425  -  ptg000068l  16516  17539  669  1024  0
12583429/ccs  12445  11494  12425  -  ptg000068l  17154  18348  707  1195  0
12583429/ccs  12445  11494  12425  -  ptg000041l  10093  11115  662  1023  0
12583429/ccs  12445  11494  12425  -  ptg000041l  10184  11297  693  1113  0
12583429/ccs  12445  11494  12425  -  ptg000041l  10457  11559  703  1114  0
12583429/ccs  12445  11494  12425  -  ptg000041l  4688   5745   690  1080  0
12583429/ccs  12445  11494  12425  -  ptg000041l  5815   6838   657  1025  0
12583429/ccs  12445  11494  12425  -  ptg000041l  7273   8385   685  1114  0
# super repet
148243936/ccs  14408  12685  14352  -  ptg000068l  16163  17810  924   1674   0
148243936/ccs  14408  12685  14390  -  ptg000068l  17935  19837  988   1969   0
148243936/ccs  14408  12727  14390  -  ptg000068l  18197  19979  959   1848   0
148243936/ccs  14408  12685  14390  -  ptg000041l  9510   11557  990   2071   0
# both suggest overlap between terminal ~2kb of 68 and first 5kb of 41
# chew back, fuse and take consensus

J=2
mapJunc ptg000041l_rc.fa ptg000107l.fa $READS
# 0 candidates
mapJunc ptg000041l_rc.fa ptg000107l_rc.fa $READS
# 0 candidates

J=3
# too repetitive for mapping
# mapJunc ptg000107l.fa ptg000035l.fa $READS

# check homology
module load blast+/2.11.0
B=10000
echo -e "ptg000107l\t1\t$B" | bedtools getfasta -fi ptg000107l.fa -bed - | fold > L107.fa
echo -e "ptg000107l\t1\t$B" | bedtools getfasta -fi ptg000107l_rc.fa -bed - | fold > R107.fa
echo -e "ptg000041l\t1\t$B" | bedtools getfasta -fi ptg000041l.fa -bed - | fold > L41.fa
echo -e "ptg000041l\t1\t$B" | bedtools getfasta -fi ptg000041l_rc.fa -bed - | fold > R41.fa
echo -e "ptg000035l\t1\t$B" | bedtools getfasta -fi ptg000035l.fa -bed - | fold > L35.fa
echo -e "ptg000035l\t1\t$B" | bedtools getfasta -fi ptg000035l_rc.fa -bed - | fold > R35.fa
for i in 35 41; do 
    blastn -query L$i.fa -db L107.fa -outfmt 6 -dust no > Q$i.LL.blast
    blastn -query L$i.fa -db R107.fa -outfmt 6 -dust no > Q$i.LR.blast
    blastn -query R$i.fa -db L107.fa -outfmt 6 -dust no > Q$i.RL.blast
    blastn -query R$i.fa -db R107.fa -outfmt 6 -dust no > Q$i.RR.blast
done

# for both 41 and 35, RR hits to 107 only
# but 41 is confidently - so take 41(-) 107(+) 35(-)

# pseudochromosome
# make 68-41 join, patch all others
~/bin/seqtk/seqtk seq -r ptg000028l.fa | grep -v ">" > tmp
echo "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" >> tmp
echo -e "ptg000068l\t1\t1601494" | bedtools getfasta -fi ptg000068l.fa -bed - |\
    grep -v ">" >> tmp
echo -e "ptg000041l\t5000\t2263744" | bedtools getfasta -fi ptg000041l_rc.fa -bed - |\
    grep -v ">" >> tmp    
echo "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" >> tmp
grep -v ">" ptg000107l.fa >> tmp
grep -v ">" ptg000035l_rc.fa >> tmp
echo ">LG$LG" > lg$LG.pchrom.fa
cat tmp | tr -d "\n" | fold >> lg$LG.pchrom.fa

~/bin/mummer-4.0.0rc1/nucmer --maxmatch -l 100 -c 500 -p stitch2v1 lg4.pchrom.fa ../../stitch1//lg4/lg4.pchrom.fa
gzip stitch2v1.delta


LG=1
# dropped ~100kb in stitch1 LG1 
# Hifiasm: 17(-) <1cM> 30(-) <0> 54(?) <0> 6(+)
# hifiasm is expanded in a few places WRT scaffold_8, 
# scaffold_8 links hifiasm 6/54 - but collapsed repeats
# hifiasm HiC 8 links 17 and 30
# DIY
J=1
# 17=1360669, 30=8592421, sum=9953090
# merge=9953590
cat ptg000017l_rc.fa ptg000030l_rc.fa > $J.hifiasm.fa
merge_wrapper.py -l 5000 -pre LG${LG}_J${J} HiC_scaffold_8.fa $J.hifiasm.fa 

J=2
# 26 candidate reads
mapJunc merged_LG1_J1.fasta ptg000054l_rc.fa $READS
doAln merged_LG1_J1.fasta ptg000054l_rc.fa
# taking only spanning Q60 alignments > 2kb
# ~clean spanning alignments
# but highly variable in where they align!
-5  150276163/ccs  15342  16     6952   -  HiC_scaffold_8  8089   15006  6916   6937   60
-12 150276163/ccs  15342  6951   12749  -  ptg000054l      11861  17652  5791   5798   60
-3  165872725/ccs  14591  5410   14565  +  HiC_scaffold_8  8089   17234  9136   9163   60
-12 165872725/ccs  14591  15     5411   +  ptg000054l      12260  17652  5387   5400   60
-6  90114898/ccs   12514  16     5755   -  HiC_scaffold_8  8089   13827  5735   5742   60
-12 90114898/ccs   12514  5754   11547  -  ptg000054l      11861  17652  5789   5795   60
-8  134547784/ccs  15260  11018  15243  +  HiC_scaffold_8  8089   12312  4221   4227   60
-12 134547784/ccs  15260  5223   11019  +  ptg000054l      11861  17652  5791   5796   60

-0  125043886/ccs  17349  2856   14773  +  HiC_scaffold_8  8089   20000  11911  11917  60
-15 125043886/ccs  17349  16     2857   +  ptg000054l      14810  17652  2840   2843   60
-10 143198580/ccs  14677  12297  14660  +  HiC_scaffold_8  8089   10452  2363   2363   60
-12 143198580/ccs  14677  6507   12298  +  ptg000054l      11861  17652  5791   5791   60
-8  43843896/ccs   12904  17     3447   -  HiC_scaffold_8  8089   11518  3429   3430   60
-12 43843896/ccs   12904  3446   9236   -  ptg000054l      11861  17652  5790   5791   60
-7  39321679/ccs   12258  17     4495   -  HiC_scaffold_8  8089   12567  4478   4478   60
-12 39321679/ccs   12258  4494   10285  -  ptg000054l      11861  17652  5791   5791   60
-4  150275482/ccs  12060  3860   12043  +  HiC_scaffold_8  8106   16288  8165   8199   60
-14 150275482/ccs  12060  17     3843   +  ptg000054l      13839  17652  3810   3827   60

# abort 
echo -e "ptg000054l\t11861\t489597" | bedtools getfasta -fi ptg000054l_rc.fa -bed - > tmp.fa
echo -e "HiC_scaffold_8\t1\t9379608" | bedtools getfasta -fi merged_LG1_J1.fasta -bed - > tmp2.fa

mapJunc tmp2.fa tmp.fa $READS
doAln tmp2.fa tmp.fa

ALN=HiC_scaffold_8_ptg000054l.alns
ALN=
sort -k2,2gr -k6,6 $ALN | awk '($11>0 && $10 > 5000)' | cut -f1 | uniq -c | grep -w 2 | trst | cut -f3 > LG$LG.J$J.stitchr; wc -l LG$LG.J$J.stitchr

grep -A3 -f LG$LG.J$J.stitchr $FQ | grep -v "\--" > LG$LG.J$J.stitchr.fq
echo -e "ptg000054l\t1\t1601494" | bedtools getfasta -fi ptg000068l.fa -bed - |\
    grep -v ">" >> tmp


# 54=489597, 6=4773269, sum=5262866
mapJunc ptg000054l_rc.fa ptg000006l.fa $READS
# 17 candidates, 8 mapped
doAln ptg000054l_rc.fa ptg000006l.fa
# from tview, 3 good spanning reads, deleted for 268bp of 54, 202bp of 6
# stitch with 
126223382
118294992
28508883
93325937
minimap2 -cx asm5 --cs ptg000054l_ptg000006l.j.fa ptg000054l_ptg000006l.stitchr.fa
# nope nope nope inverted repeat. abandon. 

# pseudochromosome
# ignore J1 for now
grep -v ">" merged_LG1_J1.fasta > tmp
echo "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" >> tmp
grep -v ">" ptg000054l_rc.fa >> tmp
echo "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN" >> tmp
grep -v ">" ptg000006l.fa >> tmp
echo ">LG$LG" > lg$LG.pchrom.fa
cat tmp | tr -d "\n" | fold >> lg$LG.pchrom.fa


# merge, flip all to + strand
cd ../
~/bin/seqtk/seqtk seq -r lg1/*pchrom.fa | fold > lg1.pchrom_rc.fa
~/bin/seqtk/seqtk seq -r lg4/*pchrom.fa | fold > lg4.pchrom_rc.fa
~/bin/seqtk/seqtk seq -r ../stitch1/lg5/*pchrom.fa | fold  > lg5.pchrom_rc.fa
~/bin/seqtk/seqtk seq -r ../stitch1/lg6/*pchrom.fa | fold > lg6.pchrom_rc.fa
ln -s ../stitch1/lg2/*pchrom.fa .
ln -s ../stitch1/lg3/*pchrom.fa .

cat lg1*fa lg2*fa lg3*fa lg4*fa lg5*fa lg6*fa | sed 's/>/\n>/' | tail -n+2 > stitch2.fa

# OK
# polish and go
REF=stitch2.fa
LGconsensus ${REF}_1.pchrom $REF
for j in {1..3}; do LGconsensus ${REF}_$((j+1)).pchrom ${REF}_${j}.pchrom.cons.fa; done

# take iter 6 on error rate, bases mapped
ln -s stitch2.fa_6.pchrom.cons.fa polished.fa

# rename and reorder by elegans homologs
minimap2 -t 10 -x asm20 stitch2.fa /scratch/lmn3/genome/WS220/c_elegans.WS220.genomic.fa > elegans.paf
for i in {1..6}; do 
    echo $i
    grep "LG${i}" elegans.paf | awk '($12==60)' |\
        awk '{a[$1]+=$11} END {for (i in a) print i, a[i]}' |\
        sort -k2,2gr | awk '($2 > 100000)'
done

# 1
# CHROMOSOME_V 325481
# 2
# CHROMOSOME_II 281240
# 3
# CHROMOSOME_I 231519
# 4
# CHROMOSOME_X 1030023
# 5
# CHROMOSOME_III 196656
# 6
# CHROMOSOME_IV 209037

fastaSubsetSingle.py "LG1" polished.fa | sed 's/LG1/V/' > out1.fa
fastaSubsetSingle.py "LG2" polished.fa | sed 's/LG2/II/' > out2.fa
fastaSubsetSingle.py "LG3" polished.fa | sed 's/LG3/I/' > out3.fa
fastaSubsetSingle.py "LG4" polished.fa | sed 's/LG4/X/' > out4.fa
fastaSubsetSingle.py "LG5" polished.fa | sed 's/LG5/III/' > out5.fa
fastaSubsetSingle.py "LG6" polished.fa | sed 's/LG6/IV/' > out6.fa
cat out*.fa > QG2082_genome_15022021.fa


