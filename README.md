# YoungColorectalCancer
Script to process 111 young colorectal cancer





### 210 retrieved


## Filter by GATK4 ( NO efficient memory )
```shell
for i in *.mutect.vcf.gz
do
  nohup ./filter.sh $i >log.`basename $i .mutect.vcf.gz` &
done
```

## Retain the passed variants
```shell
for i in *oncefiltered.vcf.gz
do
  zcat $i |awk '/^#/ || $7~/PASS/' - >`basename $i .somatic_oncefiltered.vcf.gz`.clean.vcf
done
```
## Convert the vcf to bed
```shell
for i in *.clean.vcf
do
  vcf2bed < $i  >`basename $i .vcf`.bed
done
```
## Combine the bed files
```shell
if [ -s merged.bed ];then rm merged.bed;fi
for i in *.clean.bed
do
  name=`basename $i .clean.bed`
  awk 'BEGIN{OFS="\t"}{print $1,$2-length($6)+2,$3,$6,$7,"'$name'"}' $i >>merged.bed
done
sort -V merged.bed >merged.sorted.bed
```

## Annotation
```shell
annotate_variation.pl -buildver hg19 -outfile merged --otherinfo merged.sorted.bed  /database/annotation/annovar/humandb/
```
