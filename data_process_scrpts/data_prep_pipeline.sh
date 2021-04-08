vcf=''
filter=''
filter_file='None'
verbose='false'

print_usage() {
  printf "Usage: ..."
}

while getopts 'i:h:f:v:' opt; do
  case $opt in
    i) vcf="$OPTARG" ;;
    h) filter="$OPTARG" ;;
    f) filter_file="$OPTARG" ;;
    v) verbose='true' ;;
    *) print_usage
       exit 1 ;;
  esac
done
echo $vcf
echo $filter
echo $filter_file
if [ $filter_file == "None" ]
then
vcftools --vcf $vcf --mac 2 --max-missing 0.8 --recode --out ${filter}_${vcf/.vcf}
else
vcftools --vcf $vcf --exclude $filter_file --mac 2 --max-missing 0.8 --recode --out ${filter}_${vcf/.vcf/}
fi
tail -n +16 ${filter}_${vcf/.vcf/}.recode.vcf| awk '{print $3}' > loci_list.txt
shuf -n 4000 loci_list.txt > loci_sublist_r1.txt
vcftools --vcf ${filter}_${vcf/.vcf/}.recode.vcf --snps loci_sublist_r1.txt --recode  --out ${filter}_${vcf/.vcf/}_subrep1
shuf -n 4000 loci_list.txt > loci_sublist_r2.txt
vcftools --vcf ${filter}_${vcf/.vcf/}.recode.vcf --snps loci_sublist_r2.txt --recode  --out ${filter}_${vcf/.vcf/}_subrep2
shuf -n 4000 loci_list.txt > loci_sublist_r3.txt
vcftools --vcf ${filter}_${vcf/.vcf/}.recode.vcf --snps loci_sublist_r3.txt --recode  --out ${filter}_${vcf/.vcf/}_subrep3
shuf -n 4000 loci_list.txt > loci_sublist_r4.txt
vcftools --vcf ${filter}_${vcf/.vcf/}.recode.vcf --snps loci_sublist_r4.txt --recode  --out ${filter}_${vcf/.vcf/}_subrep4
shuf -n 4000 loci_list.txt > loci_sublist_r5.txt
vcftools --vcf ${filter}_${vcf/.vcf/}.recode.vcf --snps loci_sublist_r5.txt --recode  --out ${filter}_${vcf/.vcf/}_subrep5
shuf -n 4000 loci_list.txt > loci_sublist_r6.txt
vcftools --vcf ${filter}_${vcf/.vcf/}.recode.vcf --snps loci_sublist_r6.txt --recode  --out ${filter}_${vcf/.vcf/}_subrep6
shuf -n 4000 loci_list.txt > loci_sublist_r7.txt
vcftools --vcf ${filter}_${vcf/.vcf/}.recode.vcf --snps loci_sublist_r7.txt --recode  --out ${filter}_${vcf/.vcf/}_subrep7
shuf -n 4000 loci_list.txt > loci_sublist_r8.txt
vcftools --vcf ${filter}_${vcf/.vcf/}.recode.vcf --snps loci_sublist_r8.txt --recode  --out ${filter}_${vcf/.vcf/}_subrep8
shuf -n 4000 loci_list.txt > loci_sublist_r9.txt
vcftools --vcf ${filter}_${vcf/.vcf/}.recode.vcf --snps loci_sublist_r9.txt --recode  --out ${filter}_${vcf/.vcf/}_subrep9
shuf -n 4000 loci_list.txt > loci_sublist_r10.txt
vcftools --vcf ${filter}_${vcf/.vcf/}.recode.vcf --snps loci_sublist_r10.txt --recode  --out ${filter}_${vcf/.vcf/}_subrep10

for FILE in *recode.vcf
do
plink --vcf $FILE --recode structure --out ${FILE/.recode.vcf/}.str
done

# mv ${filter}_${vcf/.vcf/}_rep1.str.recode.strct_in ${filter}_${vcf/.vcf/}_rep1.str

for FILE in *strct_in
do
mv $FILE ${FILE/.recode.strct_in/}
done
rm *nosex
rm *log

#plink1.9 --vcf ${filter}_${vcf/.vcf/}_rep1.recode.vcf --recode structure --out ${filter}_${vcf/.vcf/}_rep1.str
