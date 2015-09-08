#!/bin/csh


# vcf file
set vcf = $argv[1]
# name of the column that contains tumor data 
set name = $argv[2]
# name of the output file 
set out = $argv[3]
# only extract het positions supported by that many reads in the normal
set mindepth = $argv[4]

set field = `awk '{if($1=="#CHROM"){for( i=10;i<=11;i++){if($i=="'$name'"){print i} };exit}}' $vcf`


if( $field == 10 ) then 
echo "Found sample name" $name "in column" $field 
echo CHR POS REF_COUNT VAR_COUNT > $out 
grep PASS $vcf |\
 awk '{nf=split($9,format,":"); ng1=split($11,g1,":"); ng2=split($10,g2,":")}\
      {for(i=1;i<=nf;i++){ g1r[format[i]]=g1[i]; g2r[format[i]]=g2[i]} }\
      g1r["GT"]=="0/1" && ng1==nf && ng2==nf  {print $1,$2,g1r["AD"],g2r["AD"]}'  |\
 sed 's/,/ /g' | awk 'NF==6 && $3+$4>='$mindepth' {print $1,$2,$5,$6}'  >>  $out 
else if( $field == 11 ) then 
echo "Found sample name" $name "in column" $field 
echo CHR POS REF_COUNT VAR_COUNT > $out 
grep PASS $vcf |\
 awk '{nf=split($9,format,":"); ng1=split($10,g1,":"); ng2=split($11,g2,":")}\
      {for(i=1;i<=nf;i++){ g1r[format[i]]=g1[i]; g2r[format[i]]=g2[i]} }\
      g1r["GT"]=="0/1" && ng1==nf && ng2==nf  {print $1,$2,g1r["AD"],g2r["AD"]}'  |\
 sed 's/,/ /g' | awk 'NF==6 && $3+$4>='$mindepth'  {print $1,$2,$5,$6}'  >>  $out 
else 
 echo "Could not find" $name "in vcf's header, make sure it is found in column 10 or 11"
endif







