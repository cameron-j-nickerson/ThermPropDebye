#!/bin/bash

for i in */
do

cd ${i}

grep -A6 "Elastic constant" elastics.out > tmp
tail -n +2 tmp > tmp2
cat tmp2 > ${i%/}.clmz
rm tmp
rm tmp2
echo "" >> ${i%/}.clmz

grep "a:" elastics.out | sed 's/a://' > tmp
grep "b:" elastics.out | sed 's/b://' >> tmp
grep "c:" elastics.out | sed 's/c://' >> tmp
awk '{$1*=0.529177249;$2*=0.529177249;$3*=0.529177249;print}' tmp > tmp2
cat tmp2 >> ${i%/}.clmz
rm tmp
rm tmp2
echo "" >> ${i%/}.clmz

grep "Molar mass" elastics.out | awk '{print $4}' >> ${i%/}.clmz
echo "" >> ${i%/}.clmz

echo "crystal geometry.in" > ding
critic2 < ding > dong
grep "molecules per cell" dong | awk '{print $8}' >> ${i%/}.clmz
rm ding
rm dong

cd ..

done

