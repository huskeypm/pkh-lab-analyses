for i in 3atp 2atp 1atp 2atp.2mg 1atp.1mg

do

for j in o1g o2g o3g o3b o1b o2b o3a o1a o2a o5\' o2\' o3\' o4\'

do

paste <(cat mg-$i-dist.noe.dat | grep "STATISTICS") <(cat mg-$i-dist.noe.dat | grep "AVERAGE") | awk '{$1=""; $3=""}1' | awk '{split($1,a,"-"); print a[1], a[2], $2, $3, $4}' | awk -v i="$i" -v j="$j" -F " " '$2==j{sum3+=$3; count ++; mean=sum3/count;} END {printf "%-5s %-5s %.2f %-5s %.2f\n",  i, j, sum3, count, mean}' >> summary.txt

done
done
