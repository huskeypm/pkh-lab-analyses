rm mean-sd.txt 

for k in 3atp 2atp.2mg 1atp.1mg 2atp 1atp

do
cat $k-*.txt | awk '{$1=""}1' | xargs -n1 | awk -v k="$k" '{sum+=$1; mean=sum/NR; sumsq+=$1*$1} END {print k, mean, sqrt((sumsq/NR)-(sum/NR)^2)}' >> mean-sd.txt
done
