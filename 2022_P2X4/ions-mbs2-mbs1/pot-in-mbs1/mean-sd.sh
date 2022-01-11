for k in 1atp 2atp 3atp 1atp.1mg 2atp.2mg

do

paste $k-1.txt $k-2.txt $k-3.txt | awk 'NR<990{print $2, $4, $6}' | awk '{sum=0; sumsq=0; for (i=1; i<=NF; i++) {sum+=$i; sumsq+=$i*$i}}{mean=sum/NF}; {print NR, $1, $2, $3, mean, sqrt((sumsq/NF)-(sum/NF)^2);}' > $k-mean-sd-potct.txt

done
