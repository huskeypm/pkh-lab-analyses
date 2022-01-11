for k in 1atp 2atp 3atp 2atp.2mg 1atp.1mg

do
paste $k-1-wrap.txt $k-2-wrap.txt $k-3-wrap.txt | head -990 | awk '{print $2, $5, $8}' | awk '{sum=0; sumsq=0; for (i=1; i<=NF; i++) {sum+=$i; sumsq+=$i*$i}}{mean=sum/NF}; {print NR, $1, $2, $3, mean, sqrt((sumsq/NF)-(sum/NF)^2);}' > $k-mean-sd.txt
done

