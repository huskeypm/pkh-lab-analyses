#!/bin/bash


# dimention of hd
HD=32
LD=2


FILEHD='landmark'
FILELD='lowd'

# high dimension sigma, a, b [e.g. 6.0 2 6 ]
SIGMAHD=18
AHD=8
BHD=6

# low dimension sigma, a, b 
SIGMALD=18
ALD=1
BLD=2


GW=100
#Step1 is the 'distance matching optimization' 
# No sigmoid Function 
# Purpose: get a relative initial start conformation

grep -v "#" $FILEHD > t
mv t $FILEHD
rm -rf initial_error
dimred -vv -D $HD -d $LD -w -preopt 20 -grid $GW,20,200 -gopt 5 < $FILEHD > tmp 2>>initial_error


#Step2: mixing 'distance matching optimization' and 'sketch-map optimization' using the
# '-imix' option. the value of 'imix' should be slowly annealed from 1.0 to 0.0
NERR=`awk '/Error/{print $(NF)}'  tmp | tail -n 1`
SMERR=`awk '/Error/{print $(NF)}'  tmp | tail -n 1`



IMIX=1.0
MAXITER=50
for ((ITER=1; ITER<=$MAXITER; ITER++)); do
   MDERR=$NERR
   echo "Mixing in $IMIX"
   if [ ! -e $FILELD.gmds_$ITER ]; then
       
	grep -v "#" tmp > t
        mv t tmp	
      echo "Now running dimred -D $HD -d $LD -w -preopt 50 -grid $GW,20,200 -fun-hd $SIGMAHD,$AHD,$BHD -fun-ld $SIGMALD,$ALD,$BLD -init tmp -gopt 5 -imix $IMIX < $FILEHD > $FILELD.gmds_$ITER 2>>log"
      grep -v \#  $FILEHD | dimred -vv -D $HD -d $LD -w -preopt 50 -grid 100,20,200 -fun-hd $SIGMAHD,$AHD,$BHD -fun-ld $SIGMALD,$ALD,$BLD -init tmp -gopt 5 -imix $IMIX < $FILEHD > $FILELD.gmds_$ITER 2>>log
   fi
   grep -v "#" $FILELD.gmds_$ITER | awk '{print $1, $2}' > tmp
   GW=`awk 'BEGIN{maxr=0} !/#/{r=sqrt($1*$1+$2*$2); if (r>maxr) maxr=r} END{print maxr*1.2}' $FILELD.gmds_$ITER`;
   NERR=`awk '/Error/{print $(NF)}' $FILELD.gmds_$ITER  | tail -n 1 `
   echo "Residual error is $NERR"
   IMIX=`echo "$IMIX  $SMERR  $NERR" | awk '{new=$2/($2+$3); if (new<0.1) new=0.1; if (new>0.5) new=0.5; print new*$1 }'`
   echo "DEBUG  $MDERR $NERR"
   if [ ` echo $MDERR $NERR | awk -v i=$ITER '{ if (i>1 && (($1-$2)/$2)*(($1-$2)/$2)<1e-4) print "done"; else print "nope";}' ` = "done" ]; then ((ITER++)); break; fi;
done

echo "Doing final fit"
((ITER--))
grep -v "#" $FILELD.gmds_$ITER | awk '{print $1, $2}' > tmp
grep -v \#  $FILEHD | dimred -vv -D $HD -d $LD -w -preopt 50 -grid $GW,20,200 -fun-hd $SIGMAHD,$AHD,$BHD -fun-ld $SIGMALD,$ALD,$BLD -init tmp -gopt 5 > $FILELD.gmds 2>>log

