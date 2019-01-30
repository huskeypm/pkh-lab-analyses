#!/bin/bash

#export LOC=/home/bch265/sources
#export MYPATH=$LOC/mypython
#export PYTHONPATH=$PYTHONPATH:/home/bch265/.local/lib/python2.7/site-packages/

#python2.7 -c "import instant"
#python2.7 -c "import modelparameters"
#python2.7 -c "import gotran"

############## Control Panel #####################
#export ODEFILEfull=microglia_revision-backup19.ode
#export ODEFILElumped=microglia_revision_lumped.ode
#export ODEFILEfull=microglia_revision.ode
export ODEFILEfull=microglia_revision_newpot.ode

## Phase 1 #########
ptxfvalid=0 # P2X4 Channel validation
ptxfvalidCa=0 # P2X4 Ca transient validation
ptxsvalid=0 # P2X7 Channel validation
## Phase 2 #########
Eq=0        ### No stimulation
caivalid=0  ## Ca transient validation against Hide's data
######## Buffer Test #######
BufferEQ=0  ##
Buffered=0  ## This must take microglia_revision-buffer.ode as an input file
NoFura=0    ## No Fura2
NoBuffer=0  ## No buffer
#########################
CaN=0       ### Special Edition for CaN Validation -> Deactivation test
CaN2=0      ### Test for CaMCN activation with respect to ATP dosage
## Phase 3 #########
leng1=0      ### (no pulsatile simulation) NFAT fitting
leng2=0      ### (no pulsatile simulation) pp38 fitting
leng3=0      ### (no pulsatile simulation) TNFa -> 3 hours
## Phase 4 #########
dura=0      ### Duration variation
freq1=0     ### Frequency variation 0.033 Hz - 0.5 Hz Stim duration = 1 sec
freq2=0     ### Frequency variation 0.5 Hz - 5 Hz Stim duration = 0.1 sec
freq3=1     ### Frequency variation 2 - 5 Hz stim duration 0.1 sec Activated vs. Rest
## Phase 5 #########
sens1=0     ### Sensitivity Analysis
sens2=0
sens3=0
sens4=0

#####################################
# for the periodicity control
# check the google spreadsheet
# https://docs.google.com/spreadsheets/d/1UpXEe90q5Tlcu7qWv07B5SrVTGP-gYaGT6FGdsQigr8/edit?usp=sharing
# test with pulse control tester
# microglia/gotran/Pulse-control-tester.ipynb in Bitbucket
# P2X Validation ####################

## Test Zone
#python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 10800000 -iters 1 -var pulse_switch 0 -var stim_amplitude 3000 -name ~/Data_storage/180min_MG_3mMATP_test
###

if [ $ptxfvalid -eq 1 ]
then
  echo "Validation of P2X4-fitted is ON" # This code needs to be fixed ( millisecond time unit )
 python2.7 daisychain.py -dt 0.1 -dSr 1000 -jit -odeName $ODEFILEfull -T 500e3 -iters 1 -var V_ptxf -0.06 -var stim_amplitude 10 -var stim_period 300e3 -var stim_gap1 270e3 -var stim_gap2 270e3 -var stim_low 1e3 -var stim_high 30e3 -name ~/Data_storage/p2x4_Toulme30st_10uMATP_total_full &
 python2.7 daisychain.py -dt 0.1 -dSr 1000 -jit -odeName $ODEFILEfull -T 500e3 -iters 1 -var V_ptxf -0.06 -var stim_amplitude 100 -var stim_period 300e3 -var stim_gap1 270e3 -var stim_gap2 270e3 -var stim_low 1e3 -var stim_high 30e3 -name ~/Data_storage/p2x4_Toulme30st_100uMATP_total_full &
 python2.7 daisychain.py -dt 0.1 -dSr 1000 -jit -odeName $ODEFILEfull -T 500e3 -iters 1 -var V_ptxf -0.06 -var stim_amplitude 1000 -var stim_period 300e3 -var stim_gap1 270e3 -var stim_gap2 270e3 -var stim_low 1e3 -var stim_high 30e3 -name ~/Data_storage/p2x4_Toulme30st_1000uMATP_total_full

fi

if [ $ptxfvalidCa -eq 1 ]
then
  echo "Validation of P2X4-Ca Transient-fitted is ON" # This code needs to be fixed ( millisecond time unit )
 python2.7 daisychain.py -dt 0.1 -dSr 1000 -jit -odeName $ODEFILEfull -T 500e3 -iters 1 -var rMiG 2.68e-6 -var expF_ptxs 0 -var V_ptxf -0.06 -var stim_amplitude 10 -var stim_period 300e3 -var stim_gap1 270e3 -var stim_gap2 270e3 -var stim_low 1e3 -var stim_high 30e3 -name ~/Data_storage/p2x4_Toulme30st_10uMATP_CaT &
 python2.7 daisychain.py -dt 0.1 -dSr 1000 -jit -odeName $ODEFILEfull -T 500e3 -iters 1 -var rMiG 2.68e-6 -var expF_ptxs 0 -var V_ptxf -0.06 -var stim_amplitude 100 -var stim_period 300e3 -var stim_gap1 270e3 -var stim_gap2 270e3 -var stim_low 1e3 -var stim_high 30e3 -name ~/Data_storage/p2x4_Toulme30st_100uMATP_CaT &
 python2.7 daisychain.py -dt 0.1 -dSr 1000 -jit -odeName $ODEFILEfull -T 500e3 -iters 1 -var rMiG 2.68e-6 -var expF_ptxs 0 -var V_ptxf -0.06 -var stim_amplitude 1000 -var stim_period 300e3 -var stim_gap1 270e3 -var stim_gap2 270e3 -var stim_low 1e3 -var stim_high 30e3 -name ~/Data_storage/p2x4_Toulme30st_1000uMATP_CaT

fi

if [ $ptxsvalid -eq 1 ]
then
  echo "Validation of P2X7-fitted is ON" # This code needs to be fixed ( millisecond time unit )
 python2.7 daisychain.py -dt 0.1 -dSr 1000 -jit -odeName $ODEFILEfull -T 60e3 -iters 1 -var V_ptxs -0.09 -var stim_amplitude 100 -var stim_period 10e3 -var stim_gap1 9e3 -var stim_gap2 9e3 -var stim_low 1e3 -var stim_high 2e3 -name ~/Data_storage/p2x7_Chessel1st_100uMATP_full &
 python2.7 daisychain.py -dt 0.1 -dSr 1000 -jit -odeName $ODEFILEfull -T 60e3 -iters 1 -var V_ptxs -0.09 -var stim_amplitude 1000 -var stim_period 10e3 -var stim_gap1 9e3 -var stim_gap2 9e3 -var stim_low 1e3 -var stim_high 2e3 -name ~/Data_storage/p2x7_Chessel1st_1000uMATP_full &
 python2.7 daisychain.py -dt 0.1 -dSr 1000 -jit -odeName $ODEFILEfull -T 60e3 -iters 1 -var V_ptxs -0.09 -var stim_amplitude 3000 -var stim_period 10e3 -var stim_gap1 9e3 -var stim_gap2 9e3 -var stim_low 1e3 -var stim_high 2e3 -name ~/Data_storage/p2x7_Chessel1st_3000uMATP_full

fi

if [ $Eq -eq 1 ]
then
  echo "No stimulation is ON"
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 20000000 -iters 1 -var pulse_switch 0 -var stim_amplitude 0 -name ~/Data_storage/MG_0mMATP
fi

if [ $caivalid -eq 1 ]
then
  echo "Validation of Cai is ON"
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 1000 -name ~/Data_storage/8min_MG_1mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 100 -name ~/Data_storage/8min_MG_0.1mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 10 -name ~/Data_storage/8min_MG_0.01mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 0 -name ~/Data_storage/8min_MG_0mMATP
fi

if [ $BufferEQ -eq 1 ]
then
  echo "No stimulation for Buffer is ON"
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var pulse_switch 0 -var stim_amplitude 0 -name ~/Data_storage/MG_buffer_0mMATP
fi

if [ $Buffered -eq 1 ]
then
  echo "All Buffer is ON"
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 1000 -name ~/Data_storage/8min_MG_buffered_1mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 100 -name ~/Data_storage/8min_MG_buffered_0.1mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 10 -name ~/Data_storage/8min_MG_buffered_0.01mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 0 -name ~/Data_storage/8min_MG_buffered_0mMATP
fi

if [ $NoFura -eq 1 ]
then
  echo "No Fura is ON"
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var Fura 0 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 1000 -name ~/Data_storage/8min_MG_NoFura_1mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var Fura 0 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 100 -name ~/Data_storage/8min_MG_NoFura_0.1mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var Fura 0 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 10 -name ~/Data_storage/8min_MG_NoFura_0.01mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var Fura 0 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 0 -name ~/Data_storage/8min_MG_NoFura_0mMATP
fi

if [ $NoBuffer -eq 1 ]
then
  echo "No Buffer is ON"
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var Buffer 0 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 1000 -name ~/Data_storage/8min_MG_NoBuffer_1mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var Buffer 0 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 100 -name ~/Data_storage/8min_MG_NoBuffer_0.1mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var Buffer 0 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 10 -name ~/Data_storage/8min_MG_NoBuffer_0.01mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 2000000 -iters 1 -var Buffer 0 -var stim_period 800e3 -var stim_gap1 400e3 -var stim_gap2 750e3 -var stim_low 1e3 -var stim_high 400e3 -var stim_amplitude 0 -name ~/Data_storage/8min_MG_NoBuffer_0mMATP
fi

if [ $CaN -eq 1 ]
then
  echo "Validation of CaN is ON"

  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 500000 -iters 1 -var stim_period 2e3 -var stim_gap1 1e3 -var stim_gap2 1e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1000 -name ~/Data_storage/1st1rsMGCaN_1mMATPnew
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 500000 -iters 1 -var stim_period 201e3 -var stim_gap1 200e3 -var stim_gap2 200e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1000 -name ~/Data_storage/1st200rsMGCaN_1mMATPnew

fi

if [ $CaN2 -eq 1 ]
then
  echo "Validation of CaN is ON"

  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 500000 -iters 1 -var stim_period 201e3 -var stim_gap1 200e3 -var stim_gap2 200e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 0 -name ~/Data_storage/1st200rsMGCaN_0mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 500000 -iters 1 -var stim_period 201e3 -var stim_gap1 200e3 -var stim_gap2 200e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 10 -name ~/Data_storage/1st200rsMGCaN_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 500000 -iters 1 -var stim_period 201e3 -var stim_gap1 200e3 -var stim_gap2 200e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 100 -name ~/Data_storage/1st200rsMGCaN_0.1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 500000 -iters 1 -var stim_period 201e3 -var stim_gap1 200e3 -var stim_gap2 200e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1000 -name ~/Data_storage/1st200rsMGCaN_1mMATPnew

fi
#####################################

#####################################

# Make sure your pre_rest is less than 500; otherwise, the calculation may crash.
# -var pre_rest 100
# Various duration with broad spacing at 1 mM ATP: 1,2,5, and 10 are only size working in this model due to periodic equation

if [ $dura -eq 1 ]
then
 echo "Duration variation is ON"

  # Various duration with broad spacing at 1 mM ATP: 1,2,5, and 10 are only size working in this model due to periodic equation

  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 19e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1000 -name ~/Data_storage/1st19rsMG_1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 15e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 5e3 -var stim_amplitude 1000 -name ~/Data_storage/5st15rsMG_1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 10e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 10e3 -var stim_amplitude 1000 -name ~/Data_storage/10st10rsMG_1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 5e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 15e3 -var stim_amplitude 1000 -name ~/Data_storage/15st5rsMG_1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 1e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 19e3 -var stim_amplitude 1000 -name ~/Data_storage/19st1rsMG_1mMATPnew


  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 19e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 100 -name ~/Data_storage/1st19rsMG_0.1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 15e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 5e3 -var stim_amplitude 100 -name ~/Data_storage/5st15rsMG_0.1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 10e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 10e3 -var stim_amplitude 100 -name ~/Data_storage/10st10rsMG_0.1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 5e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 15e3 -var stim_amplitude 100 -name ~/Data_storage/15st5rsMG_0.1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 1e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 19e3 -var stim_amplitude 100 -name ~/Data_storage/19st1rsMG_0.1mMATPnew


  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 19e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 10 -name ~/Data_storage/1st19rsMG_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 15e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 5e3 -var stim_amplitude 10 -name ~/Data_storage/5st15rsMG_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 10e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 10e3 -var stim_amplitude 10 -name ~/Data_storage/10st10rsMG_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 5e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 15e3 -var stim_amplitude 10 -name ~/Data_storage/15st5rsMG_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 1e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 19e3 -var stim_amplitude 10 -name ~/Data_storage/19st1rsMG_0.01mMATPnew


  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 19e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1 -name ~/Data_storage/1st19rsMG_0.001mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 15e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 5e3 -var stim_amplitude 1 -name ~/Data_storage/5st15rsMG_0.001mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 10e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 10e3 -var stim_amplitude 1 -name ~/Data_storage/10st10rsMG_0.001mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 5e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 15e3 -var stim_amplitude 1 -name ~/Data_storage/15st5rsMG_0.001mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 1e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 19e3 -var stim_amplitude 1 -name ~/Data_storage/19st1rsMG_0.001mMATPnew

#  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 200000 -iters 1 -var pulse_switch 0 -var stim_amplitude 0 -name ~/Data_storage/180min_MG_0mMATP

fi

if [ $freq1 -eq 1 ]
then
  echo " Low Frequency variation is ON"
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 2e3 -var stim_gap1 1e3 -var stim_gap2 1e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1000 -name ~/Data_storage/1st1rsMG_1mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 5e3 -var stim_gap1 4e3 -var stim_gap2 4e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1000 -name ~/Data_storage/1st4rsMG_1mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 10e3 -var stim_gap1 9e3 -var stim_gap2 9e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1000 -name ~/Data_storage/1st9rsMG_1mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 19e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1000 -name ~/Data_storage/1st19rsMG_1mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 30e3 -var stim_gap1 29e3 -var stim_gap2 29e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1000 -name ~/Data_storage/1st29rsMG_1mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 50e3 -var stim_gap1 49e3 -var stim_gap2 49e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1000 -name ~/Data_storage/1st49rsMG_1mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 60e3 -var stim_gap1 59e3 -var stim_gap2 59e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1000 -name ~/Data_storage/1st59rsMG_1mMATPnew

  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 2e3 -var stim_gap1 1e3 -var stim_gap2 1e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 100 -name ~/Data_storage/1st1rsMG_0.1mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 5e3 -var stim_gap1 4e3 -var stim_gap2 4e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 100 -name ~/Data_storage/1st4rsMG_0.1mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 10e3 -var stim_gap1 9e3 -var stim_gap2 9e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 100 -name ~/Data_storage/1st9rsMG_0.1mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 19e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 100 -name ~/Data_storage/1st19rsMG_0.1mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 30e3 -var stim_gap1 29e3 -var stim_gap2 29e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 100 -name ~/Data_storage/1st29rsMG_0.1mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 50e3 -var stim_gap1 49e3 -var stim_gap2 49e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 100 -name ~/Data_storage/1st49rsMG_0.1mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 60e3 -var stim_gap1 59e3 -var stim_gap2 59e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 100 -name ~/Data_storage/1st59rsMG_0.1mMATPnew

  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 2e3 -var stim_gap1 1e3 -var stim_gap2 1e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 10 -name ~/Data_storage/1st1rsMG_0.01mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 5e3 -var stim_gap1 4e3 -var stim_gap2 4e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 10 -name ~/Data_storage/1st4rsMG_0.01mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 10e3 -var stim_gap1 9e3 -var stim_gap2 9e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 10 -name ~/Data_storage/1st9rsMG_0.01mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 19e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 10 -name ~/Data_storage/1st19rsMG_0.01mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 30e3 -var stim_gap1 29e3 -var stim_gap2 29e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 10 -name ~/Data_storage/1st29rsMG_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 50e3 -var stim_gap1 49e3 -var stim_gap2 49e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 10 -name ~/Data_storage/1st49rsMG_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 60e3 -var stim_gap1 59e3 -var stim_gap2 59e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 10 -name ~/Data_storage/1st59rsMG_0.01mMATPnew

  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 2e3 -var stim_gap1 1e3 -var stim_gap2 1e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1 -name ~/Data_storage/1st1rsMG_0.001mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 5e3 -var stim_gap1 4e3 -var stim_gap2 4e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1 -name ~/Data_storage/1st4rsMG_0.001mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 10e3 -var stim_gap1 9e3 -var stim_gap2 9e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1 -name ~/Data_storage/1st9rsMG_0.001mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 20e3 -var stim_gap1 19e3 -var stim_gap2 19e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1 -name ~/Data_storage/1st19rsMG_0.001mMATPnew &
  #python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 30e3 -var stim_gap1 29e3 -var stim_gap2 29e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1 -name ~/Data_storage/1st29rsMG_0.001mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 50e3 -var stim_gap1 49e3 -var stim_gap2 49e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1 -name ~/Data_storage/1st49rsMG_0.001mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var stim_period 60e3 -var stim_gap1 59e3 -var stim_gap2 59e3 -var stim_low 1e3 -var stim_high 2e3 -var stim_amplitude 1 -name ~/Data_storage/1st59rsMG_0.001mMATPnew

fi

if [ $freq2 -eq 1 ]
then
  echo "High Frequency variation is ON"
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 2e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 1.85e3 -var stim_gap2 1.85e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 1000 -name ~/Data_storage/0.5HzMG_1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 1e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.85e3 -var stim_gap2 0.85e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 1000 -name ~/Data_storage/1HzMG_1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.5e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.35e3 -var stim_gap2 0.35e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 1000 -name ~/Data_storage/2HzMG_1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.25e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.11e3 -var stim_gap2 0.11e3 -var stim_low 0.1e3 -var stim_high 0.25e3 -var stim_amplitude 1000 -name ~/Data_storage/4HzMG_1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.215e3 -var IC1 0.2e3 -var IC2 0.14e3 -var stim_gap1 0.02e3 -var stim_gap2 0.02e3 -var stim_low 0.15e3 -var stim_high 0.25e3 -var stim_amplitude 1000 -name ~/Data_storage/5HzMG_1mMATPnew

  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 2e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 1.85e3 -var stim_gap2 1.85e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 100 -name ~/Data_storage/0.5HzMG_0.1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 1e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.85e3 -var stim_gap2 0.85e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 100 -name ~/Data_storage/1HzMG_0.1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.5e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.35e3 -var stim_gap2 0.35e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 100 -name ~/Data_storage/2HzMG_0.1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.25e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.11e3 -var stim_gap2 0.11e3 -var stim_low 0.1e3 -var stim_high 0.25e3 -var stim_amplitude 100 -name ~/Data_storage/4HzMG_0.1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.215e3 -var IC1 0.2e3 -var IC2 0.14e3 -var stim_gap1 0.02e3 -var stim_gap2 0.02e3 -var stim_low 0.15e3 -var stim_high 0.25e3 -var stim_amplitude 100 -name ~/Data_storage/5HzMG_0.1mMATPnew

  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 2e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 1.85e3 -var stim_gap2 1.85e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 10 -name ~/Data_storage/0.5HzMG_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 1e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.85e3 -var stim_gap2 0.85e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 10 -name ~/Data_storage/1HzMG_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.5e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.35e3 -var stim_gap2 0.35e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 10 -name ~/Data_storage/2HzMG_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.25e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.11e3 -var stim_gap2 0.11e3 -var stim_low 0.1e3 -var stim_high 0.25e3 -var stim_amplitude 10 -name ~/Data_storage/4HzMG_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.215e3 -var IC1 0.2e3 -var IC2 0.14e3 -var stim_gap1 0.02e3 -var stim_gap2 0.02e3 -var stim_low 0.15e3 -var stim_high 0.25e3 -var stim_amplitude 10 -name ~/Data_storage/5HzMG_0.01mMATPnew

  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 2e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 1.85e3 -var stim_gap2 1.85e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 2000 -name ~/Data_storage/0.5HzMG_2mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 1e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.85e3 -var stim_gap2 0.85e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 2000 -name ~/Data_storage/1HzMG_2mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.5e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.35e3 -var stim_gap2 0.35e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 2000 -name ~/Data_storage/2HzMG_2mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.25e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.11e3 -var stim_gap2 0.11e3 -var stim_low 0.1e3 -var stim_high 0.25e3 -var stim_amplitude 2000 -name ~/Data_storage/4HzMG_2mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.215e3 -var IC1 0.2e3 -var IC2 0.14e3 -var stim_gap1 0.02e3 -var stim_gap2 0.02e3 -var stim_low 0.15e3 -var stim_high 0.25e3 -var stim_amplitude 2000 -name ~/Data_storage/5HzMG_2mMATPnew

fi

if [ $freq3 -eq 1 ]
then
  echo "Activation vs Rest test is ON"

  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var rho4 30 -var stim_period 0.5e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.35e3 -var stim_gap2 0.35e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 100 -name ~/Data_storage/AMG_2HzMG_0.1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var rho4 30 -var stim_period 0.25e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.11e3 -var stim_gap2 0.11e3 -var stim_low 0.1e3 -var stim_high 0.25e3 -var stim_amplitude 100 -name ~/Data_storage/AMG_4HzMG_0.1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var rho4 30 -var stim_period 0.215e3 -var IC1 0.2e3 -var IC2 0.14e3 -var stim_gap1 0.02e3 -var stim_gap2 0.02e3 -var stim_low 0.15e3 -var stim_high 0.25e3 -var stim_amplitude 100 -name ~/Data_storage/AMG_5HzMG_0.1mMATPnew

  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var rho4 30 -var stim_period 0.5e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.35e3 -var stim_gap2 0.35e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 10 -name ~/Data_storage/AMG_2HzMG_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var rho4 30 -var stim_period 0.25e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.11e3 -var stim_gap2 0.11e3 -var stim_low 0.1e3 -var stim_high 0.25e3 -var stim_amplitude 10 -name ~/Data_storage/AMG_4HzMG_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var rho4 30 -var stim_period 0.215e3 -var IC1 0.2e3 -var IC2 0.14e3 -var stim_gap1 0.02e3 -var stim_gap2 0.02e3 -var stim_low 0.15e3 -var stim_high 0.25e3 -var stim_amplitude 10 -name ~/Data_storage/AMG_5HzMG_0.01mMATPnew

  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.5e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.35e3 -var stim_gap2 0.35e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 100 -name ~/Data_storage/RMG_2HzMG_0.1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.25e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.11e3 -var stim_gap2 0.11e3 -var stim_low 0.1e3 -var stim_high 0.25e3 -var stim_amplitude 100 -name ~/Data_storage/RMG_4HzMG_0.1mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.215e3 -var IC1 0.2e3 -var IC2 0.14e3 -var stim_gap1 0.02e3 -var stim_gap2 0.02e3 -var stim_low 0.15e3 -var stim_high 0.25e3 -var stim_amplitude 100 -name ~/Data_storage/RMG_5HzMG_0.1mMATPnew

  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.5e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.35e3 -var stim_gap2 0.35e3 -var stim_low 0e3 -var stim_high 0.25e3 -var stim_amplitude 10 -name ~/Data_storage/RMG_2HzMG_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.25e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.11e3 -var stim_gap2 0.11e3 -var stim_low 0.1e3 -var stim_high 0.25e3 -var stim_amplitude 10 -name ~/Data_storage/RMG_4HzMG_0.01mMATPnew &
  python2.7 daisychain.py -dt 1 -dSr 1 -jit -odeName $ODEFILEfull -T 30000 -iters 1 -var stim_period 0.215e3 -var IC1 0.2e3 -var IC2 0.14e3 -var stim_gap1 0.02e3 -var stim_gap2 0.02e3 -var stim_low 0.15e3 -var stim_high 0.25e3 -var stim_amplitude 10 -name ~/Data_storage/RMG_5HzMG_0.01mMATPnew

fi

### SENSITIVITY ANALYSIS ###
if [ $sens1 -eq 1 ]
then
array[0]=".01"
array[1]=".1"
array[2]=".2"
array[3]=".3"
array[4]=".4"
array[5]=".5"
array[6]=".6"
array[7]=".7"
array[8]=".8"
array[9]=".9"
array[10]=1
array[11]=2
array[12]=3
array[13]=4
array[14]=5
array[15]=6
array[16]=7
array[17]=8
array[18]=9
array[19]=10

  for ATP in 1000 100
    do

    for SA in SAca SAnfat SAptxs SAptxf SAcm SAnancx SAcn SApp38
      do

        for  (( i=0; i<=19; i++))
          do

          python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 300000 -iters 1 -var stim_amplitude $ATP -var $SA ${array[$i]} -var stim_period 1e3 -var IC1 0.2e3 -var IC2 0.2e3 -var stim_gap1 0.85e3 -var stim_gap2 0.85e3 -var stim_low 0e3 -var stim_high 0.25e3 -name ~/Data_storage/SA-$SA-${array[$i]}-$ATP

          done

      done

    done

fi

if [ $sens2 -eq 1 ]
then

  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 300000 -iters 1 -var stim_amplitude 1000 -var SAserca 100 -var stim_period 0.215e3 -var IC1 0.2e3 -var IC2 0.14e3 -var stim_gap1 0.02e3 -var stim_gap2 0.02e3 -var stim_low 0.15e3 -var stim_high 0.25e3 -name ~/Data_storage/SA-SAserca-100-1000
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 300000 -iters 1 -var stim_amplitude 1000 -var SAserca 1 -var stim_period 0.215e3 -var IC1 0.2e3 -var IC2 0.14e3 -var stim_gap1 0.02e3 -var stim_gap2 0.02e3 -var stim_low 0.15e3 -var stim_high 0.25e3 -name ~/Data_storage/SA-SAserca-1-1000
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 300000 -iters 1 -var stim_amplitude 1000 -var SAserca 0.01 -var stim_period 0.215e3 -var IC1 0.2e3 -var IC2 0.14e3 -var stim_gap1 0.02e3 -var stim_gap2 0.02e3 -var stim_low 0.15e3 -var stim_high 0.25e3 -name ~/Data_storage/SA-SAserca-0.01-1000

fi

### SENSITIVITY ANALYSIS ###
if [ $sens3 -eq 1 ]
then
array[0]=".01"
array[1]=".1"
array[2]=".2"
array[3]=".3"
array[4]=".4"
array[5]=".5"
array[6]=".6"
array[7]=".7"
array[8]=".8"
array[9]=".9"
array[10]=1
array[11]=2
array[12]=3
array[13]=4
array[14]=5
array[15]=6
array[16]=7
array[17]=8
array[18]=9
array[19]=10

  for ATP in 1000 100
    do

    for SA in SApp38
      do

        for  (( i=0; i<=19; i++))
          do

          python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 100000 -iters 1 -var stim_amplitude $ATP -var $SA ${array[$i]} -var stim_period 2e3 -var stim_gap1 1e3 -var stim_gap2 1e3 -var stim_low 1e3 -var stim_high 2e3 -name ~/Data_storage/SA-$SA-${array[$i]}-$ATP

          done

      done

    done

fi

### SENSITIVITY ANALYSIS ###
if [ $sens4 -eq 1 ]
then

  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 300000 -iters 1 -var stim_amplitude 1000 -var SAserca 0.0001 -var stim_period 0.215e3 -var IC1 0.2e3 -var IC2 0.14e3 -var stim_gap1 0.02e3 -var stim_gap2 0.02e3 -var stim_low 0.15e3 -var stim_high 0.25e3 -name ~/Data_storage/SA-SAserca-0.0001-1000
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 300000 -iters 1 -var stim_amplitude 1000 -var SAserca 1 -var stim_period 0.215e3 -var IC1 0.2e3 -var IC2 0.14e3 -var stim_gap1 0.02e3 -var stim_gap2 0.02e3 -var stim_low 0.15e3 -var stim_high 0.25e3 -name ~/Data_storage/SA-SAserca-1-1000
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 300000 -iters 1 -var stim_amplitude 1000 -var SAserca 10000 -var stim_period 0.215e3 -var IC1 0.2e3 -var IC2 0.14e3 -var stim_gap1 0.02e3 -var stim_gap2 0.02e3 -var stim_low 0.15e3 -var stim_high 0.25e3 -name ~/Data_storage/SA-SAserca-10000-1000

fi

### LENGTHY Stimulation ###

if [ $leng1 -eq 1 ]
then
  echo "NFAT calculation is ON"

  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 1000000 -iters 1 -var pulse_switch 0 -var stim_amplitude 0 -name ~/Data_storage/15min_MG_0mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 1000000 -iters 1 -var pulse_switch 0 -var stim_amplitude 100 -name ~/Data_storage/15min_MG_01mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 1000000 -iters 1 -var pulse_switch 0 -var stim_amplitude 2000 -name ~/Data_storage/15min_MG_2mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 1000000 -iters 1 -var pulse_switch 0 -var stim_amplitude 1000 -name ~/Data_storage/15min_MG_1mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 1000000 -iters 1 -var pulse_switch 0 -var stim_amplitude 3000 -name ~/Data_storage/15min_MG_3mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 1000000 -iters 1 -var pulse_switch 0 -var stim_amplitude 5000 -name ~/Data_storage/15min_MG_5mMATP

  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var pulse_switch 0 -var stim_amplitude 3000 -name ~/Data_storage/60min_MG_1mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 3600000 -iters 1 -var pulse_switch 0 -var stim_amplitude 3000 -name ~/Data_storage/60min_MG_3mMATP
fi

if [ $leng2 -eq 1 ]
then
  echo "pp38 calculation is ON"
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 1000000 -iters 1 -var pulse_switch 0 -var stim_amplitude 1000 -name ~/Data_storage/15min_MG_1mMATP
fi

if [ $leng3 -eq 1 ]
then
  echo "TNFa calculation is ON"
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 10800000 -iters 1 -var pulse_switch 0 -var stim_amplitude 10 -name ~/Data_storage/180min_MG_001mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 10800000 -iters 1 -var pulse_switch 0 -var stim_amplitude 100 -name ~/Data_storage/180min_MG_01mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 10800000 -iters 1 -var pulse_switch 0 -var stim_amplitude 1000 -name ~/Data_storage/180min_MG_1mMATP &
  python2.7 daisychain.py -dt 1 -dSr 1000 -jit -odeName $ODEFILEfull -T 10800000 -iters 1 -var pulse_switch 0 -var stim_amplitude 3000 -name ~/Data_storage/180min_MG_3mMATP

fi

echo "Calculation is done"
