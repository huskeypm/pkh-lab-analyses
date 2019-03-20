
#export LOC=$HOME/sources/            
export LOC=/home/AD/bdst227/sources/            
export LOC=/home/AD/pmke226/sources/            
export MYPATH=$LOC/mypython/
export PYTHONPATH=$PYTHONPATH:$MYPATH/lib/python2.7/site-packages/
export PATH=$PATH:$LOC/gotran/scripts/
python -c "import instant"
python -c "import modelparameters"
python -c "import gotran"

