#!/bin/bash
#set up environment variables needed to locate the code and the data
repo=nanopore_diffusion
export DATALOC=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
echo DATALOC=$DATALOC
export SRCLOC=$(cd "$DATALOC" && cd ../$repo/src && pwd )
echo SRCLOC=$SRCLOC
