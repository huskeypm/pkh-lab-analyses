#!/bin/bash
#Set up environment variables needed to locate the code and the data.
#This file may need modification when running on previously untested configurations.
export DATAFOLDER=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
echo DATAFOLDER=$DATAFOLDER
