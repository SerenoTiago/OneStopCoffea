#!/usr/bin/env bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc14-opt/setup.sh
set -eux 
mkdir -p .application_data/bin
g++ -O3 -std=c++20 $(root-config --cflags --libs) \
  scripts/signal_splitting/split_signals.cpp \
  -o .application_data/bin/split_signals
.application_data/bin/split_signals "$1" "$2"
xrdcp -r analysis_inputs/signals/split/ root://cmseos.fnal.gov///store/user/ckapsiak/SingleStop/raw_official_samples/
echo "JOB COMPLETED SUCCESSFULL"
