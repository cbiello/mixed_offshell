#!/bin/bash
# one line of jobs.txt: "nsc <analysis> <boson> <corr> <ch> <sec> <seed>" or "pwg <order> <boson> <scale> <seed>"
k=$1; shift
/scratch2/cbiello/mixed_runs/run_$k.sh "$@"
