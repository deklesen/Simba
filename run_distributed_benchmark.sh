#!/bin/sh
echo Starting benchmark with $1 processes...
for i in $(seq $1); do python benchmark.py $i $1 & done