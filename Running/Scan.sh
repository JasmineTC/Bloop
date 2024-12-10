#!/bin/bash
echo ../src/ThreeHiggs/BmGenerator.py --mode random --randNum 1000000

echo time OMP_NUM_THREADS=1 python3 runBenchmark.py --loopOrder 1 \
	--benchMarkFile Benchmarks/randomScan.json \
	--cores 46 \
	--bProcessMin \
	--TRangeEnd 400 \
	--TRangeStepSize 2

python3 SFOPTBM.py

python3 ../src/ThreeHiggs/BmGenerator.py --mode randomSSS

time OMP_NUM_THREADS=1 python3 runBenchmark.py --loopOrder 1 \
	--benchmarkFile Benchmarks/randomScanSSS.json \
	--cores 46 \
	--bProcessMin \
	--TRangeEnd 400 \
	--TRangeStepSize 0.1 \
	--resultsDirectory ResultsSSS
