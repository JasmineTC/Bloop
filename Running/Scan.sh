#!/bin/bash
echo ../src/ThreeHiggs/BmGenerator.py --mode random --randNum 1000000

echo time OMP_NUM_THREADS=1 python3 runBenchmark.py --loopOrder 1 \
	--benchMarkFile Benchmarks/randomScan.json \
	--bPool \
	--bNumba \
	--cores 46 \
	--bProcessMin \
	--TRangeEnd 400 \
	--TRangeStepSize 2

echo python3 SFOPTBM.py

echo python3 ../src/ThreeHiggs/BmGenerator.py --mode randomSSS

time OMP_NUM_THREADS=1 python3 runBenchmark.py --loopOrder 1 \
	--benchmarkFile Benchmarks/randomScanSSS.json \
	--bPool \
	--bNumba \
	--cores 46 \
	--bProcessMin \
	--TRangeEnd 400 \
	--TRangeStepSize 0.1 \
	--resultsDirectory ResultsSSS
