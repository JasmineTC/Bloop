#!/bin/bash

# Fine tests providing partial coverage
python3 ../src/ThreeHiggs/UnitTests.py

# Super coarse test providing full coverage
python3 runBenchmark.py -l 2 -n 1 --save --resultsDirectory TestResults --TRangeStart 50 --TRangeEnd 100 --TRangeStepSize 10 --DiagAlgo scipy --minimizationAlgo directGlobal --v1Bounds 1e-6 1e-6 --v2Bounds 1e-6 1e-6 --v3Bounds 1e-6 100
diff TestResults/*.txt ReferenceTestResults/*.txt

