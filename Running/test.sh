#!/bin/bash

# Fine tests providing partial coverage
python3 ../src/ThreeHiggs/UnitTests.py

# Super coarse test providing full coverage
python3 runBenchmark.py -l 2 -n 1 --save --resultsDirectory TestResults --TRangeStart 50 --TRangeEnd 100 --TRangeStepSize 10 --DiagAlgo scipy --minimizationAlgo directGlobal
diff TestResults/*.txt ReferenceTestResults/*.txt

