#!/bin/bash

# Fine tests providing partial coverage
python3 ../src/ThreeHiggs/UnitTests.py

# Super coarse test providing full coverage
python3 runBenchmark.py -l 2 \
                        --firstBenchmark 1 \
                        --lastBenchmark 1 \
                        --save \
                        --resultsDirectory TestResults \
                        --TRangeStart 50 \
                        --TRangeEnd 100 \
                        --TRangeStepSize 10 \
                        --DiagAlgo scipy \
                        --minimizationAlgo directGlobal \
                        --v1Bounds 1e-6 1e-6 \
                        --v2Bounds 1e-6 1e-6 \
                        --v3Bounds 1e-6 100 \
                        --absGlobalTolerance 0.1 \
                        --relGlobalTolerance 0.1 \
                        --absLocalTolerance 0.00001 \
                        --relLocalTolerance 0.00001 \
                        --bAbsMass

diff TestResults/BM_1.json ReferenceTestResults/BM_1.json

# Super coarse test providing full coverage, permit complex mass
python3 runBenchmark.py -l 2 \
                        --firstBenchmark 1 \
                        --lastBenchmark 1 \
                        --save \
                        --resultsDirectory TestComplexMassResults \
                        --TRangeStart 50 \
                        --TRangeEnd 100 \
                        --TRangeStepSize 10 \
                        --DiagAlgo scipy \
                        --minimizationAlgo directGlobal \
                        --v1Bounds 1e-6 1e-6 \
                        --v2Bounds 1e-6 1e-6 \
                        --v3Bounds 1e-6 100 \
                        --absGlobalTolerance 0.1 \
                        --relGlobalTolerance 0.1 \
                        --absLocalTolerance 0.00001 \
                        --relLocalTolerance 0.00001

diff TestComplexMassResults/BM_1.json ReferenceTestComplexMassResults/BM_1.json


