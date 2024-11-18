#!/bin/bash

# Fine tests providing partial coverage
echo Unit tests...
python3 ../src/ThreeHiggs/UnitTests.py

echo Intergration tests...

# Super coarse test providing full coverage
echo Intergration test: Running code at NNLO, abs mass mode on...
#rm TestResults/*
python3 runBenchmark.py -l 2 \
                        --firstBenchmark 1 \
                        --lastBenchmark 1 \
                        --bSave \
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
                        --bAbsMass \
                        --bPool

diff TestResults/BM_1.json ReferenceTestResults/BM_1.json

# Super coarse test providing full coverage, permit complex mass
echo Intergration test: Running code at NNLO, complex mass mode on...
#rm TestComplexMassResults/*
python3 runBenchmark.py -l 2 \
                        --firstBenchmark 1 \
                        --lastBenchmark 1 \
                        --bsave \
                        --resultsDirectory TestComplexMassResults \
                        --TRangeStart 50 \
                        --TRangeEnd 100 \
                        --TRangeStepSize 10 \
                        --diagAlgo scipy \
                        --minimizationAlgo directGlobal \
                        --v1Bounds 1e-6 1e-6 \
                        --v2Bounds 1e-6 1e-6 \
                        --v3Bounds 1e-6 100 \
                        --absGlobalTolerance 0.1 \
                        --relGlobalTolerance 0.1 \
                        --absLocalTolerance 0.00001 \
                        --relLocalTolerance 0.00001 \
                        --bPool

diff TestComplexMassResults/BM_1.json ReferenceTestComplexMassResults/BM_1.json

