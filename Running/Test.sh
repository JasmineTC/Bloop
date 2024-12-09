#!/bin/bash

# Fine tests providing partial coverage
echo Unit tests...
python3 ../src/ThreeHiggs/UnitTests.py

echo Intergration tests...

# Super coarse test providing full coverage
echo Intergration test: Running code at NNLO, abs mass mode on...
rm IntegrationTests/AbsMass/OutputResult/* 
python3 runBenchmark.py --loopOrder 2 \
                        --firstBenchmark 1 \
                        --lastBenchmark 1 \
                        --bSave \
                        --resultsDirectory IntegrationTests/AbsMass/OutputResult/  \
                        --TRangeStart 50 \
                        --TRangeEnd 100 \
                        --TRangeStepSize 10 \
                        --minimizationAlgo combo \
                        --v1Bounds 1e-6 1e-6 \
                        --v2Bounds 1e-6 1e-6 \
                        --v3Bounds 1e-6 60 \
                        --absGlobalTolerance 0.1 \
                        --relGlobalTolerance 0.1 \
                        --absLocalTolerance 0.00001 \
                        --relLocalTolerance 0.00001 \
                        --bAbsMass \

diff IntegrationTests/AbsMass/OutputResult/BM_1.json IntegrationTests/AbsMass/ReferenceResult/BM_1.json
# Super coarse test providing full coverage, permit complex mass
echo Intergration test: Running code at NNLO, complex mass mode on...
rm IntegrationTests/ComplexMass/OutputResult/* 
python3 runBenchmark.py --loopOrder 2 \
                        --firstBenchmark 1 \
                        --lastBenchmark 1 \
                        --bSave \
                        --resultsDirectory IntegrationTests/ComplexMass/OutputResult/  \
                        --TRangeStart 50 \
                        --TRangeEnd 100 \
                        --TRangeStepSize 10 \
                        --minimizationAlgo combo \
                        --v1Bounds 1e-6 1e-6 \
                        --v2Bounds 1e-6 1e-6 \
                        --v3Bounds 1e-6 60 \
                        --absGlobalTolerance 0.1 \
                        --relGlobalTolerance 0.1 \
                        --absLocalTolerance 0.00001 \
                        --relLocalTolerance 0.00001 \

diff IntegrationTests/ComplexMass/OutputResult/BM_1.json IntegrationTests/ComplexMass/ReferenceResult/BM_1.json

# Super coarse test providing full coverage, complex mass+numba
echo Intergration test: Running code at NNLO, complex mass mode on with Numba...
rm IntegrationTests/Numba/OutputResult/* 
python3 runBenchmark.py --loopOrder 2 \
                        --firstBenchmark 1 \
                        --lastBenchmark 1 \
                        --bSave \
                        --resultsDirectory IntegrationTests/Numba/OutputResult/ \
                        --TRangeStart 50 \
                        --TRangeEnd 100 \
                        --TRangeStepSize 10 \
                        --bNumba \
                        --minimizationAlgo combo \
                        --v1Bounds 1e-6 1e-6 \
                        --v2Bounds 1e-6 1e-6 \
                        --v3Bounds 1e-6 60 \
                        --absGlobalTolerance 0.1 \
                        --relGlobalTolerance 0.1 \
                        --absLocalTolerance 0.00001 \
                        --relLocalTolerance 0.00001 \

diff IntegrationTests/Numba/OutputResult/BM_1.json IntegrationTests/Numba/ReferenceResult/BM_1.json
