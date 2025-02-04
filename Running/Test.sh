#!/bin/bash

# Fine tests providing partial coverage
echo Unit tests...
python3 ../src/ThreeHiggs/UnitTests.py

echo Intergration tests...

# Super coarse test providing full coverage
echo Intergration test: Running code at NNLO, complex mass mode on...
rm -f IntegrationTests/ComplexMass/OutputResult/* 
python3 runStages.py --loopOrder 2 \
                        --firstBenchmark 1 \
                        --lastBenchmark 1 \
                        --bSave \
                        --resultsDirectory IntegrationTests/ComplexMass/OutputResult/  \
                        --TRangeStart 50 \
                        --TRangeEnd 100 \
                        --TRangeStepSize 10 \
                        --bProcessMin
diff IntegrationTests/ComplexMass/OutputResult/BM_1.json IntegrationTests/ComplexMass/ReferenceResult/BM_1.json
diff IntegrationTests/ComplexMass/OutputResult/BM_1_interp.json IntegrationTests/ComplexMass/ReferenceResult/BM_1_interp.json

# Super coarse test providing full coverage, complex mass+numba
echo Intergration test: Running code at NNLO, complex mass mode on with Numba...
rm -f IntegrationTests/Numba/OutputResult/* 
python3 runStages.py --loopOrder 2 \
                        --firstBenchmark 1 \
                        --lastBenchmark 1 \
                        --bSave \
                        --resultsDirectory IntegrationTests/Numba/OutputResult/ \
                        --TRangeStart 50 \
                        --TRangeEnd 100 \
                        --TRangeStepSize 10 \
                        --bNumba \
                        --bProcessMin
diff IntegrationTests/Numba/OutputResult/BM_1.json IntegrationTests/Numba/ReferenceResult/BM_1.json
diff IntegrationTests/Numba/OutputResult/BM_1_interp.json IntegrationTests/Numba/ReferenceResult/BM_1_interp.json

