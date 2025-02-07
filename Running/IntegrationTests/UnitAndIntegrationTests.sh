#!/bin/bash

# Fine tests providing partial coverage
echo Unit tests...
python3 ../src/ThreeHiggs/UnitTests.py

echo Intergration tests...
# Super coarse test providing full coverage
echo Running code at NLO with pool using 2 cores...
rm -f IntegrationTests/Pool/OutputResult/* 
python3 runStages.py --loopOrder 1 \
                        --firstBenchmark 0 \
                        --lastBenchmark 3 \
                        --bSave \
                        --resultsDirectory IntegrationTests/Pool/OutputResult/  \
                        --TRangeStart 100 \
                        --TRangeEnd 200 \
                        --TRangeStepSize 2 \
                        --bPool \
                        --cores 2
diff IntegrationTests/Pool/OutputResult/BM_0.json IntegrationTests/Pool/ReferenceResult/BM_0.json
diff IntegrationTests/Pool/OutputResult/BM_1.json IntegrationTests/Pool/ReferenceResult/BM_1.json
diff IntegrationTests/Pool/OutputResult/BM_2.json IntegrationTests/Pool/ReferenceResult/BM_2.json
diff IntegrationTests/Pool/OutputResult/BM_3.json IntegrationTests/Pool/ReferenceResult/BM_3.json

echo Running code at NNLO...
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


echo Running code at NNLO with Numba...
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

