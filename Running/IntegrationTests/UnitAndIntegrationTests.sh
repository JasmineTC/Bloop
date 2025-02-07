#!/bin/bash

# Fine tests providing partial coverage
echo Unit tests...
python3 ../../src/ThreeHiggs/UnitTests.py

echo Intergration tests...
# Super coarse test providing full coverage
echo Running code at NLO with pool using 2 cores...
rm -f Pool/OutputResult/* 
python3 runStages.py --loopOrder 1 \
                        --firstBenchmark 0 \
                        --lastBenchmark 3 \
                        --bSave \
                        --resultsDirectory Pool/OutputResult/  \
                        --TRangeStart 100 \
                        --TRangeEnd 200 \
                        --TRangeStepSize 2 \
                        --bPool \
                        --cores 2
diff Pool/OutputResult/BM_0.json Pool/ReferenceResult/BM_0.json
diff Pool/OutputResult/BM_1.json Pool/ReferenceResult/BM_1.json
diff Pool/OutputResult/BM_2.json Pool/ReferenceResult/BM_2.json
diff Pool/OutputResult/BM_3.json Pool/ReferenceResult/BM_3.json

echo Intergration test: Running code at NNLO...
rm -f ComplexMass/OutputResult/* 
python3 runStages.py --loopOrder 2 \
                        --firstBenchmark 1 \
                        --lastBenchmark 1 \
                        --bSave \
                        --resultsDirectory ComplexMass/OutputResult/  \
                        --TRangeStart 50 \
                        --TRangeEnd 100 \
                        --TRangeStepSize 10 \
                        --bProcessMin
diff ComplexMass/OutputResult/BM_1.json ComplexMass/ReferenceResult/BM_1.json
diff ComplexMass/OutputResult/BM_1_interp.json ComplexMass/ReferenceResult/BM_1_interp.json


echo Intergration test: Running code at NNLO with Numba...
rm -f Numba/OutputResult/* 
python3 runStages.py --loopOrder 2 \
                        --firstBenchmark 1 \
                        --lastBenchmark 1 \
                        --bSave \
                        --resultsDirectory Numba/OutputResult/ \
                        --TRangeStart 50 \
                        --TRangeEnd 100 \
                        --TRangeStepSize 10 \
                        --bNumba \
                        --bProcessMini
diff Numba/OutputResult/BM_1.json Numba/ReferenceResult/BM_1.json
diff Numba/OutputResult/BM_1_interp.json Numba/ReferenceResult/BM_1_interp.json

