#!/bin/bash

python3 runBenchmark.py -l 2 -n 1 --save --resultsDirectory TestResults --TRangeStart 50 --TRangeEnd 100 --TRangeStepSize 10

diff TestResults/*.txt ReferenceTestResults/*.txt

