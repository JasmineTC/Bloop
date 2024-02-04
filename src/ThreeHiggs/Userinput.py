import argparse
import os
import sys
import multiprocessing


class Userinput(argparse.ArgumentParser):
    
    def __init__(self):
        
        super().__init__()
        
        self.add_argument('-b', action = 'append', default = ['Benchmarks/Benchmarks_3HDM.py'], dest = 'bench_marks')
        self.add_argument('-l', action = 'store', default = 1, dest = 'loop_order', type = int, choices = [1, 2])
        self.add_argument('-c', action = 'store', default = 1, dest = 'cores', type = int, choices = list(range(1, multiprocessing.cpu_count() + 1)))
    
        
    def parse(self):
        
        args = super().parse_args()
        
        for benchmarkFile in args.bench_marks:
            if not os.path.isfile(benchmarkFile):
                print(f"File {benchmarkFile} could not be found, gg.")
                sys.exit(-1)
            
        return args