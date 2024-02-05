import argparse
##May not be needed depending on how this is done in parallel
import multiprocessing


class UserInput(argparse.ArgumentParser):
    
    def __init__(self):
        
        super().__init__()
        ##Takes user arguement to define what loop order to calculate the effective potential to
        ##TODO make safe by checking if int given actually has a bm point
        ##something like if benchMarkNumber > len(bm list) then exit
        self.add_argument('-n', '--benchMarkNumber', action = 'store', default = 0, dest = 'benchMarkNumber', type = int,
                          help = "Used to specify a particular bench mark point in the list to run")
        ##Takes user arguement to define what loop order to calculate the effective potential to
        self.add_argument('-l', 'loopOrder', action = 'store', default = 1, dest = 'loopOrder', type = int, choices = [1, 2],
                          help = "Used to specify if the effective potential should be calculated to one or two loop")
        ##Takes user bool to decide if plots should be made after saving results
        self.add_argument('-p', 'plot', action = 'store_true', default=False, dest = 'plot',  
                          help = "Used to specify if a plot of minimium vs temp should be made")
        ##Takes user arguement to define how many cores to run the benchmarks on
        ##multiprocessing.cpu_count gets from the the system how many cores are avaviable, +1 needed because of how range works
        self.add_argument('-c', '--cores', action = 'store', default = 1, dest = 'cores', type = int, choices = list(range(1, multiprocessing.cpu_count() + 1)),
                          help = "Used to specify how many cores to run the bench mark list on")
        
    ##Used to check userinputs are valid, mostly done with the choice keyword above now though
    def parse(self):
        
        args = super().parse_args()
        
        return args 