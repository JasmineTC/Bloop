'''---Makes a dir of just the SFOPT results for faster loading when doing scatter plots or something ---
Loads all the file names in the results directory 
Strips the Results/ part of the file name
Loads the benchmark information from the file
If the benchmark is strong save the file to a subsetresult dir'''
from os.path import join
from glob import glob ## I think glob lets you do the * thingy
from json import load, dump

for filePointer in glob(join('Results', '*.json')):
    fileName = filePointer.split('/')[1].lstrip().split(' ')[0]
    with open(filePointer, "r") as f:
        resultDic = load(f)
        if resultDic["strong"]:
            dump(resultDic,open(f"Results/StrongSubSetResults/{fileName}", "w"))

