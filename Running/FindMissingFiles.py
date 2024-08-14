import os
from glob import glob ## I think glob lets you do the * thingy
from json import load

count = 0
for filename in glob(os.path.join('Running/Results/ScanResult', '*.json')):
    count += 1
    if count >10:
        break
    strongPT = False
    with open(filename) as json_file:
        resultDic = load(json_file)
        if resultDic['jumpsv1']:
            for jump in resultDic['jumpsv1']:
                if abs(jump[0]) > 0.7:
                    strongPT = True
                    
        if resultDic['jumpsv2']:
            for jump in resultDic['jumpsv2']:
                if abs(jump[0]) > 0.7:
                    strongPT = True
        
        if resultDic['jumpsv3']:
            for jump in resultDic['jumpsv3']:
                if abs(jump[0]) > 0.7:
                    strongPT = True
    if strongPT:
        print(filename)