#!/usr/bin/env python3
import sys, argparse, os
import pandas as pd
from sklearn import preprocessing
from lib import combine_methods as cm
#from lib import tm_score as tm
import multiprocessing
from multiprocessing.dummy import Pool as ThreadPool
import warnings
warnings.filterwarnings("ignore")
thread=999
# Change this to false to switch off MultiProcessing. Can be useful for debugging.
MULTI = True
def run_meetdock(mypath, myoutdir, recChain = True, ligChain = True):#, electro = True, jones = True):#, thread,  depth, shape = True, electro = True, jones = True, recChain = True, ligChain = True, proba = True, foldx = True, pH = 7, dist = 8.6):
    pdbs = os.listdir(mypath)
    ouput=open('score.txt','w')
    pdbs = [i for i in pdbs if i.endswith(".pdb")]
    all_res = []
    execdir = os.path.dirname(os.path.realpath(__file__))
    if MULTI == True:
        try:
            def runpdbs(apdb):
                thepdb = mypath+"/"+apdb
                print("Computing scores for ", apdb)
                print("for combine_score////////////////")
                res = cm.combine_score(thepdb, recepChain=recChain, ligChain=ligChain,  vdwrun=jones, electrorun=electro)#, depth=depth, dist=dist)
                all_res.append(res)

            if thread == 999:
                num_cores = multiprocessing.cpu_count()

            pool = ThreadPool(num_cores)

            print("MultiProcessing available\nRunning on",num_cores, "threads")
            pool.map(runpdbs, pdbs)

        except:
            print("...")
    else :
        print("MultiProcessing turned off, running on a single thread.")
        for apdb in pdbs:
            thepdb = mypath+"/"+apdb
            print("Computing scores for ", apdb)
            res = cm.combine_score(thepdb, recepChain=recChain, ligChain=ligChain, statpotrun=proba, vdwrun=jones, electrorun=electro, shaperun=shape, pH = pH, depth=depth, dist=dist)
            all_res.append(res)
            #ouput.write(str(res))

    mydf = pd.DataFrame(all_res)
    mydf = mydf.set_index('pdb')
    print(mydf.to_string())
    ouput.write(mydf.to_string())
    #tm.tm_score(mydf, outdir = myoutdir, execdir = execdir)
    #tr=open('tmscore.txt','w')
    #tr.write(tt)



#if __name__ == "__main__":
#    mypath, myoutdir, shape, electro, jones, recChain, ligChain, proba, pH, depth, dist, thread = get_args()
#    ff=run_meetdock(mypath=mypath, myoutdir=myoutdir, ligChain=ligChain, recChain=recChain, shape = shape, electro = electro, jones =jones, proba =proba, pH =pH, dist =dist, depth =depth, thread=thread)
    #print('last')
    #for ii in ff:
        #print(ii)

    

    #################the rest is added###############
    

