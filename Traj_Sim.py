import numpy as np
import importlib
import pandas as pd
import T_SimUtilities as Utils
importlib.reload(Utils)


def runTrajSimulation(dur, use_trajstats, use_envgeo, policy):
    
    columns = ['Envtype','X', 'Y','HeadingDir',
               'AlloDis','EgoDis', 'EgoAngle','AlloAngle','ClosestPoint']
    thedf = pd.DataFrame(columns = columns)
    
    fs,b,mu,sigma = Utils.set_trajstats()

    numenv = 1 
    env_types, boundary_points = Utils.make_envirs(numenv)
        
    simx = dict()
    simy = dict()
    headingdir = dict()
    allodis = dict()
    egodis = dict()
    egoangle = dict()
    posangl =dict()
    closestpoint=dict()
     
    #Simulate trajectory   
    for curenv in env_types:   
             
        if policy == 'WallWalk':
           
           simx[curenv],simy[curenv], headingdir[curenv] = Utils.WallWalk(dur,fs,b,mu,sigma,boundary_points,curenv)

        else:
            
            print('Running Policy not known')
            
        allodis[curenv], egodis[curenv], egoangle[curenv], posangl[curenv], closestpoint[curenv]= Utils.analyzetraj(simx[curenv],simy[curenv],headingdir[curenv], boundary_points[curenv])                    
        tempdf = pd.DataFrame({

            'Envtype': curenv,
            'X' : simx[curenv],
            'Y': simy[curenv],
            'HeadingDir': headingdir[curenv],
            'AlloDis':list(allodis[curenv]),
            'EgoDis':list(egodis[curenv]),
            'EgoAngle':list(egoangle[curenv]),
            'AlloAngle':list(posangl[curenv]),
            'ClosestPoint':list(closestpoint[curenv])
            }
                                )
         
        thedf = pd.concat([thedf,tempdf])
        
    
    return thedf, fs, boundary_points 