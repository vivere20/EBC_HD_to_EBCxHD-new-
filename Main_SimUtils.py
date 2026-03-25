import numpy as np
import scipy.stats as spt
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

pi=np.pi

def vonmises_sim(covariates, Params):
    covariates = np.mod(covariates, 2*np.pi)
                    #if angle is greater than pi, subtract 2 pi
                    #to give us (-pi pi] range without flipping zero
    this_pdf = spt.vonmises.pdf(covariates, Params[1],loc=Params[0])
    maxval = spt.vonmises.pdf(Params[0],Params[1],loc=Params[0]) #for scaling purposes.  Max value is when exactly at means for gauss
    lambdas = this_pdf/maxval #parameter to put into the poisson
    
    return lambdas

def place_sim(covariates,Params):
    covariates = np.asarray(covariates).T
    mu = Params[0]
    cov = np.identity(2) * Params[1] #put in variance directly
    this_pdf = spt.multivariate_normal.pdf(covariates,mu,cov)
    maxval = spt.multivariate_normal.pdf(mu,mu,cov)
    lambdas = this_pdf/maxval
    return lambdas

def spiking(lambdas,fs,celltype):
        #allows for different peak firing rates
        #to help out of rf firing when lambda is 0 going to restrict things a bit more
        
        #remember that maxfr*fs is true firing rate since fs = 30 dt = 0.033 anything about 3 spikes in a bin is excessive
    
    if celltype == 'EBC':
            maxfr = 2 #Expected spikes per time step is 2/0.033 = 60 Hz
            lambdas *=maxfr
            spikes = np.random.poisson(lambdas)
  
    if celltype == 'ABC':
            maxfr = 2
            lambdas *=maxfr
            spikes = np.random.poisson(lambdas)

            
    if celltype == 'HD':
            
            maxfr = 2
            lambdas *=maxfr
            spikes = np.random.poisson(lambdas)

    if celltype == 'Speed':
            #2024-02-20 finally have a "fix" based on some quick testing
            
            maxfr = 1 #this is handled in the parametric rf itself
            lambdas[lambdas<0] = 0
            lambdas *=maxfr
            #lambdas = np.exp(lambdas) #do the link function properly here
            spikes = np.round(np.random.poisson(lambdas)*(1/fs) )#sim firing rates directly then get emperical spike counts in each bin by scaling by bin duration
            #only way to have proper firing rate and spike pattern (emperical spike counts only though, after round it will be back to step issue)
            #but note need to skip smooth fr for speed cells and put the firing rate in directly as smoothing obscures pattern

    if celltype == 'Place':
            
            maxfr = 2
            lambdas *=maxfr
            spikes = np.random.poisson(lambdas)

    return spikes

def rf_simactivity(fs,celltypes, nper, MasterCovlist, envlist, df_columns,n_sims):
    "Set up cell ID and params for each inputted cell type"            
 
    Params = dict()
    cell_id_list = np.array([0])
    celltypelist = []

    for c in celltypes: 
    
        if c =='None': 
            cell_ids = np.arange(1,nper*2) + np.max(cell_id_list)
            cell_id_list = np.concatenate([cell_id_list,cell_ids])
            celltypelist = celltypelist + [c] * nper*2
        
        if c == 'EBC':
            # von mises
            centers = np.hstack([np.linspace(0,np.pi*2,nper, endpoint=False)])
            widths = np.hstack([np.random.uniform(10,11,nper)])
            cell_ids = np.arange(1,nper+1) + np.max(cell_id_list)
            cell_id_list = np.concatenate([cell_id_list,cell_ids])
            celltypelist = celltypelist + [c] * nper
            together = np.reshape(np.concatenate([centers,widths]),(2,nper))
            Params['EBC'] = together
            
        if c == 'HD':
            centers = np.hstack([np.linspace(0,np.pi*2,nper, endpoint=False)])
            widths = np.hstack([np.random.uniform(6,7,nper)])
            cell_ids = np.arange(1,nper+1) + np.max(cell_id_list)
            cell_id_list = np.concatenate([cell_id_list,cell_ids])
            celltypelist = celltypelist + [c] * nper
            together = np.reshape(np.concatenate([centers,widths]),(2,nper))
            Params['HD'] = together
            
        if c == 'Place':
            #Note: can't do centers here because they have to change for each environment
            #Can do widths though
            widths = np.hstack([np.random.uniform(30,50,nper)])
            cell_ids = np.arange(1,nper+1) + np.max(cell_id_list)
            cell_id_list = np.concatenate([cell_id_list,cell_ids])
            celltypelist = celltypelist + [c] * nper
            Params['Place'] = widths
        

    row_list = [] #set up a row list of dictionaries with key that are column names for each row.  Then push into df_data
    PlaceCenters = [] #need centers of place cells for later sims
    for cn in cell_id_list:
        print('Simulating Cell #: ' + str(cn))
        envdict = dict() 
        envdict['Cell ID'] = cn
        envdict['Cell Type'] = celltypelist[cn]
        envdict['Sim Protocol'] = 'Parametric RF'
        for s in range(0,n_sims):
            print('Running Session #: ' + str(s))
            envdict['Sim Round'] = s
            Cov_session = MasterCovlist[0][s]
            Inds_session = MasterCovlist[1][s] #indicies for when the animal is within set range of the wall
                   
            for env in envlist:
                print('Current environment: '+ env)
                envdict['Environment'] = env
                #Get covariates and simulate activity according to cell type
                
                if celltypelist[cn] == 'None':
                    lambdas = np.random.uniform(0,1,len(Cov_session['HD'][env]))
                    spikes = spiking(lambdas,fs,'EBC')
                
                if celltypelist[cn] == 'EBC':
                    inds = Inds_session[env]
                    covariates = Cov_session['EBC'][env][inds]
                    templambdas = vonmises_sim(covariates, Params['EBC'][:,np.mod(cn,nper)])
                    lambdas = np.zeros(len(Cov_session['HD'][env])) #have to do this to get lambdas for whole session rather than just when by the wall
                    lambdas[inds] = np.squeeze(templambdas)
                    spikes = spiking(lambdas,fs,'EBC')
                   
                if celltypelist[cn] == 'HD':
                    covariates = Cov_session['HD'][env].values #this is right from the dataframe so have to do values to call correct parts of array
                    lambdas = vonmises_sim(covariates, Params['HD'][:,np.mod(cn,nper)])
                    spikes = spiking(lambdas,fs,'HD')
                 
                if celltypelist[cn] == 'Place':
                    covariates = Cov_session['Place'][env]
                    widths = Params['Place'][np.mod(cn,nper)] #only need the one width for each cell
                    #limit centers by where the animal actually went so that each cell has at least some spiking
                    x_low = np.min(Cov_session['Place'][env][0])
                    x_up = np.max(Cov_session['Place'][env][0])
                    y_low = np.min(Cov_session['Place'][env][1])
                    y_up = np.max(Cov_session['Place'][env][1])
                    centers =  np.concatenate([np.random.uniform(x_low,x_up,1),np.random.uniform(y_low,y_up,1)]) #only need the one center
                    PlaceCenters.append(centers)
                    place_params = [centers,widths] #for each environment 
                    lambdas = place_sim(covariates, place_params)
                    spikes = spiking(lambdas,fs,'Place')
                
                envdict['Spikes'] = spikes
                row_list.append(envdict.copy())
        
    df_data = pd.DataFrame(row_list)
                
    return df_data, PlaceCenters, Params