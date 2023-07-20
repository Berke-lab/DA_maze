"""Assortment of functions to help with Krausz et al. Maze analyses."""

__author__ = "Tim Krausz"
__email__ = "krausz.tim@gmail.com"
__status__ = "development"
__license__ = "MIT"

import math
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
scale = StandardScaler()
from statsmodels.discrete.discrete_model import Logit

fs = 250

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) \
        < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

def closestlong(array,value):
    return(array[np.abs(np.subtract(array,value)).argmin()])

def project_onto_path(posmat,path):
    #path = path.reshape((1, 2))  # A row vector
    c_values = path.dot(posmat)  # c values for scaling u
    # scale u by values to get projection
    projected = path.T.dot(c_values)

def dtg(x1,y1,xgoal,ygoal):
    '''Returns absolute distance to goal coordinates'''
    dist = np.sqrt((xgoal - x1)**2 + (ygoal - y1)**2)
    return dist
vecdtg = np.vectorize(dtg)

def distsbw(path):
    dists = []
    for p in range(len(path)-1):
        dists.append(dtg(path[p,0],path[p,1],path[p+1,0],path[p+1,1]))
    return dists

def naninterp(B):
    A = B
    ok = ~np.isnan(A)
    xp = ok.ravel().nonzero()[0]
    fp = A[~np.isnan(A)]
    x  = np.isnan(A).ravel().nonzero()[0]
    A[np.isnan(A)] = np.interp(x, xp, fp)
    return(A)

def make_lr(ports):
    '''returns list of left-right choices for each trial. 0 is right, 1 is left'''
    lr = [2]
    for p in range(1,len(ports)):
        if ports[p-1]==2 and ports[p]==0:
            lr.append(0)
        elif ports[p-1]==1 and ports[p]==2:
            lr.append(0)
        elif ports[p-1]==0 and ports[p]==1:
            lr.append(0)
        else:
            lr.append(1)
    return lr

def ramptest(data,visnum):
    y = np.array(data.iloc[visinds[visnum]-time[0]*250:visinds[visnum]])
    x = np.linspace(0,len(y),len(y))#column vector with x values. linspace(0,Fs*6,Fs*6)
    x = x.reshape(-1,1)
    model = LinearRegression()
    model.fit(x,y)
    r_sq = model.score(x, y)
    coef = model.coef_[0]
    inter = model.intercept_
    return r_sq,coef,inter

def linreg(data,v1,v2):
    '''Return linear regression r squared value, coefficient, and intercept for plotting.
    v1 and v2 must already be columns in data. Ensure no nans in v1 or v2'''
    y = np.array(data.loc[:,v2])
    x = np.array(data.loc[:,v1])
    x = x.reshape(-1,1)
    model = LinearRegression()
    model.fit(x,y)
    r_sq = model.score(x, y)
    coef = model.coef_[0]
    inter = model.intercept_
    return r_sq,coef,inter

def findpath(fr,to):
    '''Return number of path taken, given which port animal went to and from'''
    if fr==0 and to==1: #AB
        path = 0
    elif fr==1 and to==0: #BA
        path = 1
    elif fr==0 and to==2: #AC
        path = 2
    elif fr==2 and to==0: #CA
        path = 3
    elif fr==1 and to==2: #BC
        path = 4
    elif fr==2 and to==1: #CB
        path = 5
    else:
        path = np.nan
    return path

vfpath = np.vectorize(findpath)


def findpath(fr,to):
    if fr==0 and to==1: #AB
        path = 0
    elif fr==1 and to==0: #BA
        path = 1
    elif fr==0 and to==2: #AC
        path = 2
    elif fr==2 and to==0: #CA
        path = 3
    elif fr==1 and to==2: #BC
        path = 4
    elif fr==2 and to==1: #CB
        path = 5
    else:
        path = np.nan
    return path

vfpath = np.vectorize(findpath)

def perc_overlap(lastp,p):
    '''Returns percent of hex identities in current path (p) that are common
    to last taken path betwen same ports. lastp should be array with hexes of
    previous trial, p should be array with hexes of current trial.'''
    return len(np.intersect1d(lastp,p))/len(p)*100

vec_overlap = np.vectorize(perc_overlap)

def simple_rr(data,ports=[0,1,2]):
    '''Returns rolling reward rate per minute for ports of interest. Also
    returns indices of visits to port of interest. Useful for calculating rr
    by trial, etc.'''
    win = 60*250
    if len(ports)==1:
        rwds = data.loc[data.port==ports[0],'rwd'].copy()
        rwds.loc[rwds==-1]=0
        rwds = rwds.reindex(np.arange(data.index[0],data.index[-1]+1),fill_value=np.nan)
        #inds = data.loc[data.port==ports[0]].index
        return rwds.rolling(window=win,min_periods=1).sum()#,inds
    elif len(ports)==2:
        #inds = data.loc[(data.port==ports[0])|(data.port==ports[1])].index
        #make column with 1s at points of reward
        rwds = data.loc[(data.port==ports[0])|(data.port==ports[1]),'rwd'].copy()
        rwds.loc[rwds==-1]=0
        rwds = rwds.reindex(np.arange(data.index[0],data.index[-1]+1),fill_value=np.nan)
        return rwds.rolling(window=win,min_periods=1).sum()#,inds
    elif len(ports)==3:
        rwds = data.rwd.copy()
        rwds.loc[rwds==-1]=0
        #inds = data.loc[(data.rwd==0)|(data.rwd==1)].index
        #create new df/series with all rows from original, 1s in rwd indices
        return rwds.rolling(window=win,min_periods=1).sum()#,inds

#calculate number of rewards in past t trials for specified ports
def tri_rr(data,ports=[0,1,2],t=5):
    '''Return trial-by-trial estimte of previous reward history over past t trials'''
    if len(ports)==1:
        rwds = data.loc[data.port==ports[0],'rwd'].copy()
        rwds.loc[rwds==-1]=0
        rwds = rwds.rolling(window=t).sum()
        rwds = rwds.reindex(np.arange(data.index[0],data.index[-1]+1),fill_value=np.nan)
        return rwds.fillna(method='ffill')
    elif len(ports)==2:
        rwds = data.loc[(data.port==ports[0])|(data.port==ports[1]),'rwd'].copy()
        rwds.loc[rwds==-1]=0
        rwds = rwds.rolling(window=t).sum()
        rwds = rwds.reindex(np.arange(data.index[0],data.index[-1]+1),fill_value=np.nan)
        return rwds.fillna(method='ffill')
    elif len(ports)==3:
        rwds = data.loc[data.port.isnull()==False,'rwd'].copy()
        rwds.loc[rwds==-1]=0
        rwds = rwds.rolling(window=t).sum()
        rwds = rwds.reindex(np.arange(data.index[0],data.index[-1]+1),fill_value=np.nan)
        return rwds.fillna(method='ffill')


def reg_in_time(data,ttp=4,tfp=4,trace='green',factors=['pchosen','Qc_allo','Qc_ego','speed']):
    '''Calculate regression of behavioral factors to
    dLight values during approach. can add extra factors to regress'''
    dframe = data
    tridat = dframe.loc[dframe.port.isnull()==False]
    lags = np.arange(-fs*ttp,fs*tfp)
    vinds = dframe.loc[dframe.port.isnull()==False].index
    faclen = len(factors)
    rweights = np.zeros((faclen,len(lags)))
    rserrors = np.zeros((faclen,len(lags)))
    speed = data.vel.fillna(method='ffill')
    speed = speed.fillna(method='bfill')
    for n in range(len(lags)):
        y = dframe.loc[vinds+lags[n],trace]
        y = y.reset_index(drop=True)
        X = pd.DataFrame()
        for f in factors:
            if f=='pchosen':
                X[f] = data.loc[vinds,f].values/100
            elif f=='speed':
                X[f] = speed.loc[vinds+lags[n]].values
            else:
                X[f] = data.loc[vinds,f].values
        #scale factors to normalize for interpretable comparison of beta values
        X[factors] = scale.fit_transform(X[factors].as_matrix())
        mod = sm.GLS(y, X).fit()
        rweights[:,n] = mod.params.values
        rserrors[:,n] = mod.bse.values
        facnames = X.columns
    return facnames,rweights,rserrors

#create function to take factor of interest only for port(s) of interest
def get_avail(p):
    allports = [0,1,2]
    allports.pop(p)
    return allports

def avg_factor(factor,portz,sampledata):
    '''ports is a list of ports to average the factors over (one for each trial).
    factor must be either Q_ego, Q_allo, nom_rwd, or rhist'''
    visinds = sampledata.loc[sampledata.port.isnull()==False].index
    if factor == 'Q_ego':
        df = sampledata.loc[visinds,['Q_ego_a','Q_ego_b','Q_ego_c']].values
    elif factor == 'Q_allo':
        df = sampledata.loc[visinds,['Q_allo_a','Q_allo_b','Q_allo_c']].values
    elif factor == 'nom_rwd':
        df = sampledata.loc[visinds,['nom_rwd_a','nom_rwd_b','nom_rwd_c']].values
    elif factor == 'rhist':
        df = sampledata.loc[visinds,['rhist_a','rhist_b','rhist_c']].values
    else:
        print("factor not yet defined for grouping.")
        return None
    fval = []
    for i in range(len(portz)):
        fval.append(np.mean(df[i,portz[i]]))
    return fval

def factor_by_p_type(factor,sampledata,p_type='all'):
    '''p_type can be all, avail, or chosen. returns value of factor given p_type for every trial.'''
    visinds = sampledata.loc[sampledata.port.notnull()].index
    if p_type == 'chosen':
        portz = sampledata.loc[visinds,'port'].values.astype(int)
    elif p_type == 'avail':
        portz = [[0,1,2]]
        for p in sampledata.loc[visinds[:-1],'port'].values.astype(int):
            portz.append(get_avail(p))
    elif p_type == 'all':
        portz = np.tile([0,1,2],(len(visinds),1))
    elif p_type == 'other':
        portz = [0]
        for p in range(len(visinds)-1):
            avail = get_avail(int(sampledata.loc[visinds[p],'port']))
            chos = int(sampledata.loc[visinds[p+1],'port'])
            avail.remove(chos)
            portz.append(avail)
    elif p_type == 'last':
        portz = np.concatenate([[0],sampledata.loc[visinds[:-1],'port'].values.astype(int)])
    return avg_factor(factor,portz,sampledata)

def reg_by_ptype(factor,sampledata):
    addto = sampledata.loc[sampledata.port.notnull()].index
    bytype = pd.DataFrame()
    bytype[factor+'_avail'] = factor_by_p_type(factor,sampledata,'avail')
    bytype[factor+'_chosen'] = factor_by_p_type(factor,sampledata,'chosen')
    bytype[factor+'_all'] = factor_by_p_type(factor,sampledata,'all')
    bytype[factor+'_other'] = factor_by_p_type(factor,sampledata,'other')
    bytype[factor+'_last'] = factor_by_p_type(factor,sampledata,'last')
    bytype = bytype.set_index(addto)
    bytype = bytype.reindex(np.arange(sampledata.index[0],sampledata.index[-1]+1),fill_value=np.nan)
    
    sampledata[factor+'_avail'] = bytype[factor+'_avail'].fillna(method='bfill')
    sampledata[factor+'_avail'] = sampledata[factor+'_avail'].fillna(method='bfill')
    sampledata[factor+'_last'] = bytype[factor+'_last'].fillna(method='bfill')
    sampledata[factor+'_last'] = sampledata[factor+'_last'].fillna(method='bfill')
    sampledata[factor+'_other'] = bytype[factor+'_other'].fillna(method='bfill')
    sampledata[factor+'_other'] = sampledata[factor+'_other'].fillna(method='bfill')
    sampledata[factor+'_chosen'] = bytype[factor+'_chosen'].fillna(method='bfill')
    sampledata[factor+'_chosen'] = sampledata[factor+'_chosen'].fillna(method='bfill')
    sampledata[factor+'_all'] = bytype[factor+'_all'].fillna(method='bfill')
    sampledata[factor+'_all'] = sampledata[factor+'_all'].fillna(method='bfill')
    return sampledata

#make vectorized leakyrr calculator
from scipy.stats import spearmanr     
def leaky(vec,kernel,ind):
    vec[ind:ind+len(kernel)] = vec[ind:ind+len(kernel)] + kernel
    return vec
vleaky = np.vectorize(leaky) #this might not work, loop needs to be executed in order

def tauloop(n,tt,rwdentries,cors):
    rrtau = np.zeros(endtime)
    kernel = np.exp((-np.linspace(0,5*tau)/tau))
    for e in rwdentries:
        rrtau[e:e+len(kernel)] = rrtau[e:e+len(kernel)] + kernel
    if len(rrtau) != len(tt):
        rrtau = rrtau[:len(tt)]
    cor = spearmanr([rrtau,tt],axis=1,nan_policy='omit')[0]
    cors.append(cor)
    return cors
vtauloop = np.vectorize(tauloop,otypes=[list])
    
def optimize_tau(data,rwdentries,m=1,n=200,pad=3500):
    '''call to iterate through n tau values and return tau that gives maximum negative
    correlation between rr and tt. tt is trialtimes'''
    endtime = len(data) + fs*pad
    tt = data.vel_AUC.fillna(method='bfill')#behavioral parameter for optimization, velocity AUC
    tt = tt.fillna(method='ffill')
    cors = [] #store n row vector of correlation coefficients between rrs of dif tau and tt
    taus = np.arange(m,n+1)
    for tau in taus:
        rrtau = np.zeros(endtime)
        kernel = np.exp((-np.linspace(0,5*tau,fs*tau)/tau))
        for e in rwdentries:
            rrtau[e:e+len(kernel)] = rrtau[e:e+len(kernel)] + kernel
        if len(rrtau) != len(data):
            rrtau = rrtau[:len(data)]
        cor = spearmanr([rrtau,tt],axis=1,nan_policy='omit')[0]
        cors.append(cor)
    maxcor = np.max(cors)
    opt_tau = taus[cors.index(maxcor)]
    return opt_tau,maxcor,cors
    
def est_leakint_rr(data,ports,tau=0,optimize=True):
    '''Returns leaky integrator reward rate for an input vector of reward outcomes over n timesteps.
    grabs reward outcomes from data dataframe. tau = time constant. Optimize selects whether to optimize
    tau as value that results in maximum correlation between trial speed metric and rewrate. method
    specifies whether to consider all ports, available ports, or chosen port in leaky integrator'''
    pad = 3500 #not sure why this value, or what function is
    if len(ports) == 3:
        rwdinds = data.loc[data.rwd==1].index
    elif len(ports) == 2:
        rwddat = data.loc[(data.rwd==1)]
        rwdinds = rwddat.loc[(rwddat.port==ports[0])|(rwddat.port==ports[1])].index
    elif len(ports) == 1:
        rwdinds = data.loc[(data.rwd==1)&(data.port==ports[0])].index
    if tau==0 or optimize==True:
        tau,maxcor,corr = optimize_tau(data,rwdinds,1,400,pad)
    #now calculate optimal leaky integrator for reward rate
    endtime = len(data) + fs*pad
    kernel = np.exp((-np.linspace(0,5*tau,fs*tau)/tau))
    rr = np.zeros(endtime)
    for rwdind in rwdinds:
        rr[rwdind:rwdind+len(kernel)] = rr[rwdind:rwdind+len(kernel)] + kernel
    if len(rr) != len(data): #make sure length of rr vector == number of samples
        rr = rr[:len(data)]
    if optimize==False:
        tt = data.vel_AUC.fillna(method='bfill')#behavioral parameter for optimization, velocity AUC
        tt = tt.fillna(method='ffill')
        cor = spearmanr([rr,tt],axis=1,nan_policy='omit')[0]
    return rr,corr,tau


def tauloop_tri(n,tt,rwdentries,cors):
    rrtau = np.zeros(endtime)
    kernel = np.exp((-np.linspace(0,5*tau)/tau))
    for e in rwdentries:
        rrtau[e:e+len(kernel)] = rrtau[e:e+len(kernel)] + kernel
    if len(rrtau) != len(tt):
        rrtau = rrtau[:len(tt)]
    cor = spearmanr([rrtau,tt],axis=1,nan_policy='omit')[0]
    cors.append(cor)
    return cors
vtauloop_tri = np.vectorize(tauloop,otypes=[list])
    
def optimize_tau_tri(data,rwdentries,m=1,n=50):
    '''call to iterate through n tau values and return tau that gives maximum negative
    correlation between rr and tt. tt is trialtimes'''
    tdat = data.loc[data.port!=-100]
    endtime = len(tdat)
    tt = data.visinds.diff() #behavioral parameter for optimization, velocity AUC
    cors = [] #store n row vector of correlation coefficients between rrs of dif tau and tt
    taus = np.arange(m,n+1)
    for tau in taus:
        rrtau = np.zeros(endtime)
        kernel = np.exp((-np.linspace(0,5*tau,tau)/tau))
        for e in rwdentries:
            rrtau[e:e+len(kernel)] = rrtau[e:e+len(kernel)] + kernel
        if len(rrtau) != len(tdat):
            rrtau = rrtau[:len(tdat)]
        cor = spearmanr([rrtau,tt],axis=1,nan_policy='omit')[0]
        cors.append(cor)
    maxcor = np.max(cors)
    opt_tau = taus[cors.index(maxcor)]
    return opt_tau,maxcor,cors
    
def est_leakint_rr_tri(data,ports,tau=0,optimize=True):
    '''Returns leaky integrator reward rate for an input vector of reward outcomes over n timesteps.
    grabs reward outcomes from data dataframe. tau = time constant. Optimize selects whether to optimize
    tau as value that results in maximum correlation between trial speed metric and rewrate. method
    specifies whether to consider all ports, available ports, or chosen port in leaky integrator'''
    if len(ports) == 3:
        rwdinds = data.loc[data.rwd==1].index
    elif len(ports) == 2:
        rwddat = data.loc[(data.rwd==1)]
        rwdinds = rwddat.loc[(rwddat.port==ports[0])|(rwddat.port==ports[1])].index
    elif len(ports) == 1:
        rwdinds = data.loc[(data.rwd==1)&(data.port==ports[0])].index
    if tau==0 or optimize==True:
        tau,maxcor,corr = optimize_tau_tri(data,rwdinds,1,50)
    #now calculate optimal leaky integrator for reward rate
    endtime = len(data) + fs*pad
    kernel = np.exp((-np.linspace(0,5*tau,fs*tau)/tau))
    rr = np.zeros(endtime)
    for rwdind in rwdinds:
        rr[rwdind:rwdind+len(kernel)] = rr[rwdind:rwdind+len(kernel)] + kernel
    if len(rr) != len(data): #make sure length of rr vector == number of samples
        rr = rr[:len(data)]
    if optimize==False:
        tt = data.vel_AUC.fillna(method='bfill')#behavioral parameter for optimization, velocity AUC
        tt = tt.fillna(method='ffill')
        cor = spearmanr([rr,tt],axis=1,nan_policy='omit')[0]
    return rr,corr,tau

def choice_reg(factors,data):
    y = data.loc[visinds,'lrchoice']
    X = data.loc[visinds,factors]
    X = X.drop(y.loc[y==2].index,axis=0)
    y = y.drop(y.loc[y==2].index,axis=0)
    X[factors] = scale.fit_transform(X[factors].as_matrix())
    mod = Logit(y, X).fit()
    pdf = Logit(y, X).pdf(np.linspace(X[factors[0]].min(),X[factors[0]].max(),100))
    rweights = mod.params.values
    rserrors = mod.bse.values
    return rweights,rserrors,pdf


def get_latencies(sampledata,tridat):
    '''sampledata must contain apzone column'''
    visinds = sampledata.loc[sampledata.port.isnull()==False].index
    apstarts = [] #indices of approach start (used for latency, smearing, etc.)
    latencies = [] #list of latency for each trial
    for i in range(len(visinds)):
        apstarts.append(visinds[i]+np.where(np.diff(sampledata.loc[visinds[i]:,'apzone'])==1)[0][0])
        if i==0:
            latencies.append(visinds[i]/250)
        if i<len(visinds)-1:
            latencies.append((visinds[i+1]-apstarts[i])/250)
            
    tridat['latency']=latencies
    addto = visinds
    lats = pd.DataFrame(latencies)
    lats = lats.set_index(addto)
    lats = lats.reindex(np.arange(sampledata.index[0],sampledata.index[-1]+1),fill_value=np.nan)
    sampledata['latency'] = lats.fillna(method='bfill')
    lat = plt.figure()
    plt.title('Distribution of latencies')
    plt.hist(latencies,bins=100)
    plt.ylabel('# of trials')
    plt.xlabel('latency (s)')
    lat.savefig(datepath+'latency_dist.pdf')
    return sampledata,tridat
