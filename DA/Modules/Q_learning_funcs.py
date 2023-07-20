'''functions to estimate q values for each port/path based
on parameters optimized in Julia to fit choice behavior.
'''

__author__ = "Tim Krausz"
__email__ = "krausz.tim@gmail.com"
__status__ = "development"
__license__ = "MIT"

from scipy.special import erf
from scipy.stats import logistic
import numpy as np

def hyb_q_choice(data,params):
    # basic rescorla wagner learning over egocentric responses
    Qs = (np.zeros((len(np.unique(data.port)),1))+0.5).tolist() #list of nchoices lists for Q values
    Q = np.zeros((len(np.unique(data.port)),2))+0.5
    beta = params[0]
    lrmf = 0.5 + 0.5*erf(params[1]/np.sqrt(2))
    lrmb = 0.5 + 0.5*erf(params[2]/np.sqrt(2))
    
    # these are ordered in a left-right logic, unlike the ordering in the
    # allocentric version   

    s = data.port.astype(int).values[:-1]
    c = data.lrchoice.astype(int).values[1:]
    d = data.port.astype(int).values[1:]
    r = data.rwd.astype(int).values[1:]
    #r[np.where(r==0)]=-1
    
    Qs[0].append(0.5)
    Qs[1].append(0.5)
    Qs[2].append(0.5)
    seshs = data.loc[:,"session"].values
    sesh = seshs[0]
    for i in range(len(c)):
        if seshs[i]!=sesh:
            Q = np.zeros((len(np.unique(data.port)),2))+0.5
            Qs
        alts = np.setdiff1d([0,1,2],[s[i],d[i]])
        altc = 1 - c[i]
        Q[s[i],c[i]] = (1-lrmf)*Q[s[i],c[i]] + lrmf *r[i]
        Q[alts,altc] = (1-lrmb) * Q[alts,altc] + lrmb * r[i]
        if s[i] == 0:
            Qs[1].append(Q[s[i],0])
            Qs[2].append(Q[s[i],1])
            Qs[0].append(Qs[0][-1])#this is not updated, same as last trial
        elif s[i] == 1:
            Qs[2].append(Q[s[i],0])
            Qs[0].append(Q[s[i],1])
            Qs[1].append(Qs[1][-1])
        else:
            Qs[0].append(Q[s[i],0])
            Qs[1].append(Q[s[i],1])
            Qs[2].append(Qs[2][-1])
        sesh = seshs[i]
    return np.transpose(Qs)[:-1]

def port_q_choice(data,params):
    # basic rescorla wagner learning over egocentric responses
    Qs = (np.zeros((len(np.unique(data.port)),1))+0.5).tolist() #list of nchoices lists for Q values
    Q = np.zeros(3)+0.5
    beta = params[0]
    lr = 0.5 + 0.5*erf(params[1]/np.sqrt(2))
    

    c = data.port.astype(int).values[1:]
    r = data.rwd.astype(int).values[1:]
    
    Qs[0].append(0.5)
    Qs[1].append(0.5)
    Qs[2].append(0.5)
    seshs = data.loc[:,"session"].values[1:]
    sesh = seshs[0]
    for i in range(len(c)):
        if seshs[i]!=sesh:
            Q = np.zeros(3)+0.5
        Q[c[i]] = (1-lr)*Q[c[i]] + lr *r[i]
        #why all this extra shit,
        Qs[0].append(Q[0])
        Qs[1].append(Q[1])
        Qs[2].append(Q[2])
        sesh = seshs[i]

    return np.transpose(Qs)[:-1]
