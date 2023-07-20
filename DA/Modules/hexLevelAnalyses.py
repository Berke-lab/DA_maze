"""Functions to analyze and vizualize hex-level photometry data.
TODO: create combined photometry rat child class and add functions as methods."""

__author__ = "Tim Krausz"
__email__ = "krausz.tim@gmail.com"
__status__ = "development"
__license__ = "MIT"

from statsmodels.api import OLS,add_constant
from scipy.stats import wilcoxon
import numpy as np
import pandas as pd
from __main__ import *

def create_tridf(photrats):
    photrats.tridf = photrats.df.loc[photrats.df.port!=-100,:].copy()

def create_triframe(photrats):
    photrats.get_visinds()
    photrats.triframe = photrats.df.loc[photrats.visinds]
    photrats.triframe['lrchoice'] = make_lr(photrats.triframe.port.values)
    photrats.triframe['lrchoice'] = photrats.triframe.lrchoice.astype("int8")
    photrats.triframe.reset_index(inplace=True)

def add_withinBlockTriNumberToDf(hexData):
    hexData.df.loc[:,"tri"] = 0
    for s in hexData.df.session.unique():
        t = 1
        for b in hexData.df.loc[hexData.df.session==s,"block"].unique():
            hexData.df.loc[(hexData.df.session==s)&(hexData.df.block==b),"tri"] = \
            hexData.df.loc[(hexData.df.session==s)&(hexData.df.block==b),"trial"] - \
            hexData.df.loc[(hexData.df.session==s)&(hexData.df.block==b),"trial"].min()+1
    hexData.df.loc[:,"tri"]

def calc_distToPortFromHex(hexData,sesh,block,chosenPort,hexLabel,ses_type):
    portstrings = ['A','B','C']
    dists = []
    for state in phexdf.loc[phexdf.to_state==hexLabel].statecodes:
        if ses_type == 'prob':
            dists.append(hexData.sesh_hexDists['dto'+\
                portstrings[chosenPort]][sesh][state])
        else:
            try:
                dists.append(hexData.sesh_hexDists['dto'+\
                    portstrings[chosenPort]][sesh][int(block-1)][state])
            except:
                print(f"No dists for session {sesh} block {block}")
                dists.append(np.nan)
    return np.nanmin(dists)

def add_mapHexDistFromPort2Df(hexData):
    '''add a column to the dataframe indicating the shortest distance
    from the current hex to the chosen port, according to the current
    hex map.'''
    hexData.df.loc[:,"hexesFromPort"] = -1
    chosenPorts = hexData.df.port.replace(-100,method='bfill')
    hexData.df.loc[:,["hexlabel","session","block","session_type"]].values
    mapDist2port = []
    for i in hexData.df.index.values:
        sesh = hexData.df.loc[i,"session"].astype(int)
        block = hexData.df.loc[i,"block"].astype(int)
        hexLabel = hexData.df.loc[i,"hexlabel"].astype(int)
        chosenPort = chosenPorts.loc[i].astype(int)
        ses_type = hexData.df.loc[i,"session_type"]
        mapDist2port.append(calc_distToPortFromHex(hexData,
            sesh,block,chosenPort,hexLabel,ses_type))
    hexData.df.loc[:,"hexesFromPort"] = mapDist2port

def add_hexesFromPort2Df(hexData):
    hexData.get_visinds()
    hexData.df.loc[:,"hexesFromPort"] = 0
    nHexes = hexData.df.loc[:hexData.visinds[0],:].shape[0]
    hexData.df.loc[:hexData.visinds[0],"hexesFromPort"] = np.arange(-nHexes+1,1)
    for t in range(1,len(hexData.visinds)):
        nHexes = hexData.df.loc[hexData.visinds[t-1]:hexData.visinds[t],:].shape[0]
        hexData.df.loc[hexData.visinds[t-1]:hexData.visinds[t],"hexesFromPort"] = np.arange(-nHexes+1,1)
    hexData.df.loc[:,"hexesFromPort"] = hexData.df.loc[:,"hexesFromPort"]*-1
    #hexData.df.loc[hexData.df.hexlabel.isin([1,2,3]),"hexesFromPort"]=0

def calc_hexRampRatDict(hexData,stHex=10):
    try:
        hexData.df.loc[:,"hexesFromPort"]
    except:
        print("calculating distance from port")
        add_hexesFromPort2Df(hexData)
    hexData.ratRamps = {r:[] for r in hexData.df.rat.unique()}
    hexData.ratRampSems = {r:[] for r in hexData.df.rat.unique()}
    hexData.ratTris = {r:[] for r in hexData.df.rat.unique()}
    for r in hexData.ratRamps.keys():
        hexData.ratRamps[r] = np.flip(hexData.df.loc[(hexData.df.rat==r),["DA","hexesFromPort"]].\
                groupby("hexesFromPort").mean().loc[0:stHex,"DA"].values)
        hexData.ratRampSems[r] = np.flip(hexData.df.loc[(hexData.df.rat==r),["DA","hexesFromPort"]].\
                groupby("hexesFromPort").sem().loc[0:stHex,"DA"].values)
        hexData.ratTris[r] = [len(hexData.df.loc[(hexData.df.rat==r),"session"].unique()),
                              len(hexData.df.loc[(hexData.df.rat==r)&(hexData.df.port!=-100),])]

def get_hexRampRatFiberDict(hexData,stHex=10):
    try:
        hexData.df.loc[:,"hexesFromPort"]
    except:
        print("calculating distance from port")
        #add_hexesFromPort2Df(hexData)
        add_mapHexDistFromPort2Df(hexData)
    hexData.ratFiberRamps = {r:[] for r in hexData.seshInfo.ratLoc.unique()}
    hexData.ratFiberRampSems = {r:[] for r in hexData.seshInfo.ratLoc.unique()}
    hexData.ratFiberTris = {r:[] for r in hexData.seshInfo.ratLoc.unique()}
    for rloc in hexData.ratFiberRamps.keys():
        seshs = hexData.seshInfo.loc[hexData.seshInfo.ratLoc==rloc,"session"].unique()
        hexData.ratFiberRamps[rloc] = \
        np.flip(hexData.df.loc[(hexData.df.session.isin(seshs)),["DA","hexesFromPort"]].\
                groupby("hexesFromPort").mean().loc[0:stHex,"DA"].values)
        hexData.ratFiberRampSems[rloc] = \
        np.flip(hexData.df.loc[(hexData.df.session.isin(seshs)),["DA","hexesFromPort"]].\
                groupby("hexesFromPort").sem().loc[0:stHex,"DA"].values)
        hexData.ratFiberTris[rloc] = \
        [len(hexData.df.loc[(hexData.df.session.isin(seshs)),"session"].unique()),
            len(hexData.df.loc[(hexData.df.session.isin(seshs))\
                &(hexData.df.port!=-100),])][:stHex]

def plot_individualRatRamps(hexData):
    calc_hexRampRatDict(hexData,stHex=16)
    fig = plt.figure(figsize=(5,4))
    xvals = np.arange(-15,1)
    for rat in hexData.ratRamps.keys():
        plt.plot(xvals,hexData.ratRamps[rat],alpha=1,color="darkgreen",lw=4)
        plt.fill_between(xvals,hexData.ratRamps[rat]+hexData.ratRampSems[rat],
                        hexData.ratRamps[rat]-hexData.ratRampSems[rat],color="k",alpha=0.3)
        plt.xlabel("Distance from port (hexes)",fontsize=20,fontweight="bold")
        plt.ylabel("DA (z-scored)",fontsize=20,fontweight="bold")
        plt.title(str(hexData.ratTris[rat][0])+" sessions; "+str(hexData.ratTris[rat][1])+" trials.")
        plt.tight_layout()
        fig.savefig(hexData.directory_prefix+"avgRampRat_"+rat+".pdf")
        plt.clf()

def plot_individualRatRampsByFiber(hexData,saveFigs=True):
    get_hexRampRatFiberDict(hexData,stHex=16)
    fig = plt.figure(figsize=(5,4))
    xvals = np.arange(-16,1)
    lastrat = None
    ratNums = {rl:'' for rl in hexData.ratFiberRamps.keys()}
    for rloc in hexData.ratFiberRamps.keys():
        rat = rloc.split(":")[0]
        plt.plot(xvals,hexData.ratFiberRamps[rloc],alpha=1,lw=4)#color="darkgreen"
        plt.fill_between(xvals,hexData.ratFiberRamps[rloc]+hexData.ratFiberRampSems[rloc],
                        hexData.ratFiberRamps[rloc]-hexData.ratFiberRampSems[rloc],color="k",alpha=0.3)
        plt.xlabel("Distance from port (hexes)",fontsize=20)
        plt.ylabel("DA (z-scored)",fontsize=20)
        ratNums[rloc] = str(hexData.ratFiberTris[rloc][0])+" sessions; "+str(hexData.ratFiberTris[rloc][1])+" trials."
        #plt.title(str(hexData.ratFiberTris[rloc][0])+" sessions; "+str(hexData.ratFiberTris[rloc][1])+" trials.")
        plt.tight_layout()
        if (lastrat == rat or rat=="IM-1532") and saveFigs==True:
            fig.savefig(hexData.directory_prefix+f"dataByRat/avgRampRatFiber_{rat}_mapDists.pdf")
            plt.clf()
        lastrat = rat
    with open(hexData.directory_prefix+"rat_and_fiber_numbers.txt", 'w') as f:
        f.write(str(ratNums))

def calc_distFromOptPathLen(hexData):
    dfopt = []
    for s in hexData.df.session.unique():
        sdat = hexData.df.loc[hexData.df.session==s,:].copy()
        svisinds = sdat.loc[sdat.port!=-100].index
        dropInds = svisinds[:-1][sdat.loc[svisinds[:-1]+1,'hexlabel'].isin([1,2,3]).values]+1
        if len(dropInds)>0:
            sdat.drop(dropInds,axis=0,inplace=True)
            sdat.reset_index(inplace=True)
            svisinds = sdat.loc[sdat.port!=-100].index
        optlengths = sdat.loc[svisinds,'dtop'].values
        if sdat.index.min()==svisinds[0]:
            takenlengths = np.diff(np.concatenate([[svisidns[0]-1],svisinds]))
        else:
            takenlengths = np.diff(np.concatenate([[sdat.index.min()],svisinds]))
        dist_from_optlen = takenlengths - optlengths
        dist_from_optlen[dist_from_optlen<0]=0
        dfopt = np.concatenate([dfopt,dist_from_optlen])
    hexData.df.loc[:,'d_from_opt_len']=-100
    hexData.df.loc[hexData.df.port!=-100,'d_from_opt_len']=dfopt
    hexData.df.loc[:,'d_from_opt_len'] = hexData.df.d_from_opt_len.replace(-100,method='bfill')
    hexData.df.loc[:,'d_from_opt_len'] = hexData.df.d_from_opt_len.astype("int8")

def add_rewardLagColumns2Df(photrats,lags=5):
    for l in range(1,lags+1):
        photrats.tridf.loc[:,"rt-"+str(l)] = np.nan
        photrats.tridf.loc[:,"samePath_t-"+str(l)] = 0
        photrats.df.loc[:,"rt-"+str(l)] = np.nan
        photrats.df.loc[:,"samePath_t-"+str(l)] = np.nan
        for p in range(3):
            lagInds = photrats.tridf.loc[photrats.tridf.port==p,:].index
            photrats.tridf.loc[lagInds,"rt-"+str(l)] = photrats.tridf.loc[lagInds,'rwd'].shift(l).values
            photrats.tridf.loc[lagInds,"samePath_t-"+str(l)] = (photrats.tridf.loc[lagInds,\
                'fromport']==photrats.tridf.loc[lagInds,'fromport'].shift(l)).astype("int8")
            photrats.df.loc[lagInds,"rt-"+str(l)] = photrats.tridf.loc[lagInds,"rt-"+str(l)]
            photrats.df.loc[lagInds,"samePath_t-"+str(l)] = photrats.tridf.loc[lagInds,"samePath_t-"+str(l)]
        #block out trials before first occurrence of lab
        firstInd = photrats.df.loc[photrats.df["rt-"+str(l)].notnull(),:].index.min()-1
        photrats.df.loc[firstInd,"rt-"+str(l)] = -1
        firstInd = photrats.df.loc[photrats.df["samePath_t-"+str(l)].notnull(),:].index.min()-1
        photrats.df.loc[firstInd,"samePath_t-"+str(l)] = -1
        photrats.df.loc[:,"rt-"+str(l)] = photrats.df.loc[\
            :,"rt-"+str(l)].fillna(method='bfill').astype("float16")
        photrats.df.loc[:,"samePath_t-"+str(l)] = photrats.df.loc[\
            :,"samePath_t-"+str(l)].\
            fillna(method='bfill').astype("float16")
        photrats.df.loc[photrats.df["samePath_t-"+str(l)]==-1,"samePath_t-"+str(l)] = np.nan
        photrats.df.loc[photrats.df["rt-"+str(l)]==-1,"rt-"+str(l)] = np.nan

def add_pathSpecificRwdLagColumns2Df(photrats,lags=5):
    for l in range(1,lags+1):
        photrats.tridf.loc[:,"path_rt-"+str(l)] = np.nan
        photrats.df.loc[:,"path_rt-"+str(l)] = np.nan
        for path in [[0,1],[1,0],[0,2],[2,0],[1,2],[2,1]]:
            lagInds = photrats.tridf.loc[(photrats.tridf.fromport==path[0])&(photrats.tridf.port==path[1]),:].index
            photrats.tridf.loc[lagInds,"path_rt-"+str(l)] = photrats.tridf.loc[lagInds,'rwd'].shift(l).values
            photrats.df.loc[lagInds,"path_rt-"+str(l)] = photrats.tridf.loc[lagInds,"path_rt-"+str(l)]
        #block out trials before first occurrence of lab
        firstInd = photrats.df.loc[photrats.df["path_rt-"+str(l)].notnull(),:].index.min()-1
        photrats.df.loc[firstInd,"path_rt-"+str(l)] = -1
        photrats.df.loc[:,"path_rt-"+str(l)] = photrats.df.loc[:,"path_rt-"+str(l)].fillna(method='bfill').astype("float16")
        photrats.df.loc[photrats.df["path_rt-"+str(l)]==-1,"path_rt-"+str(l)] = np.nan

def add_altPathSpecificRwdLagColumns2Df(photrats,lags=5):
    for l in range(1,lags+1):
        photrats.tridf.loc[:,"altPath_rt-"+str(l)] = np.nan
        photrats.df.loc[:,"altPath_rt-"+str(l)] = np.nan
        paths = [[0,1],[1,0],[0,2],[2,0],[1,2],[2,1]]
        altPaths = [[2,1],[2,0],[1,2],[1,0],[0,2],[0,1]]
        for p in range(len(paths)):
            altpathRwds = []
            path = paths[p]
            altPath = altPaths[p]
            pathLagInds = photrats.tridf.loc[(photrats.tridf.fromport==path[0])&(photrats.tridf.port==path[1]),:].index
            altLagInds = photrats.tridf.loc[(photrats.tridf.fromport==altPath[0])&(photrats.tridf.port==altPath[1]),:].index
            #need to find latest altLagInd that occurred before each pathLagInd
            for i in pathLagInds:
                prevAltInds = altLagInds[altLagInds<i]
                if len(prevAltInds)==0:
                    altpathRwds.append(np.nan)
                else:
                    altpathRwds.append(photrats.tridf.loc[prevAltInds[-l],"rwd"])
            photrats.tridf.loc[pathLagInds,"altPath_rt-"+str(l)] = altpathRwds
            photrats.df.loc[pathLagInds,"altPath_rt-"+str(l)] = photrats.tridf.loc[pathLagInds,"altPath_rt-"+str(l)]
        #block out trials before first occurrence of lab
        firstInd = photrats.df.loc[photrats.df["altPath_rt-"+str(l)].notnull(),:].index.min()-1
        photrats.df.loc[firstInd,"altPath_rt-"+str(l)] = -1
        photrats.df.loc[:,"altPath_rt-"+str(l)] = photrats.df.loc[:,"altPath_rt-"+str(l)].fillna(method='bfill').astype("float16")
        photrats.df.loc[photrats.df["altPath_rt-"+str(l)]==-1,"altPath_rt-"+str(l)] = np.nan

def add_otherPortPriorRewardOutcome(photrats,l=1):
    photrats.tridf.loc[:,"otherPort_rt-1"] = np.nan
    photrats.df.loc[:,"otherPort_rt-1"] = np.nan
    photrats.tridf.loc[:,"otherPort"] =[np.setdiff1d([0,1,2],[photrats.tridf.loc[:,"port"].values[i],\
                        photrats.tridf.loc[:,"port"].shift(1).values[i]])[0] for i in range(len(photrats.tridf))]
    for p in range(3):
        otherRwds = photrats.tridf.loc[photrats.tridf.port==p,"rwd"]
        lagInds = photrats.tridf.loc[photrats.tridf.otherPort==p].index
        for i in lagInds:
            try:
                photrats.tridf.loc[i,"otherPort_rt-1"] = otherRwds.loc[otherRwds.index<i].values[-1]
                photrats.df.loc[i,"otherPort_rt-1"] = photrats.tridf.loc[i,"otherPort_rt-1"]
            except:
                continue
    #block out trials before first occurrence of lab
    firstInd = photrats.df.loc[photrats.df["otherPort_rt-1"].notnull(),:].index.min()-1
    photrats.df.loc[firstInd,"otherPort_rt-1"] = -1
    photrats.df.loc[:,"otherPort_rt-1"] = photrats.df.loc[:,"otherPort_rt-1"].fillna(method='bfill').astype("float16")
    photrats.df.loc[photrats.df["otherPort_rt-1"]==-1,"otherPort_rt-1"] = np.nan
    
def add_lastRewardColumn(photrats):
    photrats.df.loc[:,"lastRwd"] = -100
    photrats.df.loc[photrats.df.port!=-100,"lastRwd"] = photrats.df.loc[photrats.df.port!=-100,"rwd"]
    photrats.df.loc[:,"lastRwd"] = photrats.df.loc[:,"lastRwd"].replace(-100,method="ffill")

def add_pathSimilarity2df(photrats):
    photrats.get_visinds()
    photrats.df.loc[photrats.visinds,"same2"]=photrats.tridf.same2.values
    photrats.df.loc[:,"same2"] = photrats.df.loc[:,"same2"].fillna(method='bfill')
    photrats.df.loc[photrats.visinds,"fromport"]=photrats.tridf.fromport.values
    photrats.df.loc[:,"fromport"] = photrats.df.loc[:,"fromport"].fillna(method='bfill')
    photrats.df.loc[photrats.visinds,"fromport2"]=photrats.tridf.fromport2.values
    photrats.df.loc[:,"fromport2"] = photrats.df.loc[:,"fromport2"].fillna(method='bfill')

def get_back2samePortTrials(photrats):
    photrats.tridf.loc[:,"same2"] = 0
    photrats.tridf.loc[photrats.tridf.port==photrats.tridf.fromport2,"same2"] = 1
    
def get_portHist(photrats):
    photrats.tridf.loc[:,"fromport"] = photrats.tridf.loc[:,"port"].shift(1).values
    photrats.tridf.loc[:,"fromport2"] = photrats.tridf.loc[:,"port"].shift(2).values
    
def add_nHexesTraversed2port(photrats):
    photrats.df.loc[:,"hexesFromPort"] = 0
    nHexes = photrats.df.loc[:photrats.visinds[0],:].shape[0]
    photrats.df.loc[:photrats.visinds[0],"hexesFromPort"] = np.arange(-nHexes+1,1)
    for t in range(1,len(photrats.visinds)):
        nHexes = photrats.df.loc[photrats.visinds[t-1]:photrats.visinds[t],:].shape[0]
        photrats.df.loc[photrats.visinds[t-1]:photrats.visinds[t],"hexesFromPort"] = np.arange(-nHexes+1,1)
    photrats.df.loc[:,"hexesFromPort"] = photrats.df.loc[:,"hexesFromPort"]*-1


def get_scaledRampsByRwdSeq(hexDf,rwd_pattern,stHex = 15,alongPath=False,getPriorTri=True):
    rwd_pattern_inds,rat_IDs = getRwdPatternInds(hexDf,rwd_pattern,alongPath,getPriorTri)
    
    scaledRamps = [[] for _ in range(len(rwd_pattern)+int(getPriorTri))]
    for patternedEntryInds,rat in zip(rwd_pattern_inds,rat_IDs):
        for e in range(len(scaledRamps)):
            trace = hexDf.loc[patternedEntryInds[e]-stHex:\
                                        patternedEntryInds[e],"DA"].values
            gain,icept = calcRampGain(hexData,rat,trace)
            scaledRamps[e].append(hexData.ratRamps[rat]*gain+icept)
    return scaledRamps

def plot_rampGainByRwdSeq(hexData,rwd_pattern,scaledRamps,stHex = 15,
    fsize=[4.5,8],ylim=[-.14,0.4],linwid=1,pltError=False):
    '''Plot scaled ramp overlaid with actual DA'''
    xvals = np.arange(-stHex,1)
    fig  = plt.figure(figsize=(fsize[0],fsize[1]))
    ax=plt.subplot(len(scaledRamps)-1,1,1)
    for i in range(1,len(scaledRamps)):
        ax=plt.subplot(len(scaledRamps),1,i,sharey=ax,sharex=ax)
        toplot = np.mean(scaledRamps[i],axis=0)
        plt.plot(xvals,toplot,color="k")
        plt.plot(xvals,np.mean(scaledRamps[i-1],axis=0),color="grey")
        if pltError:
            plt.fill_between(xvals,toplot+sem(scaledRamps[i]),\
                        toplot-sem(scaledRamps[i]),color='darkgrey',alpha=.2)
        plt.grid(visible=True,lw=.5)
        plt.xticks(xvals[::3])
        plt.ylabel("Fitted Ramp",fontsize='xx-large')
        plt.ylim(ylim[0],ylim[1])
        if i < len(scaledRamps)-1:
            ax.set_xticklabels('')
    plt.xlabel("Distance to port (hexes)",fontsize='xx-large')
    plt.tight_layout()

    
def get_ResidualsByRwdSeq(hexData,rwd_pattern,stHex = 15,alongPath=False,getPriorTri=True):
    '''Fit ramp to individusl traces in each reward pattern, subtract
    fitted ramps, and return residuals.'''
    rwd_pattern_inds,rat_IDs = getRwdPatternInds(hexData.df,rwd_pattern,alongPath,getPriorTri)
    
    residuals = [[] for _ in range(len(rwd_pattern)+int(getPriorTri))]
    for patternedEntryInds,rat in zip(rwd_pattern_inds,rat_IDs):
        for e in range(len(residuals)):
            trace = hexData.df.loc[patternedEntryInds[e]-stHex:\
                                        patternedEntryInds[e],"DA"].values
            residuals[e].append(calcRampResiduals(hexData,rat,trace))
    return residuals

def getRwdPatternInds(hexDf,rwd_pattern,alongPath=False,altPath=False,getPriorTri=True):
    patLen = len(rwd_pattern)
    rwd_pattern_inds = []
    rat_IDs = []

    for s in hexDf.session.unique():
        for b in hexDf.loc[hexDf.session==s,"block"].unique():
            dat = hexDf.loc[(hexDf.port!=-100)&(hexDf.session==s)\
                                  &(hexDf.block==b),["port","rwd","samePath_t-1","rt-1"]].copy()
            dat.loc[:,"path"] = dat.port.astype(str)+dat.port.shift(1).fillna(-1).astype(str)
            if alongPath:
                dat.loc[:,"path"] = dat.port.astype(str)+dat.port.shift(1).fillna(-1).astype(str)
                ps = dat.loc[(dat.path!='0.0-1.0')&(dat.path!='2-1.0')&(dat.path!='0-1.0')&(dat.path!='1-1.0'),"path"].unique()
            else:
                ps = [0,1,2]
            for p in ps:
                rwdseq = dat.loc[(dat.path==p)&(dat["samePath_t-1"]==1),"rt-1"] if alongPath else dat.loc[dat.port==p,"rt-1"]
                rwdInds = rwdseq.index.values
                if altPath:
                    inds2check = dat.loc[(dat.port==p)&(dat["samePath_t-1"]==0),].reset_index().index.values-patLen
                    inds2check = inds2check[inds2check>=0]
                else:
                    inds2check = range(len(rwdseq)-patLen)
                patternInds = [(rwdInds[i+1-int(getPriorTri)], rwdInds[i+patLen])
                               for i in inds2check\
                 if np.all(rwdseq.values[i+1:i+patLen+1] == rwd_pattern)]
                if len(patternInds)>0:
                    for r in patternInds:
                        rwd_pattern_inds.append(rwdseq.loc[r[0]:r[1]].index.values)
                    rat_IDs.append(hexDf.loc[hexDf.session==s,"rat"].values[0])
    return rwd_pattern_inds,rat_IDs

def getRwdPatternInds_longSeq(hexDf,rwd_pattern,alongPath=False,getPriorTri=True):
    rwd_pattern_inds = []
    rat_IDs = []

    for s in hexDf.session.unique():
        for b in hexDf.loc[hexDf.session==s,"block"].unique():
            dat = hexDf.loc[(hexDf.port!=-100)&(hexDf.session==s)\
                                  &(hexDf.block==b),["port","rwd"]].copy()
            dat.loc[:,"path"] = dat.port.astype(str)+dat.port.shift(1).fillna(-1).astype(str)
            if alongPath:
                dat.loc[:,"path"] = dat.port.astype(str)+dat.port.shift(1).fillna(-1).astype(str)
                ps = dat.loc[(dat.path!='0.0-1.0')&(dat.path!='2-1.0')&(dat.path!='0-1.0')&(dat.path!='1-1.0'),"path"].unique()
            else:
                ps = [0,1,2]
            for p in ps:
                rwdseq = dat.loc[dat.path==p,"rwd"] if alongPath else dat.loc[dat.port==p,"rwd"]
                rwdInds = rwdseq.index.values
                patternInds = [(rwdInds[i+1-int(getPriorTri)], rwdInds[i+len(rwd_pattern)])
                               for i in range(len(rwdseq)-len(rwd_pattern))\
                 if np.all(rwdseq.values[i:i+len(rwd_pattern)] == rwd_pattern)]
                if len(patternInds)>0:
                    for r in patternInds:
                        rwd_pattern_inds.append(rwdseq.loc[r[0]:r[1]].index.values)
                    rat_IDs.append(hexDf.loc[hexDf.session==s,"rat"].values[0])
    return rwd_pattern_inds,rat_IDs

def get_NullResidualsByRwdSeq(hexData,rwd_pattern,stHex = 10,getPriorTri=True):
    '''Fit ramp to individusl traces in each reward pattern, subtract
    fitted ramps, and return residuals.'''
    shuff_pattern_inds,rat_IDs = getShuffledRwdPatternInds(hexData.df.copy(),rwd_pattern,getPriorTri)
    
    Null_residuals = [[] for _ in range(len(rwd_pattern)+int(getPriorTri))]
    for patternedEntryInds,rat in zip(shuff_pattern_inds,rat_IDs):
        for e in range(len(Null_residuals)):
            trace = hexData.df.loc[patternedEntryInds[e]-stHex:\
                                        patternedEntryInds[e],"DA"].values
            Null_residuals[e].append(calcRampResiduals(hexData,rat,trace))
    return np.mean(Null_residuals,axis=1)

def getShuffledRwdPatternInds(hexDf,rwd_pattern,alongPath=False,altPath=False,getPriorTri=True):
    shuffled_inds = []
    rat_IDs = []
    patLen = len(rwd_pattern)
    for s in hexDf.session.unique():
        for b in hexDf.loc[hexDf.session==s,"block"].unique():
            dat = hexDf.loc[(hexDf.port!=-100)&(hexDf.session==s)\
                                  &(hexDf.block==b),["port","rwd","samePath_t-1","rt-1"]].copy()
            dat.loc[:,"path"] = dat.port.astype(str)+dat.port.shift(1).fillna(-1).astype(str)
            if alongPath:
                dat.loc[:,"path"] = dat.port.astype(str)+dat.port.shift(1).fillna(-1).astype(str)
                ps = dat.loc[(dat.path!='0.0-1.0')&(dat.path!='2-1.0')&(dat.path!='0-1.0')\
                             &(dat.path!='1-1.0'),"path"].unique()
            else:
                ps = [0,1,2]
            for p in ps:
                rwdseq = dat.loc[(dat.path==p)&(dat["samePath_t-1"]==1),"rt-1"] if alongPath\
                    else dat.loc[dat.port==p,"rt-1"]
                if len(rwdseq)<patLen+int(getPriorTri):
                    continue
                rwdInds = rwdseq.index.values
                inds2check = range(len(rwdseq)-patLen)
                pInds = [(rwdInds[i+1-int(getPriorTri)], rwdInds[i+patLen])
                               for i in inds2check\
                 if np.all(rwdseq.values[i+1:i+patLen+1] == rwd_pattern)]
                if len(pInds)==0:
                    continue
                for i in range(len(pInds)):
                    shufInds = rwdInds[np.random.choice(np.arange(0,len(rwdInds)),\
                                                        patLen+int(getPriorTri),replace=False)]
                    shuffled_inds.append(shufInds)
                    rat_IDs.append(hexDf.loc[hexDf.session==s,"rat"].values[0])
    return shuffled_inds,rat_IDs

def plotResidualsByRwdSeq(hexData,rwd_pattern,residuals,null_residuals,
    stHex = 10,fsize=[5.5,7],ylim=[-.11,0.12],linwid=1):
    xvals = np.arange(-stHex,1)
    fig  = plt.figure(figsize=(fsize[0],fsize[1]))
    ax=plt.subplot(len(residuals)-1,1,1)
    for i in range(1,len(residuals)):
        ax=plt.subplot(len(residuals),1,i,sharey=ax)#sharex=ax)
        toplot = np.mean(residuals[i],axis=0)
        plt.plot(xvals,np.mean(nullResids,axis=0)[i]+np.percentile(nullResids[:,i,:],97.5,axis=0),\
                 label="97.5% Null",color='grey',ls='--',lw=1,alpha=0.5)
        plt.plot(xvals,np.mean(nullResids,axis=0)[i]+np.percentile(nullResids[:,i,:],2.5,axis=0),\
                 label="2.5% Null",color='grey',ls=':',lw=1,alpha=0.5)
        plt.scatter(xvals[toplot>=0],toplot[toplot>=0],color='darkred',alpha=1,s=25*linwid,marker="o")
        plt.vlines(xvals[toplot>=0],ymin=np.zeros(len(toplot[toplot>=0])),\
                   ymax=toplot[toplot>=0],color='darkred',lw=linwid)
        plt.scatter(xvals[toplot<0],toplot[toplot<0],color='darkblue',alpha=1,s=25*linwid,marker="o")
        plt.vlines(xvals[toplot<0],ymin=toplot[toplot<0],ymax=np.zeros(len(toplot[toplot<0])),\
                   color='darkblue',lw=linwid)
        plt.xticks(xvals)
        plt.ylabel("")#"DA residual",fontsize='xx-large')
        plt.ylim(ylim[0],ylim[1])
        if i < len(residuals)-1:
            ax.set_xticklabels('')
        else:
            ax.set_xticklabels(xvals)
    plt.xlabel("Distance to port (hexes)",fontsize='xx-large')
    plt.legend()
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.00,bottom=0.00)
    return fig

def calcRampResiduals(hexData,rat,daTrace):
    gain,icept = calcRampGain(hexData,rat,daTrace)
    return daTrace - (hexData.ratRamps[rat]*gain+icept)
    
def calcRampGain(hexData,rat,daTrace):
    ramp = hexData.ratRamps[rat]
    mod = LR(fit_intercept=True,normalize=False).fit(ramp.reshape(-1,1),daTrace.reshape(-1,1))
    return mod.coef_[0][0],mod.intercept_

def get_entriesByRwdSeq(hexDf,rwd_pattern = [0,1,0,0],stHex = 15,alongPath=False,altPath=False,getPriorTri=True):
    rwd_pattern_inds,rat_IDs = getRwdPatternInds(hexDf,rwd_pattern,alongPath,altPath,getPriorTri)
    
    entries = [[] for _ in range(len(rwd_pattern)+int(getPriorTri))]
    for patternedEntryInds in rwd_pattern_inds:
        for e in range(len(entries)):
            entries[e].append(hexDf.loc[patternedEntryInds[e]-stHex:\
                                        patternedEntryInds[e],"DA"].values)
    return entries


def get_entriesByRwdSeq_maxPathGaps(hexDf,rwd_pattern = [0,1,0,0],stHex = 15,alongPath=True,getPriorTri=True):
    rwd_pattern_inds = []
    for s in hexDf.session.unique():
        for b in hexDf.loc[hexDf.session==s,"block"].unique():
            dat = hexDf.loc[(hexDf.port!=-100)&(hexDf.session==s)\
                                  &(hexDf.block==b),["port","rwd"]].copy()
            if alongPath:
                dat.loc[:,"path"] = dat.port.astype(str)+dat.port.shift(1).fillna(-1).astype(str)
                ps = dat.loc[(dat.path!='0.0-1.0')&(dat.path!='2-1.0'),"path"].unique()
            else:
                ps = [0,1,2]
            for p in ps:
                gaps_bw_visits = []
                rwdseq = dat.loc[dat.path==p,"rwd"] if alongPath else dat.loc[dat.port==p,"rwd"]
                rwdInds = rwdseq.index.values
                patternInds = [(rwdInds[i+1-int(getPriorTri)], rwdInds[i+len(rwd_pattern)])
                               for i in range(len(rwdseq)-len(rwd_pattern))\
                 if np.all(rwdseq.values[i:i+len(rwd_pattern)] == rwd_pattern)]
                if len(patternInds)>0:
                    allIndsBwPattern = dat.loc[patternInds[0][0]:patternInds[0][1],"port"].index.values
                    for r in patternInds:
                        rpatInds = rwdseq.loc[r[0]:r[1]].index.values
                    for i in range(1,len(rwd_pattern)):
                        gaps_bw_visits.append(len(allIndsBwPattern[(allIndsBwPattern<rpatInds[i])&\
                               (allIndsBwPattern>rpatInds[i-1])]))
                    if len(gaps_bw_visits) ==0 or np.max(gaps_bw_visits)<=2:
                        rwd_pattern_inds.append(rpatInds)
    
    entries = [[] for _ in range(len(rwd_pattern)+int(getPriorTri))]
    for patternedEntryInds in rwd_pattern_inds:
        for e in range(len(entries)):
            entries[e].append(hexDf.loc[patternedEntryInds[e]-stHex:\
                                        patternedEntryInds[e],"DA"].values)
    return entries

def plot_traceByRwdSeq(rwd_pattern,entries,stHex = 15,pltError=True,
    fsize=[4.5,8],ylim=[-.14,0.45]):
    xvals = np.arange(-stHex,1)
    fig  = plt.figure(figsize=(fsize[0],fsize[1]))
    ax=plt.subplot(len(entries)-1,1,1) if len(entries)>1 else plt.subplot(1,1,1)
    for i in range(1,len(entries)):
        ax=plt.subplot(len(entries)-1,1,i,sharey=ax) if len(entries)>1 else plt.subplot(1,1,1)
        toplot = np.mean(entries[i-1],axis=0)
        plt.plot(xvals,toplot,color='darkgrey',alpha=1)
        if pltError:
            plt.fill_between(xvals,toplot+sem(entries[i-1]),\
                        toplot-sem(entries[i-1]),color='darkgrey',alpha=.2)
        toplot = np.mean(entries[i],axis=0)
        plt.plot(xvals,toplot,color="k")
        if pltError:
            plt.fill_between(xvals,toplot+sem(entries[i]),\
                            toplot-sem(entries[i]),color='k',alpha=.3)
        plt.grid(visible=True,lw=.5)
        plt.xticks(xvals)
        plt.ylabel("DA z-scored",fontsize='xx-large')
        plt.ylim(ylim[0],ylim[1])
        if i < len(entries)-1:
            ax.set_xticklabels('')
    plt.xlabel("Distance to port (hexes)",fontsize='xx-large')
    plt.tight_layout()
    return fig

def plotDifInTracesByRwdSeq(rwd_pattern,entries,stHex = 15,remove_mean=False,pltError=True,linwid=1):    
    xvals = np.arange(-stHex,1)
    fig = plt.figure()
    ax = plt.gca()
    for i in range(0,len(entries)-1):
        ax=plt.subplot(len(entries)-1,1,i+1,sharey=ax,sharex=ax)
        if remove_mean:
            toplot = np.mean(np.subtract(np.subtract(entries[i+1],np.mean(entries[i+1])),\
                                    np.subtract(entries[i],np.mean(entries[i]))),axis=0)
        else:
            toplot = np.mean(np.subtract(entries[i+1],entries[i]),axis=0)
        plt.scatter(xvals[toplot>=0],toplot[toplot>=0],color='darkred',alpha=1,s=35*linwid,marker="o")#"^")
        plt.vlines(xvals[toplot>=0],ymin=np.zeros(len(toplot[toplot>=0])),ymax=toplot[toplot>=0],color='darkred',lw=linwid)
        plt.scatter(xvals[toplot<0],toplot[toplot<0],color='darkblue',alpha=1,s=35*linwid,marker="o")#"v")
        plt.vlines(xvals[toplot<0],ymin=toplot[toplot<0],ymax=np.zeros(len(toplot[toplot<0])),color='darkblue',lw=linwid)
        if pltError:
            plt.fill_between(xvals,toplot+sem(np.subtract(entries[i+1],entries[i])),\
                        toplot-sem(np.subtract(entries[i+1],entries[i])),color='darkgrey',alpha=.2)
        plt.grid(visible=True,lw=.5)
        plt.xticks(xvals[::3])
        plt.ylabel("$\Delta$ DA",fontsize='xx-large')
    plt.xlabel("hexes from port entry",fontsize='xx-large')
    plt.tight_layout()
    return fig

def plot_traceByRwdSeqRatAvg(rwd_pattern,ratDict,stHex = 15,remove_mean = False):    
    xvals = np.arange(-stHex,1)
    plt.figure(figsize=(4.5,8))
    ax=plt.subplot(len(rwd_pattern)+1,1,1)
    for i in range(1,len(rwd_pattern)+1):
        ax=plt.subplot(len(rwd_pattern)+1,1,i,sharey=ax)#,sharex=ax)
        meanTrace,semTrace = get_ratAvgEntry4Plot(ratDict,i-1)
        toplot = meanTrace - np.mean(meanTrace) if remove_mean else meanTrace
        plt.plot(xvals,toplot,color='darkgrey',alpha=1)
        plt.fill_between(xvals,toplot+semTrace,\
                        toplot-semTrace,color='darkgrey',alpha=.2)
        meanTrace,semTrace = get_ratAvgEntry4Plot(ratDict,i)
        toplot = meanTrace - np.mean(meanTrace) if remove_mean else meanTrace
        plt.plot(xvals,toplot,color="k")
        plt.fill_between(xvals,toplot+semTrace,\
                        toplot-semTrace,color='k',alpha=.3)
        plt.grid(visible=True,lw=.5)
        plt.xticks(xvals[::3])
        plt.ylabel("DA z-scored",fontsize='xx-large')
        if i < len(rwd_pattern):
            ax.set_xticklabels('')
    plt.xlabel("hexes from port",fontsize='xx-large')
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1)

def calcSigDifVs0Vector_usingWSR(sampleArray,greater=False,less=False):
    '''for each column in the sample array, return significance level using WSR of distribution vs 0.
    Shape must be (nTests x nSamples)'''
    if greater:
        alternativeString = "greater"
    elif less:
        alternativeString = "less"
    else:
        alternativeString = "two-sided"
    if np.ndim(sampleArray)==1:
        return wilcoxon(sampleArray,alternative=alternativeString)[1]
    sigVals = []
    for sub in sampleArray:
        sigVals.append(wilcoxon(sub,alternative=alternativeString)[1])
    return np.array(sigVals)

def get_sigRats_fromMeanList(meanList,altString="two-sided"):
    pval = wilcoxon(meanList,alternative=altString)[1]
    print("p-value = ",pval)
    return pval<0.001,pval<0.01,pval<0.05

def get_sigRatsPaired_from2samples(A,B,altString="two-sided"):
    pval = wilcoxon(A,B,alternative=altString)[1]
    print("p-value = ",pval)
    return pval<0.001,pval<0.01,pval<0.05

def calcDifInTracesByRwdSeq(entries,stHex = 15,remove_mean=False,pltError=True,linwid=1):
    '''Returns the difference between one run and the next for each pair of runs.
    Entries must be shape (nRunTypes x nRuns x nHexes)'''
    hexDifs = []
    for i in range(0,len(entries)-1):
        hexDifs.append(np.subtract(entries[i+1],entries[i]))
    return hexDifs
    
def get_ratAvgEntry4Plot(ratDict,entryNum):
    ratmeans = []
    for rat in ratDict:
        ratmeans.append(np.nanmean(ratDict[rat][entryNum],axis=0))
    return np.nanmean(ratmeans,axis=0),sem(ratmeans)

def plot_hexRampByPrwd(hexData,stHex=15,ylim=[-.14,.45],var='DA'):
    highmeans = []
    midmeans = []
    lowmeans = []
    midTri=25
    for r in hexData.df.rat.unique():
        highmeans.append(hexData.df.loc[(hexData.df.rat==r)&(hexData.df.tri>midTri)&\
                                         (hexData.df.nom_rwd_chosen>79),[var,"hexesFromPort"]].\
            groupby("hexesFromPort").mean().loc[0:stHex,var].values)
        if len(highmeans[-1]) == 0:
            highmeans = highmeans[:-1]
        midmeans.append(hexData.df.loc[(hexData.df.rat==r)&(hexData.df.tri>midTri)&\
                                        (hexData.df.nom_rwd_chosen==50),\
            [var,"hexesFromPort"]].groupby("hexesFromPort").mean().loc[0:stHex,var].values)
        if len(midmeans[-1]) == 0:
            midmeans = midmeans[:-1]
        lowmeans.append(hexData.df.loc[(hexData.df.rat==r)&(hexData.df.tri>midTri)&\
                                        (hexData.df.nom_rwd_chosen<=20)&\
                    (hexData.df.nom_rwd_chosen>=0),[var,"hexesFromPort"]].\
            groupby("hexesFromPort").mean().loc[0:stHex,var].values)
        if len(lowmeans[-1]) == 0:
            lowmeans = lowmeans[:-1]
    
    xvals = np.arange(-stHex,1)
    fig = plt.figure(figsize=(6,3.7))
    toplt = np.flip(np.mean(highmeans,axis=0))
    plt.plot(xvals,toplt,color='lime',lw=3,label="high pRwd")
    plt.fill_between(xvals,toplt+np.flip(sem(highmeans)),toplt-np.flip(sem(highmeans)),color='forestgreen',alpha=0.5)
    toplt = np.flip(np.mean(midmeans,axis=0))
    plt.plot(xvals,toplt,color='green',lw=3,label="mid pRwd")
    plt.fill_between(xvals,toplt+np.flip(sem(midmeans)),toplt-np.flip(sem(midmeans)),color='green',alpha=0.5)
    toplt = np.flip(np.mean(lowmeans,axis=0))
    plt.plot(xvals,toplt,color='darkgreen',lw=3,label="low pRwd")
    plt.fill_between(xvals,toplt+np.flip(sem(lowmeans)),toplt-np.flip(sem(lowmeans)),color='darkgreen',alpha=0.5)
    plt.xticks(xvals)
    plt.grid(visible=True,lw=.5)
    plt.ylabel(f"{var}",fontsize=20)
    plt.xlabel("Distance to port (hexes)",fontsize=20)
    #plt.title("high traces = "+str(len(highmeans))+"\nmid traces = "+str(len(midmeans))+"\nlow traces = "+str(len(lowmeans)))
    plt.legend()
    plt.ylim(ylim)
    plt.tight_layout()
    return fig

def plot_hexRamp(hexData,stHex=14,rampcolor='k',ylim=[-.17,0.47]):
    tracemeans = []
    velTracemeans = []
    for r in hexData.df.rat.unique():
        tracemeans.append(hexData.df.loc[(hexData.df.rat==r),["DA","hexesFromPort"]].\
            groupby("hexesFromPort").mean().loc[0:stHex,"DA"].values)
        if len(tracemeans[-1]) == 0:
            tracemeans = tracemeans[:-1]
        velTracemeans.append(hexData.df.loc[(hexData.df.rat==r),["vel","hexesFromPort"]].\
            groupby("hexesFromPort").mean().loc[0:stHex,"vel"].values)
        if len(tracemeans[-1]) == 0:
            velTracemeans = tracemeans[:-1]
    
    xvals = np.arange(-stHex,1)
    fig = plt.figure(figsize=(6,3.7))
    toplt = np.flip(np.mean(tracemeans,axis=0))
    plt.plot(xvals,toplt,color=rampcolor,lw=3)
    plt.fill_between(xvals,toplt+np.flip(sem(tracemeans)),\
        toplt-np.flip(sem(tracemeans)),color='k',alpha=0.5)
    plt.grid(visible=True,lw=.5)
    plt.xticks(xvals[::3])
    plt.ylabel("DA (z-scored)",fontsize=20)
    plt.xlabel("Distance to port (hexes)",fontsize=20)
    plt.ylim(ylim)
    ax = plt.gca()
    ax1 = ax.twinx()
    toplt = np.flip(np.mean(velTracemeans,axis=0))
    ax1.plot(xvals,toplt,color='grey',lw=2,ls='--',label='vel')
    #ax1.fill_between(xvals,toplt+np.flip(sem(velTracemeans)),\
    #     toplt-np.flip(sem(velTracemeans)),color='k',alpha=0.5)
    ax1.set_ylabel("speed (cm/s)",fontsize=20)
    plt.tight_layout()
    plt.legend()
    return fig

def calcValRegByRatAndSesh(photrats):
    ratpRwdRegs = {r:[] for r in photrats.df.loc[:,"rat"].unique()}
    ratNs = {r:[] for r in photrats.df.loc[:,"rat"].unique()}
    for r in tqdm(range(len(photrats.df.loc[:,"rat"].unique()))):
        rat = photrats.df.loc[:,"rat"].unique()[r]
        for s in photrats.df.loc[photrats.df.rat==rat,"session"].unique():
            dat = photrats.df.loc[(photrats.df.session==s)&(photrats.df.tri>25),:]
            pRwds = dat.loc[dat.port!=-100,"nom_rwd_chosen"].values/100
            rweights = calcValLagRegWeights(photrats,dat,pRwds)
            ratpRwdRegs[rat].append(rweights[0])
            ratNs[rat].append(len(pRwds))
    return ratpRwdRegs,ratNs

def calcValLagRegWeights(photrats,dat,pRwds):
    dat_vinds = dat.loc[(dat.port!=-100),:].index
    lags = np.arange(-16,0)
    rweights = np.zeros((1,len(lags)))
    for n in range(len(lags)):
        preLagInds = dat_vinds+lags[n]
        lagInds = preLagInds[np.isin(preLagInds,dat.index)]
        y = dat.loc[lagInds,"DA"]
        y.reset_index(drop=True,inplace=True)
        X = pd.DataFrame({"pRwd":pRwds[np.isin(preLagInds,dat.index)]})
        X = X.drop(y.loc[y.isnull()].index,axis=0)
        y = y.drop(y.loc[y.isnull()].index,axis=0)
        modRwd = LR(fit_intercept=True,normalize=False).fit(X,y)
        rweights[:,n] = modRwd.coef_
    return rweights

def plot_distOfPeakRwdEffect(hexData,hbin):
    regWeights = calcDaRhistRegInSpace(hexData,hexbin=hbin,rwdsBack=[1,5],alongPath=True)
    fig = plotDaRhistRegInSpce(hexData,regWeights,hexbin=hbin,rwdsBack=[1,5],alongPath=True)
    fig.savefig(hexData.directory_prefix+f"/rHistRegWeightsInSpaceSamePath_{hbin}hexBin_mapDists.pdf")
    
    xvals=['R(t-1)','R(t-2)','R(t-3)','R(t-4)','R(t-5)']
    fig = plt.figure(figsize=(5,7))
    plt.bar(xvals,np.mean(np.argmax(regWeights,axis=1),axis=0)*hbin,alpha=.5,color='maroon')
    plt.errorbar(xvals,np.mean(np.argmax(regWeights,axis=1)*hbin,axis=0),\
                 yerr=sem(np.argmax(regWeights,axis=1)*hbin,axis=0),fmt='o',color='k')
    plt.ylabel("Distance of peak effect on DA (hexes)",fontsize=20)
    plt.xlabel("Prior reward at port",fontsize=20)
    plt.tight_layout()
    fig.savefig(hexData.directory_prefix+f"/distOfPeakRwdEffect_{hbin}hexBin_mapDists.pdf")
    
    cors = []
    for r in range(hexData.df.loc[hexData.df.rat.notnull(),"rat"].unique().shape[0]):
        cors.append(spearmanr(np.arange(1,6),np.argmax(regWeights,axis=1)[r])[0])
        
    fig = plt.figure(figsize=(3,6))
    plt.bar(0,np.mean(cors),color='k',alpha=0.5)
    sns.stripplot(np.zeros(hexData.df.loc[hexData.df.rat.notnull(),"rat"].unique().shape[0]),\
                  cors,size=15,marker='D',color='maroon',jitter=1)
    plt.ylabel("correlation coefficient")
    plt.xlabel("DA vs prior reward")
    plt.xticks([])
    plt.title("p = "+str(wilcoxon(cors,np.zeros(hexData.df.loc[hexData.df.rat.notnull(),"rat"].unique().shape[0]))[1]))
    plt.tight_layout()
    fig.savefig(hexData.directory_prefix+f"distOfPeakRwdEffect_CorPlot_{hbin}hexBin_mapDists.pdf")

def plot_ratMeans(xvals,ratDict,pltColor='darkred',pltLabel="",plot_ratTraces=False):
    ratmeans = []
    for rat in ratDict:
        ratmeans.append(np.nanmean(ratDict[rat],axis=0))
        if plot_ratTraces:
            plt.plot(xvals,ratmeans[-1],color=pltColor,alpha=0.2,lw=.5)
    plt.plot(xvals,np.nanmean(ratmeans,axis=0),color=pltColor,lw=3,label=pltLabel)
    plt.fill_between(xvals,np.nanmean(ratmeans,axis=0)-sem(ratmeans),\
          np.mean(ratmeans,axis=0)+sem(ratmeans),color=pltColor,lw=3,alpha=0.3)

def plot_sigPoints(xvals,ratDict,pltColor='darkred',plot99=False):
    rat95Errors = []
    rat99Errors = []
    for rat in ratDict:
        er95 = sem(ratDict[rat],axis=0)*1.96
        er99 = sem(ratDict[rat],axis=0)*2.58
        meanTrace = np.nanmean(ratDict[rat],axis=0)
        rat95Errors.append(((meanTrace-er95)>0).astype(int))
        rat99Errors.append(((meanTrace-er99)>0).astype(int))
    if plot99:
        plt.plot(xvals,np.sum(rat99Errors,axis=0)/10,color=pltColor)
    else:
        plt.plot(xvals,np.sum(rat95Errors,axis=0)/10,color=pltColor)

def calcDaRhistRegInSpace(hexData,hexbin = 2,rwdsBack=[1,4],alongPath=True):
    rDf = hexData.df.loc[(hexData.df.notnull().all(axis=1))&(hexData.df.hexesFromPort>=0),:]
    
    reg_weights = []
    distbins = [np.arange(i,i+hexbin) for i in np.arange(0,16,hexbin)]
    rlags = np.arange(rwdsBack[0],rwdsBack[1]+1)
    for r in hexData.df.rat.unique():
        regDf = rDf.loc[(rDf.rat==r),:]
        if len(regDf)==0:
            continue
        r_weights = []
        for dbin in distbins:
            rhistString = "path_rt-" if alongPath else "rt-"
            X = regDf.loc[regDf.hexesFromPort.isin(dbin),\
                          [rhistString+str(n) for n in rlags]].values
            if len(X)==0:
                break
            y = regDf.loc[regDf.hexesFromPort.isin(dbin),"DA"]
            mod = LR(fit_intercept=True,normalize=False).fit(X,y)
            r_weights.append(mod.coef_) #regression at specified distance
        reg_weights.append(r_weights)
    return reg_weights

def plotDaRhistRegInSpce(hexData,reg_weights,hexbin = 2,rwdsBack=[1,4],alongPath=True):
    reg_sems = sem(reg_weights)
    reg_weights = np.mean(reg_weights,axis=0)
    distbins = [np.arange(i,i+hexbin) for i in np.arange(1,17,hexbin)]
    xvals = [str(d[0])+'-'+str(d[0]+hexbin) for d in distbins]
    fig = plt.figure(figsize=(5,3.5))
    for d in range(len(reg_weights.T)):
        plt.plot(xvals,reg_weights.T[d],marker='o',markersize=9,\
            label="t-"+str(d),color='k',alpha=1-((d)/len(reg_weights.T)))
        plt.errorbar(xvals,reg_weights.T[d],yerr=np.transpose(reg_sems)[d],\
            ecolor='k',alpha=.4,lw=.1,elinewidth=2)
    ax=plt.gca()
    ax.invert_xaxis()
    plt.xticks(xvals)
    #plt.grid("on")
    plt.xlabel("Hexes from port",fontsize='xx-large')
    plt.ylabel("Prior Reward "+r"$\beta $",fontsize='xx-large')
    titleString = 'along same path' if alongPath else "along any path"
    plt.title(titleString,fontsize='x-large')
    plt.tight_layout()
    return fig


def calc_rampSlopes(hexData,dat,stHex=16):
    daRamps = []
    datVinds = dat.loc[dat.port!=-100].index
    for i in range(0,len(datVinds)-1):
        tdat = hexData.df.loc[datVinds[i]:datVinds[i+1]]
        tdat = tdat.loc[tdat.hexesFromPort<=stHex,:]
        if len(tdat)<2:
            continue
        X = tdat.index.values - tdat.index.min()
        y = tdat.loc[:,"DA"].values
        modo = LR(fit_intercept=True,normalize=False).fit(X.reshape(-1,1),y)
        daRamps.append(modo.coef_[0])
    return np.mean(daRamps)

def get_rampslopesByRatAndSesh(photrats,stHex=16):
    ratSlopes = {r:[] for r in photrats.df.loc[:,"rat"].unique()}
    ratNs = {r:[] for r in photrats.df.loc[:,"rat"].unique()}
    for r in tqdm(range(len(photrats.df.loc[:,"rat"].unique()))):
        rat = photrats.df.loc[:,"rat"].unique()[r]
        for s in photrats.df.loc[photrats.df.rat==rat,"session"].unique():
            dat = photrats.df.loc[(photrats.df.session==s),:]
            ratSlopes[rat].append(calc_rampSlopes(photrats,dat,stHex))
            ratNs[rat].append(len(dat.loc[dat.port!=-100,:]))
    return ratSlopes,ratNs

def get_rampslopesByRatSeshAndHem(photrats,stHex=16):
    ratSlopes_bar = {r:{fl:[] for fl in photrats.seshInfo.loc[(photrats.seshInfo.rat==r),'fiberloc'].unique()} for r in photrats.df.loc[:,"rat"].unique()}
    ratNs_bar = {r:[] for r in photrats.df.loc[:,"rat"].unique()}
    ratSlopes_prob = {r:{fl:[] for fl in photrats.seshInfo.loc[(photrats.seshInfo.rat==r),'fiberloc'].unique()} for r in photrats.df.loc[:,"rat"].unique()}
    ratNs_prob = {r:[] for r in photrats.df.loc[:,"rat"].unique()}
    ratSlopes_all = {r:{fl:[] for fl in photrats.seshInfo.loc[(photrats.seshInfo.rat==r),'fiberloc'].unique()} for r in photrats.df.loc[:,"rat"].unique()}
    for r in tqdm(range(len(photrats.df.loc[:,"rat"].unique()))):
        rat = photrats.df.loc[:,"rat"].unique()[r]
        for s in photrats.df.loc[(photrats.df.rat==rat)&(photrats.df.session_type=="barrier"),"session"].unique():
            dat = photrats.df.loc[(photrats.df.session==s),:]
            floc = photrats.seshInfo.loc[photrats.seshInfo.session==s,"fiberloc"].values[0]
            ratSlopes_bar[rat][floc].append(calc_rampSlopes(photrats,dat,stHex))
            ratNs_bar[rat].append(len(dat.loc[dat.port!=-100,:]))
        for s in photrats.df.loc[(photrats.df.rat==rat)&(photrats.df.session_type=="prob"),"session"].unique():
            dat = photrats.df.loc[(photrats.df.session==s),:]
            floc = photrats.seshInfo.loc[photrats.seshInfo.session==s,"fiberloc"].values[0]
            ratSlopes_prob[rat][floc].append(calc_rampSlopes(photrats,dat,stHex))
            ratNs_prob[rat].append(len(dat.loc[dat.port!=-100,:]))
        for s in photrats.df.loc[(photrats.df.rat==rat),"session"].unique():
            dat = photrats.df.loc[(photrats.df.session==s),:]
            floc = photrats.seshInfo.loc[photrats.seshInfo.session==s,"fiberloc"].values[0]
            ratSlopes_all[rat][floc].append(calc_rampSlopes(photrats,dat,stHex))
    return ratSlopes_all,ratSlopes_bar,ratNs_bar,ratSlopes_prob,ratNs_prob

def plot_ratMeans_box4EachRat(ratDict,pltColor='darkred'):
    '''Plots distribution of session values for each rat in ratDict
    in form of box plot.'''
    ratDF = pd.DataFrame({ key:pd.Series(value) for key, value in ratDict.items() })
    sns.barplot(data=ratDF,color='lightgrey',ci=None)
    sns.stripplot(data=ratDF,color="k")
    plt.axhline(0,ls=':',color='k')
    #plot astericks over significant rats
    sigrats = np.transpose(get_sigRats(ratDict))
    for r in range(len(ratDict)):
        if sigrats[r][0]==1:
            plt.text(x=r-.25, y=ratDF.max().max()+ratDF.std().max()/1.5, 
                s='***',fontweight='bold',fontsize='x-large')
        elif sigrats[r][1]==1:
            plt.text(x=r-.125, y=ratDF.max().max()+ratDF.std().max()/1.5, 
                s='**',fontweight='bold',fontsize='x-large')
        elif sigrats[r][2]==1:
            plt.text(x=r, y=ratDF.max().max()+ratDF.std().max()/1.5, 
                s='*',fontweight='bold',fontsize='x-large')
        else:
            plt.text(x=r-.25, y=ratDF.max().max()+ratDF.std().max()/1.5, 
                s='ns',fontweight='bold',fontsize='x-large')
            
def plot_ratMeans_box4EachRat_byHem(ratDictall,ratDictbar,ratDictprob,pltColor='darkred'):
    '''Plots distribution of session values for each rat in ratDict
    in form of box plot.'''
    ratDF = pd.DataFrame({ key:pd.Series(value) for key, value in ratDictall.items() })
    ratDFbar = pd.DataFrame({ key:pd.Series(value) for key, value in ratDictbar.items() })
    ratDFprob = pd.DataFrame({ key:pd.Series(value) for key, value in ratDictprob.items() })
    sns.barplot(data=ratDF.loc['NAcc-Core-Left ',],color="#1f77b4",ci=None,alpha=0.4)
    for rat in ratDFbar.columns:
        if not np.any(np.isnan(ratDFbar.loc['NAcc-Core-Left ',rat])):
            plt.scatter(np.tile(rat,len(ratDFbar.loc['NAcc-Core-Left ',rat])),\
                        ratDFbar.loc['NAcc-Core-Left ',rat],color="#1f77b4",marker='X',s=40,edgecolor='k',lw=1)
        plt.scatter(np.tile(rat,len(ratDFbar.loc['NAcc-Core-Right',rat])),\
                    ratDFbar.loc['NAcc-Core-Right',rat],color="#ff7f0e",marker='X',s=40,edgecolor='k',lw=1)
    for rat in ratDFprob.columns:
        if not np.any(np.isnan(ratDFprob.loc['NAcc-Core-Left ',rat])):
            plt.scatter(np.tile(rat,len(ratDFprob.loc['NAcc-Core-Left ',rat])),\
                        ratDFprob.loc['NAcc-Core-Left ',rat],color="#1f77b4",marker='o',s=25,edgecolor='k',lw=1)
        plt.scatter(np.tile(rat,len(ratDFprob.loc['NAcc-Core-Right',rat])),\
                    ratDFprob.loc['NAcc-Core-Right',rat],color="#ff7f0e",marker='o',s=25,edgecolor='k',lw=1)
    sns.barplot(data=ratDF.loc['NAcc-Core-Right',],color="#ff7f0e",ci=None,alpha=0.4)
    plt.axhline(0,ls=':',color='k')

def get_sigRats(ratDict):
    rat95Errors = []
    rat99Errors = []
    rat999Errors = []
    for rat in ratDict:
        er95 = sem(ratDict[rat])*1.96
        er99 = sem(ratDict[rat])*2.58
        er999 = sem(ratDict[rat])*3.291
        meanTrace = np.nanmean(ratDict[rat])
        rat95Errors.append(((meanTrace-er95)>0).astype(int))
        rat99Errors.append(((meanTrace-er99)>0).astype(int))
        rat999Errors.append(((meanTrace-er999)>0).astype(int))
    return rat999Errors,rat99Errors,rat95Errors

def createSegmentLevelDf(hexData):
    segRows = [[] for _ in range(7)]
    cpInds = hexData.df.loc[(hexData.df.choicePoint==1)\
                             |(hexData.df.port!=-100)].index
    for i in range(1,len(cpInds)):
        if hexData.df.loc[cpInds[i],"port"]!=-100:
            segVals = hexData.df.loc[cpInds[i],["DA","port","rwd","block","trial","session","rat"]].values
            for i in range(7):
                segRows[i].append(segVals[i])
            continue
        if hexData.df.loc[cpInds[i],"hexlabels"] in [0,1,2] or hexData.df.loc[cpInds[i-1],"hexlabels"] in [0,1,2]:
            continue
        #segRows.append(hexData.df.loc[cpInds[i-1]:cpInds[i]-1,["DA","port","rwd","block","trial","session","rat"]].mean(axis=0).values)
        segRows[0].append(hexData.df.loc[cpInds[i-1]:cpInds[i]-1,"DA"].mean())
        segRows[1].append(-100)
        segRows[2].append(0)
        segRows[3].append(hexData.df.loc[cpInds[i-1]:cpInds[i]-1,"block"].values[0])
        segRows[4].append(hexData.df.loc[cpInds[i-1]:cpInds[i]-1,"tri"].values[0])
        segRows[5].append(hexData.df.loc[cpInds[i-1]:cpInds[i]-1,"session"].values[0])
        segRows[6].append(hexData.df.loc[cpInds[i-1]:cpInds[i]-1,"rat"].values[0])
    segDf = pd.DataFrame(np.transpose(segRows),columns=["DA","port","rwd","block","trial","session","rat"])
    dtype4col = ["float32","int8","uint8","uint8","uint16","uint8","object"]
    for c in range(len(segDf.columns)):
        segDf.loc[:,segDf.columns[c]] = segDf.loc[:,segDf.columns[c]].astype(dtype4col[c])  

def calcChosVotherVlastRegWeights(hexData,hexbin=1):
    rDf = hexData.df.loc[(hexData.df.hexesFromPort>=0),
    ['DA','port', 'rwd', 'session', 'block', 'trial', 'rat',
       'date', 'session_type', 'nom_rwd_chosen',
       'pairedHexState', 
                'tri', 'hexesFromPort', 'rt-1',
       'samePath_t-1','otherPort_rt-1', 'lastRwd']]
    rDf = rDf.loc[(rDf.notnull().all(axis=1)),:]
    
    reg_weights = []
    pvalues = []
    distbins = [np.arange(i,i+hexbin) for i in np.arange(0,15,hexbin)]
    for r in hexData.df.rat.unique():
        regDf = rDf.loc[(rDf.rat==r),:]
        if len(regDf)==0:
            continue
        r_weights = []
        pvals = []
        for dbin in distbins:
            X = regDf.loc[regDf.hexesFromPort.isin(dbin),\
                          ["rt-1","otherPort_rt-1","lastRwd"]].values
            if len(X)==0:
                break
            X = add_constant(X)
            y = regDf.loc[regDf.hexesFromPort.isin(dbin),"DA"]
            #mod = LR(fit_intercept=True,normalize=False).fit(X,y)
            mod = OLS(y.astype(float),X).fit()
            #r_weights.append(mod.coef_) #regression at specified distance
            r_weights.append(mod.params.values[1:])
            pvals.append(mod.pvalues.values[1:])
        reg_weights.append(r_weights)
        pvalues.append(pvals)
        reg_sems = sem(reg_weights)    
    return np.array(reg_weights),reg_sems,np.array(pvalues)

def plotDaChosenVotherVlastRegInSpace(hexData,hexbin=1,alphaVal=0.05,regStyles=['-','--',':']):
    
    reg_weights_all,reg_sems,pvalues = calcChosVotherVlastRegWeights(hexData,hexbin)
    reg_weights = np.mean(reg_weights_all,axis=0)
    propSig = np.mean((pvalues<alphaVal).astype(int),axis=0)
    distbins = [np.arange(i,i+hexbin) for i in np.arange(0,15,hexbin)]
    xvals = np.array(distbins)[:,0]
    fig = plt.figure(figsize=(6,5.5))
    ax1 = plt.subplot2grid((7,4),(2,0),colspan = 4, rowspan =5)
    sigmarktyles = ['-',':','--']
    labels = ["chosen","other","previous"]
    for d in range(len(reg_weights.T)):
        plt.plot(xvals,reg_weights.T[d],\
            color='k',alpha=1,label=labels[d],ls=regStyles[d])
        plt.fill_between(xvals,reg_weights.T[d]+reg_sems.T[d],\
            reg_weights.T[d]-reg_sems.T[d],color='k',alpha=.4)
    ax=plt.gca()
    ax.invert_xaxis()
    plt.xticks(xvals)
    plt.xlabel("Hexes from port",fontsize='xx-large')
    plt.ylabel("Prior Reward ß",fontsize='xx-large')
    plt.grid(visible=True,lw=0.5)
    plt.legend()
    
    ax2 = plt.subplot2grid((7,4),(0,0),colspan = 4, rowspan =2,sharex=ax1)
    for i in range(3):
        plt.plot(xvals,propSig[:,i],color='k',ls=regStyles[i],lw=2)
    ax2.tick_params('x', labelbottom=False)
    plt.ylabel("Fraction\nsignificant")
    plt.grid(visible=True,lw=0.5)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.2)
    return fig

def get_rDf(hexData):
    return hexData.df.loc[(hexData.df.hexesFromPort>0),
    ['DA','port', 'rwd', 'session', 'block', 'trial', 'rat',
       'date', 'session_type', 'nom_rwd_chosen', 'vel',
       'pairedHexState', 'critCP', 'hexesFromCp', 'postCP', 'hexesFromFirstCp',
                'tri', 'hexesFromPort', 'rt-1',
       'samePath_t-1', 'rt-2', 'samePath_t-2', 'rt-3', 'samePath_t-3', 'rt-4',
       'samePath_t-4', 'rt-5', 'samePath_t-5', 'path_rt-1', 'path_rt-2',
       'path_rt-3', 'path_rt-4', 'path_rt-5', 'same2', 'fromport', 'fromport2',
    'otherPort_rt-1', 'lastRwd']]

def calcDaRhistRegInSpace(hexData,hexbin = 2,rwdsBack=[1,4],alongPath=True):
    rDf = hexData.df.loc[(hexData.df.notnull().all(axis=1))&(hexData.df.hexesFromPort>0),:]
    
    reg_weights = []
    distbins = [np.arange(i,i+hexbin) for i in np.arange(1,16,hexbin)]
    rlags = np.arange(rwdsBack[0],rwdsBack[1]+1)
    for r in hexData.df.rat.unique():
        regDf = rDf.loc[(rDf.rat==r),:]
        if len(regDf)==0:
            continue
        r_weights = []
        for dbin in distbins:
            rhistString = "path_rt-" if alongPath else "rt-"
            X = regDf.loc[regDf.hexesFromPort.isin(dbin),\
                          [rhistString+str(n) for n in rlags]].values
            if len(X)==0:
                break
            y = regDf.loc[regDf.hexesFromPort.isin(dbin),"DA"]
            mod = LR(fit_intercept=True,normalize=False).fit(X,y)
            r_weights.append(mod.coef_) #regression at specified distance
        reg_weights.append(r_weights)
    return reg_weights

def plotDaRhistRegInSpace(hexData,hexbin = 2,rwdsBack=[1,4],alongPath=True):
    #rDf = hexData.df.loc[(hexData.df.notnull().all(axis=1))&(hexData.df.hexesFromPort>0),:]
    rDf = get_rDf(hexData)
    rDf = rDf.loc[(rDf.notnull().all(axis=1)),:]
    
    reg_weights = []
    distbins = [np.arange(i,i+hexbin) for i in np.arange(1,16,hexbin)]
    rlags = np.arange(rwdsBack[0],rwdsBack[1]+1)
    for r in hexData.df.rat.unique():
        regDf = rDf.loc[(rDf.rat==r),:]
        if len(regDf)==0:
            continue
        r_weights = []
        for dbin in distbins:
            rhistString = "path_rt-" if alongPath else "rt-"
            X = regDf.loc[regDf.hexesFromPort.isin(dbin),\
                          [rhistString+str(n) for n in rlags]].values
            if len(X)==0:
                break
            y = regDf.loc[regDf.hexesFromPort.isin(dbin),"DA"]
            mod = LR(fit_intercept=True,normalize=False).fit(X,y)
            r_weights.append(mod.coef_) #regression at specified distance
        reg_weights.append(r_weights)
    reg_sems = sem(reg_weights)    
    reg_weights = np.mean(reg_weights,axis=0)
    
    xvals = np.array(distbins)[:,0]#[str(d[0])+'-'+str(d[0]+hexbin) for d in distbins]
    fig = plt.figure(figsize=(5,3.5))
    for d in range(len(reg_weights.T)):
        plt.plot(xvals,reg_weights.T[d],marker='o',markersize=9,\
            label="t-"+str(d),color='k',alpha=1-((d)/len(reg_weights.T)))
        plt.errorbar(xvals,reg_weights.T[d],yerr=np.transpose(reg_sems)[d],\
            ecolor='k',alpha=.4,lw=.1,elinewidth=2)
    ax=plt.gca()
    ax.invert_xaxis()
    #plt.xticks(xvals,["0-2","3-5","6-8","8-10","11-13"])
    plt.xticks(xvals)
    #plt.grid("on")
    plt.xlabel("Hexes from port",fontsize='xx-large',fontweight='bold')
    plt.ylabel("Prior Reward "+r"$\beta $",fontsize='xx-large',fontweight='bold')
    titleString = 'along same path' if alongPath else "along any path"
    plt.title(titleString,fontsize='x-large')
    plt.tight_layout()
    return fig

def plot_hexRampByPortVal(photrats,stHex=15,secondHalfOnly=True,
                          poolFactor="nom_rwd_chosen",useRatGroupLevel=True):
    photrats.fs=1
    photrats.set_plot_window([-stHex,-1])
    photrats.set_pool_factor(poolFactor)
    photrats.set_plot_trace("DA")
    fig = plt.figure(figsize = (7,5))
    rwdhigh,omhigh,rwdmid,ommid,rwdlow,omlow = photrats.getSessionTercMeans(\
        secondHalf=secondHalfOnly,useRat=useRatGroupLevel)
    # to plot grouping by hexesFromPort, need to rewrite function
        
    xvals = np.arange(-stHex+1,1)
    ax1 = plt.gca()
    plot_avgWithSem(photrats,ax1,xvals,np.vstack((rwdhigh,omhigh)),'lightgrey','-',\
                    [None,None],"high")
    plot_avgWithSem(photrats,ax1,xvals,np.vstack((rwdmid,ommid)),'darkgrey','-',\
                    [None,None],"mid")
    ax1.set_xlabel('time (s) from port entry')
    plot_avgWithSem(photrats,ax1,xvals,np.vstack((rwdlow,omlow)),'dimgrey','-',\
                    [None,None],"low")
    ax1.axvline(x=0.0,ymin=-.1,ymax=1.0,color='k',linestyle='--')
    ax1.legend()
    plt.xlabel("Distance to port (hexes)",fontsize='xx-large')
    plt.ylabel("DA (z-scored)",fontsize='xx-large')
    plt.xticks(xvals)
    plt.grid()
    plt.tight_layout()
    return fig

def plot_avgWithSem(photratz,ax,xvals,plot_trace,colorString,linstyle='-',subset=[None,None],traceLabel=None):
    ax.plot(xvals[subset[0]:subset[1]],np.nanmean\
        (plot_trace,axis=0)[subset[0]:subset[1]],color=colorString,ls=linstyle,label=traceLabel)
    ax.fill_between(xvals[subset[0]:subset[1]],(np.nanmean\
        (plot_trace,axis=0)+sem(plot_trace,nan_policy="omit"))\
        [subset[0]:subset[1]],(np.nanmean\
        (plot_trace,axis=0)-sem(plot_trace,nan_policy="omit"))\
        [subset[0]:subset[1]],color=colorString,alpha=.5)

def plot_deviationFromOptLen(hexData,binSize=5,binByMedian=True,laterProb=False,laterBar=False):
    bins = [i for i in range(0,71,binSize)]#[0,10,20,30,40,50,60,70]
    dfopt_means = pd.DataFrame(index=[str(i)+':'+str(i+binSize) for i in bins[:-1]])
    for s in hexData.df.session.unique():
        if laterProb:
            titleString = "Blocks > 1, prob variants only."
            toplt = hexData.df.loc[(hexData.df.session==s)&(hexData.df.port!=-100)&(hexData.df.session_type=="prob")&\
                              (hexData.df.block!=1)&(hexData.df.d_from_opt_len>=0)].groupby('tri').d_from_opt_len.mean()
        elif laterBar:
            titleString = "Blocks > 1, bar variants only."
            toplt = hexData.df.loc[(hexData.df.session==s)&(hexData.df.port!=-100)&(hexData.df.session_type=="barrier")&\
                              (hexData.df.block!=1)&(hexData.df.d_from_opt_len>=0)].groupby('tri').d_from_opt_len.mean()
        else:
            titleString = "First block only, all sessions"
            toplt = hexData.df.loc[(hexData.df.session==s)&(hexData.df.port!=-100)&\
                              (hexData.df.block==1)&(hexData.df.d_from_opt_len>=0)].groupby('tri').d_from_opt_len.mean()
        toplt = pd.DataFrame(toplt)
        toplt['bin'] = pd.cut(toplt.index,bins=bins)
        if binByMedian:
            dfopt_means[s]=toplt.groupby('bin').median().values.T[0]
        else:
            dfopt_means[s]=toplt.groupby('bin').mean().values.T[0]
    fig = plt.figure()
    sns.pointplot(data=dfopt_means.T,ci=95,estimator=np.median)
    #sns.stripplot(data=dfopt_means.T,color='darkred')
    plt.ylabel('Deviation from optimal\npath length (hexes)',fontsize='xx-large')
    plt.xlabel('Trial in block',fontsize='xx-large')
    plt.title('Median and individual points from '+str(len(hexData.df.session.unique()))+' sessions.\n'+titleString\
             ,fontsize='large')
    plt.tight_layout()
    return fig

def create_sameValtPreCpDADict(hexData):
    ratSameValtPreCpDA = {r:{'DARwdVsOmSame':[],'DARwdVsOmAlt':[]} 
    for r in hexData.df.loc[:,"rat"].unique()}
    ratNsPreCpDa = {r:{'DARwdVsOmSameRwd':[],'DARwdVsOmAltRwd':[],\
                'DARwdVsOmSameOm':[],'DARwdVsOmAltOm':[]} 
                for r in hexData.df.loc[:,"rat"].unique()}
    for rat in hexData.df.rat.unique():
        for s in hexData.df.loc[hexData.df.rat==rat,"session"].unique():
            dat = hexData.df.loc[(hexData.df.session==s)\
                                 &(hexData.df.same2==1)&(hexData.df.preCP==1),:]
            ratSameValtPreCpDA[rat]['DARwdVsOmSame'].append(\
               dat.loc[(dat["rt-1"]==1)&(dat["samePath_t-1"]==1),"DA"].mean()-\
                 dat.loc[(dat["rt-1"]==0)&(dat["samePath_t-1"]==1),"DA"].mean())
            ratSameValtPreCpDA[rat]['DARwdVsOmAlt'].append(\
               dat.loc[(dat["rt-1"]==1)&(dat["samePath_t-1"]==0),"DA"].mean()-\
                 dat.loc[(dat["rt-1"]==0)&(dat["samePath_t-1"]==0),"DA"].mean())
            ratNsPreCpDa[rat]['DARwdVsOmSameRwd'].append(dat.loc[(dat["rt-1"]==1)&\
                (dat["samePath_t-1"]==1),"DA"].shape[0])
            ratNsPreCpDa[rat]['DARwdVsOmSameOm'].append(dat.loc[(dat["rt-1"]==0)&\
                (dat["samePath_t-1"]==1),"DA"].shape[0])
            ratNsPreCpDa[rat]['DARwdVsOmAltRwd'].append(dat.loc[(dat["rt-1"]==1)&\
                (dat["samePath_t-1"]==0),"DA"].shape[0])
            ratNsPreCpDa[rat]['DARwdVsOmAltOm'].append(dat.loc[(dat["rt-1"]==0)&\
                (dat["samePath_t-1"]==0),"DA"].shape[0])
    return ratSameValtPreCpDA,ratNsPreCpDa

def create_sameValtPreCpDARatMeans(hexData):
    ratSameValtPreCpDA = {r:{'DARwdVsOmSame':[],'DARwdVsOmAlt':[]}
     for r in hexData.df.loc[:,"rat"].unique()}
    for rat in hexData.df.rat.unique():
        dat = hexData.df.loc[(hexData.df.rat==rat)\
                             &(hexData.df.same2==1)&(hexData.df.preCP==1),:]
        ratSameValtPreCpDA[rat]['DARwdVsOmSame'] = \
           dat.loc[(dat["rt-1"]==1)&(dat["samePath_t-1"]==1),"DA"].mean()-\
            dat.loc[(dat["rt-1"]==0)&(dat["samePath_t-1"]==1),"DA"].mean()
        ratSameValtPreCpDA[rat]['DARwdVsOmAlt'] = \
           dat.loc[(dat["rt-1"]==1)&(dat["samePath_t-1"]==0),"DA"].mean()-\
            dat.loc[(dat["rt-1"]==0)&(dat["samePath_t-1"]==0),"DA"].mean()
    return ratSameValtPreCpDA


def plot_rampSequenceRatAvg(hexData,rPat=[0,1],altPath=False,alongPath=False)
    alphaVal = 0.05
    entriesByRat = np.zeros((len(rPat),9,16))
    totEntries = 0
    for r in range(len(hexData.df.loc[:,"rat"].unique())):
        rat = hexData.df.loc[:,"rat"].unique()[r]
        entries = get_entriesByRwdSeq(hexData.df.loc[hexData.df.rat==rat,:],\
                                      rwd_pattern=rPat,alongPath=alongPath,altPath=altPath,getPriorTri=False)
        totEntries += np.shape(entries)[1]
        entriesByRat[:,r,:] = np.mean(entries,axis=1)
    fig = plot_traceByRwdSeq(rPat,entriesByRat,pltError=True,fsize=[6.5,5],ylim=[-.197,.51])
    entryDifs = calcDifInTracesByRwdSeq(entriesByRat)
    if rPat[-1]>0:
        sigInds = np.arange(-15,1)[calcSigDifVs0Vector_usingWSR(entryDifs[0].T,greater=True)<alphaVal]
    else:
        sigInds = np.arange(-15,1)[calcSigDifVs0Vector_usingWSR(entryDifs[0].T,less=True)<alphaVal]
    for i in sigInds:
        plt.text(x=i-.1, y=0.35, s='*',fontweight='bold',fontsize='x-large')
    if altPath:
        pathSpecString = "altPathOnly"
    elif alongPath:
        pathSpecString = "samePathOnly"
    else:
        pathSpecString = "anyPath"
    plt.suptitle(str(rPat)+" sequence averaged over rats\n"+pathSpecString+"; n = "+str(totEntries)+" events")
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0)
    print(f"n = {totEntries}")
    return fig

def plot_DAdifSameValt(hexData,ratSameValtPreCpDA):
    fig = plt.figure(figsize=(4.5,5))
    ratmeans_samePath = []
    ratmeans_altPath = []
    for rat in ratSameValtPreCpDA:
        ratmeans_samePath.append(np.mean(ratSameValtPreCpDA[rat]["DARwdVsOmSame"]))
        ratmeans_altPath.append(np.mean(ratSameValtPreCpDA[rat]["DARwdVsOmAlt"]))
        plt.scatter(x=["Same\nPath", "Alternative\nPath"],\
            y=[ratmeans_samePath[-1],ratmeans_altPath[-1]],color='k',s=30,marker='D',alpha=1)
    sns.barplot(x=["Same\nPath", "Alternative\nPath"],\
        y=[np.mean(ratmeans_samePath),np.mean(ratmeans_altPath)],\
        palette=["#ff7f0e","#2ca02c"],alpha=.5)
    plt.ylabel("$\Delta$ DA before\nchoice point",fontsize=20)
    plt.suptitle("DA change on next run\nfollowing reward - omission",fontsize='xx-large')
    plt.tight_layout()
    plt.xticks(fontsize=20)
    sameSigRats = get_sigRats_fromMeanList(ratmeans_samePath,altString="greater")
    plot_sigMarkers(sameSigRats,0,yval=0.17)
    altSigRats = get_sigRats_fromMeanList(ratmeans_altPath,altString="greater")
    plot_sigMarkers(altSigRats,.95,yval=0.17)
    plt.axhline(y=0,ls='--',color='k',lw=1)
    plt.ylim(-.23,0.23)
    return fig
    
def get_sigRats_fromMeanList(meanList,altString="two-sided"):
    pval = wilcoxon(meanList,alternative=altString)[1]
    print(pval)
    return pval<0.001,pval<0.01,pval<0.05

def add_prevRwdTakenAndAlt(hexData):
    '''add columns indicating whether rat was rewarded
    when previously taking the same (taken) or alternative path
    to the currently chosen port.'''
    hexData.df.loc[:,["prev_rwd_taken","prev_rwd_alt"]] = np.nan
    hexData.df.loc[(hexData.df.port!=-100),["prev_rwd_taken","prev_rwd_alt"]] = -100
    hexData.df.loc[(hexData.df["samePath_t-1"]==0)&(hexData.df.port!=-100),"prev_rwd_alt"]\
     = hexData.tridf.loc[(hexData.tridf["samePath_t-1"]==0),"rwd"].values
    hexData.df.loc[(hexData.df["samePath_t-1"]==1)&(hexData.df.port!=-100),"prev_rwd_taken"]\
     = hexData.tridf.loc[(hexData.tridf["samePath_t-1"]==1),"rwd"].values
    hexData.df.loc[:,"prev_rwd_alt"] = hexData.df.loc[:,"prev_rwd_alt"].fillna(method='bfill')
    hexData.df.loc[:,"prev_rwd_taken"] = hexData.df.loc[:,"prev_rwd_taken"].fillna(method='bfill')
    hexData.df.loc[hexData.df.prev_rwd_alt==-100,"prev_rwd_alt"] = np.nan
    hexData.df.loc[hexData.df.prev_rwd_taken==-100,"prev_rwd_taken"] = np.nan


def plot_portQevolutionOneSession(hexData,s):
    tridat = hexData.df.loc[(hexData.df.session==s)&(hexData.df.port!=-100),:].copy()
    tridat.reset_index(inplace=True)
    seshQs = tridat.loc[:,["Q_a","Q_b","Q_c"]].values
    fig = plt.figure(figsize=(10,5))
    x1 = np.arange(0,len(tridat))
    A1 = tridat.rwd.loc[tridat.port==0] + 1.5
    yA1 = np.zeros(len(tridat))
    yA1[A1.index.values] = A1
    B1 = tridat.rwd.loc[tridat.port==1] + 1.5
    yB1 = np.zeros(len(tridat))
    yB1[B1.index.values] = B1
    C1 = tridat.rwd.loc[tridat.port==2] + 1.5
    yC1 = np.zeros(len(tridat))
    yC1[C1.index.values] = C1
    ax1 = plt.subplot2grid((10,1),(0,0),colspan = 1, rowspan = 1)
    ax1.bar(x1,yA1,color = '#1f77b4')
    ax1.axis('off')
    ax2 = plt.subplot2grid((10,1),(1,0),colspan = 1, rowspan = 1,sharex=ax1)
    ax2.bar(x1,yB1,color = '#ff7f0e')
    ax2.axis('off')
    ax3 = plt.subplot2grid((10,1),(2,0),colspan = 1, rowspan = 1,sharex=ax1)
    ax3.bar(x1,yC1,color = '#2ca02c')
    ax3.axis('off')
    ax4 = plt.subplot2grid((10,1),(3,0),colspan = 1, rowspan = 7,sharex=ax1)
    ax4.plot(x1,seshQs[:,0],label = "port A",color="#1f77b4")
    ax4.plot(x1,seshQs[:,1],label = "port B",color="#ff7f0e")
    ax4.plot(x1,seshQs[:,2],label = "port C",color="#2ca02c")
    plt.xlabel("Trial",fontsize=20,fontweight="bold")
    plt.ylabel("Q Value",fontsize=20,fontweight="bold")
    plt.legend()
    plt.tight_layout()
    return fig