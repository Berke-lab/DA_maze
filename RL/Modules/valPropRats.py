#!/usr/bin/env python3

"""Functions to plot value and DA across maze.
Class with methods to run value iteration and optimize
gamma values to maximize fit to DA."""

__author__ = "Tim Krausz"
__email__ = "krausz.tim@gmail.com"
__status__ = "development"

import scipy.optimize as so
import matplotlib as mpl
from __main__ import *

norm = mpl.colors.Normalize(vmin=-1, vmax=2.5)
DAmap = cm.ScalarMappable(norm=norm,cmap="viridis")
def plot_hex_outline(photrats,sesh,block,ax,size='sm'):
    sz = 300 if size=='sm' else 1000
    if photrats.df.loc[photrats.df.session==sesh,'session_type'].unique()[0]=='barrier':
        bardf = centdf.drop(photrats.sesh_barIDs[sesh][block-1]+1,axis=0)
    else:
        bardf = centdf.drop(photrats.sesh_barIDs[sesh]+1,axis=0)
    ax.scatter(bardf[0].values,bardf[1].values,marker='H',color='grey',\
        edgecolors= "darkgrey",alpha=0.3,s=sz)
    
def plot_arrowMapFromStates(statelist,arrowvals,plotDA=False):
    if plotDA:
        plot_arrowMap(arrowdf.loc[np.isin(arrowdf.statecodes,statelist)],arrowvals,mapToUse=DAmap.to_rgba)
    else:
        plot_arrowMap(arrowdf.loc[np.isin(arrowdf.statecodes,statelist)],arrowvals)
        
def viz_value_map_pairedState(photrats,valmap,i):
    plt.clf()
    plt.suptitle('trial '+str(i)+' value map',fontsize='xx-large',fontweight='bold')
    plot_arrowMapFromStates(np.arange(0,126),arrowvals=valmap)

def plot_arrowMap(arrowdf,arrowvals,mapToUse=valcols):
    #plt.figure(figsize=(10,10))
    plt.xlim(-1,15)
    plt.ylim(-1,22)
    ax = plt.gca()
    visited = []
    cind = 0
    for i in arrowdf.sort_values("statecodes").index:
        if np.isnan(arrowvals[cind]):
            cind +=1
            continue
        state = arrowdf.loc[i,["arrow_start_x","arrow_start_y","arrow_end_x","arrow_end_y"]].values#/100
        jit =  .2 if str(arrowdf.loc[i,["to_state","from_state"]].values) in visited else 0
        if arrowdf.loc[i,"to_state"] in [2,3,47,]:
            jit = -.2
        plt.arrow(x=state[0]+jit, y=state[1], dx=state[2]-state[0],\
         dy=state[3]-state[1],fc=mapToUse(arrowvals[cind]),#valcols(arrowvals[cind]),
        head_width=.6, head_length=.4,width = .25,alpha=0.9,shape='full',\
        length_includes_head=True,ec=mapToUse(arrowvals[cind]))#valcols(arrowvals[cind]))
        visited.append(str(arrowdf.loc[i,["from_state","to_state"]].values))
        cind +=1

def plot_pairedStateValsVsDA(s,b,t,onlyPointing2port=False):
    portz = ['A','B','C']
    plt.clf()
    trial = vrats.valDaDf.loc[(vrats.valDaDf.session==s)&(vrats.valDaDf.block==b)\
                              &(vrats.valDaDf.tri==t),"trial"].values[0]
    nextp = vrats.df.loc[(vrats.df.session==s)&(vrats.df.block==b)&\
                    (vrats.df.tri==t)&(vrats.df.port!=-100),"port"].values[0]
    fromp = vrats.df.loc[(vrats.df.session==s)&(vrats.df.block==b)&\
                                (vrats.df.tri==t-1),"port"].values[0]
    toPort = portz[nextp]
    if onlyPointing2port:
        if len(vrats.sesh_arrows2goal['to'+toPort][s])==126:
            arrows2include = np.where(vrats.sesh_arrows2goal['to'+toPort][s])[0]
        else:
            arrows2include = np.where(vrats.sesh_arrows2goal['to'+toPort][s][b-1])[0]
    else:
        arrows2include = np.arange(126)
    meanHexDa = vrats.valDaDf.loc[(vrats.valDaDf.session==s)&\
                 (vrats.valDaDf.block==b)&(vrats.valDaDf.tri==t)&\
                  (vrats.valDaDf.pairedHexState.isin(arrows2include)),:].groupby(\
                    "pairedHexState").DA.mean()
    
    statevals = vrats.valDaDf.loc[(vrats.valDaDf.session==s)&\
                 (vrats.valDaDf.block==b)&(vrats.valDaDf.tri==t)&\
                  (vrats.valDaDf.pairedHexState.isin(arrows2include)),:].groupby(\
                    "pairedHexState").Value.mean()
    
    statelist = meanHexDa.index.values
    stateDa = meanHexDa.values
    plt.subplot(121)
    plt.title("value map")
    plot_hex_outline(vrats,sesh=s,block=b,ax=plt.gca(),size='sm')
    plot_arrowMapFromStates(statelist,statevals.values)
    plt.axis("off")
    
    plt.subplot(122)
    plt.title("DA")
    plot_hex_outline(vrats,sesh=s,block=b,ax=plt.gca(),size='sm')
    plot_arrowMapFromStates(statelist,stateDa,plotDA=True)
    print("v(chosen) = ",vrats.df.loc[(vrats.df.session==s)&\
                                (vrats.df.tri==t)&(vrats.df.block==b)&\
                                (vrats.df.port!=-100),"nom_rwd_chosen"].values[0])
    plt.axis("off")

class ValIterRats(PhotRats):
    
    def __init__(self,seshdict,use_hexByDirection=True,use_offset=False):
        super().__init__(seshdict)
        self.use_hexByDirection = use_hexByDirection
        self.use_offset = use_offset
        if self.use_hexByDirection:
            self.nstates = 126
        
    def get_optValIterGams(self):
        self.opt_gams = []
        seshs = self.df.session.unique()
        for s in tqdm(range(len(seshs))):
            sesh = seshs[s]
            seshGam = self.calc_optimalValIterGamma4sesh(sesh)
            self.opt_gams.append(0.5 + 0.5*erf(seshGam/np.sqrt(2)))
    
    def calc_optimalValIterGamma4sesh(self,sesh,startparam=1):
        '''for given session, maximize p(DA|gamma)'''
        lossfunc = lambda gam: self.calc_gammaLoss4ValIterOpt(sesh,gam)
        est = so.minimize_scalar(lossfunc)#so.minimize(lossfunc,startparam,method = "L-BFGS-B")
        return est.x
    
    def calc_gammaLoss4ValIterOpt(self,sesh,param,portValType="chosen"):
        #squeeze param bw 0 and 1
        gam = 0.5 + 0.5*erf(param/np.sqrt(2))
        valsByTri = self.calc_valIterValsByTri(sesh,gam,portValType,"nan")
        nll = self.calc_nllBwValAndDA(valsByTri,sesh)
        return nll
    
    def calc_nllBwValAndDA(self,valdat,sesh):
        seshDA = np.array([])
        seshVals = np.array([])
        vinds = self.df.loc[(self.df.session==sesh)&(self.df.port!=-100)].index
        for t in range(1,len(vinds)):
            tdat = self.df.loc[vinds[t-1]+1:vinds[t]]
            if len(tdat) < 2:
                continue
            seshDA = np.concatenate([seshDA,tdat.loc[tdat.pairedHexState!=-1,"DA"].values])
            stateInds = tdat.loc[tdat.pairedHexState!=-1,"pairedHexState"].values
            seshVals = np.concatenate([seshVals,valdat[t,stateInds]])
        X = np.ones((len(seshVals),2))
        X[:,1] = seshVals
        Y = seshDA[~np.isnan(X[:,1])]
        X = X[~np.isnan(X[:,1]),:]
        betas = np.dot(np.linalg.inv(np.dot(X.T,X)),np.dot(X.T,Y))
        std = np.sqrt(np.mean((Y-(betas[0]+betas[1]*X[:,1]))**2))
        loglik = -len(X[:,1])*np.log(np.sqrt(2*np.pi)*std) - len(X[:,1])/2
        return -loglik
    
    def calc_valIterValsByTri(self,sesh,gam,portValType="chosen",handleBars="zero"):
        for b in self.df.loc[self.df.session==sesh,'block'].unique():
            if self.df.loc[self.df.session==sesh,'session_type'].values[0]=='prob':
                self.tmat = self.sesh_tmats[sesh]
                self.bars = self.sesh_barIDs[sesh]
            else:
                try:
                    self.tmat = self.sesh_tmats[sesh][int(b-1)]
                except:
                    continue
                self.bars = self.sesh_barIDs[sesh][int(b-1)]
            self.hexdata = self.df.loc[(self.df.session==sesh)&(self.df.block==b)]
            vals = self.calcValIter(self.hexdata.loc[self.hexdata.port!=-100],\
                gam,portValType=portValType,handleBars=handleBars)
            if b>1:
                valmap = np.vstack([valmap,vals])
            else:
                valmap = vals
        return valmap
    
    def calcValIter(self,data,gam,portValType="chosen",handleBars="zero"):
        actions,V,portQs,portz = self.get_vars4ValIter(data)
        porthexes = [0,1,2]
        self.add_bars2tmatrix_pairedState()
        for t in range(1,len(data)):
            rwdvec = np.zeros(126)
            if portValType=="chosen":
                avail = np.array([portz[data.iloc[t].port]])
                availPortHexes = np.array([porthexes[data.iloc[t].port]])
                avail = np.array([portz[availPortHexes[0]]])
            else:
                avail = np.delete(portz,data.iloc[t-1].port)
                availPortHexes = np.setdiff1d(porthexes,data.iloc[t-1].port)        
            rwdvec[avail] = portQs[t,availPortHexes]
            i = 0
            oldV = V[t,:]
            while True:
                allVs = []
                for a in actions:
                    allVs.append(gam*np.dot(oldV,self.bar_tmatrix[:,a].T)+rwdvec)
                newV = np.amax(allVs,axis=0)
                if np.max(newV-oldV)<10e-4:
                    V[t] = newV
                    V[t,self.paired_barstates]=0 if handleBars=="zero" else np.nan
                    break
                oldV = newV
                i += 1
        return V
    
    def get_vars4ValIter(self,data):
        actions = [0,1]
        V = np.full((len(data),126),0.0,dtype="float16")
        portQs = data.loc[:,['Q_a', 'Q_b', 'Q_c']].values
        portz = [phexdf.loc[phexdf.to_state==1,'statecodes'].values[0],
        phexdf.loc[phexdf.to_state==2,'statecodes'].values[0],
         phexdf.loc[phexdf.to_state==3,'statecodes'].values[0]]
        return actions,V,portQs,portz

    def create_optValDf(self,portValType="chosen"):
        seshs = self.df.session.unique()
        valDaArray = np.array([[]for i in range(9)]).T
        for s in range(len(self.opt_gams)):
            sesh = seshs[s]
            gam = self.opt_gams[s]
            valsByTri = self.calc_valIterValsByTri(sesh,gam,portValType,handleBars="nan")
            valDaArray = np.vstack([valDaArray,self.get_seshValVsDa(valsByTri,sesh)])
        self.valDaDf = pd.DataFrame(valDaArray,columns = ["DA",\
            "vel","pairedHexState","tri","trial","block","session","Value","rat"])
        self.valDaDf.loc[:,"session"] = self.valDaDf.session.astype(float).astype("uint8")
        self.valDaDf.loc[:,"block"] = self.valDaDf.block.astype(float).astype("uint8")
        self.valDaDf.loc[:,"tri"] = self.valDaDf.tri.astype(float).astype("uint8")
        self.valDaDf.loc[:,"trial"] = self.valDaDf.trial.astype(float).astype("int16")
        self.valDaDf.loc[:,"DA"] = self.valDaDf.DA.astype(float).astype("float16")
        self.valDaDf.loc[:,"Value"] = self.valDaDf.Value.astype(float).astype("float16")
        self.valDaDf.loc[:,"vel"] = self.valDaDf.vel.astype(float).astype("float16")
        self.valDaDf.loc[:,"pairedHexState"] = \
        self.valDaDf.pairedHexState.astype(float).astype("uint8")

    def get_seshValVsDa(self,valdat,sesh):
        seshVals = np.array([])
        seshInfo = np.array([[]for i in range(6)])
        vinds = self.df.loc[(self.df.session==sesh)&(self.df.port!=-100)].index
        rat = self.df.loc[(self.df.session==sesh),"rat"].values[0]
        for t in range(1,len(vinds)):
            tdat = self.df.loc[vinds[t-1]+1:vinds[t]]
            if len(tdat) < 2:
                continue
            seshInfo = np.hstack([seshInfo,\
                        tdat.loc[tdat.pairedHexState!=-1,\
                        ["DA","vel","pairedHexState","tri","trial","block"]].values.T])
            stateInds = tdat.loc[tdat.pairedHexState!=-1,"pairedHexState"].values
            seshVals = np.concatenate([seshVals,valdat[t,stateInds]])
        seshInfo = np.vstack([seshInfo,np.full((np.shape(seshInfo)[1]),sesh)])
        seshInfo = np.vstack([seshInfo,seshVals])
        valDaArray = np.vstack([seshInfo,np.full((np.shape(seshInfo)[1]),rat)])
        return valDaArray.T

    def plot_nllCurve(self,s):
        plt.figure()
        for i in np.linspace(-1,3):
            plt.scatter(i,self.calc_gammaLoss4ValIterOpt(s,i),color='k')
        plt.title("session"+str(s))
        plt.xticks(np.linspace(-1,3,20),np.round(0.5 + 0.5*erf(np.linspace(-1,3,5)/np.sqrt(2)),decimals=2))
        plt.xlabel("gamma",fontsize='x-large',fontweight='bold')
        plt.ylabel("-log likelihood",fontsize='x-large',fontweight='bold')
        plt.tight_layout()