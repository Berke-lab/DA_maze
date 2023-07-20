__author__ = "Tim Krausz"
__email__ = "krausz.tim@gmail.com"
__status__ = "development"
__license__ = "MIT"

from sklearn.linear_model import LogisticRegression as LogReg
import statsmodels.api as sm
from scipy import stats
from __main__ import *

def add_prevRwdTakenAndAlt(photrats):
    '''add columns indicating whether rat was rewarded
    when previously taking the same (taken) or alternative path
    to the currently chosen port.'''
    photrats.triframe.loc[:,["prev_rwd_taken","prev_rwd_alt"]] = np.nan
    photrats.triframe.loc[:,["prev_rwd_taken","prev_rwd_alt"]] = -100
    photrats.triframe.loc[(photrats.triframe["samePath_t-1"]==0),"prev_rwd_alt"] = \
        photrats.triframe.loc[(photrats.triframe["samePath_t-1"]==0),"rwd"].values
    photrats.triframe.loc[(photrats.triframe["samePath_t-1"]==1),"prev_rwd_taken"] = \
        photrats.triframe.loc[(photrats.triframe["samePath_t-1"]==1),"rwd"].values
    photrats.triframe.loc[photrats.triframe.prev_rwd_alt==-100,"prev_rwd_alt"] = np.nan
    photrats.triframe.loc[photrats.triframe.prev_rwd_taken==-100,"prev_rwd_taken"] = np.nan

def get_back2samePortTrials(photrats):
    photrats.triframe.loc[:,"same2"] = 0
    photrats.triframe.loc[photrats.triframe.port==photrats.triframe.fromport2,"same2"] = 1
    
def get_portHist(photrats):
    photrats.triframe.loc[:,"fromport"] = photrats.triframe.loc[:,"port"].shift(1).values
    photrats.triframe.loc[:,"fromport2"] = photrats.triframe.loc[:,"port"].shift(2).values
    
def add_rewardLagColumns2triframe(photrats,lags=1):
    for l in range(1,lags+1):
        photrats.triframe.loc[:,"rt-"+str(l)] = np.nan
        photrats.triframe.loc[:,"samePath_t-"+str(l)] = 0
        for p in range(3):
            lagInds = photrats.triframe.loc[photrats.triframe.port==p,:].index
            photrats.triframe.loc[lagInds,"rt-"+str(l)] = photrats.triframe.loc[lagInds,'rwd'].shift(l).values
            photrats.triframe.loc[lagInds,"samePath_t-"+str(l)] = (photrats.triframe.loc[lagInds,\
                'fromport']==photrats.triframe.loc[lagInds,'fromport'].shift(l)).astype("int8")

def add_prevRwdTakenAndAlt(photrats):
    '''add columns indicating whether rat was rewarded
    when previously taking the same (taken) or alternative path
    to the currently chosen port.'''
    photrats.triframe.loc[:,["prev_rwd_taken","prev_rwd_alt"]] = np.nan
    photrats.triframe.loc[:,["prev_rwd_taken","prev_rwd_alt"]] = -100
    photrats.triframe.loc[(photrats.triframe["samePath_t-1"]==0),"prev_rwd_alt"] = \
        photrats.triframe.loc[(photrats.triframe["samePath_t-1"]==0),"rwd"].values
    photrats.triframe.loc[(photrats.triframe["samePath_t-1"]==1),"prev_rwd_taken"] = \
        photrats.triframe.loc[(photrats.triframe["samePath_t-1"]==1),"rwd"].values

def add_rewardLagColumns2triframeByPort(photrats,lags=5):
    portStrings = ["a","b","c"]
    for l in range(1,lags+1):
        photrats.triframe.loc[:,["rt-"+str(l)+"_"+s for s in portStrings]] = np.nan
        for p in range(3):
            lagInds = photrats.triframe.loc[photrats.triframe.port==p,:].index
            photrats.triframe.loc[lagInds,"rt-"+str(l)+"_"+portStrings[p]] = photrats.triframe.loc[lagInds,'rwd'].shift(l).values
        for sesh in photrats.triframe.session.unique():
            photrats.triframe.loc[photrats.triframe.session==sesh,["rt-"+str(l)+"_"+s for s in portStrings]] = \
            photrats.triframe.loc[photrats.triframe.session==sesh,["rt-"+str(l)+"_"+s for s in portStrings]].fillna(method='ffill')
            
def add_samePathLagColumns2triframeByPort(photrats,lags=1):
    portStrings = ["a","b","c"]
    for l in range(1,lags+1):
        photrats.triframe.loc[:,["samePath_t-"+str(l)+"_"+s for s in portStrings]] = np.nan
        for p in range(3):
            lagInds = photrats.triframe.loc[photrats.triframe.port==p,:].index
            photrats.triframe.loc[lagInds,"samePath_t-"+str(l)+"_"+portStrings[p]] = (photrats.triframe.loc[lagInds,\
                'fromport']==photrats.triframe.loc[lagInds,'fromport'].shift(l)).astype("int8")
        for sesh in photrats.triframe.session.unique():
            photrats.triframe.loc[photrats.triframe.session==sesh,["samePath_t-"+str(l)+"_"+s for s in portStrings]] = \
            photrats.triframe.loc[photrats.triframe.session==sesh,["samePath_t-"+str(l)+"_"+s for s in portStrings]].fillna(method='ffill')
            
def add_leftLastColumn(photrats):
    '''adds column indicating the left port was visited last
    (two port entries prior to choice of interest). Make sure
    triframe index has been reset'''
    photrats.triframe.loc[:,"leftLast"] = np.nan
    # for each trial, identify left port, and ask whether prior entry == left
    leftChoices = [2,0,1]
    for sesh in photrats.triframe.session.unique():
        dat = photrats.triframe.loc[photrats.triframe.session==sesh,:]
        for i in dat.index[2:]:
            p = int(dat.loc[:,"port"].shift(1)[i])
            if dat.loc[:,"port"].shift(2)[i]==leftChoices[p]:
                photrats.triframe.loc[i,"leftLast"] = 1
            else:
                photrats.triframe.loc[i,"leftLast"] = 0

def create_sameValtChoiceDict(photrats):
    ratSameValtChoiceProbs = {r:{'pDiffRwdVsOmSame':[],'pDiffRwdVsOmAlt':[]}\
     for r in photrats.triframe.loc[:,"rat"].unique()}
    ratNs = {r:{'pDiffRwdVsOmSameRwd':[],'pDiffRwdVsOmAltRwd':[],\
                'pDiffRwdVsOmSameOm':[],'pDiffRwdVsOmAltOm':[]}\
                 for r in photrats.triframe.loc[:,"rat"].unique()}
    for rat in photrats.triframe.rat.unique():
        for s in photrats.triframe.loc[photrats.triframe.rat==rat,"session"].unique():
            dat = photrats.triframe.loc[(photrats.triframe.rat==rat)&\
                (photrats.triframe.session==s)&(photrats.triframe.leftLast==1),:]
            ratSameValtChoiceProbs[rat]['pDiffRwdVsOmSame'].append(\
               dat.loc[(dat["samePath_t-1_left"]==1)&(dat["rt-1_left"]==1),"lrchoice"].mean()\
               -dat.loc[(dat["samePath_t-1_left"]==1)&(dat["rt-1_left"]==0),"lrchoice"].mean())
            ratSameValtChoiceProbs[rat]['pDiffRwdVsOmAlt'].append(\
               dat.loc[(dat["samePath_t-1_left"]==0)&(dat["rt-1_left"]==1),"lrchoice"].mean()\
               -dat.loc[(dat["samePath_t-1_left"]==0)&(dat["rt-1_left"]==0),"lrchoice"].mean())
            ratNs[rat]['pDiffRwdVsOmSameRwd'].append(dat.loc[\
                (dat["samePath_t-1_left"]==1)&(dat["rt-1_left"]==1),].shape[0])
            ratNs[rat]['pDiffRwdVsOmSameOm'].append(dat.loc[\
                (dat["samePath_t-1_left"]==1)&(dat["rt-1_left"]==0),].shape[0])
            ratNs[rat]['pDiffRwdVsOmAltRwd'].append(dat.loc[\
                (dat["samePath_t-1_left"]==0)&(dat["rt-1_left"]==1),].shape[0])
            ratNs[rat]['pDiffRwdVsOmAltOm'].append(dat.loc[\
                (dat["samePath_t-1_left"]==0)&(dat["rt-1_left"]==0),].shape[0])
    return ratSameValtChoiceProbs,ratNs

def plotChoosePortSameVAlt(photrats,ratSameValtChoiceProbs):
    fig = plt.figure(figsize=(4,5))
    plt.axhline(y=0,ls='--',color='k',lw=1)
    ratmeans_samePath = []
    ratmeans_altPath = []
    for rat in ratSameValtChoiceProbs:
        ratmeans_samePath.append(np.nanmean(ratSameValtChoiceProbs[rat]["pDiffRwdVsOmSame"]))
        ratmeans_altPath.append(np.nanmean(ratSameValtChoiceProbs[rat]["pDiffRwdVsOmAlt"]))
        plt.scatter(x=["Same\nPath", "Alternative\nPath"],y=[ratmeans_samePath[-1],
                  ratmeans_altPath[-1]],color='k',s=50,marker='D')
    sns.barplot(x=["Same\nPath", "Alternative\nPath"],y=[np.mean(ratmeans_samePath),
                     np.mean(ratmeans_altPath)],palette=["#ff7f0e","#2ca02c"],alpha=.5)
    plt.ylabel("$\Delta$ P(choose port)",fontsize='xx-large')
    plt.suptitle("Port choice increase\nfollowing reward - omission",fontsize='xx-large')
    plt.tight_layout()
    plt.xticks()
    sameSigRats = get_sigRats_fromMeanList(ratmeans_samePath,altString="greater")
    plot_sigMarkers(sameSigRats,-.1,yval=0.31)
    altSigRats = get_sigRats_fromMeanList(ratmeans_altPath,altString="greater")
    plot_sigMarkers(altSigRats,.95,yval=0.31)
    return fig

def get_sigRats_fromMeanList(meanList,altString="two-sided"):
    stat,pval = wilcoxon(meanList,alternative=altString)
    print(stat,pval)
    return pval<0.001,pval<0.01,pval<0.05

def plot_choiceRegWeightsByRat_MF(photrats):
    regWeightsByRat = calc_choiceRegWeightsByRat(photrats)
    regWeightsByRat.loc[:,"rat"] = np.array(photrats.triframe.rat.unique())
    xvals = ["p(Reward)","Distance\ncost"]#["left bias","relative\np(reward)","relative\ndistance"]
    females = ["IM-1291","IM-1292","IM-1531","IM-1532","IM-1533"]
    fig = plt.figure(figsize=(3.55,6))
    plt.bar([0,1],np.mean(regWeightsByRat,axis=0)[1:],color='grey',alpha=.5)
    sns.stripplot(data=regWeightsByRat.loc[~regWeightsByRat.rat.isin(females),\
                                           ["relative p(R)","relative distance"]],color='k',size=10,marker='D',alpha=.9)
    sns.stripplot(data=regWeightsByRat.loc[regWeightsByRat.rat.isin(females),\
                                           ["relative p(R)","relative distance"]],color='deeppink',size=10,marker='D',alpha=.9)
    plt.axhline(y=0,ls='--',color='k',lw=1)
    plt.xticks([0,1],xvals,fontsize="large",fontstyle="italic")
    #plt.ylim(-10,2.2)
    plt.ylabel("Port Choice ß",fontsize='xx-large',fontweight='bold')
    plt.xlabel("Choice Feature",fontsize='xx-large',fontweight='bold')
    plt.tight_layout()
    fig.savefig(photrats.directory_prefix+"pChooseLregCoefsIndividualRats_1.pdf")

from scipy import stats
def run_sklearnLogRegWithPval(X,y):
    lm = LogReg(fit_intercept=True).fit(X,y)
    params = np.append(lm.intercept_,lm.coef_)
    predictions = lm.predict(X)
    new_X = np.append(np.ones((len(X),1)), X, axis=1)
    M_S_E = (sum((y-predictions)**2))/(len(new_X)-len(new_X[0]))
    v_b = M_S_E*(np.linalg.inv(np.dot(new_X.T,new_X)).diagonal())
    s_b = np.sqrt(v_b)
    t_b = params/ s_b
    p_val =[2*(1-stats.t.cdf(np.abs(i),(len(new_X)-len(new_X[0])))) for i in t_b]
    return np.concatenate(params),p_val

def run_smLogRegWithPval(X,y):
    X = np.hstack([np.ones(len(X)).reshape(-1,1),X])
    mod = sm.Logit(y, X).fit()
    coefs = mod.summary2().tables[1]['Coef.'].values
    pvals = mod.summary2().tables[1]["P>|z|"].values
    return coefs,pvals
    
def calc_choiceRegWeightsByRat(photrats):
    try:
        photrats.regdf
    except:
        print("choose left regDf not yet saved. Getting now...")
        get_log_pchoos_v_costNben(photrats)
        add_scaledVarsToRegDf(photrats)
    regWeights = np.zeros((len(photrats.triframe.rat.unique()),3))
    ratSigLevels = np.zeros((len(photrats.triframe.rat.unique()),3))
    for r in range(len(photrats.regdf.rat.unique())):
        regDf = photrats.regdf.loc[\
            (photrats.regdf.rat==photrats.regdf.rat.unique()[r])&(photrats.regdf.tri>25),\
            ["session","choose_L","rwdDifScaled","lengthDifScaled"]]
        regDf = regDf.loc[regDf.notnull().all(axis=1),:]
        y = regDf.loc[:,"choose_L"]
        X = regDf.loc[:,["rwdDifScaled","lengthDifScaled"]]
        betas,pvals = run_smLogRegWithPval(X,y)
        regWeights[r,:] = betas
        ratSigLevels[r,:] = pvals
        ratRegDf = pd.DataFrame(regWeights,columns=["intercept","relative p(R)","relative distance"])
        ratRegDf.loc[:,["intercept p value","relative p(R) p value","relative distance p value"]] = ratSigLevels
    return ratRegDf


def create_PortChoiceSameValtDf(wsdf_all):
    omSamePath = []
    rwdSamePath = []
    omAltPath = []
    rwdAltPath = []
    omSameSem = []
    rwdSameSem = []
    omAltSem = []
    rwdAltSem = []
    for r in wsdf_all.rat.unique():
        dat = wsdf_all.loc[(wsdf_all.rat==r)&(wsdf_all.same2==1),:]
        omSamePath.append(dat.loc[(dat["prev_rwd_taken"]==0),"stay"].mean())
        rwdSamePath.append(dat.loc[(dat["prev_rwd_taken"]==1),"stay"].mean())
        omAltPath.append(dat.loc[(dat["prev_rwd_alt"]==0),"stay"].mean())
        rwdAltPath.append(dat.loc[(dat["prev_rwd_alt"]==1),"stay"].mean())
        omSameSem.append(dat.loc[(dat["prev_rwd_taken"]==1),"stay"].sem())
        rwdSameSem.append(dat.loc[(dat["prev_rwd_taken"]==1),"stay"].sem())
        omAltSem.append(dat.loc[(dat["prev_rwd_alt"]==0),"stay"].sem())
        rwdAltSem.append(dat.loc[(dat["prev_rwd_alt"]==1),"stay"].sem())
    choiceSameValtDf = pd.DataFrame({"omSame":omSamePath,"rwdSame":rwdSamePath,"omAlt":omAltPath,"rwdAlt":rwdAltPath,\
                             "omSameSem":omSameSem,"rwdSameSem":rwdSameSem,"omAltSem":omAltSem,"rwdAltSem":rwdAltSem})
    return choiceSameValtDf

#todo: add marker difference for sig vs ns
def plot_choiceRegWeightsByRat(photrats):
    regWeightsByRat = calc_choiceRegWeightsByRat(photrats)
    xvals = ["Left\nbias","Reward","Dist"]#["left bias","relative\np(reward)","relative\ndistance"]
    fig = plt.figure(figsize=(3.55,4))
    plt.bar([0,1,2],np.mean(regWeightsByRat,axis=0),color='grey',alpha=.5)
    sigpoints = (regWeightsByRat.values[:,3:]<0.05).astype(int)
    for i in range(len(regWeightsByRat)):
        #markers = 
        plt.scatter(regWeightsByRat.columns[:3],regWeightsByRat)
    sns.stripplot(data=regWeightsByRat,color='k',size=5,marker='D',alpha=.9)
    plt.axhline(y=0,ls='--',color='k',lw=1)
    plt.xticks([0,1,2],xvals,fontsize="large",fontstyle="italic")
    #plt.ylim(-10,2.2)
    plt.ylabel("Port Choice ß",fontsize='xx-large',fontweight='bold')
    plt.xlabel("Choice Feature",fontsize='xx-large',fontweight='bold')
    plt.tight_layout()
    fig.savefig(photrats.directory_prefix+"pChooseLregCoefsIndividualRats.pdf")

def plot_choiceSigmoidsByRelativeCostsAndBen(photrats,saveFigs=True):
    fig = plot_pChooseVrelativeVar_byRat(photrats,use_scaled=False,use_prob=True)
    fig1 = plot_pChooseVrelativeVar_byRat(photrats,use_scaled=True,use_prob=True)
    fig2 = plot_pChooseVrelativeVar_byRat(photrats,use_scaled=False,use_prob=False)
    fig3 = plot_pChooseVrelativeVar_byRat(photrats,use_scaled=True,use_prob=False)
    ax = fig.axes[0]
    ax.set_xticks([-80,-40,0,40,80])
    ax.set_ylim(0.15,0.76)
    ax.set_yticks(np.linspace(.2,.7,6))
    ax = fig1.axes[0]
    ax.set_xticks([0,.5,1])
    ax.set_ylim(0.15,0.76)
    ax.set_yticks(np.linspace(.2,.7,6))
    ax = fig3.axes[0]
    ax.set_xticks([0,0.5,1])
    ax.set_yticks(np.linspace(0,1,5))
    ax = fig2.axes[0]
    ax.set_xticks([-10,-5,0,5,10])
    ax.set_yticks(np.linspace(0,1,5))
    if saveFigs:
        fig.savefig(photrats.directory_prefix+"rat_logPchooseVpRwd.pdf")
        fig1.savefig(photrats.directory_prefix+"rat_logPchooseVscaledPrwdDif.pdf")
        fig2.savefig(photrats.directory_prefix+"rat_logPchooseVlengthDif.pdf")
        fig3.savefig(photrats.directory_prefix+"rat_logPchooseVscaledLengthDif.pdf")

def get_sessionProbAndDistChangeInds(photrats):
    probseshs = photrats.triframe.loc[photrats.triframe.session_type=='prob','session'].unique()
    portz = ['a','b','c']
    high2lowProb = {}
    low2highProb = {}
    
    for s in probseshs:
        tdat = photrats.triframe.loc[photrats.triframe.session==s]
        if len(tdat.block.unique())<2:
            continue
        high2lowProb[s] = {b:np.array([0,0,0]) for b in tdat.block.unique()}
        low2highProb[s] = {b:np.array([0,0,0]) for b in tdat.block.unique()}    
        probChanges = np.diff(tdat.loc[tdat.block.diff()!=0,\
            ["nom_rwd_a","nom_rwd_b","nom_rwd_c"]].values,axis=0)
        probIncreasesByBlock = np.array(probChanges>0).astype(int)
        probDecreasesByBlock = np.array(probChanges<0).astype(int)
        for i in range(len(probDecreasesByBlock)):
            low2highProb[s][i+1] = probIncreasesByBlock[i]
            high2lowProb[s][i+1] = probDecreasesByBlock[i]
    
    barseshs = photrats.triframe.loc[photrats.triframe.session_type=='barrier','session'].unique()
    
    pathz = ["AB","AC","BC"]
    high2lowDist = {}
    low2highDist = {}
    for s in barseshs:
        tdat = photrats.triframe.loc[photrats.triframe.session==s]
        if len(tdat.block.unique())<2:
            continue
        high2lowDist[s] = {b:np.array([0,0,0]) for b in tdat.block.unique()}
        low2highDist[s] = {b:np.array([0,0,0]) for b in tdat.block.unique()}    
        distChanges = np.diff(tdat.loc[tdat.block.diff()!=0,\
            ["lenAC","lenBC","lenAB"]].values,axis=0)
        distIncreasesByBlock = np.array(distChanges>0).astype(int)
        distDecreasesByBlock = np.array(distChanges<0).astype(int)
        for i in range(len(distDecreasesByBlock)):
            low2highDist[s][i+1] = distIncreasesByBlock[i]
            high2lowDist[s][i+1] = distDecreasesByBlock[i]
    return high2lowProb,high2lowDist,low2highProb,low2highDist

def get_ProbDecreasePVisit(photrats,rat,high2lowProb):
    high2lowProbPVisit = []
    for s in high2lowProb:
        if photrats.triframe.loc[photrats.triframe.session==s,'rat'].values[0] != rat:
            continue
        tdat = photrats.triframe.loc[photrats.triframe.session==s,:]
        for b in high2lowProb[s]:
            if 1 not in high2lowProb[s][b]:
                continue        
            bports = np.where(high2lowProb[s][b]==1)[0]
            for bport in bports:
                bind = tdat.loc[tdat.block.diff()!=0,:].index[b]
                dat = photrats.triframe.loc[bind-49:bind+50,:].copy()
                dat.loc[:,'choosePort'] = 0
                dat.loc[dat.port==bport,'choosePort']=1
                #mark trials where port of interest is available for choice
                dat.loc[:,"portAvail"] = np.array(dat.port.shift(1)!=bport).astype(int)
                availdat = dat.loc[dat.portAvail==1,:]
                try:
                    blkTransInd = availdat.loc[availdat.block.diff()==1,:].index[0]
                except:
                    continue
                preTransSum = availdat.loc[:blkTransInd,"choosePort"].values
                postTransSum = availdat.loc[blkTransInd:,"choosePort"].values
                pvisitPort = pd.DataFrame({"chosePort":np.full(100,np.nan)})
                try:
                    pvisitPort.loc[blkTransInd-dat.index.min()-len(preTransSum):\
                               blkTransInd-dat.index.min()-1,"chosePort"] = preTransSum
                except:
                    print("didn't work for session "+str(s)+" block "+str(b))
                    continue
                pvisitPort.loc[blkTransInd-dat.index.min():blkTransInd-\
                               dat.index.min()+len(postTransSum)-1,"chosePort"] = postTransSum
                high2lowProbPVisit.append(pvisitPort.values.T[0])
    return np.mean(high2lowProbPVisit,axis=0),len(high2lowProbPVisit)

def get_ProbIncreasePVisit(photrats,rat,low2highProb):
    low2highProbPVisit = []
    for s in low2highProb:
        if photrats.triframe.loc[photrats.triframe.session==s,'rat'].values[0] != rat:
            continue
        tdat = photrats.triframe.loc[photrats.triframe.session==s,:]
        for b in low2highProb[s]:
            if 1 not in low2highProb[s][b]:
                continue        
            bports = np.where(low2highProb[s][b]==1)[0]
            for bport in bports:
                bind = tdat.loc[tdat.block.diff()!=0,:].index[b]
                dat = photrats.triframe.loc[bind-49:bind+50,:].copy()
                dat.loc[:,'choosePort'] = 0
                dat.loc[dat.port==bport,'choosePort']=1
                dat.loc[:,"portAvail"] = np.array(dat.port.shift(1)!=bport).astype(int)
                availdat = dat.loc[dat.portAvail==1,:]
                try:
                    blkTransInd = availdat.loc[availdat.block.diff()==1,:].index[0]
                except:
                    continue
                preTransSum = availdat.loc[:blkTransInd,"choosePort"].values
                postTransSum = availdat.loc[blkTransInd:,"choosePort"].values
                pvisitPort = pd.DataFrame({"chosePort":np.full(100,np.nan)})
                try:
                    pvisitPort.loc[blkTransInd-dat.index.min()-len(preTransSum):\
                               blkTransInd-dat.index.min()-1,"chosePort"] = preTransSum
                except:
                    print("didn't work for session "+str(s)+" block "+str(b))
                    continue
                pvisitPort.loc[blkTransInd-dat.index.min():blkTransInd-\
                               dat.index.min()+len(postTransSum)-1,"chosePort"] = postTransSum
                low2highProbPVisit.append(pvisitPort.values.T[0])
    return np.mean(low2highProbPVisit,axis=0),len(low2highProbPVisit)


pathPortPairs = [["01","10"],["02","20"],["21","12"]]
def get_DistIncreasePVisit(photrats,rat,high2lowDist):
    high2lowDistPVisit = []
    for s in high2lowDist:
        if s not in photrats.triframe.loc[photrats.triframe.rat==rat,"session"].unique():
            continue
        tdat = photrats.triframe.loc[photrats.triframe.session==s,:]
        tdat.loc[:,"path"]= tdat.port.shift(1).fillna(-1).astype("uint8").astype(str).values\
        + tdat.port.astype("uint8").astype(str).values
        for b in high2lowDist[s]:
            if 1 not in high2lowDist[s][b]:
                continue
            bpaths = pathPortPairs[np.where(np.array(high2lowDist[s][b])==1)[0][0]]
            for i in range(len(bpaths)):
                pathports = [int(i) for i in bpaths[i]]
                bind = tdat.loc[tdat.block.diff()!=0].index[b]
                dat = tdat.loc[bind-49:bind+50].copy()
                if len(dat)<100:
                    continue
                dat.loc[:,'choosePath'] = 0
                dat.loc[np.isin(dat.path,bpaths),'choosePath']=1
                dat.loc[:,"pathAvail"] = np.isin(dat.port.shift(1),pathports).astype(int)
                availdat = dat.loc[dat.pathAvail==1,:]
                blkTransInd = availdat.loc[availdat.block.diff()==1,:].index[0]
                preTransSum = availdat.loc[:blkTransInd,"choosePath"].values
                postTransSum = availdat.loc[blkTransInd:,"choosePath"].values
                pTakePath = pd.DataFrame({"takeSum":np.full(100,np.nan)})
                pTakePath.loc[blkTransInd-dat.index.min()-len(preTransSum)+1:\
                               blkTransInd-dat.index.min(),"takeSum"] = preTransSum
                pTakePath.loc[blkTransInd-dat.index.min():blkTransInd-\
                               dat.index.min()+len(postTransSum)-1,"takeSum"] = postTransSum
                high2lowDistPVisit.append(pTakePath.values.T[0])
    return np.mean(high2lowDistPVisit,axis=0),len(high2lowDistPVisit)

def get_DistDecreasePVisit(photrats,rat,low2highDist):
    low2highDistPVisit = []
    for s in low2highDist:
        if s not in photrats.triframe.loc[photrats.triframe.rat==rat,"session"].unique():
            continue
        tdat = photrats.triframe.loc[photrats.triframe.session==s,:]
        tdat.loc[:,"path"]= tdat.port.shift(1).fillna(-1).astype("uint8").astype(str).values\
        + tdat.port.astype("uint8").astype(str).values
        for b in low2highDist[s]:
            if 1 not in low2highDist[s][b]:
                continue
            bpaths = pathPortPairs[np.where(np.array(low2highDist[s][b])==1)[0][0]]
            for i in range(len(bpaths)):
                pathports = [int(i) for i in bpaths[i]]
                bind = tdat.loc[tdat.block.diff()!=0].index[b]
                dat = tdat.loc[bind-49:bind+50].copy()
                if len(dat)<100:
                    continue
                dat.loc[:,'choosePath'] = 0
                dat.loc[np.isin(dat.path,bpaths),'choosePath']=1
                dat.loc[:,"pathAvail"] = np.isin(dat.port.shift(1),pathports).astype(int)
                availdat = dat.loc[dat.pathAvail==1,:]
                blkTransInd = availdat.loc[availdat.block.diff()==1,:].index[0]
                preTransSum = availdat.loc[:blkTransInd,"choosePath"].values
                postTransSum = availdat.loc[blkTransInd:,"choosePath"].values
                pTakePath = pd.DataFrame({"takeSum":np.full(100,np.nan)})
                pTakePath.loc[blkTransInd-dat.index.min()-len(preTransSum)+1:\
                               blkTransInd-dat.index.min(),"takeSum"] = preTransSum
                pTakePath.loc[blkTransInd-dat.index.min():blkTransInd-\
                               dat.index.min()+len(postTransSum)-1,"takeSum"] = postTransSum
                low2highDistPVisit.append(pTakePath.values.T[0])
    return np.mean(low2highDistPVisit,axis=0),len(low2highDistPVisit)

def calc_pChoosePortAtProbBlkTrans(photrats,high2lowProb,low2highProb):
    lows = []
    highs = []
    lowcount = 0
    highcount = 0
    for rat in photrats.triframe.rat.unique().astype(list):
        if rat =="IM-1322":
            continue
        lowVisInfo = get_ProbDecreasePVisit(photrats,rat,high2lowProb)
        highVisInfo = get_ProbIncreasePVisit(photrats,rat,low2highProb)
        lows.append(lowVisInfo[0])
        lowcount += lowVisInfo[1]
        highs.append(highVisInfo[0])
        highcount += highVisInfo[1]
    return lows,highs,lowcount,highcount

def calc_pChoosePathAtBarBlkTrans(photrats,high2lowDist,low2highDist):
    lows = []
    highs = []
    lowcount = 0
    highcount = 0
    for rat in photrats.triframe.rat.unique().astype(list):
        lowVisInfo = get_DistDecreasePVisit(photrats,rat,high2lowDist)
        highVisInfo = get_DistIncreasePVisit(photrats,rat,low2highDist)
        lows.append(lowVisInfo[0])
        lowcount += lowVisInfo[1]
        highs.append(highVisInfo[0])
        highcount += highVisInfo[1]
    return lows,highs,lowcount,highcount

def plot_portChoiceAtProbBlkTrans(photrats,lows,highs,
    smoothWin=5,blineStart=4,legend_on=True):
    xvals = np.arange(-50,50)
    fig = plt.figure(figsize=(5,3.5))
    toplt = np.mean(highs,axis=0)
    toplt = toplt-np.mean(toplt[int(len(toplt)/2)-blineStart:int(len(toplt)/2)+1])
    toplt = pd.Series(toplt).rolling(window=smoothWin).mean()
    plt.plot(xvals,toplt,label = 'p(rwd) Increase',lw=3,color='k')
    plt.fill_between(xvals,toplt+pd.Series(sem(highs)).rolling(smoothWin).mean(),
                     toplt-pd.Series(sem(highs)).rolling(smoothWin).mean(),color='k',alpha=.2)
    toplt = np.mean(lows,axis=0)
    toplt = toplt-np.mean(toplt[int(len(toplt)/2)-blineStart:int(len(toplt)/2)+1])
    toplt = pd.Series(toplt).rolling(window=smoothWin).mean()
    plt.plot(xvals,toplt,label = 'p(rwd) Decrease',lw=3,color='k',ls='--')
    plt.fill_between(xvals,toplt+pd.Series(sem(lows)).rolling(smoothWin).mean(),
                     toplt-pd.Series(sem(lows)).rolling(smoothWin).mean(),color='k',alpha=.2)
    plt.axvline(x=0,color='k')
    plt.ylabel("\u0394 p(visit port)",fontsize='xx-large')
    plt.xlabel("trials from p(rwd) change",fontsize='xx-large')
    ax = plt.gca()
    ax.tick_params(axis="x", direction="inout")
    ax.tick_params(axis="y", direction="inout")
    plt.xlim(0,20)
    if legend_on:
        plt.legend()
    plt.tight_layout()
    return fig


def plot_portChoiceAtBarBlkTrans(photrats,lows,highs,
    smoothWin=5,blineStart=10,legend_on=True):
    xvals = np.arange(-48,52)
    fig = plt.figure(figsize=(4.8,3.1))
    toplt = np.mean(lows,axis=0)
    toplt = toplt-np.mean(toplt[int(len(toplt)/2)-blineStart:int(len(toplt)/2)])
    toplt = pd.Series(toplt).rolling(window=smoothWin).mean()
    plt.plot(xvals,toplt,label = 'Path length decrease',lw=3,color='k')
    plt.fill_between(xvals,toplt+pd.Series(sem(highs)).rolling(smoothWin).mean(),
                     toplt-pd.Series(sem(highs)).rolling(smoothWin).mean(),color='k',alpha=.2)
    toplt = np.mean(highs,axis=0)
    toplt = toplt-np.mean(toplt[int(len(toplt)/2)-blineStart:int(len(toplt)/2)])
    toplt = pd.Series(toplt).rolling(window=smoothWin).mean()
    plt.plot(xvals,toplt,label = 'Path length increase',lw=3,color='k',ls='--')
    plt.fill_between(xvals,toplt+pd.Series(sem(lows)).rolling(smoothWin).mean(),
                     toplt-pd.Series(sem(lows)).rolling(smoothWin).mean(),color='k',alpha=.2)
    plt.axvline(x=0,color='k')
    plt.ylabel("\u0394 p(take path)",fontsize='xx-large')
    plt.xlabel("trials from barrier change",fontsize='xx-large')
    ax = plt.gca()
    ax.tick_params(axis="x", direction="inout")
    ax.tick_params(axis="y", direction="inout")
    plt.xlim(0,20)
    if legend_on:
        plt.legend()
    plt.tight_layout()
    return fig

def plot_logChoosPortVsPrwdDif(photrats):
    plt.figure(figsize=(3.2,3.5))
    for r in photrats.regdf.rat.unique():
        sns.regplot(x="pRwdDif",y="choose_L",data=photrats.regdf.loc[
            (photrats.regdf.rat==r)&(photrats.regdf.tri>25)],logistic=True,
                   scatter_kws={"color":'lightblue','alpha':0.00,"marker":''},line_kws={"color":'k',"lw":2,"alpha":0.4},ci=None)
    sns.regplot(x="pRwdDif",y="choose_L",data=photrats.regdf.loc[
            (photrats.regdf.tri>25)],logistic=True,
                   scatter_kws={"color":'lightblue','alpha':0.00,"marker":''},line_kws={"color":'k',"lw":3,"alpha":1})
    plt.ylabel("P(choose port)",fontsize="xx-large")
    plt.xlabel("relative P(rwd)",fontsize="xx-large")
    plt.ylim(.1,.83)#0.18,0.72)
    plt.xlim(-100,100)
    plt.xticks([-80,-40,0,40,80])
    plt.tight_layout()

def plot_logChoosPortVsDistDif(photrats):
    plt.figure(figsize=(3.2,3.5))
    for r in photrats.regdf.rat.unique():
        sns.regplot(x="ldif",y="choose_L",data=photrats.regdf.loc[photrats.regdf.rat==r],logistic=True,
                   scatter_kws={"color":'lightblue','alpha':0.00,"marker":''},line_kws={"color":'k',"lw":2,"alpha":0.4},ci=None)
    sns.regplot(x="ldif",y="choose_L",data=photrats.regdf,logistic=True,
                   scatter_kws={"color":'lightblue','alpha':0.00,"marker":''},line_kws={"color":'k',"lw":3,"alpha":1.0})
    plt.ylabel("P(choose port)",fontsize="xx-large")
    plt.xlabel("relative path length\n(hexes)",fontsize="xx-large")
    plt.ylim(0,1)
    plt.xticks([-10,-5,0,5,10])
    plt.tight_layout()


def add_scaledVarsToRegDf(photrats):
    photrats.regdf.loc[:,"rwdDifScaled"] = valscale.fit_transform(\
        photrats.regdf.loc[:,"pRwdDif"].values.reshape(-1,1))
    photrats.regdf.loc[:,"lengthDifScaled"] = valscale.fit_transform(\
        photrats.regdf.loc[:,"ldif"].values.reshape(-1,1))

def plot_pChooseVrelativeVar_byRat(photrats,use_scaled=False,use_prob=True):
    if use_prob:
        varString = "rwdDifScaled" if use_scaled else "pRwdDif"
        xstring = "relative value\n(scaled p(reward))" if use_scaled\
        else "relative value\n(reward probability)"
    else:
        varString = "lengthDifScaled" if use_scaled else "ldif"
        xstring = "relative path length\n(transformed hexes)" if use_scaled\
        else "relative path length\n(hexes)"
    fig = plt.figure(figsize=(3.2,3.5))
    for r in photrats.regdf.rat.unique():
        sns.regplot(x=varString,y="choose_L",data=photrats.regdf.loc[photrats.regdf.rat==r],logistic=True,
                   scatter_kws={"color":'lightblue','alpha':0.00,"marker":''},line_kws={"color":'k',"lw":2,"alpha":0.4},ci=None)
    sns.regplot(x=varString,y="choose_L",data=photrats.regdf,logistic=True,
                   scatter_kws={"color":'lightblue','alpha':0.00,"marker":''},line_kws={"color":'k',"lw":3,"alpha":1.0})
    plt.ylabel("P(choose port)",fontsize="xx-large")
    plt.xlabel(xstring,fontsize="xx-large")
    plt.ylim(0,1)
    #plt.xticks([-10,-5,0,5,10])
    plt.tight_layout()
    return fig