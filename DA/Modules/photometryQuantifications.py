"""Functions to analyze and vizualize photometry data.
Specifically analyses for DA maze paper main figure plots.
TODO: break into multiple PhotRats child classes"""

__author__ = "Tim Krausz"
__email__ = "krausz.tim@gmail.com"
__status__ = "development"
__license__ = "MIT"

from __main__ import *


def add_pairedHexStates2df(photrats):
    photrats.df.loc[:,"pairedHexStates"] = -1
    photrats.df.loc[photrats.df.hexlabels.diff()!=0,\
    "pairedHexStates"] = photrats.convert2pairedState(
           photrats.df.loc[photrats.df.hexlabels.diff()!=0,\
           "hexlabels"].values)
    photrats.df.loc[:,"pairedHexStates"] = photrats.df.loc[:,\
    "pairedHexStates"].replace(-1,method="ffill")

def plot_portAlignedDaInTime(photrats,secondHalfOnly=True,
    poolFactor="nom_rwd_chosen",useRatGroupLevel=True):
    photrats.set_pool_factor(poolFactor)
    photrats.set_plot_trace("green_z_scored")
    fig = plt.figure(figsize = (7,5))
    rwdhigh,omhigh,rwdmid,ommid,rwdlow,omlow = photrats.getSessionTercMeans(\
        secondHalf=secondHalfOnly,useRat=useRatGroupLevel)
    
    high_color = "red"#"indianred"
    mid_color = "firebrick"
    low_color = "maroon"
    high_colorOm = "dodgerblue"
    mid_colorOm = "blue"
    low_colorOm ="darkblue"
    
    xvals = np.arange(photrats.fs*photrats.plot_window[0],\
        photrats.fs*photrats.plot_window[1]+1)/photrats.fs
    ax1 = plt.gca()
    plot_avgWithSem(photrats,ax1,xvals,np.vstack((rwdhigh,omhigh)),'lightgrey','-',\
                    [None,-photrats.plot_window[0]*photrats.fs],"high p(Reward)")
    plot_avgWithSem(photrats,ax1,xvals,rwdhigh,high_color,'-',\
        [-photrats.plot_window[0]*photrats.fs,None])
    plot_avgWithSem(photrats,ax1,xvals,omhigh,high_colorOm,':',\
        [-photrats.plot_window[0]*photrats.fs,None])
    plot_avgWithSem(photrats,ax1,xvals,np.vstack((rwdmid,ommid)),'darkgrey','-',\
                    [None,-photrats.plot_window[0]*photrats.fs],"medium p(Reward)")
    plot_avgWithSem(photrats,ax1,xvals,rwdmid,mid_color,'-',\
        [-photrats.plot_window[0]*photrats.fs,None])
    plot_avgWithSem(photrats,ax1,xvals,ommid,mid_colorOm,':',\
        [-photrats.plot_window[0]*photrats.fs,None])
    ax1.axvline(x=0.0,ymin=-.1,ymax=1.0,color='k',linestyle='--')
    ax1.set_xlabel('time (s) from port entry')
    plot_avgWithSem(photrats,ax1,xvals,np.vstack((rwdlow,omlow)),'dimgrey','-',\
                    [None,-photrats.plot_window[0]*photrats.fs],"low p(Reward)")
    plot_avgWithSem(photrats,ax1,xvals,rwdlow,low_color,'-',\
        [-photrats.plot_window[0]*photrats.fs,None])
    plot_avgWithSem(photrats,ax1,xvals,omlow,low_colorOm,':',\
        [-photrats.plot_window[0]*photrats.fs,None])
    ax1.axvline(x=0.0,ymin=-.1,ymax=1.0,color='k',linestyle='--')
    ax1.set_xlabel('time (s) from port entry')
    ax1.legend()
    plt.xlabel("time from port entry (s)",fontsize='xx-large')
    plt.ylabel("DA (z-scored)",fontsize='xx-large')
    plt.tight_layout()
    return fig


def plot_portAlignedDaInTime_byQ(photrats):
    photrats.set_pool_factor("Q_chosen")
    photrats.set_plot_trace("green_z_scored")
    fig = plt.figure(figsize = (7,5))
    rwdhigh,omhigh,rwdmid,ommid,rwdlow,omlow = photrats.getSessionTercMeans(secondHalf=True)
    
    high_color = "red"#"indianred"
    mid_color = "firebrick"
    low_color = "maroon"
    high_colorOm = "dodgerblue"
    mid_colorOm = "blue"
    low_colorOm ="darkblue"
    
    xvals = np.arange(photrats.fs*photrats.plot_window[0],photrats.fs*photrats.plot_window[1]+1)/photrats.fs
    ax1 = plt.gca()
    plot_avgWithSem(photrats,ax1,xvals,np.vstack((rwdhigh,omhigh)),'lightgrey','-',\
                    [None,-photrats.plot_window[0]*photrats.fs],"high Q")
    plot_avgWithSem(photrats,ax1,xvals,rwdhigh,high_color,'-',[-photrats.plot_window[0]*photrats.fs,None])
    plot_avgWithSem(photrats,ax1,xvals,omhigh,high_colorOm,':',[-photrats.plot_window[0]*photrats.fs,None])
    plot_avgWithSem(photrats,ax1,xvals,np.vstack((rwdmid,ommid)),'darkgrey','-',\
                    [None,-photrats.plot_window[0]*photrats.fs],"medium Q")
    plot_avgWithSem(photrats,ax1,xvals,rwdmid,mid_color,'-',[-photrats.plot_window[0]*photrats.fs,None])
    plot_avgWithSem(photrats,ax1,xvals,ommid,mid_colorOm,':',[-photrats.plot_window[0]*photrats.fs,None])
    ax1.axvline(x=0.0,ymin=-.1,ymax=1.0,color='k',linestyle='--')
    ax1.set_xlabel('time (s) from port entry')
    plot_avgWithSem(photrats,ax1,xvals,np.vstack((rwdlow,omlow)),'dimgrey','-',\
                    [None,-photrats.plot_window[0]*photrats.fs],"low Q")
    plot_avgWithSem(photrats,ax1,xvals,rwdlow,low_color,'-',[-photrats.plot_window[0]*photrats.fs,None])
    plot_avgWithSem(photrats,ax1,xvals,omlow,low_colorOm,':',[-photrats.plot_window[0]*photrats.fs,None])
    ax1.axvline(x=0.0,ymin=-.1,ymax=1.0,color='k',linestyle='--')
    ax1.set_xlabel('time (s) from port entry')
    ax1.legend()
    plt.xlabel("time from port entry (s)",fontsize='xx-large')
    plt.ylabel("DA (z-scored)",fontsize='xx-large')
    plt.tight_layout()

def plot_avgWithSem(photratz,ax,xvals,plot_trace,colorString,linstyle='-',subset=[None,None],traceLabel=None):
    ax.plot(xvals[subset[0]:subset[1]],np.nanmean\
        (plot_trace,axis=0)[subset[0]:subset[1]],color=colorString,ls=linstyle,label=traceLabel)
    ax.fill_between(xvals[subset[0]:subset[1]],(np.nanmean\
        (plot_trace,axis=0)+sem(plot_trace,nan_policy="omit"))\
        [subset[0]:subset[1]],(np.nanmean\
        (plot_trace,axis=0)-sem(plot_trace,nan_policy="omit"))\
        [subset[0]:subset[1]],color=colorString,alpha=.5)


def plot_meanRatDAafterPortEntry(photrats,highInds,midInds,lowInds,pltCol="firebrick"):
    highMeans = calc_DaPeakTroughAfterInds(photrats,highInds)
    midMeans = calc_DaPeakTroughAfterInds(photrats,midInds)
    lowMeans = calc_DaPeakTroughAfterInds(photrats,lowInds)
    highMeans = [np.mean(rm) for rm in highMeans if len(rm)>0]
    midMeans = [np.mean(rm) for rm in midMeans if len(rm)>0]
    lowMeans = [np.mean(rm) for rm in lowMeans if len(rm)>0]
    plt.bar([0,1,2],[np.mean(highMeans),\
            np.mean(midMeans),np.mean(lowMeans)],color=pltCol,alpha=0.3)
    plt.ylabel("mean $\Delta$DA",fontsize='xx-large',fontweight='bold')
    plt.xlabel("p(reward)",fontsize='xx-large',fontweight='bold')
    for r in range(len(highMeans)):
        plt.scatter(x=np.add([0,1,2],np.random.randn(1)/10),\
            y=[np.mean(highMeans[r]),np.mean(midMeans[r]),np.mean(lowMeans[r])],\
                    c=pltCol,edgecolors='k',lw=2,s=45)
        plt.plot([0,1,2],[np.mean(highMeans[r]),\
            np.mean(midMeans[r]),np.mean(lowMeans[r])],color='k',alpha=0.5,lw=1)
    plt.xticks([0,1,2],["High","Med","Low"],fontsize="x-large")
    plt.tight_layout()

def plot_peakTroughDAafterPortEntry_barWithRats(photrats,highInds,midInds,lowInds,peak=True,pltCol="firebrick"):
    highMeans = calc_DaPeakTroughAfterInds(photrats,highInds,peak)
    midMeans = calc_DaPeakTroughAfterInds(photrats,midInds,peak)
    lowMeans = calc_DaPeakTroughAfterInds(photrats,lowInds,peak)
    highMeans = [np.mean(rm) for rm in highMeans if len(rm)>0]
    midMeans = [np.mean(rm) for rm in midMeans if len(rm)>0]
    lowMeans = [np.mean(rm) for rm in lowMeans if len(rm)>0]
    plt.bar([0,1,2],[np.mean(highMeans),\
            np.mean(midMeans),np.mean(lowMeans)],color=pltCol,alpha=0.3)
    plt.ylabel("mean $\Delta$DA",fontsize='xx-large',fontweight='bold')
    plt.xlabel("p(reward)",fontsize='xx-large',fontweight='bold')
    for r in range(len(highMeans)):
        plt.scatter(x=np.add([0,1,2],np.random.randn(1)/10),\
            y=[np.mean(highMeans[r]),np.mean(midMeans[r]),np.mean(lowMeans[r])],\
                    c=pltCol,edgecolors='k',lw=2,s=45)
        plt.plot([0,1,2],[np.mean(highMeans[r]),\
            np.mean(midMeans[r]),np.mean(lowMeans[r])],color='k',alpha=0.5,lw=1)
    plt.xticks([0,1,2],["High","Med","Low"],fontsize="x-large")
    plt.tight_layout()

def calc_DaChangeAtInds(photrats,indices):
    photrats.set_plot_window([0,1])
    tracesPost = get_TracesAroundIndex(photrats,indices)
    photrats.set_plot_window([-1,0])
    tracesPre = get_TracesAroundIndex(photrats,indices)
    traceChangeRats = photrats.df.loc[indices,"rat"].astype(str).values
    traceChangeRatMeans = []
    for rat in photrats.df.rat.unique():
        traceChangeRatMeans.append(np.mean(tracesPost[np.where(traceChangeRats==rat)[0]],axis=1)-\
                            np.mean(tracesPre[np.where(traceChangeRats==rat)[0]],axis=1))
    return traceChangeRatMeans

def calc_DaPeakTroughAfterInds(photrats,indices,peak=True):
    photrats.set_plot_window([0,1])
    tracesPost = get_TracesAroundIndex(photrats,indices)
    tracePeakRats = photrats.df.loc[indices,"rat"].astype(str).values
    tracePeakRatMeans = []
    for rat in photrats.df.rat.unique():
        if peak:
            tracePeakRatMeans.append(np.max(tracesPost[np.where(tracePeakRats==rat)[0]],axis=1))
        else:
            tracePeakRatMeans.append(np.min(tracesPost[np.where(tracePeakRats==rat)[0]],axis=1))
    return tracePeakRatMeans


def plot_peakTroughDaDifAfterPortEntry_barWithRats(photrats,highInds,
    midInds,lowInds,peak=True,pltCol="firebrick"):
    highMeans = calc_DaPeakTroughDiffAfterPortInds(photrats,highInds,peak)
    midMeans = calc_DaPeakTroughDiffAfterPortInds(photrats,midInds,peak)
    lowMeans = calc_DaPeakTroughDiffAfterPortInds(photrats,lowInds,peak)
    plt.bar([0,1,2],[np.mean(highMeans),\
            np.mean(midMeans),np.mean(lowMeans)],color=pltCol,alpha=0.3)
    plt.ylabel("mean $\Delta$DA",fontsize='xx-large',fontweight='bold')
    plt.xlabel("p(reward)",fontsize='xx-large',fontweight='bold')
    for r in range(len(highMeans)):
        plt.scatter(x=np.add([0,1,2],np.random.randn(1)/10),\
            y=[np.mean(highMeans[r]),np.mean(midMeans[r]),np.mean(lowMeans[r])],\
                    c=pltCol,edgecolors='k',lw=2,s=55)
        plt.plot([0,1,2],[np.mean(highMeans[r]),\
            np.mean(midMeans[r]),np.mean(lowMeans[r])],color='k',alpha=0.5,lw=1)
    plt.xticks([0,1,2],["High","Med","Low"],fontsize="x-large")
    plt.tight_layout()

def calc_DaPeakTroughDiffAfterPortInds(photrats,indices,peak=True):
    winMax = 0.5 if peak else 1.0
    photrats.set_plot_window([0,winMax])
    tracePeakRats = photrats.df.loc[indices,"rat"].astype(str).values
    tracePeakRatMeans = []
    for rat in photrats.df.rat.unique():
        tracesPost = get_TracesAroundIndex(photrats,
            indices[np.isin(indices,photrats.df.loc[photrats.df.rat==rat].index)])
        tracePost = np.mean(tracesPost,axis=0)
        if peak:
            tracePeakRatMeans.append(np.max(tracePost)-tracePost[0])
        else:
            tracePeakRatMeans.append(np.min(tracePost)-tracePost[0])
    return tracePeakRatMeans

def get_TracesAroundIndex(photrats,indices):
    traces = []
    for i in indices:
        if np.isnan(i) or i == -1:
            continue
        traces.append(photrats.df.loc[i+photrats.fs*photrats.plot_window[0]:\
                i+photrats.fs*photrats.plot_window[1],photrats.plot_trace].values)
    return np.array(traces)

def calc_DaChangeVprobCors(photrats):
    rwdCors = []
    omCors = []
    for rat in photrats.df.rat.unique():
        photrats.dat = photrats.df.loc[(photrats.df.rat==rat)\
                                       &(photrats.df.tri>25)&(photrats.df.rwd==1),]
        photrats.dat_visinds = photrats.dat.loc[photrats.dat.port!=-100].index
        lowInds,midInds,highInds = photrats.getTriIndsByTerc()
        daChanges = np.concatenate([calc_DaChangeAtIndsOneRat(photrats,highInds,peak=True),
                     calc_DaChangeAtIndsOneRat(photrats,midInds,peak=True),
                     calc_DaChangeAtIndsOneRat(photrats,lowInds,peak=True)])
        probs = np.concatenate([photrats.dat.loc[highInds,"nom_rwd_chosen"].values/100,\
        photrats.dat.loc[midInds,"nom_rwd_chosen"].values/100,\
        photrats.dat.loc[lowInds,"nom_rwd_chosen"].values/100])
        rwdCors.append(pearsonr(probs,daChanges))
        photrats.dat = photrats.df.loc[(photrats.df.rat==rat)\
                                       &(photrats.df.tri>25)&(photrats.df.rwd==0),]
        photrats.dat_visinds = photrats.dat.loc[photrats.dat.port!=-100].index
        lowInds,midInds,highInds = photrats.getTriIndsByTerc()
        daChanges = np.concatenate([calc_DaChangeAtIndsOneRat(photrats,highInds,peak=False),
                     calc_DaChangeAtIndsOneRat(photrats,midInds,peak=False),
                     calc_DaChangeAtIndsOneRat(photrats,lowInds,peak=False)])
        probs = np.concatenate([photrats.dat.loc[highInds,"nom_rwd_chosen"].values/100,\
        photrats.dat.loc[midInds,"nom_rwd_chosen"].values/100,\
        photrats.dat.loc[lowInds,"nom_rwd_chosen"].values/100])
        omCors.append(pearsonr(probs,daChanges))
        pd.DataFrame(rwdCors,columns=["coef","p-val"]).to_csv(photrats.directory_prefix+"pearsonR_result_DaVsRpe_rwd.csv")
        pd.DataFrame(omCors,columns=["coef","p-val"]).to_csv(photrats.directory_prefix+"pearsonR_result_DaVsRpe_om.csv")
    return rwdCors,omCors

def calc_DaChangeAtIndsOneRat(photrats,indices,peak=True):
    winMax = 0.5 if peak else 1.0
    photrats.set_plot_window([0,winMax])
    tracesPost = get_TracesAroundIndex(photrats,indices)
    daChange = np.max(tracesPost,axis=1) - tracesPost[:,0]\
     if peak else np.min(tracesPost,axis=1) - tracesPost[:,0]
    return daChange

def calc_DaPeakTroughAfterIndsOneRat(photrats,indices,rat,peak=True):
    photrats.set_plot_window([0,1])
    tracesPost = get_TracesAroundIndex(photrats,indices)
    tracePeakRats = photrats.df.loc[indices,"rat"].astype(str).values
    tracePeakRatMeans = []
    if peak:
        tracePeakRatMeans = np.max(tracesPost[np.where(tracePeakRats==rat)[0]],axis=1)
    else:
        tracePeakRatMeans = np.min(tracesPost[np.where(tracePeakRats==rat)[0]],axis=1)
    return tracePeakRatMeans

def calcRpeRegByRatAndSesh(photrats,useQ=False):
    ratRwdRpes = {r:[] for r in photrats.df.loc[:,"rat"].unique()}
    ratOmRpes = {r:[] for r in photrats.df.loc[:,"rat"].unique()}
    ratRwdNs = {r:[] for r in photrats.df.loc[:,"rat"].unique()}
    ratOmNs = {r:[] for r in photrats.df.loc[:,"rat"].unique()}
    for r in tqdm(range(len(photrats.df.loc[:,"rat"].unique()))):
        rat = photrats.df.loc[:,"rat"].unique()[r]
        for s in photrats.df.loc[photrats.df.rat==rat,"session"].unique():
            if useQ:
                dat = photrats.df.loc[(photrats.df.session==s),:]
                qs = dat.loc[dat.port!=-100,"Q_chosen"].values
            else:
                dat = photrats.df.loc[(photrats.df.session==s)&(photrats.df.tri>25),:]
                pRwds = dat.loc[dat.port!=-100,"nom_rwd_chosen"].values
            rwds = dat.loc[dat.port!=-100,"rwd"].values
            rpes = rwds-qs if useQ else rwds-(pRwds/100)
            rweightsRwd,rweightsOm = calcRpeLagRegWeightsBinned(photrats,dat,rpes,rwds,100)
            ratRwdRpes[rat].append(rweightsRwd[0])
            ratOmRpes[rat].append(rweightsOm[0])
            ratRwdNs[rat].append(len(rwds[rwds==1]))
            ratOmNs[rat].append(len(rwds[rwds==0]))
    return ratRwdRpes,ratOmRpes,ratRwdNs,ratOmNs

def calcRpeLagRegWeights(photrats,dat,pRwdRpes,rwds):
    dat_vinds = dat.loc[(dat.port!=-100),:].index
    lags = np.arange(0,photrats.fs*2.5)
    rweightsRwd = np.zeros((1,len(lags)))
    rweightsOm = np.zeros((1,len(lags)))
    for n in range(len(lags)):
        preLagInds = dat_vinds[np.where(rwds==1)]+lags[n]
        lagInds = preLagInds[np.isin(preLagInds,dat.index)]
        yRwd = dat.loc[lagInds,"green_z_scored"]
        yRwd.reset_index(drop=True,inplace=True)
        XRwd = pd.DataFrame({"posRpe":pRwdRpes[np.where(rwds==1)][np.isin(preLagInds,dat.index)]})
        XRwd = XRwd.drop(yRwd.loc[yRwd.isnull()].index,axis=0)
        yRwd = yRwd.drop(yRwd.loc[yRwd.isnull()].index,axis=0)
        modRwd = LR(fit_intercept=True,normalize=False).fit(XRwd,yRwd)
        rweightsRwd[:,n] = modRwd.coef_
        preLagInds = dat_vinds[np.where(rwds==0)]+lags[n]
        lagInds = preLagInds[np.isin(preLagInds,dat.index)]
        yOm = dat.loc[lagInds,"green_z_scored"]
        yOm.reset_index(drop=True,inplace=True)
        XOm = pd.DataFrame({"posRpe":pRwdRpes[np.where(rwds==0)][np.isin(preLagInds,dat.index)]})
        XOm = XOm.drop(yOm.loc[yOm.isnull()].index,axis=0)
        yOm = yOm.drop(yOm.loc[yOm.isnull()].index,axis=0)
        modOm = LR(fit_intercept=True,normalize=False).fit(XOm,yOm)
        rweightsOm[:,n] = modOm.coef_
    return rweightsRwd,rweightsOm

def calcRpeLagRegWeightsBinned(photrats,dat,pRwdRpes,rwds,binsize=50):
    dat_vinds = dat.loc[(dat.port!=-100),:].index
    lags = np.arange(0,photrats.fs*2.5,int(photrats.fs*(binsize/1000)))
    rweightsRwd = np.zeros((1,len(lags)))
    rweightsOm = np.zeros((1,len(lags)))
    for n in range(0,len(lags)-1):
        preLagInds = dat_vinds[np.where(rwds==1)]+lags[n]
        lagInds = preLagInds[np.isin(preLagInds,dat.index)]
        das = []
        for l in lagInds:
            das.append(dat.loc[l:l+lags[n+1],"green_z_scored"].mean())
        yRwd = pd.Series(das)#dat.loc[lagInds,"green_z_scored"]
        XRwd = pd.DataFrame({"posRpe":pRwdRpes[np.where(rwds==1)][np.isin(preLagInds,dat.index)]})
        XRwd = XRwd.drop(yRwd.loc[yRwd.isnull()].index,axis=0)
        yRwd = yRwd.drop(yRwd.loc[yRwd.isnull()].index,axis=0)
        modRwd = LR(fit_intercept=True,normalize=False).fit(XRwd,yRwd)
        rweightsRwd[:,n] = modRwd.coef_
        preLagInds = dat_vinds[np.where(rwds==0)]+lags[n]
        lagInds = preLagInds[np.isin(preLagInds,dat.index)]
        das = []
        for l in lagInds:
            das.append(dat.loc[l:l+lags[n+1],"green_z_scored"].mean())
        yOm = pd.Series(das)
        XOm = pd.DataFrame({"posRpe":pRwdRpes[np.where(rwds==0)][np.isin(preLagInds,dat.index)]})
        XOm = XOm.drop(yOm.loc[yOm.isnull()].index,axis=0)
        yOm = yOm.drop(yOm.loc[yOm.isnull()].index,axis=0)
        modOm = LR(fit_intercept=True,normalize=False).fit(XOm,yOm)
        rweightsOm[:,n] = modOm.coef_
    return rweightsRwd,rweightsOm

def plot_rpeLagRegCoefs(photrats,binsize=100):
    xvals = np.arange(0,photrats.fs*2,250/(1000/binsize))/photrats.fs#np.arange(0,photrats.fs*2)/photrats.fs
    fig = plt.figure(figsize=(4.5,5))
    ax1 = plt.subplot2grid((6,4),(2,0),colspan = 4, rowspan =4)
    plot_ratMeans(xvals,ratRwdRpes,'darkred',pltLabel="reward")
    plot_ratMeans(xvals,ratOmRpes,'darkblue',pltLabel="omission")
    plt.xlabel("time from port entry (s)",fontsize='xx-large')
    plt.axhline(0,ls=':',color='k')
    plt.ylabel("RPE ß",fontsize='xx-large')
    plt.xlim(0,2)
    plt.legend()
    ax2 = plt.subplot2grid((6,4),(0,0),colspan = 4, rowspan =2,sharex=ax1)
    plot_sigPoints(xvals,ratRwdRpes,'darkred',plot99=False)
    plot_sigPoints(xvals,ratOmRpes,'darkblue',plot99=False)
    ax2.tick_params('x', labelbottom=False)
    plt.ylabel("Fraction\nsignificant")
    ax2.set_ylim(0,1)
    plt.tight_layout()
    fig.savefig(photrats.directory_prefix+"rpeLagReg_binned.pdf")

def plot_ratMeans(xvals,ratDict,pltColor='darkred',pltLabel=""):
    ratmeans = []
    for rat in ratDict:
        ratmeans.append(np.nanmean(ratDict[rat],axis=0))
        #plt.plot(xvals,ratmeans[-1],color=pltColor,alpha=0.2,lw=.5)
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
        
def save_RpeNumbers(photrats,ratRwdNs,ratOmNs):
    numbers = 'rewarded:\n'
    for rat in ratRwdNs:
        numbers += rat +" had " + str(len(ratRwdNs[rat])) +\
        " sessions with " + str(ratRwdNs[rat]) + " = " + str(np.sum(ratRwdNs[rat])) + " RPEs.\n"
    numbers += "omissions:\n"
    for rat in ratOmNs:
        numbers += rat +" had " + str(len(ratOmNs[rat])) +\
        " sessions  with " + str(ratOmNs[rat]) + " = " + str(np.sum(ratOmNs[rat])) + " RPEs.\n"
    
    with open(photrats.directory_prefix+"rpeLagRegNumbers.txt", 'w') as f:
        f.write(numbers)

def removeErrantBlock1Assignments(photrats):
    '''Identify and remove data appended to end of certain sessions
    marked as "block 1" incorrectly'''
    start_extra_block1_inds = photrats.df.loc[\
    (photrats.df.block.diff()<0)&(photrats.df.session.diff()==0),:].index
    seshStartInds = photrats.df.loc[(photrats.df.session.diff()!=0)].index
    end_extraBlock1_inds = []
    for i in start_extra_block1_inds[:-1]:
        end_extraBlock1_inds.append(seshStartInds[np.where(seshStartInds>i)[0][0]])
    end_extraBlock1_inds.append(photrats.df.index.max())
    for i in range(len(start_extra_block1_inds)):
        photrats.df.drop(index=np.arange(\
            start_extra_block1_inds[i],end_extraBlock1_inds[i]),axis=0,inplace=True)
    photrats.df.drop(index=photrats.df.index.max(),axis=0,inplace=True)
    photrats.df.reset_index(inplace=True)
    photrats.get_visinds()

def get_newPathTracesByDistToPort(adjHexIndsSortedByDist):
    smoothWin = int(photrats.fs/4)
    shortTrace = []
    midTrace = []
    longTrace = []
    terc_cutoff = int(len(adjHexIndsSortedByDist)/3)
    for i in range(len(adjHexIndsSortedByDist)):
        adjInd = adjHexIndsSortedByDist[i]
        trace = photrats.df.loc[adjInd+photrats.plot_window[0]*photrats.fs:\
                            adjInd+photrats.plot_window[1]*photrats.fs,photrats.plot_trace].rolling(smoothWin).mean().values
        if i<=terc_cutoff:
            shortTrace.append(trace)
        elif i<=terc_cutoff*2:
            midTrace.append(trace)
        else:
            longTrace.append(trace)
    return shortTrace,midTrace,longTrace

def get_newPathTracesByDistToPort_absoluteDist(newHexAdjInds,distsToPort,dist_cutoff=7):
    smoothWin = int(photrats.fs/4)
    shortTrace = []
    midTrace = []
    longTrace = []
    for i in range(len(newHexAdjInds)):
        adjInd = newHexAdjInds[i]
        dist = distsToPort[i]
        trace = photrats.df.loc[adjInd+photrats.plot_window[0]*photrats.fs:\
                adjInd+photrats.plot_window[1]*photrats.fs,photrats.plot_trace].rolling(smoothWin).mean().values
        if dist<=dist_cutoff:
            shortTrace.append(trace)
        elif dist<=dist_cutoff*2:
            midTrace.append(trace)
        else:
            longTrace.append(trace)
    return shortTrace,midTrace,longTrace

def find_newHexAdjInds(photrats):
    newHexAdjInds = []
    newHexAdjStates = []
    enteredHex = []
    enteredHexSoon = []
    allEntryInds = photrats.df.loc[photrats.df.hexlabels.diff()!=0,:].index
    for s in photrats.sesh_newlyAvailHexes:
        for b in range(len(photrats.sesh_newlyAvailHexes[s])):
            #get indices of first entry into an adjacent-to-newly-available hex
            newHexAdjInds.append(photrats.df.loc[(photrats.df.session==s)\
                &(photrats.df.block==b+2)&(photrats.df.adj2newlyAvail.diff()==1)].index.min())
            if np.isnan(newHexAdjInds[-1]):
                print("session ",str(s)," block ",b+1)
                newHexAdjStates.append(-1)
                enteredHex.append(-1)
                enteredHexSoon.append(-1)
            else:
                nextHexInd = allEntryInds[np.where(allEntryInds==newHexAdjInds[-1])[0][0]+1]
                nextHex = photrats.df.loc[nextHexInd,"hexlabels"]#newly avail are in hexID format
                nextHexes = photrats.df.loc[allEntryInds[np.where(allEntryInds==newHexAdjInds[-1])[0][0]+1:\
                         np.where(allEntryInds==newHexAdjInds[-1])[0][0]+4],"hexlabels"].values
                enteredHex.append((nextHex==photrats.sesh_newlyAvailHexes[s][b]).astype(int)[0])
                enteredHexSoon.append(int(photrats.sesh_newlyAvailHexes[s][b][0] in nextHexes))
                newHexAdjStates.append(photrats.df.loc[newHexAdjInds[-1],"pairedHexStates"])
    return newHexAdjInds,newHexAdjStates,enteredHex,enteredHexSoon

def find_enteredVignoredNewlyAvailInds(newHexAdjInds,enteredHex):
    enteredInds = [n for n,z in zip(newHexAdjInds,enteredHex) if z==1]
    ignoredInds = [n for n,z in zip(newHexAdjInds,enteredHex) if z==0]
    return enteredInds,ignoredInds

def find_blockedHexAdjInds(photrats):
    blockedHexAdjInds = []
    blockedHexAdjStates = []
    for s in photrats.sesh_newlyBlockedHexes:
        for b in range(len(photrats.sesh_newlyBlockedHexes[s])):
            blockedHexAdjInds.append(photrats.df.loc[(photrats.df.session==s)\
                &(photrats.df.block==b+2)&(photrats.df.adj2newlyBlocked.diff()==1)].index.min())
            if np.isnan(blockedHexAdjInds[-1]):
                print("session ",str(s)," block ",b+1)
                blockedHexAdjStates.append(-1)
            else:
                blockedHexAdjStates.append(photrats.df.loc[blockedHexAdjInds[-1],"pairedHexStates"])
    return blockedHexAdjInds,blockedHexAdjStates

# can plot aligned to newly avail entry, or
def find_newHexEntryAndPriorHexInds(photrats):
    newHexInds = []
    newHexStates = []
    adjHexStates = [] #previously entered hex
    adjHexInds = []
    allEntryInds = photrats.df.loc[photrats.df.pairedHexStates.diff()!=0,:].index
    for s in photrats.sesh_newlyAvailHexes:
        for b in range(len(photrats.sesh_newlyAvailHexes[s])):
            newHexInds.append(photrats.df.loc[(photrats.df.session==s)\
                &(photrats.df.block==b+2)&(photrats.df.newlyAvailHex.diff()==1)].index.min())
            if np.isnan(newHexInds[-1]):
                print("session ",str(s)," block ",b+1)
                newHexStates.append(-1)
                adjHexInds.append(-1)
                adjHexStates.append(-1)
            else:
                newHexStates.append(photrats.df.loc[newHexInds[-1],"pairedHexStates"])
                try:
                    adjHexInds.append(allEntryInds[np.where(allEntryInds==newHexInds[-1])[0][0]-1])
                    adjHexStates.append(photrats.df.loc[adjHexInds[-1],"pairedHexStates"])
                except:
                    adjHexInds.append(-1)
                    adjHexStates.append(-1)
                    print("No previous adjacent hex entry detected for session ",str(s)," block ",b+1)
    return newHexInds,newHexStates,adjHexStates,adjHexInds

def plotFirstEntryHexChangeMeanOverRats(photrats,availratmeans,blockedratmeans,legend_on=False,pltCol1='deeppink',pltCol2='k',ls2='-'):
    photrats.set_plot_trace("green_z_scored")
    photrats.set_plot_window([-5,5])
    fig = plt.figure(figsize=(7,5))#(4.8,5))
    xvals = np.arange(photrats.plot_window[0]*photrats.fs,photrats.plot_window[1]*photrats.fs+1)/photrats.fs
    smoothWin = int(photrats.fs/10)
    toplt = pd.Series(np.mean(availratmeans,axis=0)).rolling(smoothWin).mean().values
    topltSem = pd.Series(sem(availratmeans,axis=0)).rolling(smoothWin).mean().values
    plt.plot(xvals,toplt,label="Newly available",color=pltCol1,lw=3)
    plt.fill_between(xvals,toplt-topltSem,toplt+topltSem,color=pltCol1,alpha=.3)
    toplt = pd.Series(np.mean(blockedratmeans,axis=0)).rolling(smoothWin).mean().values
    topltSem = pd.Series(sem(blockedratmeans,axis=0)).rolling(smoothWin).mean().values
    plt.plot(xvals,toplt,label="Newly blocked",color=pltCol2,ls=ls2,lw=3)
    plt.fill_between(xvals,toplt-topltSem,toplt+topltSem,color=pltCol2,alpha=.3)
    plt.xlabel("time from hex entry (s)",fontsize="xx-large")
    plt.ylabel("Mean z-scored DA",fontsize="xx-large")
    plt.axvline(x=0,ls='--',color='k',alpha=.8,lw=1)
    plt.xticks(np.arange(-5,6))
    plt.tight_layout()
    plt.ylim(-.9,2.8)
    if legend_on:
        plt.legend()
    plt.tight_layout()
    return fig

def calc_DaPeakIndividualDiffsAfterNewPathInds(photrats,indices):
    photrats.set_plot_window([-1,0.25])
    tracePeakRats = photrats.df.loc[indices,"rat"].astype(str).values
    tracePeakChanges = []
    missingInds = []
    for rat in photrats.df.rat.unique():
        if rat not in tracePeakRats:
            missingInds.append(np.where(photrats.df.rat.unique()==rat)[0][0])
            continue
        tracesPost = get_TracesAroundIndex(photrats,indices[tracePeakRats==rat])
        bline = tracesPost[:,0]#np.mean(tracesPost,axis=0)[0]
        #tracePost = np.mean(tracesPost,axis=0)[photrats.fs*1:]
        daChanges = np.max(tracesPost[:,photrats.fs*1:],axis=1)-bline
        tracePeakChanges += list(daChanges)
    return tracePeakChanges

def calc_DaPeakDiffAfterNewPathInds(photrats,indices):
    photrats.set_plot_window([-1,0.25])
    tracePeakRats = photrats.df.loc[indices,"rat"].astype(str).values
    tracePeakRatMeans = []
    missingInds = []
    for rat in photrats.df.rat.unique():
        if rat not in tracePeakRats:
            missingInds.append(np.where(photrats.df.rat.unique()==rat)[0][0])
            continue
        tracesPost = get_TracesAroundIndex(photrats,indices[tracePeakRats==rat])
        bline = tracesPost[:,0]#np.mean(tracesPost,axis=0)[0]
        #tracePost = np.mean(tracesPost,axis=0)[photrats.fs*1:]
        daChanges = np.max(tracesPost[:,photrats.fs*1:],axis=1)-bline
        tracePeakRatMeans.append(np.mean(daChanges))#max(tracePost)-bline)
    return tracePeakRatMeans,missingInds

def get_ratTracesAtHexChangeDiscovery(photrats,xvals,indices,plotTraces=True):
    tracePeakRats = photrats.df.loc[indices,"rat"].astype(str).values
    ratmeans = []
    n_PerRat = {r:[] for r in photrats.df.rat.unique()}
    for rat in photrats.df.rat.unique():
        tracesPost = get_TracesAroundIndex(photrats,indices[tracePeakRats==rat])
        tracePost = np.mean(tracesPost,axis=0)
        n_PerRat[rat] = len(tracesPost)
        if len(tracesPost)>=3:
            ratmeans.append(tracePost)
            if plotTraces:
                plt.plot(xvals,tracePost,color='k',alpha=0.3,lw=1)
    return ratmeans,n_PerRat

def plot_ratTracesAtHexChangeDiscovery(photrats,inds,pltCol='k'):
    '''Plot individual rat averages at discovery of newly available and newly blocked paths.
    Return average of rat average traces.'''
    xvals = np.arange(photrats.plot_window[0]*photrats.fs,\
        photrats.plot_window[1]*photrats.fs+1)/photrats.fs
    fig = plt.figure()
    tracePeakRats = photrats.df.loc[inds,"rat"].astype(str).values
    ratmeans,n_PerRat = get_ratTracesAtHexChangeDiscovery(photrats,xvals,inds)
    plt.plot(xvals,np.mean(ratmeans,axis=0),color=pltCol)
    #plt.ylim(-1.5,4.9)
    plt.ylabel("mean DA")
    plt.xlabel("time from port entry (s)")
    plt.tight_layout()
    return fig,ratmeans,n_PerRat

def plot_FirstEntryHexChange(photrats,adjHexInds,blockedHexAdjInds,legend_on=False):
    photrats.set_plot_window([-5,5])
    smoothWin = int(photrats.fs/4)
    fig = plt.figure(figsize=(7,5))#(4.8,5))
    xvals = np.arange(photrats.plot_window[0]*photrats.fs,photrats.plot_window[1]*photrats.fs+1)/photrats.fs
    
    avail_traces,blocked_traces = get_availAndBlockedTraces(photrats,adjHexInds,blockedHexAdjInds)
    
    toplt = pd.Series(np.mean(avail_traces,axis=0)).rolling(smoothWin).mean().values
    topltSem = pd.Series(sem(avail_traces,axis=0)).rolling(smoothWin).mean().values
    plt.plot(xvals,toplt,label="Newly available",color='deeppink',lw=3)
    plt.fill_between(xvals,toplt-topltSem,toplt+topltSem,color='deeppink',alpha=.3)
    toplt = pd.Series(np.mean(blocked_traces,axis=0)).rolling(smoothWin).mean().values
    topltSem = pd.Series(sem(blocked_traces,axis=0)).rolling(smoothWin).mean().values
    plt.plot(xvals,toplt,label="Newly blocked",color='k',ls=':',lw=3)
    plt.fill_between(xvals,toplt-topltSem,toplt+topltSem,color='k',ls=':',alpha=.3)
    plt.xlabel("time from hex entry (s)",fontsize="xx-large")
    plt.ylabel("Mean z-scored DA",fontsize="xx-large")
    plt.axvline(x=0,ls='--',color='k',alpha=.8,lw=1)
    plt.xticks(np.arange(-5,6))
    if legend_on:
        plt.legend()
    plt.tight_layout()
    return figayout()
    return fig

def plot_FirstAdjEntryByEnteredVsIgnored(photrats,enteredInds,ignoredInds,legend_on=False,pltCol1='deeppink',pltCol2='k',ls2='-'):
    photrats.set_plot_trace("green_z_scored")
    photrats.set_plot_window([-5,5])
    smoothWin = int(photrats.fs/4)
    fig = plt.figure(figsize=(7,5))
    xvals = np.arange(photrats.plot_window[0]*photrats.fs,photrats.plot_window[1]*photrats.fs+1)/photrats.fs
    
    avail_traces,blocked_traces = get_availAndBlockedTraces(photrats,enteredInds,ignoredInds)
    
    toplt = pd.Series(np.mean(avail_traces,axis=0)).rolling(smoothWin).mean().values
    topltSem = pd.Series(sem(avail_traces,axis=0)).rolling(smoothWin).mean().values
    plt.plot(xvals,toplt,label="entered",color=pltCol1,lw=3)
    plt.fill_between(xvals,toplt-topltSem,toplt+topltSem,color=pltCol1,alpha=.2)
    toplt = pd.Series(np.mean(blocked_traces,axis=0)).rolling(smoothWin).mean().values
    topltSem = pd.Series(sem(blocked_traces,axis=0)).rolling(smoothWin).mean().values
    plt.plot(xvals,toplt,label="ignored",ls=ls2,color=pltCol2,lw=3)
    plt.fill_between(xvals,toplt-topltSem,toplt+topltSem,color=pltCol2,alpha=.2)
    plt.xlabel("Time from changed-hex discovery (s)",fontsize="xx-large")
    plt.ylabel("Mean z-scored DA",fontsize="xx-large")
    plt.axvline(x=0,ls='--',color='k',alpha=.9,lw=1)
    plt.xticks(np.arange(-5,6))
    if legend_on:
        plt.legend()
    plt.tight_layout()
    return fig

def get_availAndBlockedTraces(photrats,adjHexInds,blockedHexAdjInds):
    avail_traces = []
    blocked_traces = []
    for i in adjHexInds:
        if np.isnan(i) or i == -1:
            continue
        avail_traces.append(photrats.df.loc[i+photrats.fs*photrats.plot_window[0]:\
                i+photrats.fs*photrats.plot_window[1],photrats.plot_trace].values)
    for i in blockedHexAdjInds:
        if np.isnan(i):
            continue
        blocked_traces.append(photrats.df.loc[i+photrats.fs*photrats.plot_window[0]:\
                i+photrats.fs*photrats.plot_window[1],photrats.plot_trace].values)
    return np.array(avail_traces),np.array(blocked_traces)

def plot_meanRatDAafterHexEntry(photrats,adjHexInds,blockedHexAdjInds,pltCol1="#27aeef",pltCol2= "#b33dc6"):
    availRatMeans,blockedRatMeans = calc_DaChangeAtHexEntry(photrats,adjHexInds,blockedHexAdjInds)
    blockedMeans = [np.mean(rm) for rm in blockedRatMeans if len(rm)>0]
    availMeans = [np.mean(rm) for rm in availRatMeans if len(rm)>0]
    fig = plt.figure(figsize=(4,5.5))
    plt.bar([0,1],[np.mean(availMeans),\
            np.mean(blockedMeans)],color='k',alpha=0.3)
    plt.ylabel("mean $\Delta$DA",fontsize='xx-large',fontweight='bold')
    plt.xlabel("hex type",fontsize='xx-large',fontweight='bold')
    for r in range(len(availMeans)):
        plt.scatter(x=0,y=np.mean(availMeans[r]),color=pltCol1,marker='o')
        try:
            plt.scatter(x=1,y=np.mean(blockedMeans[r]),color=pltCol2,marker='o')
            plt.plot([0,1],[np.mean(availMeans[r]),np.mean(blockedMeans[r])],color='k',alpha=0.5,lw=1)
        except:
            continue
    print("blocked hex: ")
    sigBlocked = get_sigRats_fromMeanList(blockedMeans)
    print("avail hex: ")
    sigAvail = get_sigRats_fromMeanList(availMeans)
    missingInd = [i for i in range(len(blockedRatMeans)) if len(blockedRatMeans[i])==0]
    if len(missingInd)>0:
        [availMeans.pop(i) for i in missingInd]
    sigPaired = get_sigRatsPaired_from2samples(availMeans,blockedMeans,"greater")
    plot_sigMarkers(sigPaired,0.5,2.4)
    plot_sigMarkers(sigAvail,-.05,2.1)
    plot_sigMarkers(sigBlocked,1,2.1)
    plt.tight_layout()
    return fig

def plot_meanRatDaChangeAfterHexEntry(photrats,adjHexInds,blockedHexAdjInds,pltCol1="#27aeef",pltCol2= "#b33dc6"):
    #availRatMeans,blockedRatMeans = calc_DaChangeAtHexEntry(photrats,adjHexInds,blockedHexAdjInds)
    availMeans,_ = calc_DaPeakDiffAfterNewPathInds(photrats,adjHexInds)
    blockedMeans,missingInd = calc_DaPeakDiffAfterNewPathInds(photrats,blockedHexAdjInds)
    fig = plt.figure(figsize=(4,5.5))
    plt.bar([0,1],[np.mean(availMeans),\
            np.mean(blockedMeans)],color='k',alpha=0.3)
    plt.ylabel("mean $\Delta$DA",fontsize='xx-large',fontweight='bold')
    plt.xlabel("hex type",fontsize='xx-large',fontweight='bold')
    for r in range(len(availMeans)):
        plt.scatter(x=0,y=np.mean(availMeans[r]),color=pltCol1,marker='o')
        try:
            plt.scatter(x=1,y=np.mean(blockedMeans[r]),color=pltCol2,marker='o')
            plt.plot([0,1],[np.mean(availMeans[r]),np.mean(blockedMeans[r])],color='k',alpha=0.5,lw=1)
        except:
            continue
    #for r in [np.mean(rm) for rm in blockedRatMeans]:
    #    plt.scatter(x=1,y=np.mean(r),color=pltCol2,marker='o')
    print("blocked hex: ")
    sigBlocked = get_sigRats_fromMeanList(blockedMeans)
    print("avail hex: ")
    sigAvail = get_sigRats_fromMeanList(availMeans)
    #missingInd = [i for i in range(len(blockedRatMeans)) if len(blockedRatMeans[i])==0]
    if len(missingInd)>0:
        availMeans = np.delete(availMeans,missingInd)
    print("paired test")
    sigPaired = get_sigRatsPaired_from2samples(availMeans,blockedMeans,"greater")
    plot_sigMarkers(sigPaired,0.5,2.4)
    plot_sigMarkers(sigAvail,-.05,2.1)
    plot_sigMarkers(sigBlocked,1,2.1)
    plt.tight_layout()
    return fig

def calc_DaChangeAtHexEntry(photrats,adjHexInds,blockedHexAdjInds):
    photrats.set_plot_window([0,1])
    avail_tracesPost,blocked_tracesPost = get_availAndBlockedTraces(photrats,adjHexInds,blockedHexAdjInds)
    photrats.set_plot_window([-1,0])
    avail_tracesPre,blocked_tracesPre = get_availAndBlockedTraces(photrats,adjHexInds,blockedHexAdjInds)
    availRats = photrats.df.loc[adjHexInds,"rat"].astype(str).values
    blockedRats = photrats.df.loc[blockedHexAdjInds,"rat"].astype(str).values
    availRatMeans = []
    blockedRatMeans = []
    for rat in photrats.df.rat.unique():
        availRatMeans.append(np.mean(avail_tracesPost[np.where(availRats==rat)[0]],axis=1)-\
                            np.mean(avail_tracesPre[np.where(availRats==rat)[0]],axis=1))
        blockedRatMeans.append(np.mean(blocked_tracesPost[np.where(blockedRats==rat)[0]],axis=1)-\
                              np.mean(blocked_tracesPre[np.where(blockedRats==rat)[0]],axis=1))
    return availRatMeans,blockedRatMeans

def create_triframe(photrats):
    photrats.get_visinds()
    photrats.triframe = photrats.df.loc[photrats.visinds]
    photrats.triframe['lrchoice'] = make_lr(photrats.triframe.port.values)
    photrats.triframe['lrchoice'] = photrats.triframe.lrchoice.astype("int8")
    photrats.triframe.drop(['x','y','vel',"green_z_scored","beamA",'beamB',
                            'beamC'],axis=1,inplace=True)
    photrats.triframe.reset_index(inplace=True)
    
def plotRegSigLevel(regSum,x,y,regInd):
    if "Pr(>|z|)" in regSum.columns:
        if regSum.loc[regInd,"Pr(>|z|)"]<0.001:
            plt.text(x=x-.1,y=y,s="***",fontweight='bold',fontsize='xx-large')
        elif regSum.loc[regInd,"Pr(>|z|)"]<0.01:
            plt.text(x=x-.05,y=y,s="**",fontweight='bold',fontsize='xx-large')
        elif regSum.loc[regInd,"Pr(>|z|)"]<0.05:
            plt.text(x=x,y=y,s="*",fontweight='bold',fontsize='xx-large')
        else:
            plt.text(x=x,y=y,s="ns",fontweight='bold',fontsize='xx-large')
    else:
        if regSum.loc[regInd,"t value"]>3.291:
            plt.text(x=x-.1,y=y,s="***",fontweight='bold',fontsize='xx-large')
        elif regSum.loc[regInd,"t value"]>2.58:
            plt.text(x=x-.05,y=y,s="**",fontweight='bold',fontsize='xx-large')
        elif regSum.loc[regInd,"t value"]>1.96:
            plt.text(x=x,y=y,s="*",fontweight='bold',fontsize=25)
        else:
            plt.text(x=x,y=y,s="ns",fontweight='bold',fontsize='xx-large')

def plot_sigMarkers(sigrats,r,yval):
    if sigrats[0]==1:
        plt.text(x=r, y=yval, s='***',fontweight='bold',fontsize='xx-large')
    elif sigrats[1]==1:
        plt.text(x=r, y=yval, s='**',fontweight='bold',fontsize='xx-large')
    elif sigrats[2]==1:
        plt.text(x=r, y=yval, s='*',fontweight='bold',fontsize='xx-large')
    else:
        plt.text(x=r-.25, y=yval, s='ns',fontweight='bold',fontsize='x-large')
        
