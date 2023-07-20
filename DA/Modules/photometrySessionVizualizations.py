"""Functions to vizualize single-session photometry data.
Specifically analyses for DA maze paper main figure plots."""

__author__ = "Tim Krausz"
__email__ = "krausz.tim@gmail.com"
__status__ = "development"

import matplotlib.colors as mc
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter
from __main__ import *

hexlist = [2,47,46,45,44,43,3,\
49,42,41,40,39,48,\
38,37,36,35,34,33,\
32,31,30,29,28,\
27,26,25,24,23,\
22,21,20,19,\
18,17,16,15,\
14,13,12,\
11,10,9,\
8,7,\
6,5,\
4,\
1]
#hexlist = np.subtract(hexlist,1) #convert to index-based states
#scaleFactor = [(460)/13,350/18]
scaleFactor = [(550-70)/13,380/18]
coords = []
cols = [7,6,6,5,5,4,4,3,3,2,2,1,1]
maxrows = 13
r = 0
x = 1
y = 18
startr = 1
while r < maxrows:
    maxcols = cols[r]
    c = 0
    if r%2!=0:
        startr+=1
    x=startr
    while c < maxcols:
        #coords.append([x*scaleFactor[0]+80,y*scaleFactor[1]+55])
        coords.append([x*scaleFactor[0]+70,y*scaleFactor[1]+40])
        x += 2
        c += 1
    if r%2!=0:
        y -= 2
    else:
        y-=1
    r += 1
cents = {h: c for h,c in zip(hexlist,coords)}
centdf = pd.DataFrame(cents)
centdf = centdf.T

def plot_hex_outline(barriers):
    bardf = centdf.drop(barriers,axis=0)
    plt.scatter(bardf.loc[:,0].values,bardf.loc[:,1].values,c=\
            'darkgrey',marker='H',s=1000,edgecolors="k",alpha=1,lw=2)

def plot_hex_barriers(barriers):
    plt.scatter(centdf.loc[barriers,0].values,centdf.loc[barriers,1].values,c=\
            'darkred',marker='o',edgecolor='k',s=750,alpha=.4)
    plt.scatter(centdf.loc[barriers,0].values,centdf.loc[barriers,1].values,c=\
            '',marker='o',edgecolor='darkred',lw=3,s=750,alpha=1)
    
def plot_sesh_pathChoices_ticks(dat,blocks=[],startrows=0):
    if len(blocks)!=0:
        dat = dat.loc[dat.block.isin(blocks)]
    visinds = dat.loc[dat.port>=0].index.values
    x1 = np.linspace(0,len(dat),len(dat))/photrats.fs/60
    #bdpath = plt.figure(figsize = (11.3,3.9))#(22,8))
    blks = blocks
    ax1 = plt.subplot2grid((6+startrows,len(blks)),(2+startrows,0),\
        colspan = len(blks), rowspan =1)
    ax2 = plt.subplot2grid((6+startrows,len(blks)),(1+startrows,0),\
    colspan = len(blks), rowspan =1,sharey=ax1)#,sharex=ax1)
    ax3 = plt.subplot2grid((6+startrows,len(blks)),(0+startrows,0),\
    colspan = len(blks), rowspan =1,sharey=ax1)#,sharex=ax1)
    linw = 5
    blineLength = 12
    denomLength = 3
    
    A1 = (dat.loc[(dat.port==0),"dtop"]-blineLength)/denomLength
    yA1 = np.zeros(len(dat))
    yA1[A1.index.values-dat.index.min()] = A1
    B1 = (dat.loc[(dat.port==1),"dtop"]-blineLength)/denomLength
    yB1 = np.zeros(len(dat))
    yB1[B1.index.values-dat.index.min()] = B1
    C1 = (dat.loc[(dat.port==2),"dtop"]-blineLength)/denomLength
    yC1 = np.zeros(len(dat))
    yC1[C1.index.values-dat.index.min()] = C1
    ax1.plot(x1,yA1,color="darkblue",lw=5)
    ax2.plot(x1,yB1,color="darkorange",lw=5)
    ax3.plot(x1,yC1,color="darkgreen",lw=5)
    
    A1 = (dat.loc[(dat.port==0)&(dat.rwd==0),"dtop"]-blineLength)/denomLength
    yA1 = np.zeros(len(dat))
    yA1[A1.index.values-dat.index.min()] = A1
    B1 = (dat.loc[(dat.port==1)&(dat.rwd==0),"dtop"]-blineLength)/denomLength
    yB1 = np.zeros(len(dat))
    yB1[B1.index.values-dat.index.min()] = B1
    C1 = (dat.loc[(dat.port==2)&(dat.rwd==0),"dtop"]-blineLength)/denomLength
    yC1 = np.zeros(len(dat))
    yC1[C1.index.values-dat.index.min()] = C1
    ax1.plot(x1,yA1,color="white",lw=3)
    ax3.set_yticklabels([''])
    ax2.plot(x1,yB1,color="white",lw=3)
    ax2.set_xticklabels('')
    ax2.set_yticklabels([''])
    ax3.plot(x1,yC1,color="white",lw=3)
    ax3.set_xticklabels('')
    ax3.set_yticklabels([''])
    
    
    ax1.set_xlabel('Time (min)',fontsize=25,fontweight='bold')
    
    for b in range(dat.block.min(),int(dat.block.max())+1):
        xstart = int(min(dat.loc[dat.block==b].index-dat.index.min())/photrats.fs/60)
        ind = max(dat.loc[dat.block==b].index)
        xmin = int(min(dat.loc[dat.block==b].index-dat.index.min())/photrats.fs/60)
        xmid = int(xmin+(xstart-xmin)/2)+5
        if b==1:
            ax1.axvline(x=xstart,color ='r',linestyle='--', label = 'Block Change')
        ax1.axvline(x=xstart,color ='red',linestyle='--',lw=5)
        ax2.axvline(x=xstart,color ='red',linestyle='--',lw=5)
        ax3.axvline(x=xstart,color ='red',linestyle='--',lw=5)
    ax1.set_xlim(0,x1[-1])
    ax2.set_xlim(0,x1[-1])
    ax3.set_xlim(0,x1[-1])
    plt.ylim(0.2,4.2)
    plt.tight_layout()

def plot_posOverlayAndTickPlot(photrats,s,blks = [1,2],
        posColor='cyan',saveFig=True,secondHalf=False,plot_ticks=False,
        edgCol='k', plotOverlay=True,plotProbs=True,trans=0.2,density=False,vmax=60,vmin=5):
    dat = photrats.df.loc[photrats.df.session==s].copy()
    halfString = '_only2ndHalfBlkPos' if secondHalf else ''
    densString = "_density" if density else ""
    fig = plt.figure(figsize=(15.5,9))
    startrows = 10
    if plot_ticks:
        plot_sesh_pathChoices_ticks(dat,blocks=blks,startrows=startrows)
    plt.subplots_adjust(wspace=0, hspace=0,top=0.95)
    for blk in blks:
        prwds = dat.loc[dat.block==blk,["nom_rwd_a","nom_rwd_b",\
        "nom_rwd_c"]].values[0].astype(int)
        if dat.session_type.values[0] == "barrier":
            bars = np.add(photrats.sesh_barIDs[s][blk-1],1)
        else:
            bars = np.add(photrats.sesh_barIDs[s],1)
        ax = plt.subplot2grid((6+startrows,len(blks)*10),(0,(blk-np.min(blks))*10),\
            colspan = 8, rowspan =startrows)
        minTri = 25 if secondHalf else 0 
        if plotOverlay:
            x = dat.loc[(dat.block==blk)&(dat.tri>minTri)&(dat.x.notnull()),'x'].values
            y = dat.loc[(dat.block==blk)&(dat.tri>minTri)&(dat.x.notnull()),'y'].values
            if density:
                H, xedges, yedges = np.histogram2d(x,y,bins=40)
                X, Y = np.meshgrid(xedges, yedges)
                plt.pcolormesh(X,Y,H.T,cmap="Greys",norm=mc.Normalize(vmin=vmin,vmax=vmax))
                plot_hex_outline(bars)
                plot_hex_barriers(bars)
                plt.ylim(0,522)
                plt.xlim(55,590)
                ax = plt.gca()
                ax.invert_yaxis()
                plt.tight_layout()
                plt.xticks([])
                plt.yticks([])
            else:
                plot_hex_outline(bars)
                plot_hex_barriers(bars)
                plt.ylim(0,522)
                plt.xlim(55,590)
                ax = plt.gca()
                ax.invert_yaxis()
                plt.tight_layout()
                plt.xticks([])
                plt.yticks([])
                xy = np.vstack([x,y])
                z = gaussian_kde(xy)(xy)
                plt.scatter(x,y,c=z,s=20,cmap="viridis")

        if plotProbs:
            plt.text(x=300,y=0,s=str(prwds[0])+"%",fontsize=30,\
                fontweight='bold',backgroundcolor="k",color="white")
            plt.text(x=30,y=470,s=str(prwds[1])+"%",fontsize=30,\
                fontweight='bold',backgroundcolor="k",color="white")
            plt.text(x=570,y=470,s=str(prwds[2])+"%",fontsize=30,\
                fontweight='bold',backgroundcolor="k",color="white")
        plt.axis("off")
    if saveFig:
        fig.savefig(photrats.directory_prefix+"sesh"+\
            str(s)+"blocks"+str(blks)+"positionOverlayPlot"+halfString+densString+".pdf")
        fig.savefig(photrats.directory_prefix+"sesh"+\
            str(s)+"blocks"+str(blks)+"positionOverlayPlot"+halfString+densString+".png")
def plot_10tri_posOverlay(photrats,s,blks = [1,2],
        posColor='cyan',saveFig=True,groupOfTen=1,
        edgCol='k', plotOverlay=True,plotProbs=True):
    dat = photrats.df.loc[photrats.df.session==s].copy()
    triString = str(groupOfTen)
    fig = plt.figure(figsize=(15.5,9))
    startrows = 10
    plt.subplots_adjust(wspace=0, hspace=0,top=0.95)
    #plt.figure(figsize=(17.4,6))
    for blk in blks:#dat.block.unique():
        prwds = dat.loc[dat.block==blk,["nom_rwd_a","nom_rwd_b",\
        "nom_rwd_c"]].values[0].astype(int)
        if dat.session_type.values[0] == "barrier":
            bars = np.add(photrats.sesh_barIDs[s][blk-1],1)
        else:
            bars = np.add(photrats.sesh_barIDs[s],1)
        
        ax = plt.subplot2grid((6+startrows,len(blks)*10),(0,(blk-np.min(blks))*10),\
            colspan = 8, rowspan =startrows)
        plot_hex_outline(bars)
        plot_hex_barriers(bars)
        plt.ylim(0,522)
        plt.xlim(55,590)
        ax = plt.gca()
        ax.invert_yaxis()
        plt.tight_layout()
        plt.xticks([])
        plt.yticks([])
        minTri = (groupOfTen-1)*10
        maxTri = groupOfTen*10
        if plotOverlay:
            x = dat.loc[(dat.block==blk)&(dat.tri>minTri)&(dat.x.notnull()),'x'].values
            y = dat.loc[(dat.block==blk)&(dat.tri>minTri)&(dat.x.notnull()),'y'].values
            xy = np.vstack([x,y])
            z = gaussian_kde(xy)(xy)
            plt.scatter(x,y,c=z,s=20,cmap="viridis")
        if plotProbs:
            plt.text(x=300,y=0,s=str(prwds[0])+"%",fontsize=30,\
                fontweight='bold',backgroundcolor="k",color="white")
            plt.text(x=30,y=470,s=str(prwds[1])+"%",fontsize=30,\
                fontweight='bold',backgroundcolor="k",color="white")
            plt.text(x=570,y=470,s=str(prwds[2])+"%",fontsize=30,\
                fontweight='bold',backgroundcolor="k",color="white")
        plt.axis("off")
    if saveFig:
        fig.savefig(photrats.directory_prefix+"sesh"+\
            str(s)+"blocks"+str(blks)+"positionOverlayPlot"+str(triString)+"_10tris.pdf")
        fig.savefig(photrats.directory_prefix+"sesh"+\
            str(s)+"blocks"+str(blks)+"positionOverlayPlot"+str(triString)+"_10tris.png")

def plot_individualSeshTrace_RwdAndNewHexLabeled(photrats,s,
    plot_newlyAvail=False,plot_ref=False):
    if plot_ref:
        rat = photrats.df.loc[photrats.df.session==s,"rat"].values[0]
        date = photrats.df.loc[photrats.df.session==s,"date"].values[0]
        date = '0'+str(date) if int(str(date)[0]) not in [0,1] else str(date)
        df4ref = pd.read_csv(photrats.directory_prefix+rat+"/"+date+"/sampleframe.csv")
        reftrace = df4ref.ref.rolling(window=int(photrats.fs/5)).mean().values
        del df4ref
    trace = photrats.df.loc[photrats.df.session==s,"green_z_scored"].rolling(window=int(photrats.fs/5)).mean().values
    veltrace = photrats.df.loc[photrats.df.session==s,"vel"].rolling(window=int(photrats.fs/5)).mean().values
    startInd = photrats.df.loc[photrats.df.session==s,].index.min()
    if plot_newlyAvail:
        seshEnterInds = np.array(enteredInds)[(np.array(enteredInds)>photrats.df.loc[\
                    photrats.df.session==s,"green_z_scored"].index.min())&(np.array(enteredInds)<\
                        photrats.df.loc[photrats.df.session==s,"green_z_scored"].index.max())]-startInd
    rwds = photrats.df.loc[(photrats.df.session==s)&(photrats.df.port!=-100),"rwd"]
    fig = plt.figure(figsize=(18,8))
    ax1 = plt.subplot(311)
    plt.plot(trace)
    for r in range(len(rwds)):
        linsty = '-' if rwds.iloc[r]==1 else ':'
        if r>0:
            plt.axvline(rwds.index[r]-startInd,ls=linsty,color='k')
        else:
            plt.axvline(rwds.index[r]-startInd,color="k",label="rwd")
            plt.axvline(rwds.index[r]-startInd,color="k",ls=":",label="omission")
    if plot_newlyAvail: 
        for i in seshEnterInds:
            plt.axvline(i,color="deeppink")
    ax1=plt.gca()
    scalebar1 = AnchoredSizeBar(ax1.transData,size_vertical=2,label_top=True,size=2,label="2Z",loc="upper right")
    scalebar2 = AnchoredSizeBar(ax1.transData,size=photrats.fs*2,label="2s",loc="center right")
    ax1.add_artist(scalebar1)
    ax1.add_artist(scalebar2)
    if plot_newlyAvail:
        plt.axvline(i,color="deeppink",ls="-",label="entered newly avail hex")
    plt.legend()
    ax2 = plt.subplot(313,sharex=ax1)
    ax2.plot(veltrace,color="darkgrey")
    ax2=plt.gca()
    scalebar3 = AnchoredSizeBar(ax2.transData,\
        size_vertical=20,label_top=True,size=2,label="20 cm/s",loc="upper right")
    scalebar4 = AnchoredSizeBar(ax2.transData,\
        size=photrats.fs*2,label="2s",loc="center right")
    ax2.add_artist(scalebar3)
    ax2.add_artist(scalebar4)
    plt.axis("off")
    if plot_ref:
        ax3 = plt.subplot(312,sharex=ax1)
        ax3.plot(reftrace,color="darkviolet")
        ax3=plt.gca()
        scalebar5 = AnchoredSizeBar(ax3.transData,\
            size_vertical=2,label_top=True,size=2,label="2Z",loc="upper right")
        scalebar6 = AnchoredSizeBar(ax3.transData,\
            size=photrats.fs*2,label="2s",loc="center right")
        ax3.add_artist(scalebar5)
        ax2.add_artist(scalebar6)
    plt.axis("off")
    ax1.axis("off")
    plt.tight_layout()
    return fig,s
