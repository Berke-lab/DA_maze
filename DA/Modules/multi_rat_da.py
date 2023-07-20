"""Object to aggregate and analyze photometry maze
 data across sessions and rats"""

__author__ = "Tim Krausz"
__email__ = "krausz.tim@gmail.com"
__status__ = "development"
__license__ = "MIT"


from single_rat_da import *
from portVal import *
from bar_changes import *
from tmat_ops import *
from tqdm import tqdm

class PhotRats(TmatOperations,BarChanges,PortValAnalyses,Photrat):

    cols2load = ['green_z_scored',"ref",'port','rwd','x','y','nom_rwd_a','nom_rwd_b',\
              'beamA', 'beamB', 'beamC','vel','acc','tri','block',\
             'nom_rwd_c','hexlabels','rat','date','lenAC', 'lenBC',\
             'lenAB', 'dtop','fiberloc','session_type']
    plotvals = False
    use_nom_rwds = False
    bin_size = int(250/10)
    phot_directory_prefix = "/Volumes/Tim K/Photometry Data/Triangle Maze/"
    poolByTerc = True

    def __init__(self,seshdict):
        self.seshdict = seshdict

    def set_savepath(self,rat):
        self.savepath = self.directory_prefix+rat+'/'

    def add_df(self,df):
        try:
            self.df = self.df.append(df,ignore_index=True)
        except:
            self.df = df

    def load_dfs(self):
        '''iterate through IDs and dates to load all dfs and compile into one df'''
        sesh = 1
        for rat in self.seshdict:
            for date in self.seshdict[rat]:
                newdf = self.attempt_load(rat,date)
                newdf['session'] = sesh
                self.add_df(newdf)
                sesh += 1

    def attempt_load(self,rat,date):
        try:
            newdf = reduce_mem_usage(pd.read_csv(self.directory_prefix+\
                rat+'/'+date+'/'+rat+'_'+date+'_h_sampleframe.csv',\
                usecols=self.cols2load))
            return newdf
        except:
            print(f'could not load df for {rat}; {date}')

    def load_and_add_df(self,ID,date):
        newdf = reduce_mem_usage(pd.read_csv(self.directory_prefix+\
            ID+'/'+date+'/'+ID+'_'+date+'_h_sampleframe.csv',usecols=self.cols2load))
        newdf['session']=self.df.session.max()+1
        self.df = self.df.append(newdf,ignore_index=True)

    def get_barIDs(self):
        self.sesh_barIDs = {s:[] for s in self.df.session.unique()}
        for sesh in self.df.session.unique():
            self.sesh_barIDs[sesh] = self.attempt_get_barIDs(sesh)

    def attempt_get_barIDs(self,sesh):
        barIds = []
        try:
            if self.df.loc[self.df.session==sesh,'session_type'].values[0]=='prob':
                tmat = self.sesh_tmats[sesh]
                barIds = np.where(np.mean(tmat,axis=1)==0)[0][1:]-1
            else:
                for b in self.df.loc[self.df.session==sesh,'block'].unique():
                    tmat = self.sesh_tmats[sesh][int(b-1)]
                    barIds.append(np.where(np.mean(tmat,axis=1)==0)[0][1:]-1)
        except:
            print(f'unable to find barriers for session {str(sesh)} block {str(b)}')
        return barIds

    def add_hexDistToPort(self):
        '''Add a column to df that reports how many subsequent
        hexes the rat traversed, following the current time point.
        NOT map distance to port.'''
        self.df['hexDistToPort'] = np.nan
        for s in self.df.session.unique():
            print(s)
            dat = self.df.loc[self.df.session==s]
            vinds = dat.loc[dat.port!=-100].index
            hexchanges = dat.loc[dat.hexlabels.diff()!=0].index
            hexDists = []
            for i in range(len(vinds)):
                if i == 0:
                    triHexChanges = hexchanges[(hexchanges>=0)&(hexchanges<vinds[i])]
                else:
                    triHexChanges = hexchanges[(hexchanges>vinds[i-1])&(hexchanges<vinds[i])]
                hexDists = np.append(hexDists,np.flip(np.arange(1,len(triHexChanges)+1)))
            self.df.loc[hexchanges[hexchanges<vinds[-1]],'hexDistToPort']=hexDists
        self.df.loc[self.visinds,'hexDistToPort']=0
        self.df.hexDistToPort.fillna(method='bfill',inplace=True)

    def get_last_leave(self):
        '''Identify the moment when the rat last left a reward port
        before entering the next port. Values of 0 mark before this event.
        1 marks moment rat left port, 2 marks times after rat left port.'''
        leaveinds = []
        next_inds = []
        for s in self.df.session.unique():
            sdat = self.df.loc[self.df.session==s]
            leave,nextind,sdat = get_leave_inds(sdat)
            leaveinds = leaveinds + list(leave)
            next_inds = next_inds + list(nextind)
        #
        lastapproach = np.zeros(len(self.df))
        lastleave = np.zeros(len(self.df))
        for l,n in zip(leaveinds,next_inds):
            lastapproach[l+1:n+1] = 1
            lastleave[l] = 1
        self.df['lastleave']=lastleave
        self.df.loc[np.where(lastapproach==1)[0],'lastleave']=2
        self.df['lastleave']=self.df.lastleave.astype("int8")

    def check_df_exists(self):
        try:
            self.df
        except:
            raise Exception("Must define dataframe before executing function")
    
    def convertDfHexlabelsToPairedState(self):
        self.df.loc[:,'pairedHexStates'] = -1
        seshs = self.df.session.unique()
        for s in tqdm(range(len(seshs))):
            dat = self.df.loc[self.df.session==s]
            datinds = dat.loc[dat.port!=-100].index
            for i in range(1,len(datinds)):
                tdat = dat.loc[datinds[i-1]:datinds[i]]
                hexInOrder = tdat.loc[tdat.hexlabels.diff()!=0,'hexlabels']
                hexinds = hexInOrder.index
                self.df.loc[hexinds,'pairedHexStates']=convert2pairedState(self,hexInOrder.values)
        self.df.loc[:,"pairedHexStates"] = self.df.loc[:,"pairedHexStates"].astype(\
            object).replace(-1,method='ffill')
        self.df.loc[:,"pairedHexStates"] = self.df.loc[:,"pairedHexStates"].astype(\
            object).replace(-1,method='bfill')

    def convert2pairedState(self,orderedHexes):
        pairedStates = []
        for i in range(1,len(orderedHexes)):
            fr,to = orderedHexes[i-1],orderedHexes[i]
            newstate = phexdf.loc[(phexdf.from_state==fr)&
            (phexdf.to_state==to),'statecodes'].values
            if len(newstate)<1:
                pairedStates.append(-1)
            else:
                pairedStates.append(newstate[0])
        return np.append([-1],pairedStates)

    def plot_hex_outline(self,sesh,block,ax,size='sm'):
        sz = 300 if size=='sm' else 2000
        if self.df.loc[self.df.session==sesh,'session_type'].unique()[0]=='barrier':
            bardf = centdf.drop(self.sesh_barIDs[sesh][block-1]+1,axis=0)
        else:
            bardf = centdf.drop(self.sesh_barIDs[sesh]+1,axis=0)
        ax.scatter(bardf[0].values,bardf[1].values,marker='H',color='grey',\
            edgecolors= "darkgrey",alpha=0.4,s=sz)

    def getTriIndsByTerc(self,rwdtype=None):
        if rwdtype=='rwd':
            vinds = self.dat.loc[(self.dat.port!=-100)&(self.dat.rwd==1)].index
        elif rwdtype=='om':
            vinds = self.dat.loc[(self.dat.port!=-100)&(self.dat.rwd!=1)].index
        else:
            vinds = self.dat_visinds
        if 'nom_rwd_chosen' in self.pool_factor:
            low = self.dat.loc[(self.dat.port!=-100)&((self.dat.nom_rwd_chosen==20)|(self.dat.nom_rwd_chosen==10))].index
            mid = self.dat.loc[(self.dat.port!=-100)&(self.dat.nom_rwd_chosen==50)].index
            high = self.dat.loc[(self.dat.port!=-100)&((self.dat.nom_rwd_chosen==80)|(self.dat.nom_rwd_chosen==90))].index
            return low,mid,high
        if not self.poolByTerc:
            #get indices of absolute threshold groups
            poolmin = self.df.loc[self.visinds,self.pool_factor].min()
            thirds = (self.df.loc[self.visinds,self.pool_factor].max()-poolmin)/3
            thirds = [thirds+poolmin,thirds*2+poolmin,self.dat.loc[:,self.pool_factor].max()]
            low = self.dat.loc[(self.dat.port!=-100)&(self.dat[self.pool_factor]<=thirds[0])].index
            mid = self.dat.loc[(self.dat.port!=-100)&(self.dat[self.pool_factor]<=thirds[1])&
                (self.dat[self.pool_factor]>thirds[0])].index
            high = self.dat.loc[(self.dat.port!=-100)&(self.dat[self.pool_factor]>thirds[1])].index
            return low,mid,high
        approach_dat = []
        for i in range(len(vinds)):
            approach_dat.append([self.dat.loc[vinds[i],self.pool_factor],vinds[i]])
        approach_dat = np.array(approach_dat)
        sorteddat = approach_dat[np.argsort(approach_dat[:,0],axis=0)]
        sorted_tris = sorteddat[:,1]
        low = sorted_tris[:int(len(sorteddat)/3)].astype(int)
        mid = sorted_tris[int(len(sorteddat)/3):int(2*len(sorteddat)/3)].astype(int)
        high = sorted_tris[int(2*len(sorteddat)/3):].astype(int)
        return low,mid,high

    def getTracesByTerc(self):
        low,mid,high = self.getTriIndsByTerc()
        rwdhigh = []
        omhigh = []
        rwdmid = []
        ommid = []
        rwdlow = []
        omlow = []

        for i in low:
            tritrace = self.df.loc[i+self.fs*self.plot_window[0]:i+self.plot_window[1]\
            *self.fs,self.plot_trace].values
            if self.df.loc[i,'rwd'] == 1:
                rwdlow.append(tritrace)
            else:
                omlow.append(tritrace)
        for i in mid:
            tritrace = self.df.loc[i+self.fs*self.plot_window[0]:i+self.plot_window[1]\
            *self.fs,self.plot_trace].values
            if self.df.loc[i,'rwd'] == 1:
                rwdmid.append(tritrace)
            else:
                ommid.append(tritrace)
        for i in high:
            tritrace = self.df.loc[i+self.fs*self.plot_window[0]:i+self.plot_window[1]\
            *self.fs,self.plot_trace].values
            if self.df.loc[i,'rwd'] == 1:
                rwdhigh.append(tritrace)
            else:
                omhigh.append(tritrace)
        return rwdhigh,omhigh,rwdmid,ommid,rwdlow,omlow

    def getSessionTercMeans(self,secondHalf=False,useRat=True):
        rwdhigh_means = []
        omhigh_means = []
        rwdmid_means = []
        ommid_means = []
        rwdlow_means = []
        omlow_means = []
        groupLevel = "rat" if useRat else "session"
        for s in self.df.loc[:,groupLevel].unique():
            if secondHalf:
                self.dat = self.df.loc[(self.df.loc[:,groupLevel]==s)&(self.df.tri>25)]
            else:
                self.dat = self.df.loc[self.df.loc[:,groupLevel]==s]
            self.dat_visinds = self.dat.loc[self.dat.port!=-100].index
            rwdhigh,omhigh,rwdmid,ommid,rwdlow,omlow = self.getTracesByTerc()
            rs = [rwdhigh,omhigh,rwdmid,ommid,rwdlow,omlow]
            for rvec in range(len(rs)):
                if len(rs[rvec])<1:
                    rs[rvec] = np.full((10,(self.plot_window[1]-self.plot_window[0])*self.fs+1),np.nan)
                    print("rat didn't visit one port in session ",str(s),\
                        " or probability was not offered")
            rwdhigh,omhigh,rwdmid,ommid,rwdlow,omlow = rs
            rwdhigh_means.append(pd.Series(np.nanmean(rwdhigh,axis=0)).rolling(window=self.bin_size).mean())
            omhigh_means.append(pd.Series(np.nanmean(omhigh,axis=0)).rolling(window=self.bin_size).mean())
            rwdmid_means.append(pd.Series(np.nanmean(rwdmid,axis=0)).rolling(window=self.bin_size).mean())
            ommid_means.append(pd.Series(np.nanmean(ommid,axis=0)).rolling(window=self.bin_size).mean())
            rwdlow_means.append(pd.Series(np.nanmean(rwdlow,axis=0)).rolling(window=self.bin_size).mean())
            omlow_means.append(pd.Series(np.nanmean(omlow,axis=0)).rolling(window=self.bin_size).mean())
        return rwdhigh_means,omhigh_means,rwdmid_means,ommid_means,rwdlow_means,omlow_means

    def set_regressionFeatures(self,featureList):
        self.regFeatures = featureList

    def plotSeshRweightsInTime(self,rweights):
        fig = plt.figure(figsize=(12,8))
        plt.title('Regression to dLight in Time from: \n'+str(self.dat.rat.unique())+\
            ", "+str(self.dat.fiberloc.unique())+\
            ', and '+str(len(self.dat.session.unique()))+"sessions"\
              ,fontsize='xx-large',fontweight='bold')
        toplt = rweights
        for r in range(len(self.regFeatures)):
            plt.plot(self.regLags/self.fs,toplt[r,:],label=self.regFeatures[r])
        plt.ylabel('Regression Weight',fontsize='large',fontweight='bold')
        plt.axhline(y=0, color="k", linestyle=":")
        plt.axvline(x=0, color="k", linestyle="--",alpha=.5)
        ax = plt.gca()
        plt.legend()
        return fig

    def calcRegWeightsInTimeFromPort(self):
        faclen = len(self.regFeatures)
        rweights = np.zeros((faclen,len(self.regLags)))
        for n in range(len(self.regLags)):
            y = self.dat.loc[self.dat_visinds+self.regLags[n],self.plot_trace]
            y = y.reset_index(drop=True)
            X = pd.DataFrame()
            for f in self.regFeatures:
                if self.reg_to_portval:
                    X[f] = self.dat.loc[self.dat_visinds,f].values
                else:
                    X[f] = self.dat.loc[self.dat_visinds+self.regLags[n],f].values
            X[self.regFeatures] = scale.fit_transform(X[self.regFeatures].values)
            X = X.drop(y.loc[y.isnull()].index,axis=0)
            y = y.drop(y.loc[y.isnull()].index,axis=0)
            try:
                y = y.drop(X.loc[X.isnull().values].index,axis=0)
                X = X.drop(X.loc[X.isnull().values].index,axis=0)
            except:
                None #print("no null values in feature space")
            modo = LR(fit_intercept=True,normalize=False).fit(X,y)
            rweights[:,n] = modo.coef_
        return rweights

    def plotRegWeightsInTimeFromPort(self):
        self.regLags = np.arange(self.fs*self.plot_window[0],self.fs*self.plot_window[1])
        rweights = []
        for s in self.df.session.unique():
            self.dat = self.df.loc[self.df.session==s]
            self.dat_visinds = self.dat.loc[self.dat.port!=-100].index
            rweights.append(self.calcRegWeightsInTimeFromPort())
        rweights = np.array(rweights)
        toplt = np.mean(rweights,axis=0)
        fig = plt.figure(figsize=(12,8))
        plt.title('Regression to dLight in Time. '+str(len(self.df.rat.unique()))+\
            ' animals, '+str(len(self.df.session.unique()))+' sessions.'\
              ,fontsize='xx-large',fontweight='bold')
        for r in range(len(self.regFeatures)):
            plt.plot(self.regLags/self.fs,toplt[r,:],label=self.regFeatures[r])
        plt.ylabel('Regression Weight',fontsize='large',fontweight='bold')
        plt.axhline(y=0, color="k", linestyle=":")
        plt.axvline(x=0, color="k", linestyle="--",alpha=.5)
        for i in range(len(self.regFeatures)):
            plt.fill_between(self.regLags/self.fs, toplt[i,:]-sem(rweights)[i,:],\
                toplt[i,:]+sem(rweights)[i,:],color="gray",alpha=0.2)
        plt.xlabel('Time from port entry (s)',fontsize='x-large',fontweight='bold')
        plt.legend()
        fig.savefig(self.directory_prefix+'portDaRegPlot.pdf')

    def create_hexDf(self):
        '''creates a hex-level dataframe, where each row is a
        traversed hex in the maze, between port entries.'''
        hexArray = [[] for _ in range(13)]
        Qs = np.array([0,0,0])
        nomRwds = np.array([0,0,0])
        lens = np.array([0,0,0])
        dstop = np.array([0])
        for sesh in self.df.session.unique():
            r = self.df.loc[self.df.session==sesh,'rat'].values[0]
            vinds = np.concatenate([[self.df.loc[(self.df.session==sesh)].index[0]],\
                self.df.loc[(self.df.session==sesh)&(self.df.port!=-100)].index])
            lastHex = None
            seshType =  self.df.loc[self.df.session==sesh,'session_type'].values[0]
            date =  self.df.loc[self.df.session==sesh,'date'].values[0]
            for t in range(1,len(vinds)):
                tdat = self.df.loc[vinds[t-1]+1:vinds[t]]
                tdat = tdat.loc[tdat.lastleave==2]
                if len(tdat) <2:
                    continue
                try:
                    tdat.loc[tdat.port!=-100,"port"].values[0]
                except:
                    continue
                b = tdat.block.unique()[0]
                triQs = tdat.loc[:,["Q_a","Q_b","Q_c"]].values[0]
                triNomRwds = tdat.loc[:,["nom_rwd_a","nom_rwd_b","nom_rwd_c"]].values[0]
                triPathLens = tdat.loc[:,["lenAB","lenBC","lenAC"]].values[0]
                dToP = tdat.loc[:,"dtop"].values[0]
                sesh = tdat.session.unique()[0]
                rwd = tdat.loc[tdat.port!=-100,"rwd"].values[0]
                tinds = tdat.loc[tdat.hexlabels.diff()!=0].index
                prwd = tdat.nom_rwd_chosen.unique()[0]
                hexes = tdat.loc[tinds,'hexlabels'].values
                if hexes[0] ==lastHex:#in [1,2,3]:
                    tinds = tinds[1:]
                    hexes = hexes[1:]
                if t == 1:
                    tinds = tinds[1:]
                    hexes = hexes[1:]
                lastHex = hexes[-1]
                #create new numerical column counting number of hex transitions
                tdat.loc[:,"newHex"] = -1
                cnt = 0
                for h in tinds:
                    tdat.loc[h,"newHex"] = cnt
                    cnt += 1 
                tdat.loc[:,"newHex"] = tdat.loc[:,"newHex"].replace(-1,method="ffill")
                dvals = tdat.loc[tdat.newHex!=-1,:].groupby("newHex").mean().green_z_scored.values
                velvals = tdat.loc[tdat.newHex!=-1,:].groupby("newHex").mean().vel.values
                accvals = tdat.loc[tdat.newHex!=-1,:].groupby("newHex").mean().acc.values
                hexArray[0] = hexArray[0] + list(hexes)
                hexArray[1] = hexArray[1] + list(np.full(len(hexes),-100,dtype="int8"))
                hexArray[2] = hexArray[2] + list(np.full(len(hexes),0,dtype="int8"))
                hexArray[3] = hexArray[3] + list(np.full(len(dvals),sesh,dtype="int8"))
                hexArray[4] = hexArray[4] + list(dvals)
                hexArray[5] = hexArray[5] + list(np.full(len(hexes),b,dtype="int8"))
                hexArray[6] = hexArray[6] + list(np.full(len(hexes),t-1,dtype="int16"))
                hexArray[7] = hexArray[7] + list(np.full(len(hexes),r))
                hexArray[8] = hexArray[8] + list(np.full(len(hexes),date))
                hexArray[9] = hexArray[9] + list(np.full(len(hexes),seshType))
                hexArray[10] = hexArray[10] + list(np.full(len(hexes),prwd))
                hexArray[11] = hexArray[11] + list(velvals)
                hexArray[12] = hexArray[12] + list(accvals)
                Qs = np.vstack([Qs,np.full((len(hexes),3),triQs)])
                nomRwds = np.vstack([nomRwds,np.full((len(hexes),3),triNomRwds)])
                lens = np.vstack([lens,np.full((len(hexes),3),triPathLens)])
                dstop = np.vstack([dstop,np.full((len(hexes),1),dToP)])
                hexArray[1][-1] = tdat.loc[tdat.port!=-100,"port"].values[0]
                hexArray[2][-1] = tdat.loc[tdat.port!=-100,"rwd"].values[0]
        hexArray = np.transpose(hexArray)
        hexDf = pd.DataFrame(hexArray,columns=["hexlabel","port","rwd","session","DA","block",
                               "trial","rat","date","session_type","nom_rwd_chosen","vel","acc"])
        hexDf.loc[:,"hexlabel"] = hexDf.loc[:,"hexlabel"].astype(float).astype("int8")
        hexDf.loc[:,"port"] = hexDf.loc[:,"port"].astype(float).astype("int8")
        hexDf.loc[:,"rwd"] = hexDf.loc[:,"rwd"].astype(float).astype("int8")
        hexDf.loc[:,"session"] = hexDf.loc[:,"session"].astype(float).astype("int8")
        hexDf.loc[:,"block"] = hexDf.loc[:,"block"].astype(float).astype("int8")
        hexDf.loc[:,"trial"] = hexDf.loc[:,"trial"].astype(float).astype("int16")
        hexDf.loc[:,"date"] = hexDf.loc[:,"date"].astype(float).astype(int)
        hexDf.loc[:,"DA"] = hexDf.loc[:,"DA"].astype(float)
        hexDf.loc[:,["Q_a","Q_b","Q_c"]] = Qs[1:]
        hexDf.loc[:,["lenAB","lenBC","lenAC"]] = lens[1:]
        hexDf.loc[:,["nom_rwd_a","nom_rwd_b","nom_rwd_c"]] = nomRwds[1:]
        hexDf.loc[:,"dtop"] = dstop[1:]
        print("converting hexlabels to hex + direction (pairedHexState)")
        hexDf.loc[:,"pairedHexState"] = self.convert2pairedState(hexDf.hexlabel.values)
        self.hexDf = hexDf

    def create_hexDf_startEndRepeatOK(self):
        '''edited function to allow starts and ends of trials to be the same corner-hex value'''
        hexArray = [[] for _ in range(13)]
        Qs = np.array([0,0,0])
        nomRwds = np.array([0,0,0])
        lens = np.array([0,0,0])
        dstop = np.array([0])
        for sesh in self.df.session.unique():
            r = self.df.loc[self.df.session==sesh,'rat'].values[0]
            vinds = np.concatenate([[self.df.loc[(self.df.session==sesh)].index[0]],\
                self.df.loc[(self.df.session==sesh)&(self.df.port!=-100)].index])
            lastPortHex = -1
            seshType =  self.df.loc[self.df.session==sesh,'session_type'].values[0]
            date =  self.df.loc[self.df.session==sesh,'date'].values[0]
            for t in range(1,len(vinds)):
                tdat = self.df.loc[vinds[t-1]+1:vinds[t]]
                tdat = tdat.loc[tdat.lastleave==2]
                if len(tdat) <2:
                    continue
                try:
                    tdat.loc[tdat.port!=-100,"port"].values[0]
                except:
                    continue
                b = tdat.block.unique()[0]
                triQs = tdat.loc[:,["Q_a","Q_b","Q_c"]].values[0]
                triNomRwds = tdat.loc[:,["nom_rwd_a","nom_rwd_b","nom_rwd_c"]].values[0]
                triPathLens = tdat.loc[:,["lenAB","lenBC","lenAC"]].values[0]
                dToP = tdat.loc[:,"dtop"].values[0]
                sesh = tdat.session.unique()[0]
                rwd = tdat.loc[tdat.port!=-100,"rwd"].values[0]
                tinds = tdat.loc[tdat.hexlabels.diff()!=0].index
                prwd = tdat.nom_rwd_chosen.unique()[0]
                hexes = tdat.loc[tinds,'hexlabels'].values
                tdat.loc[:,"newHex"] = -1
                cnt = 0
                if t == 1:
                    tinds = tinds[1:]
                    hexes = hexes[1:]
                #create new numerical column counting number of hex transitions
                for h in tinds:
                    tdat.loc[h,"newHex"] = cnt
                    cnt += 1
                if hexes[-1] not in [1,2,3]:
                    portHex = tdat.loc[tdat.port!=-100,"port"].values[0] + 1
                    tdat.loc[tdat.index.max()-17,"newHex"] = cnt
                    cnt += 1
                    tinds = np.concatenate([tinds,[tdat.index.max()-17]])
                    hexes = np.concatenate([hexes,[portHex]])
                tdat.loc[:,"newHex"] = tdat.loc[:,"newHex"].replace(-1,method="ffill")
                dvals = tdat.loc[tdat.newHex!=-1,:].groupby("newHex").mean().green_z_scored.values
                velvals = tdat.loc[tdat.newHex!=-1,:].groupby("newHex").mean().vel.values
                accvals = tdat.loc[tdat.newHex!=-1,:].groupby("newHex").mean().acc.values
                hexArray[0] = hexArray[0] + list(hexes)
                hexArray[1] = hexArray[1] + list(np.full(len(hexes),-100,dtype="int8"))
                hexArray[2] = hexArray[2] + list(np.full(len(hexes),0,dtype="int8"))
                hexArray[3] = hexArray[3] + list(np.full(len(dvals),sesh,dtype="int8"))
                hexArray[4] = hexArray[4] + list(dvals)
                hexArray[5] = hexArray[5] + list(np.full(len(hexes),b,dtype="int8"))
                hexArray[6] = hexArray[6] + list(np.full(len(hexes),t-1,dtype="int16"))
                hexArray[7] = hexArray[7] + list(np.full(len(hexes),r))
                hexArray[8] = hexArray[8] + list(np.full(len(hexes),date))
                hexArray[9] = hexArray[9] + list(np.full(len(hexes),seshType))
                hexArray[10] = hexArray[10] + list(np.full(len(hexes),prwd))
                hexArray[11] = hexArray[11] + list(velvals)
                hexArray[12] = hexArray[12] + list(accvals)
                Qs = np.vstack([Qs,np.full((len(hexes),3),triQs)])
                nomRwds = np.vstack([nomRwds,np.full((len(hexes),3),triNomRwds)])
                lens = np.vstack([lens,np.full((len(hexes),3),triPathLens)])
                dstop = np.vstack([dstop,np.full((len(hexes),1),dToP)])
                hexArray[1][-1] = tdat.loc[tdat.port!=-100,"port"].values[0]
                hexArray[2][-1] = tdat.loc[tdat.port!=-100,"rwd"].values[0]
        hexArray = np.transpose(hexArray)
        hexDf = pd.DataFrame(hexArray,columns=["hexlabel","port","rwd","session","DA","block",
                               "trial","rat","date","session_type","nom_rwd_chosen","vel","acc"])
        hexDf.loc[:,"hexlabel"] = hexDf.loc[:,"hexlabel"].astype(float).astype("int8")
        hexDf.loc[:,"port"] = hexDf.loc[:,"port"].astype(float).astype("int8")
        hexDf.loc[:,"rwd"] = hexDf.loc[:,"rwd"].astype(float).astype("int8")
        hexDf.loc[:,"session"] = hexDf.loc[:,"session"].astype(float).astype("int8")
        hexDf.loc[:,"block"] = hexDf.loc[:,"block"].astype(float).astype("int8")
        hexDf.loc[:,"trial"] = hexDf.loc[:,"trial"].astype(float).astype("int16")
        hexDf.loc[:,"date"] = hexDf.loc[:,"date"].astype(float).astype(int)
        hexDf.loc[:,"DA"] = hexDf.loc[:,"DA"].astype(float)
        hexDf.loc[:,["Q_a","Q_b","Q_c"]] = Qs[1:]
        hexDf.loc[:,["lenAB","lenBC","lenAC"]] = lens[1:]
        hexDf.loc[:,["nom_rwd_a","nom_rwd_b","nom_rwd_c"]] = nomRwds[1:]
        hexDf.loc[:,"dtop"] = dstop[1:]
        print("converting hexlabels to hex + direction (pairedHexState)")
        hexDf.loc[:,"pairedHexState"] = self.convert2pairedState(hexDf.hexlabel.values)
        self.hexDf = hexDf  

