"""Object to help with port-value analysis"""

__author__ = "Tim Krausz"
__email__ = "krausz.tim@gmail.com"
__status__ = "development"
__license__ = "MIT"


from single_rat_da import *
from sklearn.preprocessing import StandardScaler,MinMaxScaler
scale = StandardScaler()
valscale = MinMaxScaler()
scaler = valscale.fit(np.array([-1,0,1]).reshape(-1,1))
from Q_learning_funcs import *

class PortValAnalyses(Photrat):

    def get_portQvals(self,qtype='hybrid',level="session"):
        '''Using the parameters optimized in Julia, estimate values for each port,
        on each trial. Specify which type of q value from ['hybrid','mb','mf','port'].
        Specify whether parameters were drawn from the level of rat or session.'''
        self.df.loc[:,'Q_a'] = np.nan
        self.df.loc[:,'Q_b'] = np.nan
        self.df.loc[:,'Q_c'] = np.nan
        self.load_q_params(qtype)
        if level=="session":
            seshs = self.triframe.session.unique()
            for s in range(len(seshs)):
                tsesh = self.triframe.loc[self.triframe.session==seshs[s]]
                if qtype=="hybrid":
                    Q = hyb_q_choice(tsesh,self.hyb_params.loc[s].values)
                elif qtype=="mb":
                    Q = mb_q_choice(tsesh,self.mb_params.loc[s].values)
                elif qtype=="mf":
                    Q = mf_q_choice(tsesh,self.mf_params.loc[s].values)
                elif qtype == "port":
                    Q = port_q_choice(tsesh,self.q_params.loc[s].values)
                self.addQs2dfSubset(Q,seshs[s],useSesh=True)
                
        elif level=="rat":
            rats = self.triframe.rat.unique()
            for r in range(len(rats)):
                #should maybe still do this by sesh but keep same parameters
                for s in self.triframe.loc[self.triframe.rat==rats[r],"session"].unique():
                    tsesh = self.triframe.loc[self.triframe.session==s]#rat==rats[r]]
                    if qtype=="hybrid":
                        Q = hyb_q_choice(tsesh,self.hyb_params.loc[r].values)
                    elif qtype=="mb":
                        Q = mb_q_choice(tsesh,self.mb_params.loc[r].values)
                    elif qtype=="mf":
                        Q = mf_q_choice(tsesh,self.mf_params.loc[r].values)
                    elif qtype == "port":
                        Q = port_q_choice(tsesh,self.q_params.loc[r].values)
                    #self.addQs2dfSubset(Q,rats[r],useSesh=False)
                    self.addQs2dfSubset(Q,s,useSesh=True)
        self.df.loc[:,'Q_a'] = self.df.loc[:,'Q_a'].fillna(method='bfill').astype("float16")
        self.df.loc[:,'Q_b'] = self.df.loc[:,'Q_b'].fillna(method='bfill').astype("float16")
        self.df.loc[:,'Q_c'] = self.df.loc[:,'Q_c'].fillna(method='bfill').astype("float16")

    def addQs2dfSubset(self,Qs2add,levelIndex,useSesh=True):
        if useSesh:
            poolString = "session"
        else:
            poolString = "rat"
        self.df.loc[(self.df[poolString]==levelIndex)&(self.df.port!=-100),'Q_a']=Qs2add[:,0]
        self.df.loc[(self.df[poolString]==levelIndex)&(self.df.port!=-100),'Q_b']=Qs2add[:,1]
        self.df.loc[(self.df[poolString]==levelIndex)&(self.df.port!=-100),'Q_c']=Qs2add[:,2]

    def load_q_params(self,qtype='hybrid'):
        if qtype == 'hybrid':
            self.hyb_params = pd.read_csv(self.directory_prefix+"tri_hybrid_params.csv")
        elif qtype == 'port':
            self.q_params = pd.read_csv(self.directory_prefix+"tri_q3port_params.csv")
        
    def get_vals_byChosenEtc(self,chosen_only=False):
        for s in self.df.session.unique():
            self.dat = self.df.loc[self.df.session==s,:].copy()
            if chosen_only:
                self.get_valOfChosenPort()
            else:
                self.get_vals_by_portType()
            if s==self.df.session[0]:
                newdf = self.dat
            else:
                newdf = newdf.append(self.dat,ignore_index=True)
        self.df = newdf

    def get_qRPE(self):
        r = self.triframe.loc[:,'rwd'].values
        portQs = self.df.loc[self.visinds,'Q_chosen'].values
        self.triframe.loc[:,"q_rpe"] = r - scaler.transform(portQs.reshape(-1,1)).T[0]

    def get_avail(self,p):
        allports = [0,1,2]
        allports.pop(p)
        return allports

    def get_lr_dif_val(self,data,factor):
        '''returns list of factor values for left choice for each trial.'''
        ports = data.loc[data.port!=-100,'port'].values
        vinds = data.loc[data.port!=-100].index
        fac = [0]
        if factor=='dist':
            for p in range(1,len(ports)): #append values of going left
                if ports[p-1]==2:
                    fac.append(data.loc[vinds[p],'lenBC']-data.loc[vinds[p],'lenAC'])
                elif ports[p-1]==1:
                    fac.append(data.loc[vinds[p],'lenAB']-data.loc[vinds[p],'lenBC'])
                elif ports[p-1]==0:
                    fac.append(data.loc[vinds[p],'lenAC']-data.loc[vinds[p],'lenAB'])
        else:
            for p in range(1,len(ports)): #append values of going left
                if ports[p-1]==2:
                    fac.append(data.loc[vinds[p],factor+'_b']-data.loc[vinds[p],factor+'_a'])
                elif ports[p-1]==1:
                    fac.append(data.loc[vinds[p],factor+'_a']-data.loc[vinds[p],factor+'_c'])
                elif ports[p-1]==0:
                    fac.append(data.loc[vinds[p],factor+'_c']-data.loc[vinds[p],factor+'_b'])
        return fac

    def get_left_val(self,data,factor):
        '''returns list of factor values for left choice for each trial.'''
        ports = data.loc[data.port!=-100,'port'].values
        vinds = data.loc[data.port!=-100].index
        fac = [0]
        if factor=='dist':
            for p in range(1,len(ports)): #append values of going left
                if ports[p-1]==2:
                    fac.append(data.loc[vinds[p],'lenBC'])
                elif ports[p-1]==1:
                    fac.append(data.loc[vinds[p],'lenAB'])
                elif ports[p-1]==0:
                    fac.append(data.loc[vinds[p],'lenAC'])
        else:
            for p in range(1,len(ports)): #append values of going left
                if ports[p-1]==2:
                    fac.append(data.loc[vinds[p],factor+'_b'])
                elif ports[p-1]==1:
                    fac.append(data.loc[vinds[p],factor+'_a'])
                elif ports[p-1]==0:
                    fac.append(data.loc[vinds[p],factor+'_c'])
        return fac
    
    def tri_rr(self,data,ports=[0,1,2],t=5):
        '''Return trial-by-trial estimte of previous reward history over past t trials'''
        if len(ports)==1:
            rwds = data.loc[data.port==ports[0],'rwd'].copy()
            rwds.loc[rwds==-1]=0
            rwds = rwds.rolling(window=t).sum()
            rwds = rwds.reindex(np.arange(data.index[0],data.index[-1]+1),fill_value=-100)
            return rwds.replace(-100,method='ffill')
        elif len(ports)==2:
            rwds = data.loc[(data.port==ports[0])|(data.port==ports[1]),'rwd'].copy()
            rwds.loc[rwds==-1]=0
            rwds = rwds.rolling(window=t).sum()
            rwds = rwds.reindex(np.arange(data.index[0],data.index[-1]+1),fill_value=-100)
            return rwds.replace(-100,method='ffill')
        elif len(ports)==3:
            rwds = data.loc[data.port!=-100,'rwd'].copy()
            rwds.loc[rwds==-1]=0
            rwds = rwds.rolling(window=t).sum()
            rwds = rwds.reindex(np.arange(data.index[0],data.index[-1]+1),fill_value=-100)
            return rwds.replace(-100,method='ffill')
    
    def get_log_pchoos_v_costNben(self):
        df = self.triframe.copy()
        df.rat = df.rat.astype('category')
        df['ratcodes'] = df.rat.cat.codes
        seshs=df.session.unique()
        for s in range(len(seshs)):
            sdf = df.loc[(df.session==seshs[s])]
            rdf = pd.DataFrame({'rdif':self.get_lr_dif_val(sdf,'nom_rwd'),\
                'ldif':self.get_lr_dif_val(sdf,'dist')})
            rdf['rhist_dif'] = get_lr_dif_val(sdf,'rhist')
            rdf['choose_L'] = sdf.lrchoice.values
            rdf['session']=s
            rdf['rat'] = sdf.ratcodes.values
            rdf['tri'] = sdf.tri.values
            rdf['block'] = sdf.block.values
            if s == 0:
                self.regdf = rdf
            else:
                self.regdf = self.regdf.append(rdf,ignore_index=True)
        self.regdf.loc[regdf.choose_L==2,'choose_L']=np.nan
        
    def plot_log_pchoos_v_costNben(self):
        plt.figure()
        plt.subplot(121)
        plt.title('p(choose L) vs len dif',fontsize='x-large',fontweight='bold')
        sns.regplot(x='ldif',y='choose_L',data=regdf,logistic=True)#,scatter=False)
        plt.ylim(.1,.9)
        plt.xlabel('distance dif',fontsize='large',fontweight='bold')
        plt.ylabel('choose L',fontsize='large',fontweight='bold')
        plt.subplot(122)
        #plt.twinx(plt.gca())
        plt.title('p(choose L) vs rwd dif',fontsize='x-large',fontweight='bold')
        sns.regplot(x='rhist_dif',y='choose_L',data=regdf.loc[(regdf.rhist_dif!=-100)&(regdf.rhist_dif!=100)],\
                    logistic=True,color='orange')#,scatter=False)
        plt.xlabel('rwd hist dif',fontsize='large',fontweight='bold')
        plt.ylabel('')
        plt.ylim(.1,0.9)
        

    def avg_factor(self,portz,visinds):
        #isinds = self.dat.loc[self.dat.port!=-100].index
        df = self.dat.loc[visinds,[self.factor+'_a',self.factor\
        +'_b',self.factor+'_c']].values
        fval = []
        for i in range(len(portz)):
            fval.append(np.mean(df[i,portz[i]]))
        return fval

    def factor_by_p_type(self,p_type='all'):
        '''p_type can be all, avail, or chosen. returns value of factor given p_type for every trial.'''
        visinds = self.dat.loc[self.dat.port!=-100].index
        if p_type == 'chosen':
            portz = self.dat.loc[visinds,'port'].values.astype(int)
        elif p_type == 'avail':
            portz = [[0,1,2]]
            for p in self.dat.loc[visinds[:-1],'port'].values.astype(int):
                portz.append(self.get_avail(p))
            visinds = np.concatenate([visinds[1:],[visinds[-1]]])
        elif p_type == 'all':
            portz = np.tile([0,1,2],(len(visinds),1))
        elif p_type == 'other':
            portz = [0]
            for p in range(len(visinds)-1):
                avail = self.get_avail(int(self.dat.loc[visinds[p],'port']))
                chos = int(self.dat.loc[visinds[p+1],'port'])
                avail.remove(chos)
                portz.append(avail)
            visinds = np.concatenate([visinds[1:],[visinds[-1]]])
        return self.avg_factor(portz,visinds)

    def get_vals_by_portType(self,dattype='float16'):
        addto = self.dat.loc[self.dat.port!=-100].index
        bytype = pd.DataFrame()
        bytype[self.factor+'_avail'] = self.factor_by_p_type('avail')
        bytype[self.factor+'_chosen'] = self.factor_by_p_type('chosen')
        bytype[self.factor+'_all'] = self.factor_by_p_type('all')
        bytype[self.factor+'_other'] = self.factor_by_p_type('other')
        bytype = bytype.set_index(addto)
        bytype = bytype.reindex(np.arange(self.dat.index[0],self.dat.index[-1]+1),fill_value=np.nan)
        
        self.dat[self.factor+'_avail'] = bytype.loc[:,self.factor+'_avail'].\
        fillna(method='ffill').astype(dattype)
        self.dat[self.factor+'_other'] = bytype.loc[:,self.factor+'_other'].\
        fillna(method='ffill').astype(dattype)
        self.dat[self.factor+'_chosen'] = bytype.loc[:,self.factor+'_chosen'].\
        fillna(method='bfill').astype(dattype)
        self.dat[self.factor+'_all'] = bytype.loc[:,self.factor+'_all'].\
        fillna(method='bfill').astype(dattype)

    def get_valOfChosenPort(self,dattype='float16'):
        addto = self.dat.loc[self.dat.port!=-100].index
        bytype = pd.DataFrame()
        bytype[self.factor+'_chosen'] = self.factor_by_p_type('chosen')
        bytype = bytype.set_index(addto)
        bytype = bytype.reindex(np.arange(self.dat.index[0],self.dat.index[-1]+1),fill_value=np.nan)
        self.dat[self.factor+'_chosen'] = bytype.loc[:,self.factor+'_chosen'].\
        fillna(method='bfill').astype(dattype)