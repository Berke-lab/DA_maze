"""Object with methods for transition matrix operations"""

__author__ = "Tim Krausz"
__email__ = "krausz.tim@gmail.com"
__status__ = "development"
__license__ = "MIT"


from single_rat_da import *
from math import log as mlog

class TmatOperations(Photrat):

    def get_critChoicePoints(self):
        #depreciated, better to use manual identification and entry
        self.check_seshTmatsExist()
        self.crit_cps = {s:[] for s in self.df.session.unique()}
        self.sesh_deadEnds = {s:[] for s in self.df.session.unique()}
        self.deadEndCps = {s:[] for s in self.df.session.unique()}
        self.find_cps4allSessions()
        self.add_cps2Df()

    def get_allChoicePoints(self):
        self.check_seshTmatsExist()
        self.all_cps = {s:[] for s in self.df.session.unique()}
        self.find_allCpsBySesh()
        self.add_allCps2Df()

    def check_seshTmatsExist(self):
        try:
            self.sesh_tmats
        except:
            print("transition matrices not yet loaded. Loading now...")
            self.load_tmats()

    def load_tmats(self):
        self.sesh_tmats = {s:[] for s in self.df.session.unique()}
        for sesh in self.df.session.unique():
            self.attempt_loadTmats(sesh)

    def attempt_loadTmats(self,sesh):
        rat = self.df.loc[self.df.session==sesh,'rat'].unique()[0]
        date = str(self.df.loc[self.df.session==sesh,'date'].unique()[0])
        if int(date[0])!=1:
            date = '0'+date
        if self.df.loc[self.df.session==sesh,'session_type'].values[0]=='prob':
            self.sesh_tmats[sesh] = np.load(self.directory_prefix+rat+'/'+\
                    date+'/'+'tmat.npy')
        else:
            for b in self.df.loc[self.df.session==sesh,'block'].unique():
                try:
                    self.sesh_tmats[sesh].append(np.load(self.directory_prefix+rat+'/'+\
                            date+'/'+'tmat_block_'+str(b)+'.0.npy'))
                except:
                    print('no tmat saved for sesh '+str(sesh)+' block '+str(b))


    def get_simpleStateTransitions(self,s,b):
        if self.df.loc[self.df.session==s,'session_type'].values[0]=='prob':
            self.bars = self.sesh_barIDs[s]
        else:
            self.tmat = self.tmats[s][int(b-1)]
            self.bars = self.sesh_barIDs[s][int(b-1)]
        self.add_bars2tmatrix50()
        self.to_state = [[0,0,0,0,0,0]]+[np.argmax(a,axis=1) \
            for a in self.bar_tmatrix]

    def add_bars2tmatrix(self):
        self.bar_tmatrix = tmatrix.copy()
        self.bar_tmatrix[self.bars] = np.tile(emptp,(6,1))

    def add_bars2tmatrix_pairedState(self):
        self.bar_tmatrix = paired_tmatrix.copy()
        self.paired_barstates = phexdf.loc[(np.isin(phexdf.from_state,self.bars+1))|
        (np.isin(phexdf.to_state,self.bars+1)),'statecodes'].values
        for b in self.paired_barstates:
            barloc = np.where(adj_matx==b)
            if np.shape(barloc)[1]==1:
                self.bar_tmatrix[np.array(barloc)[0,0],np.array(barloc)[1,0]]=np.zeros(126)
            else:
                for bind in range(len(np.array(barloc))):
                   self.bar_tmatrix[np.array(barloc)[0,bind],np.array(barloc)[1,bind]]=np.zeros(126)

    def get_availstates(self):
        self.availstates = np.setdiff1d(np.arange(0,49),self.bars)
        self.availstates = np.setdiff1d(self.availstates,[0,1,2])

    def find_critCps4allSessions(self):
        for s in self.df.session.unique():
            self.dat = self.df.loc[self.df.session==s]
            if self.dat.session_type.unique()[0]=='prob':
                self.tmat = self.sesh_tmats[s]
                self.filter_extraneousCPs()
                self.crit_cps[s].append(self.cps)
                self.deadEndCps[s].append(self.dEndCps)
                self.sesh_deadEnds[s].append(self.deadEnds)
            else:
                for b in self.dat.block.unique():
                    try:
                        self.tmat = self.sesh_tmats[s][b]
                    except:
                        continue
                    self.filter_extraneousCPs()
                    self.crit_cps[s].append(self.cps)
                    self.deadEndCps[s].append(self.dEndCps)
                    self.sesh_deadEnds[s].append(self.deadEnds)

    def find_allCpsBySesh(self):
        for s in self.df.session.unique():
            self.dat = self.df.loc[self.df.session==s]
            if self.dat.session_type.unique()[0]=='prob':
                self.tmat = self.sesh_tmats[s]
                self.find_allCPs()
                self.all_cps[s].append(self.cps)
            else:
                for b in self.dat.block.unique():
                    try:
                        self.tmat = self.sesh_tmats[s][b-1]
                    except:
                        continue
                    self.find_allCPs()
                    self.all_cps[s].append(self.cps)

    def find_allCPs(self):
        self.cps = np.unique(np.where((self.tmat*10).astype("int8")==3)[0])
        
    def find_allDeadEnds(self):
        self.deadEnds = np.unique(np.where(self.tmat==1)[0])
        self.deadEnds = np.delete(self.deadEnds,np.isin(self.deadEnds,[1,2,3]))

    def get_deadEndCPs(self):
        self.dEndCps = []
        for h in self.cps:
            if self.leads2deadEnd(h,[]):
                self.dEndCps.append(h)
        self.cps = np.delete(self.cps,np.isin(self.cps,self.dEndCps))

    def leads2deadEnd(self,h,visited):
        visited.append(h)
        nexthexes = np.where(self.tmat[h]>0)[0]
        toDeadEnd = []
        for hprime in nexthexes:
            if hprime in visited or hprime in self.cps or hprime in [1,2,3]:
                toDeadEnd.append(False)
            elif hprime in self.deadEnds:
                toDeadEnd.append(True)
            else:
                toDeadEnd.append(self.leads2deadEnd(hprime,visited))
        return any(toDeadEnd)

    def add_allCps2Df(self):
        '''adds a choicePoint column to the df that marks
        any moment the rat is in an intersection between hexes.
        Different from critical choice points, marked by the experimenter
        based on the hex configuration'''
        self.df.loc[:,'choicePoint'] = 0
        for s in self.df.session.unique():
            dat = self.df.loc[self.df.session==s]
            if self.df.loc[self.df.session==s,'session_type'].values[0]=='prob':
                incp = dat.loc[dat.loc[:,'hexlabels'].isin(self.all_cps[s][0])].index
                self.df.loc[incp,"choicePoint"] = 1
            else:
                for b in range(1,len(self.all_cps[s])+1):#dat.block.unique():
                    incp = np.where(dat.loc[dat.block==b,'hexlabels'].isin(\
                        self.all_cps[s][int(b-1)])==1)[0]+dat.loc[dat.block==b,:].index.min()
                    self.df.loc[incp,"choicePoint"] = 1

    def get_distanceToPort(self,getDistsFromDeadEnds=True):
        ports = [0,1,2]
        portStrs = ["A","B","C"]
        self.sesh_hexDists = {"dto"+portStrs[p]:{s:[] for s in 
            self.df.session.unique()} for p in ports}
        self.sesh_arrows2goal = {"to"+portStrs[p]:{s:[] for s in 
            self.df.session.unique()} for p in ports}
        for p in ports:
            for s in self.df.session.unique():
                dat = self.df.loc[(self.df.session==s)&(self.df.port!=-100)]
                if self.df.loc[self.df.session==s,'session_type'].values[0]=='prob':
                    self.tmat = self.sesh_tmats[s]
                    self.bars = self.sesh_barIDs[s]
                    dmap = self.compute_distanceToPort(dat,p,getDistsFromDeadEnds)
                    self.sesh_hexDists["dto"+portStrs[p]][s] = dmap
                    togoal = self.find_arrows2goal(dmap,p)
                    self.sesh_arrows2goal["to"+portStrs[p]][s] = togoal
                else:
                    for b in dat.block.unique():
                        try:
                            self.tmat = self.sesh_tmats[s][int(b-1)]
                        except:
                            continue
                        self.bars = self.sesh_barIDs[s][int(b-1)]
                        dmap = self.compute_distanceToPort(dat.loc[dat.block==b],\
                            p,getDistsFromDeadEnds)
                        self.sesh_hexDists["dto"+portStrs[p]][s].append(dmap)
                        togoal = self.find_arrows2goal(dmap,p)
                        self.sesh_arrows2goal["to"+portStrs[p]][s].append(togoal)

    def compute_distanceToPort(self,data,port):
        actions = [0,1]
        distmap = np.full(126,0,dtype="float16")
        portz = [phexdf.loc[phexdf.to_state==1,'statecodes'].values[0],
        phexdf.loc[phexdf.to_state==2,'statecodes'].values[0],
         phexdf.loc[phexdf.to_state==3,'statecodes'].values[0]]
        self.add_bars2tmatrix_pairedState()
        t=1
        avail = np.array([portz[port]])
        i = 0
        oldD = distmap
        portvec = distmap
        portvec[avail]=1
        dfact=0.9
        while True:
            allDs = []
            for a in actions:
                allDs.append(dfact*np.dot(oldD,self.bar_tmatrix[:,a].T)+portvec)
            newD = np.amax(allDs,axis=0)
            if np.max(newD-oldD)<10e-4:
                distmap = np.array([np.nan if newD[i]==0 else\
                    np.round(mlog(newD[i],dfact)).astype(int) for i in range(len(newD))])
                #set barriers to nan
                distmap[self.paired_barstates]=np.nan
                break
            oldD = newD
            i += 1
        return distmap

    def find_arrows2goal(self,distmap,port):
        arrows2goal = np.copy(distmap)
        for i in arrowdf.from_state.unique():
            options = list(arrowdf.loc[arrowdf.from_state==i,"statecodes"].values)
            if all(np.isnan(distmap[options])):
                continue
            #get state with shortest distance from port
            minarrow = np.nanargmin(distmap[options])
            onlyshort = options[minarrow]
            if i-1 == port:
                arrows2goal[onlyshort]=np.nan
            if not np.all(np.isnan(distmap[options])) and \
            np.nanmin(distmap[options])==distmap[onlyshort]:
                options.pop(minarrow)
            arrows2goal[options]=np.nan
        towardsgoal = np.array([not np.isnan(i) for i in arrows2goal])
        return towardsgoal
