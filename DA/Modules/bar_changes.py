"""Object to analyze data around barrier changes in
the maze from Krausz et al."""

__author__ = "Tim Krausz"
__email__ = "krausz.tim@gmail.com"
__status__ = "development"
__license__ = "MIT"


from single_rat_da import *

class BarChanges(Photrat):

    def get_newlyAvailHexesBySesh(self):
        self.check_seshTmatsExist()
        barseshs = self.df.loc[self.df.session_type=='barrier','session'].unique()
        self.sesh_newlyAvailHexes = {s:[] for s in barseshs}
        self.find_newlyAvailHexes()
        self.add_adjacent2newlyAvail()

    def get_newlyBlockedHexesBySesh(self):
        self.check_seshTmatsExist()
        barseshs = self.df.loc[self.df.session_type=='barrier','session'].unique()
        self.sesh_newlyBlockedHexes = {s:[] for s in barseshs}
        self.find_newlyBlockedHexes()
        self.add_adjacent2newlyBlocked()

    def find_newlyAvailHexes(self):
        for s in self.sesh_newlyAvailHexes:
            if len(self.sesh_tmats[s])==1:
                continue
            for b in range(len(self.sesh_tmats[s])-1):
                testmat = self.sesh_tmats[s][b]
                testmat1 = self.sesh_tmats[s][b+1]
                self.sesh_newlyAvailHexes[s].append(np.where((\
                    np.sum(testmat,axis=0)==0)&(np.sum(testmat1,axis=0)!=0))[0])

    def find_newlyBlockedHexes(self):
        for s in self.sesh_newlyBlockedHexes:
            for b in range(len(self.sesh_tmats[s])-1):
                testmat = self.sesh_tmats[s][b]
                testmat1 = self.sesh_tmats[s][b+1]
                self.sesh_newlyBlockedHexes[s].append(np.where((\
                    np.sum(testmat,axis=0)!=0)&(np.sum(testmat1,axis=0)==0))[0])

    def add_newlyAvailHexesToDf(self):
        self.df.loc[:,"newlyAvailHex"]=0
        self.df.loc[:,"newlyAvailHex"] = self.df.loc[:,"newlyAvailHex"].astype("int8")
        for s in self.sesh_newlyAvailHexes:
            for b in range(len(self.sesh_newlyAvailHexes[s])):
                newHex = self.sesh_newlyAvailHexes[s][b]
                dat = self.df.loc[(self.df.session==s)&(self.df.block==b+2)]
                inNewHex = np.where(dat.loc[:,'hexlabels'].isin(\
                    newHex)==1)[0]+dat.index.min()
                self.df.loc[inNewHex,"newlyAvailHex"] = 1

    def add_adjacent2newlyAvail(self):
        self.df.loc[:,"adj2newlyAvail"]=0
        self.df.loc[:,"adj2newlyAvail"] = self.df.loc[:,"adj2newlyAvail"].astype("int8")
        for s in self.sesh_newlyAvailHexes:
            for b in range(len(self.sesh_newlyAvailHexes[s])):
                newHex = self.sesh_newlyAvailHexes[s][b]
                if len(newHex)<1:
                    continue
                adjHexes = []
                for h in newHex:
                    adjHexes = np.concatenate([adjHexes,np.where(self.sesh_tmats[s][b+1][h]>0)[0]])
                dat = self.df.loc[(self.df.session==s)&(self.df.block==b+2)]
                inAdjHex = np.where(dat.loc[:,'hexlabels'].isin(\
                    adjHexes)==1)[0]+dat.index.min()
                self.df.loc[inAdjHex,"adj2newlyAvail"] = 1

    def add_newlyBlockedHexesToDf(self):
        self.df.loc[:,"newlyBlockedHex_Adj"]=0
        self.df.loc[:,"newlyBlockedHex_Adj"] = self.df.loc[:,"newlyBlockedHex_Adj"].astype("int8")
        for s in self.sesh_newlyBlockedHexes:
            if len(self.sesh_newlyBlockedHexes[s])==1:
                continue
            for b in range(len(self.sesh_newlyBlockedHexes[s])):
                newHex = self.sesh_newlyBlockedHexes[s][b]
                dat = self.df.loc[(self.df.session==s)&(self.df.block==b+2)]
                inNewHex = np.where(dat.loc[:,'hexlabels'].isin(\
                    newHex)==1)[0]+dat.index.min()
                self.df.loc[inNewHex,"newlyBlockedHex_Adj"] = 1

    def add_adjacent2newlyBlocked(self):
        self.df.loc[:,"adj2newlyBlocked"]=0
        self.df.loc[:,"adj2newlyBlocked"] = self.df.loc[:,"adj2newlyBlocked"].astype("int8")
        for s in self.sesh_newlyBlockedHexes:
            for b in range(len(self.sesh_newlyBlockedHexes[s])):
                newHex = self.sesh_newlyBlockedHexes[s][b]
                if len(newHex)<1:
                    continue
                adjHexes = []
                for h in newHex:
                    adjHexes = np.concatenate([adjHexes,np.where(self.sesh_tmats[s][b][h]>0)[0]])
                dat = self.df.loc[(self.df.session==s)&(self.df.block==b+2)]
                inAdjHex = np.where(dat.loc[:,'hexlabels'].isin(\
                    adjHexes)==1)[0]+dat.index.min()
                self.df.loc[inAdjHex,"adj2newlyBlocked"] = 1