
"""Class for processing and analysis of single photometry session.
Also includes code for hex and paired direction-hex state transition matrices"""

__author__ = "Tim Krausz"
__email__ = "krausz.tim@gmail.com"
__status__ = "development"
__license__ = "MIT"

import numpy as np
import pandas as pd
from TriFuncs import make_lr
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression as LR

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
coords = []
cols = [7,6,6,5,5,4,4,3,3,2,2,1,1]
maxrows = 13
r = 0
x = 1
y = 1
startr = 1
while r < maxrows:
    maxcols = cols[r]
    c = 0
    if r%2!=0:
        startr+=1
    x=startr
    while c < maxcols:
        coords.append([x,y])
        x += 2
        c += 1
    if r%2!=0:
        y += 2
    else:
        y+=1
    r += 1
cents = {h: c for h,c in zip(hexlist,coords)}
centdf = pd.DataFrame(cents)
centdf = centdf.T

emptp = np.zeros(49)
tprobs = np.identity(49)
tmatrix = np.array([[emptp,emptp,emptp,emptp,emptp,emptp],[emptp,emptp,emptp,emptp,emptp,emptp],\
        [emptp,emptp,emptp,emptp,emptp,emptp],[tprobs[0],emptp,tprobs[4],emptp,tprobs[5],emptp],\
        [emptp,emptp,emptp,tprobs[6],emptp,tprobs[3]],[emptp,tprobs[3],emptp,tprobs[7],emptp,emptp],\
        [tprobs[4],emptp,tprobs[8],emptp,tprobs[9],emptp],[tprobs[5],emptp,tprobs[9],emptp,tprobs[10],emptp],\
        [emptp,emptp,emptp,tprobs[11],emptp,tprobs[6]],[emptp,tprobs[6],emptp,tprobs[12],emptp,tprobs[7]],\
        [emptp,tprobs[7],emptp,tprobs[13],emptp,emptp],[tprobs[8],emptp,tprobs[14],emptp,tprobs[15],emptp],\
        [tprobs[9],emptp,tprobs[15],emptp,tprobs[16],emptp],[tprobs[10],emptp,tprobs[16],emptp,tprobs[17],emptp],\
        [emptp,emptp,emptp,tprobs[18],emptp,tprobs[11]],\
        [emptp,tprobs[11],emptp,tprobs[19],emptp,tprobs[12]],[emptp,tprobs[12],emptp,tprobs[20],emptp,tprobs[13]],\
        [emptp,tprobs[13],emptp,tprobs[21],emptp,emptp],[tprobs[14],emptp,tprobs[22],emptp,tprobs[23],emptp],\
        [tprobs[15],emptp,tprobs[23],emptp,tprobs[24],emptp],[tprobs[16],emptp,tprobs[24],emptp,tprobs[25],emptp],\
        [tprobs[17],emptp,tprobs[25],emptp,tprobs[26],emptp],[emptp,emptp,emptp,tprobs[27],emptp,tprobs[18]],\
        [emptp,tprobs[18],emptp,tprobs[28],emptp,tprobs[19]],[emptp,tprobs[19],emptp,tprobs[29],emptp,tprobs[20]],\
        [emptp,tprobs[20],emptp,tprobs[30],emptp,tprobs[21]],[emptp,tprobs[21],emptp,tprobs[31],emptp,emptp],\
        [tprobs[22],emptp,tprobs[32],emptp,tprobs[33],emptp],[tprobs[23],emptp,tprobs[33],emptp,tprobs[34],emptp],\
        [tprobs[24],emptp,tprobs[34],emptp,tprobs[35],emptp],[tprobs[25],emptp,tprobs[35],emptp,tprobs[36],emptp],\
        [tprobs[26],emptp,tprobs[36],emptp,tprobs[37],emptp],[emptp,emptp,emptp,tprobs[47],emptp,tprobs[27]],\
        [emptp,tprobs[27],emptp,tprobs[38],emptp,tprobs[28]],[emptp,tprobs[28],emptp,tprobs[39],emptp,tprobs[29]],\
        [emptp,tprobs[29],emptp,tprobs[40],emptp,tprobs[30]],[emptp,tprobs[30],emptp,tprobs[41],emptp,tprobs[31]],\
        [emptp,tprobs[31],emptp,tprobs[48],emptp,emptp],[tprobs[33],emptp,tprobs[42],emptp,tprobs[43],emptp],\
        [tprobs[34],emptp,tprobs[43],emptp,tprobs[44],emptp],[tprobs[35],emptp,tprobs[44],emptp,tprobs[45],emptp],\
        [tprobs[36],emptp,tprobs[45],emptp,tprobs[46],emptp],[emptp,tprobs[47],emptp,emptp,emptp,tprobs[38]],\
        [emptp,tprobs[38],emptp,emptp,emptp,tprobs[39]],[emptp,tprobs[39],emptp,emptp,emptp,tprobs[40]],\
        [emptp,tprobs[40],emptp,emptp,emptp,tprobs[41]],[emptp,tprobs[41],emptp,emptp,emptp,tprobs[48]],\
        [tprobs[32],emptp,tprobs[2],emptp,tprobs[42],emptp],[tprobs[37],emptp,tprobs[46],emptp,tprobs[1],emptp]])

pairedHex_adjacentDict = {1.04:[4.05,4.06],5.04:[4.06,4.01],6.04:[4.01,4.05],4.01:[np.nan,np.nan],
4.06:[6.08,np.nan],4.05:[np.nan,5.07],5.07:[7.09,7.10],6.08:[8.10,8.11],8.06:[np.nan,6.04],7.05:[5.04,np.nan],
10.08:[8.11,8.06],11.08:[8.06,8.10],11.14:[14.17,14.18],14.11:[np.nan,11.08],10.13:
[13.16,13.17],13.10:[10.08,10.07],8.10:[10.07,10.13],7.10:[10.13,10.08],8.11:[11.14,np.nan],
7.09:[np.nan,9.12],9.07:[7.10,7.05],9.12:[12.15,12.16],12.09:[9.07,np.nan],14.18:[18.22,np.nan],
18.14:[14.11,14.17],14.17:[17.13,17.21],17.14:[14.18,14.11],13.17:[17.21,17.14],17.13:
[13.10,13.16],13.16:[16.12,16.20],16.13:[13.17,13.10],12.16:[16.20,16.13],16.12:[12.09,12.15],
12.15:[np.nan,15.19],18.22:[22.26,22.27],22.18:[np.nan,18.14],17.21:[21.25,21.26],
21.17:[17.14,17.13],16.20:[20.24,20.25],20.16:[16.13,16.12],15.19:[19.23,19.24],
19.15:[15.12,np.nan],22.27:[27.32,np.nan],27.22:[22.18,22.26],22.26:[26.21,26.31],
26.22:[22.27,22.18],21.26:[26.31,26.22],26.21:[21.17,21.25],15.12:[12.16,12.09],
21.25:[25.20,25.30],25.21:[21.26,21.17],20.25:[25.30,25.21],25.20:[20.16,20.24],20.24:
[24.19,24.29],24.20:[20.25,20.16],19.24:[24.29,24.20],24.19:[19.15,19.23],19.23:[np.nan,23.28],
23.19:[19.24,19.15],27.32:[32.37,32.38],32.27:[np.nan,27.22],26.31:[31.36,31.37],
31.26:[26.22,26.21],25.30:[30.35,30.36],30.25:[25.21,25.20],24.29:[29.34,29.35],
29.24:[24.20,24.19],23.28:[28.33,28.34],28.23:[23.19,np.nan],32.38:[38.49,np.nan],
38.32:[32.27,32.37],32.37:[37.31,37.42],37.32:[32.38,32.27],31.37:[37.42,37.32],
37.31:[31.26,31.36],31.36:[36.30,36.41],36.31:[31.37,31.26],30.36:[36.41,36.31],
36.30:[30.25,30.35],30.35:[35.29,35.40],35.30:[30.36,30.25],29.35:[35.40,35.30],
35.29:[29.24,29.34],29.34:[34.28,34.39],34.29:[29.35,29.24],28.34:[34.39,34.29],
34.28:[28.23,28.33],28.33:[np.nan,33.48],33.28:[28.34,28.23],38.49:[49.47,49.02],
49.38:[np.nan,38.32],37.42:[42.46,42.47],42.37:[37.32,37.31],36.41:[41.45,41.46],
41.36:[36.31,36.30],35.40:[40.44,40.45],40.35:[35.30,35.29],34.39:[39.43,39.44],
39.34:[34.29,34.28],33.48:[48.03,48.43],49.02:[np.nan,np.nan],2.49:[49.38,49.47],
49.47:[47.42,np.nan],47.49:[49.02,49.38],42.47:[np.nan,47.49],42.46:[46.41,np.nan],
46.42:[42.47,42.37],47.42:[42.37,42.46],46.41:[41.36,41.45],41.46:[np.nan,46.42],
41.45:[45.40,np.nan],45.41:[41.46,41.36],45.40:[40.35,40.44],40.45:[np.nan,45.41],
40.44:[44.39,np.nan],44.40:[40.45,40.35],44.39:[39.34,39.43],39.44:[np.nan,44.40],
39.43:[43.48,np.nan],43.39:[39.44,39.34],43.48:[48.33,48.03],48.33:[33.28,np.nan],
48.43:[np.nan,43.39],48.03:[np.nan,np.nan],10.07:[7.05,7.09],3.48:[48.43,48.33]}

emptp = np.zeros(50)
tprobs = np.identity(50)[1:]
tmatrix50 = np.array([[emptp,emptp,emptp,tprobs[3],emptp,emptp],[emptp,tprobs[48],emptp,emptp,emptp,emptp],\
        [emptp,emptp,emptp,emptp,emptp,tprobs[47]],[tprobs[0],emptp,tprobs[4],emptp,tprobs[5],emptp],\
        [emptp,emptp,emptp,tprobs[6],emptp,tprobs[3]],[emptp,tprobs[3],emptp,tprobs[7],emptp,emptp],\
        [tprobs[4],emptp,tprobs[8],emptp,tprobs[9],emptp],[tprobs[5],emptp,tprobs[9],emptp,tprobs[10],emptp],\
        [emptp,emptp,emptp,tprobs[11],emptp,tprobs[6]],[emptp,tprobs[6],emptp,tprobs[12],emptp,tprobs[7]],\
        [emptp,tprobs[7],emptp,tprobs[13],emptp,emptp],[tprobs[8],emptp,tprobs[14],emptp,tprobs[15],emptp],\
        [tprobs[9],emptp,tprobs[15],emptp,tprobs[16],emptp],[tprobs[10],emptp,tprobs[16],emptp,tprobs[17],emptp],\
        [emptp,emptp,emptp,tprobs[18],emptp,tprobs[11]],\
        [emptp,tprobs[11],emptp,tprobs[19],emptp,tprobs[12]],[emptp,tprobs[12],emptp,tprobs[20],emptp,tprobs[13]],\
        [emptp,tprobs[13],emptp,tprobs[21],emptp,emptp],[tprobs[14],emptp,tprobs[22],emptp,tprobs[23],emptp],\
        [tprobs[15],emptp,tprobs[23],emptp,tprobs[24],emptp],[tprobs[16],emptp,tprobs[24],emptp,tprobs[25],emptp],\
        [tprobs[17],emptp,tprobs[25],emptp,tprobs[26],emptp],[emptp,emptp,emptp,tprobs[27],emptp,tprobs[18]],\
        [emptp,tprobs[18],emptp,tprobs[28],emptp,tprobs[19]],[emptp,tprobs[19],emptp,tprobs[29],emptp,tprobs[20]],\
        [emptp,tprobs[20],emptp,tprobs[30],emptp,tprobs[21]],[emptp,tprobs[21],emptp,tprobs[31],emptp,emptp],\
        [tprobs[22],emptp,tprobs[32],emptp,tprobs[33],emptp],[tprobs[23],emptp,tprobs[33],emptp,tprobs[34],emptp],\
        [tprobs[24],emptp,tprobs[34],emptp,tprobs[35],emptp],[tprobs[25],emptp,tprobs[35],emptp,tprobs[36],emptp],\
        [tprobs[26],emptp,tprobs[36],emptp,tprobs[37],emptp],[emptp,emptp,emptp,tprobs[47],emptp,tprobs[27]],\
        [emptp,tprobs[27],emptp,tprobs[38],emptp,tprobs[28]],[emptp,tprobs[28],emptp,tprobs[39],emptp,tprobs[29]],\
        [emptp,tprobs[29],emptp,tprobs[40],emptp,tprobs[30]],[emptp,tprobs[30],emptp,tprobs[41],emptp,tprobs[31]],\
        [emptp,tprobs[31],emptp,tprobs[48],emptp,emptp],[tprobs[33],emptp,tprobs[42],emptp,tprobs[43],emptp],\
        [tprobs[34],emptp,tprobs[43],emptp,tprobs[44],emptp],[tprobs[35],emptp,tprobs[44],emptp,tprobs[45],emptp],\
        [tprobs[36],emptp,tprobs[45],emptp,tprobs[46],emptp],[emptp,tprobs[47],emptp,emptp,emptp,tprobs[38]],\
        [emptp,tprobs[38],emptp,emptp,emptp,tprobs[39]],[emptp,tprobs[39],emptp,emptp,emptp,tprobs[40]],\
        [emptp,tprobs[40],emptp,emptp,emptp,tprobs[41]],[emptp,tprobs[41],emptp,emptp,emptp,tprobs[48]],\
        [tprobs[32],emptp,tprobs[2],emptp,tprobs[42],emptp],[tprobs[37],emptp,tprobs[46],emptp,tprobs[1],emptp]])

phexdf = pd.DataFrame(pairedHex_adjacentDict)
phexdf = phexdf.T
phexdf['state']=phexdf.index.astype("category")
phexdf['statecodes'] = phexdf.state.cat.codes
leftcodes = []
rightcodes = []
for i in phexdf.index:
    if np.isnan(phexdf.loc[i,0]):
        leftcodes.append(np.nan)
    else:
        leftcodes.append(phexdf.loc[phexdf.state==phexdf.loc[i,0],"statecodes"].values[0])
    if np.isnan(phexdf.loc[i,1]):
        rightcodes.append(np.nan)
    else:
        rightcodes.append(phexdf.loc[phexdf.state==phexdf.loc[i,1],"statecodes"].values[0])
phexdf['leftcodes'] = leftcodes
phexdf['rightcodes'] = rightcodes

splitstates=np.array([s.split('.') for s in phexdf.loc[:,'state'].astype(str).values])
phexdf["from_state"] = splitstates[:,0].astype("uint8")
phexdf["to_state"] = np.array([st[:2] for st in [s + ('0') for s in splitstates[:,1]]]).astype("uint8")

state_next_pairs = phexdf.loc[:,["statecodes","leftcodes","rightcodes"]].sort_values("statecodes").values
adj_matx = phexdf.loc[:,["statecodes","leftcodes","rightcodes"]].sort_values("statecodes").values[:,1:]
paired_tmat = np.zeros((len(phexdf),len(phexdf)))
paired_tmatrix = np.zeros((len(phexdf),2,len(phexdf)))
for i in state_next_pairs:
    if not np.isnan(i[1]):
        paired_tmat[int(i[0]),int(i[1])]=1
        paired_tmatrix[int(i[0]),0,int(i[1])]=1
    if not np.isnan(i[2]):
        paired_tmatrix[int(i[0]),1,int(i[2])]=1
        paired_tmat[int(i[0]),int(i[2])]=1

arrowdf = phexdf.loc[:,["statecodes","from_state","to_state"]]
for col in ["arrow_start_x","arrow_start_y","arrow_end_x","arrow_end_y"]:
    arrowdf.loc[:,col]=0
for i in arrowdf.index:
    arrowdf.loc[i,"arrow_start_x"]=centdf.loc[arrowdf.loc[i,"from_state"],0]
    arrowdf.loc[i,"arrow_start_y"]=centdf.loc[arrowdf.loc[i,"from_state"],1]
    arrowdf.loc[i,"arrow_end_x"]=centdf.loc[arrowdf.loc[i,"to_state"],0]
    arrowdf.loc[i,"arrow_end_y"]=centdf.loc[arrowdf.loc[i,"to_state"],1]

from matplotlib import cm,colors
valcols = cm.viridis
norm = colors.Normalize(vmin=0, vmax=1)

def plot_arrowMapFromStates(statelist,arrowvals):
    plot_arrowMap(arrowdf.loc[np.isin(arrowdf.statecodes,statelist)],arrowvals)

def plot_arrowMap(arrowdf,arrowvals):
    plt.figure(figsize=(10,10))
    plt.xlim(-1,15)
    plt.ylim(-1,22)
    ax = plt.gca()
    visited = []
    cind = 0
    for i in arrowdf.sort_values("from_state").index:
        state = arrowdf.loc[i,["arrow_start_x","arrow_start_y","arrow_end_x","arrow_end_y"]].values#/100
        jit = .2 if str(arrowdf.loc[i,["to_state","from_state"]].values) in visited else 0
        plt.arrow(x=state[0]+jit, y=state[1], dx=state[2]-state[0], dy=state[3]-state[1],fc=valcols(arrowvals[cind]),
        head_width=.6, head_length=.4,width = .25,alpha=0.4,shape='left',length_includes_head=True)
        visited.append(str(arrowdf.loc[i,["from_state","to_state"]].values))
        cind +=1

def plot_hex_layout(ax,bars,text=False,size='lg'):
    '''plot hex layout on the specified axis (ax)'''
    sz = 1000 if size=='lg' else 300
    bardf = centdf.drop(bars,axis=0)
    ax.scatter(bardf[0].values,bardf[1].values,marker='H',color='steelblue',alpha=0.9,s=sz)
    if text:
        for h in range(len(bardf)):
            ax.text(bardf.iloc[h][0]-.3,bardf.iloc[h][1]-.3,bardf.iloc[h].name,\
                    color='whitesmoke',fontweight='demibold')
    ax.yticks([])
    ax.xticks([])
    return mapfig

def reduce_mem_usage(df):
    start_mem = df.memory_usage().sum() / 1024**2
    print('Memory usage of dataframe is {:.2f} MB'.format(start_mem))
    
    for col in df.columns:
        col_type = df[col].dtype
        
        if col_type != object:
            c_min = df[col].min()
            c_max = df[col].max()
            if str(col_type)[:3] == 'int':
                if c_min > np.iinfo(np.int8).min and c_max < np.iinfo(np.int8).max:
                    df[col] = df[col].astype(np.int8)
                elif c_min > np.iinfo(np.int16).min and c_max < np.iinfo(np.int16).max:
                    df[col] = df[col].astype(np.int16)
                elif c_min > np.iinfo(np.int32).min and c_max < np.iinfo(np.int32).max:
                    df[col] = df[col].astype(np.int32)
                elif c_min > np.iinfo(np.int64).min and c_max < np.iinfo(np.int64).max:
                    df[col] = df[col].astype(np.int64)  
            else:
                if c_min > np.finfo(np.float16).min and c_max < np.finfo(np.float16).max:
                    df[col] = df[col].astype(np.float16)
                elif c_min > np.finfo(np.float32).min and c_max < np.finfo(np.float32).max:
                    df[col] = df[col].astype(np.float32)
                else:
                    df[col] = df[col].astype(np.float64)
        else:
            df[col] = df[col].astype('category')

    end_mem = df.memory_usage().sum() / 1024**2
    print('Memory usage after optimization is: {:.2f} MB'.format(end_mem))
    print('Decreased by {:.1f}%'.format(100 * (start_mem - end_mem) / start_mem))
    
    return df

def get_leave_inds(data):
    #get indices of last beam break before entering next
    reducedf = data.loc[:,['beamA','beamB','beamC']]
    bybeam = reducedf.fillna(0).values*np.array([1,2,3]).T
    bybeam = np.sum(bybeam,axis=1)
    bybeam = bybeam-1
    beambroken = reducedf.fillna(0).values*np.array([1,2,3]).T
    beambroken = np.sum(beambroken,axis=1)
    beambroken = beambroken-1
    beambroken = pd.DataFrame({"beam":beambroken})
    vinds = data.loc[data.port!=-100].index - data.index.min()
    if vinds[0]==0:
        justbeams = beambroken.loc[beambroken.beam!=-1,'beam']
        new_vinds = justbeams.loc[justbeams.diff()!=0].index
        beambroken['port']=-100
        beambroken.loc[new_vinds,'port']=data.loc[data.port!=-100,'port'].values
        data['port'] = beambroken.port.values
        rwds = data.loc[data.rwd!=-100,'rwd'].values
        data['rwd']=-100
        data.loc[new_vinds+data.index.min(),'rwd']=rwds
        vinds = new_vinds
    else:
        beambroken['port'] = data.port.values
    leaveinds = []
    for i in range(1,len(vinds)):
        p = beambroken.loc[vinds[i],'port']
        entries = beambroken.loc[vinds[i-1]:vinds[i]]
        inds = entries.loc[entries.beam!=-1].index.values
        sameport = True
        end = -1
        while sameport:
            try:
                lastbeam = beambroken.loc[inds[end],'beam']
            except:
                leaveinds.append(vinds[i-1])
                break
            if lastbeam != p:
                sameport = False
                leaveinds.append(inds[end])
            end-=1
    return np.add(leaveinds,data.index.min()),vinds[1:]+data.index.min(),data

def viz_value_map(valmap,i):
    plt.clf()
    plt.title('iteration '+str(i)+' value map',fontsize='x-large',fontweight='bold')
    plt.scatter(centdf[0].sort_index().values,centdf[1].sort_index().values,c=\
            valmap,marker='H',s=600,cmap='magma',edgecolors= "black",vmin=0,vmax=1)
    plt.colorbar()
    plt.pause(.1)

def viz_value_map_pairedState(valmap,i):
    plt.clf()
    plt.title('iteration '+str(i)+' value map',fontsize='x-large',fontweight='bold')
    plot_arrowMapFromStates(np.arange(0,126),arrowvals=valmap[i,:126]*100)
    plt.colorbar()
    plt.pause(.1)


class Photrat():

    fs = 250
    directory_prefix = "/Volumes/Tim K/Photometry Data/Triangle Maze/"
    pltvel = True
    plot_window = [-5,5]
    plot_heat = True

    def __init__(self,ID,date):
        self.ID = ID
        self.date = date
        self.loadpath = self.directory_prefix+ID+'/'+date+'/'

    def load_df(self):
        self.df = reduce_mem_usage(pd.read_csv(self.loadpath+self.ID\
            +'_'+self.date+'_h_sampleframe.csv'))

    def get_visinds(self):
        self.visinds = self.df.loc[self.df.port!=-100].index

    def set_plot_trace(self,trace):
        self.plot_trace = trace

    def set_pool_factor(self,factor):
        self.pool_factor = factor

    def set_plot_window(self,win):
        self.plot_window = win

    def show_vel_trace(self,velbool):
        self.pltvel = velbool

    def show_tri_heatmap(self,heatbool):
        self.plot_heat = heatbool

    def add_bars(self,bars):
        self.bars = bars

    def create_hexdf(self,s,b):
        self.hexdata = self.df.loc[(self.df.session==s)&(self.df.block==b)]

    def load_tmat(self):
        try:
            self.tmat = np.load(self.directory_prefix+self.hexdata.rat.unique()[0]+'/'+\
            str(self.hexdata.date.unique()[0])+'/'+'tmat_block_'+str(b)+'.0.npy')
        except:
            print('no tmat saved for '+str(self.hexdata.rat.unique()[0])+\
                self.hexdata.date.unique()[0])

    def get_bars(self):
        self.bars = np.where(np.mean(self.tmat,axis=1)==0)[0][1:]-1

    def add_bars2tmatrix(self):
        self.bar_tmatrix = tmatrix.copy()
        self.get_bars()
        self.bar_tmatrix[self.bars] = np.tile(emptp,(6,1))
        availstates = np.setdiff1d(np.arange(0,49),bars)
        self.availstates = np.setdiff1d(availstates,[0,1,2])
    
    def add_tot_triCount(self):
        self.df['tot_tri']=-100
        self.df.loc[self.visinds,'tot_tri']=np.arange(1,len(self.visinds)+1)
        self.df.tot_tri = self.df.tot_tri.replace(-100,method='ffill')
        self.df.tot_tri = self.df.tot_tri.replace(-100,method='bfill')

    def get_ratHemColumn(self):
        self.df['rat_fiber'] = self.df.rat.astype(str) + '_' + self.df.fiberloc.astype(str)
        self.df.loc[:,"rat_fiber"] = self.df.rat_fiber.astype('category')
    
    def get_simpleStateTransitions(self,b):
        if self.df.loc[:,'session_type'].values[0]=='prob':
            self.bars = self.barIDs[0]
        else:
            self.tmat = self.tmats[int(b-1)]
            self.bars = self.barIDs[int(b-1)]
        self.add_bars2tmatrix50()
        self.to_state = [[0,0,0,0,0,0]]+[np.argmax(a,axis=1) \
            for a in self.bar_tmatrix]

    def add_bars2tmatrix50(self):
        self.bar_tmatrix = tmatrix50.copy()
        self.bar_tmatrix[self.bars] = np.tile(emptp,(6,1))
        
    def fill_gap(self,h1,h3):
        gap = self.sharedhex(h1,h3)
        if gap!=-1:
            return [gap]
        else:
            nextgap = -1
            candidates = [h for h in self.to_state[int(h1)] if h !=0]
            for i in range(len(candidates)):
                candidateNext = candidates[i]
                nextgap = self.sharedhex(candidateNext,h3)
                if nextgap!=-1:
                    return [candidateNext,nextgap]
                else:
                    nextgapPrime = -1
                    candidatesPrime = [h for h in self.to_state[int(candidateNext)] if h !=0]
                    for i in range(len(candidatesPrime)):
                        candidateNextPrime = candidatesPrime[i]
                        nextgapPrime = self.sharedhex(candidateNextPrime,h3)
                        if nextgapPrime!=-1:
                            return [candidateNext,candidateNextPrime,nextgapPrime]
        return -1
    
    def sharedhex(self,x,y):
        comlist =  [v for v in self.to_state[int(x)] if v in self.to_state[int(y)]\
             and v!=0 and v not in self.bars]
        if len(comlist)>0:
            return comlist[np.random.randint(len(comlist))]
        else:
            return -1
    
    def fill_gaps(self,b): 
        #data = self.df.loc[self.df.hexlabels!=-1,["hexlabels","x","y"]].copy()
        data = self.df.loc[(self.df.hexlabels!=-1)&\
        (self.df.block==b),["hexlabels","x","y"]].copy()
        data = data.loc[data.hexlabels.diff()!=0,:]
        originalGapInds = data.index.values#.loc[data.hexlabels.diff()!=0,:].index
        data.reset_index(inplace=True)
        gapInds = data.index.values#.loc[data.hexlabels.diff()!=0,:].index.values
        extraInds = []
        extraHexes = []
        recordedGapInds = []
        for i,j in zip(gapInds,originalGapInds):
            if i==len(data)-1:
                continue
            prevHex = data.loc[i,"hexlabels"]
            nextHex = data.loc[i+1,"hexlabels"]
            if prevHex in self.bars: #+1
                self.df.loc[j,["hexlabels"]] = -1
                continue
            if prevHex==nextHex or prevHex in self.to_state[int(nextHex)]:
                continue
            prevInd = data.loc[i,"index"]
            nextInd = data.loc[i+1,"index"]
            filledHexes = self.fill_gap(prevHex,nextHex)
            if filledHexes == -1:
                continue
            if len(filledHexes)==1:
                extraInds.append(j+int((nextInd-prevInd)/2))
                extraHexes.append(filledHexes[0])
            elif len(filledHexes)>1:
                for g in range(len(filledHexes)):
                    #find out how much time is between detected hexes
                    extraInds.append(j+int((1+g)*(nextInd-prevInd)/(1+len(filledHexes))))
                    extraHexes.append(filledHexes[g])
        self.df.loc[extraInds,["hexlabels"]]=extraHexes
    
    
    def remove_flickers(self):
        data = self.df.loc[(self.df.hexlabels!=-1),["hexlabels","x","y"]].copy()
        data = data.loc[data.hexlabels.diff()!=0,:]
        hexlist = data.loc[:,"hexlabels"].values
        hexInds = data.index.values
        for i in range(len(hexlist)-4):
            seq = hexlist[i:i+4]
            if seq[0]==seq[2] and seq[1]==seq[3]:
                hexlist[i+1] = seq[0]
        #return hexInds,hexlist
        self.df.loc[hexInds,["hexlabels"]]=hexlist

    def get_DA_approachAUC(self):
        daAUCs = []
        for i in range(0,len(self.visinds)-1):
            dat = self.df.loc[self.visinds[i]:self.visinds[i+1]]
            daAUCs.append(dat.loc[dat.lastleave==2,'green_z_scored'].mean())
        dat = self.df.loc[self.visinds[-1]:]
        daAUCs.append(dat.loc[dat.lastleave==2,'green_z_scored'].mean())
        self.triframe.loc[:,'approachDA_auc'] = daAUCs

    def get_rampSlope(self):
        daRamps = [np.nan]
        for i in range(0,len(self.visinds)-1):
            dat = self.df.loc[self.visinds[i]:self.visinds[i+1]]
            X = (dat.index.values-dat.index.min())/250
            y = dat.green_z_scored.values
            modo = LR(fit_intercept=True,normalize=False).fit(X.reshape(-1,1),y)
            daRamps.append(modo.coef_[0])
        self.triframe.loc[:,'daRampSlope'] = daRamps
    
        
    def get_port_occupancy(self):
        newpathhexes = [47,38,43,33,6,5,32,37,42,39,34,28,8,10,7,27,31,41,44,29,23,11,13,9]
        innewpath = self.df.loc[self.df.hexlabels.isin(newpathhexes).\
        astype(int).diff()>0].index.values
        self.df['atport']=0
        for i in self.df.loc[self.df.port!=-100].index:
            lats = np.subtract(innewpath,i)
            try:
                latency = lats[lats>0][0]
                self.df.loc[i:i+latency,'atport']=1
            except:
                print('rat did not begin new trial')

    def create_triframe(self):
        self.triframe = self.df.loc[self.visinds]
        self.triframe['lrchoice'] = make_lr(self.triframe.port.values)
        self.triframe['lrchoice'] = self.triframe.lrchoice.astype("int8")
        self.triframe = self.triframe.reset_index()
        self.triframe.drop(['x','y','vel','acc'],axis=1,inplace=True)

    def save_triframe(self):
        self.triframe.to_csv(self.loadpath+"triframe.csv",columns = \
            ['port', 'rwd', 'block','tri','lrchoice',\
                'session_type','rat','date','nextp','lenAC','lenBC','lenAB',\
                'fiberloc','dtop','DA_RPE'])


