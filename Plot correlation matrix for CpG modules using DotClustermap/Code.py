import matplotlib.pylab as plt
from PyComplexHeatmap import *
from matplotlib.colors import LinearSegmentedColormap

df_corr = pd.read_csv("kycg_modules_correlations.csv",sep='\t',index_col=0)
df_ann = pd.read_csv("kycg_modules_annotations.csv",sep='\t')
betas = pd.read_csv("kycg_modules_betas.csv",sep='\t')
cpg_std=betas.std().to_dict()
cpg_mean=betas.mean().to_dict()
df_ann.set_index('CpG',inplace=True)
df_ann.Module=df_ann.Module.astype(str)

print(df_ann.Module.value_counts().head(10))

print(df_ann.loc[df_ann.Module.isin(['4','1','3','9','2','39','37'])].HM.unique())

print(df_corr.head())

print(df_ann.head())

df_ann=df_ann.loc[df_ann.Module.isin(['4','1','3','9','2','39'])]
keep_cpgs=df_ann.index.tolist()
df_corr=df_corr.loc[keep_cpgs,keep_cpgs]
df_ann['Std']=df_ann.index.to_series().map(cpg_std)
df_ann['Mean']=df_ann.index.to_series().map(cpg_mean)
data=df_corr.stack().reset_index()
data.columns=['X','Y','Correlation']
data['Module']=data.X.map(df_ann.Module.to_dict())
data['ChromHMM']=data.X.map(df_ann.ChromHMM.to_dict())
keep_hm=['H3K4me1','H3K4me3','H3K27me1','H3K27me3','H3K27me3B']
for hm in keep_hm:
    df_ann[hm]=df_ann.HM.fillna('').apply(lambda x:1 if hm in x.split(';') else 0)
    data[hm]=data.X.map(df_ann[hm].to_dict())

print(df_ann.shape)
print(df_ann.head())

print(data.shape)
print(data.head())

#Plotting the Dot clustermap

row_ha = HeatmapAnnotation(Module=anno_simple(df_ann.Module,cmap='Dark2',legend=False,height=5,
                                              add_text=True,text_kws={'color':'black','fontsize':12}),
                           axis=0,verbose=0,label_kws={'visible':False})

all_cmaps=matplotlib.pyplot.colormaps()
if 'binarize' not in all_cmaps:
    c = LinearSegmentedColormap.from_list('binarize', [(0, 'lightgray'), (1, 'black')])
    plt.register_cmap(cmap=c)

col_ha = HeatmapAnnotation(#label=anno_label(df_col.ColGroup, merge=True,rotation=45),
                           Module=anno_simple(df_ann.Module,cmap='Dark2',legend=False,height=5,
                                              add_text=True,text_kws={'color':'black','fontsize':12}),
                           ChromHMM=anno_simple(df_ann.ChromHMM,cmap='tab20'),
                           Mean=anno_barplot(df_ann.Mean,cmap='jet',linewidth=0.1),
                           # H3K4me1=anno_simple(df_ann.H3K4me1,cmap='binarize',legend=False),
                           # H3K4me3=anno_simple(df_ann.H3K4me3,cmap='binarize',legend=False),
                           # H3K27me1=anno_simple(df_ann.H3K27me1,cmap='binarize',legend=False),
                           # H3K27me3=anno_simple(df_ann.H3K27me3,cmap='binarize',legend=False),
                           # H3K27me3B=anno_simple(df_ann.H3K27me3B,cmap='binarize',legend=False),
                           verbose=0,label_side='right',label_kws={'horizontalalignment':'left'})

plt.figure(figsize=(10, 9))
cm = DotClustermapPlotter(data=data, x='X',y='Y',value='Correlation',c='Correlation',s='Correlation',
                          hue='Module', cmap='jet',#cmap={'High':'Reds','Middle':'Purples','Low':'Greens'},
                          #colors={'High':'red','Middle':'purple','Low':'green'},
                          #marker={'4':'P','1':'*','3':'D'},
                          top_annotation=col_ha,right_annotation=row_ha,
                          col_split=df_ann.Module,row_split=df_ann.Module, col_split_gap=1,row_split_gap=1,
                          row_dendrogram=True,legend_anchor="ax_heatmap",legend_hpad=7,legend_vpad=5,
                          tree_kws={'row_cmap':'Dark2'},verbose=0,legend_gap=7,alpha=2,spines=False)
# plot custom spines
for i in range(cm.heatmap_axes.shape[0]):
    for j in range(cm.heatmap_axes.shape[1]):
        if i != j:
            continue
        ax = cm.heatmap_axes[i][j]
        for side in ["top", "right", "left", "bottom"]:
            ax.spines[side].set_visible(True)
            ax.spines[side].set_color('red')
            ax.spines[side].set_linewidth(2)
plt.savefig("dotClustermap.pdf", bbox_inches='tight')
plt.show()

print(data.X.nunique())

#A smaller dot clustermap
df_ann=df_ann.loc[df_ann.Module.isin(['4','1','3'])]
keep_cpgs=df_ann.index.tolist()
df_corr=df_corr.loc[keep_cpgs,keep_cpgs]
data=df_corr.stack().reset_index()
data.columns=['X','Y','Correlation']
data['Module']=data.X.map(df_ann.Module.to_dict())
data['ChromHMM']=data.X.map(df_ann.ChromHMM.to_dict())
keep_hm=['H3K4me1','H3K4me3','H3K27me1','H3K27me3B']
for hm in keep_hm:
    df_ann[hm]=df_ann.HM.fillna('').apply(lambda x:1 if hm in x.split(';') else 0)
    data[hm]=data.X.map(df_ann[hm].to_dict())


row_ha = HeatmapAnnotation(Module=anno_simple(df_ann.Module,cmap='Dark2',legend=False,height=5,
                                              add_text=True,text_kws={'color':'black','fontsize':12}),
                           axis=0,verbose=0,label_kws={'visible':False})

all_cmaps=matplotlib.pyplot.colormaps()
if 'binarize' not in all_cmaps:
    c = LinearSegmentedColormap.from_list('binarize', [(0, 'lightgray'), (1, 'black')])
    plt.register_cmap(cmap=c)

col_ha = HeatmapAnnotation(#label=anno_label(df_col.ColGroup, merge=True,rotation=45),
                           Module=anno_simple(df_ann.Module,cmap='Dark2',legend=False,height=5,
                                              add_text=True,text_kws={'color':'black','fontsize':12}),
                           ChromHMM=anno_simple(df_ann.ChromHMM,cmap='tab20'),
                           Mean=anno_barplot(df_ann.Mean,cmap='jet',linewidth=0.1),
                           # H3K27me3=anno_simple(df_ann.H3K27me3,cmap='binarize',legend=False),
                           # H3K27me3B=anno_simple(df_ann.H3K27me3B,cmap='binarize',legend=False),
                           verbose=0,label_side='right',label_kws={'horizontalalignment':'left'})

plt.figure(figsize=(10, 9))
cm = DotClustermapPlotter(data=data, x='X',y='Y',value='Correlation',c='Correlation',s='Correlation',
                          hue='Module', cmap='jet',#cmap={'High':'Reds','Middle':'Purples','Low':'Greens'},
                          #colors={'High':'red','Middle':'purple','Low':'green'},
                          #marker={'4':'P','1':'*','3':'D'},
                          top_annotation=col_ha,right_annotation=row_ha,
                          col_split=df_ann.Module,row_split=df_ann.Module, col_split_gap=1,row_split_gap=1,
                          row_dendrogram=True,legend_anchor="ax_heatmap",legend_hpad=7,legend_vpad=5,
                          tree_kws={'row_cmap':'Dark2'},verbose=0,legend_gap=7,alpha=2,spines=True)
plt.savefig("dotClustermap2.pdf", bbox_inches='tight')
plt.show()
#https://dingwb.github.io/PyComplexHeatmap/build/html/notebooks/cpg_modules.html
