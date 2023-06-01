import os,sys
import pandas as pd
import matplotlib.pylab as plt
import pickle
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi']=300
import PyComplexHeatmap
from PyComplexHeatmap import *
import seaborn as sns

# Load the brain networks dataset, select subset, and collapse the multi-index
df = pd.read_csv("brain_networks.csv", header=[0, 1, 2], index_col=0)
used_networks = [1, 5, 6, 7, 8, 12, 13, 17]
used_columns = (df.columns
                  .get_level_values("network")
                  .astype(int)
                  .isin(used_networks))
df = df.loc[:, used_columns]

df.columns = df.columns.map("-".join)

# Compute a correlation matrix and convert to long-form
corr_mat = df.corr().stack().reset_index(name="correlation")
corr_mat['Level']=corr_mat.correlation.apply(lambda x:'High' if x>=0.7 else 'Middle' if x >= 0.3 else 'Low')
data=corr_mat.pivot(index='level_0',columns='level_1',values='correlation')
print(data.head())
print(corr_mat.Level.value_counts().index.tolist())
print(corr_mat.head())
#Plot traditional heatmap using square marker
plt.figure(figsize=(8,8))
cm = DotClustermapPlotter(data=corr_mat,x='level_0',y='level_1',value='correlation',
              c='correlation',cmap='jet',vmax=1,vmin=0,s=0.7,marker='s',spines=True,alpha=0.7)
cm.ax_heatmap.grid(which='minor',color='white',linestyle='--',alpha=0.6,linewidth=1)
plt.show()
#Simple dot heatmap using fixed dot size
plt.figure(figsize=(8,8))
cm = DotClustermapPlotter(corr_mat,x='level_0',y='level_1',value='correlation',
              c='correlation',cmap='Reds',vmax=1,vmin=0,s=0.5)
plt.show()
#Changing the size of point
plt.figure(figsize=(8,8))
cm = DotClustermapPlotter(corr_mat,x='level_0',y='level_1',value='correlation',
              c='correlation',s='correlation',cmap='Reds',vmax=1,vmin=0)
cm.ax_heatmap.grid(which='minor',color='gray',linestyle='--',alpha=0.4)
plt.show()
#Add parameter hue and use different colors for different groups
plt.figure(figsize=(8,8))
cm = DotClustermapPlotter(corr_mat,x='level_0',y='level_1',value='correlation',hue='Level',
              colors={'High':'red','Middle':'purple','Low':'green'},s='correlation',vmax=1,vmin=0)
plt.show()
#Add parameter hue and use different cmap and marker for different groups
plt.figure(figsize=(8,8))
cm = DotClustermapPlotter(corr_mat,x='level_0',y='level_1',value='correlation',hue='Level',
              cmap={'High':'Reds','Middle':'Purples','Low':'Greens'},
              colors={'High':'red','Middle':'purple','Low':'green'},
              marker={'High':'P','Middle':'*','Low':'D'},spines=True,
              vmax=1,vmin=0,alpha=0.9)
plt.show()
#Plot clustermap using seaborn brain networks dataset
print(corr_mat.head())
df_row=corr_mat['level_0'].drop_duplicates().to_frame()
df_row['RowGroup']=df_row.level_0.apply(lambda x:x.split('-')[0])
df_row.set_index('level_0',inplace=True)

df_col=corr_mat['level_1'].drop_duplicates().to_frame()
df_col['ColGroup']=df_col.level_1.apply(lambda x:x.split('-')[0])
df_col.set_index('level_1',inplace=True)

print(df_row.head())
print(df_col.head())
row_ha = HeatmapAnnotation(Row=anno_simple(df_row.RowGroup,cmap='Set1',
                                           add_text=True,text_kws={'color':'black','rotation':-90},
                                          legend=False),
                           axis=0,verbose=0,label_kws={'rotation':45,'horizontalalignment':'left'})

col_ha = HeatmapAnnotation(label=anno_label(df_col.ColGroup, merge=True,rotation=45),
                           Col=anno_simple(df_col.ColGroup,cmap='Dark2',legend=False,add_text=True),
                           verbose=0,label_side='left',label_kws={'horizontalalignment':'right'})

plt.figure(figsize=(9, 8))
cm = DotClustermapPlotter(data=corr_mat, x='level_0',y='level_1',value='correlation',
                          hue='Level', cmap={'High':'Reds','Middle':'Purples','Low':'Greens'},
                          colors={'High':'red','Middle':'purple','Low':'green'},
                          marker={'High':'P','Middle':'*','Low':'D'},
                          top_annotation=col_ha,right_annotation=row_ha,
                          col_split=2,row_split=2, col_split_gap=0.5,row_split_gap=1,
                          show_rownames=True,show_colnames=True,row_dendrogram=True,
                          tree_kws={'row_cmap': 'Set1'},verbose=0,legend_gap=7,spines=True,)
plt.show()
#Visualize up to five dimension data using DotClustermapPlotter
data=pd.read_csv("kycg_result.txt",sep='\t')
data=data.loc[data.Category.isin(['rmsk1','ChromHMM','EnsRegBuild'])]
data.SampleID.replace({'Clark2018_Argelaguet2019':'Dataset1','Luo2022':'Dataset2'},inplace=True)
max_p=np.nanmax(data['-log10(Pval)'].values)
data['-log10(Pval)'].fillna(max_p,inplace=True)
data['ID']=data.SampleID + '-' + data.CpGType
vc=data.groupby('Term').SampleID.apply(lambda x:x.nunique())

data=data.loc[data.Term.isin(vc[vc>=2].index.tolist())]
# p_max=data['-log10(Pval)'].max()
# p_min=data['-log10(Pval)'].min()
# data['-log10(Pval)']=data['-log10(Pval)'].apply(lambda x:(x-p_min)/(p_max-p_min))
df_col=data.ID.drop_duplicates().to_frame()
df_col['Dataset']=df_col.ID.apply(lambda x:x.split('-')[0])
df_col['Correlation']=df_col.ID.apply(lambda x:x.split('-')[1])
df_col.set_index('ID',inplace=True)
df_row=data.loc[:,['Term','Category']].drop_duplicates()
df_row.set_index('Term',inplace=True)
print(data.head())
print(data['-log10(Pval)'].describe())
print(data.CpGType.unique())
print(data.EnrichType.unique())
print(df_col)
print(df_row)
row_ha = HeatmapAnnotation(
                           Category=anno_simple(df_row.Category,cmap='Set1',
                                           add_text=False,legend=False),
                           label=anno_label(df_row.Category, merge=True,rotation=0),
                           axis=0,verbose=0,label_kws={'rotation':45,'horizontalalignment':'left'})

col_ha = HeatmapAnnotation(
                           Dataset=anno_simple(df_col.Dataset,cmap='Set1',legend=False,add_text=True),
                           Correlation=anno_simple(df_col.Correlation,cmap='Dark2',legend=False,add_text=True),
                           verbose=0,label_side='left',label_kws={'horizontalalignment':'right'})

plt.figure(figsize=(3.5, 5))
cm = DotClustermapPlotter(data=data, x='ID',y='Term',value='-log10(Pval)',c='-log10(Pval)',s='odds_ratio',
                          hue='EnrichType', row_cluster=False,col_cluster=False,
                          cmap={'Enrich':'RdYlGn_r','Depletion':'coolwarm_r'},
                          colors={'Enrich':'red','Depletion':'blue'},
                          #marker={'Enrich':'^','Depletion':'v'},
                          top_annotation=col_ha,right_annotation=row_ha,
                          col_split=df_col.Dataset,row_split=df_row.Category, col_split_gap=0.5,row_split_gap=1,
                          show_rownames=True,show_colnames=False,row_dendrogram=False,
                          verbose=0,legend_gap=7,alpha=0.8)
plt.savefig("dotHeatmap1.pdf",bbox_inches='tight')
plt.show()
plt.figure(figsize=(3.5, 5))
cm = DotClustermapPlotter(data=data, x='ID',y='Term',value='odds_ratio',c='-log10(Pval)',s='-log10(Pval)',
                          hue='EnrichType', row_cluster=False,cmap='jet',
                          colors={'Enrich':'red','Depletion':'blue'},
                          marker={'Enrich':'P','Depletion':'D'},value_na=25,c_na=25,
                          top_annotation=col_ha,right_annotation=row_ha,
                          col_split=df_col.Dataset,row_split=df_row.Category, col_split_gap=0.5,row_split_gap=1,
                          show_rownames=True,verbose=0,legend_gap=7,alpha=0.7)
# plt.savefig(os.path.expanduser("~/Gallery/20230227_kycg.pdf"),bbox_inches='tight')
plt.show()
print(data['-log10(Pval)'].describe())
plt.figure(figsize=(3.5, 5))
cm = DotClustermapPlotter(data=data, x='ID',y='Term',value='-log10(Pval)',c='-log10(Pval)',s='odds_ratio',
                          hue='EnrichType', row_cluster=False,col_cluster=False,cmap='jet',
                          colors={'Enrich':'red','Depletion':'blue'},
                          marker={'Enrich':'P','Depletion':'_'},value_na=25,c_na=25,
                          top_annotation=col_ha,right_annotation=row_ha,
                          col_split=df_col.Dataset,row_split=df_row.Category, col_split_gap=0.5,row_split_gap=1,
                          show_rownames=True,verbose=0,legend_gap=7,spines=True,)
plt.show()


#https://dingwb.github.io/PyComplexHeatmap/build/html/notebooks/dotHeatmap.html
