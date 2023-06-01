import matplotlib.pylab as plt
from matplotlib.colors import LinearSegmentedColormap
plt.rcParams['figure.dpi'] = 80
plt.rcParams['savefig.dpi']=300
import PyComplexHeatmap
from PyComplexHeatmap import *
import pickle

f=open("influence_of_snp_on_beta.pickle",'rb')
data=pickle.load(f)
f.close()
beta,snp,df_row,df_col,col_colors_dict,row_colors_dict=data
# beta is DNA methylation beta values matrix, df_row and df_col are row and columns annotation respectively, col_colors_dict and row_colors_dict are color for annotation
print(beta.iloc[:,list(range(5))].head(5))
print(df_row.head(5))
print(df_col.head(5))
beta=beta.sample(2000)
snp=snp.loc[beta.index.tolist()]
df_row=df_row.loc[beta.index.tolist()]

print(row_colors_dict)

row_ha = HeatmapAnnotation(Target=anno_simple(df_row.Target,colors=row_colors_dict['Target'],rasterized=True),
                           Group=anno_simple(df_row.Group,colors=row_colors_dict['Group'],rasterized=True),
                           axis=0)
col_ha= HeatmapAnnotation(label=anno_label(df_col.Strain,merge=True,rotation=15),
                          Strain=anno_simple(df_col.Strain,add_text=True),
                          Tissue=df_col.Tissue,Sex=df_col.Sex,
                          axis=1)
plt.figure(figsize=(5, 8))
cm = ClusterMapPlotter(data=beta, top_annotation=col_ha, left_annotation=row_ha,
                     show_rownames=False,show_colnames=False,
                     row_dendrogram=False,col_dendrogram=False,
                     row_split=df_row.loc[:, ['Target', 'Group']],
                     col_split=df_col['Strain'],cmap='parula',
                     rasterized=True,row_split_gap=1,legend=True,legend_anchor='ax_heatmap',legend_vpad=5)
cm.ax.set_title("Beta",y=1.03,fontdict={'fontweight':'bold'})
#plt.savefig("clustermap.pdf", bbox_inches='tight')
plt.show()


row_ha = HeatmapAnnotation(Target=anno_simple(df_row.Target, colors=row_colors_dict['Target'], rasterized=True),
                                   Group=anno_simple(df_row.Group, colors=row_colors_dict['Group'], rasterized=True),
                                   axis=0)

col_ha1 = HeatmapAnnotation(#label=anno_label(df_col.Strain, merge=True, rotation=15),,
                           Strain=anno_simple(df_col.Strain, add_text=True,
                                              text_kws={'fontweight':'bold'}),
                           Tissue=df_col.Tissue, Sex=df_col.Sex,
                           axis=1,verbose=0)  # df=df_col.loc[:,['Strain','Tissue','Sex']],

cm1 = ClusterMapPlotter(data=beta, top_annotation=col_ha1, left_annotation=None,
                       show_rownames=False, show_colnames=False,
                       row_dendrogram=False, col_dendrogram=False,
                       row_split=df_row.loc[:, ['Target', 'Group']],
                       col_split=df_col['Strain'], cmap='parula', #turbo, parula, viridis,
                       rasterized=True, row_split_gap=0.1,vmax=1,vmin=0,center=0.5,
                        plot=False,label='beta')

col_ha2 = HeatmapAnnotation(#label=anno_label(df_col.Strain, merge=True, rotation=15),,
                           Strain=anno_simple(df_col.Strain, add_text=True,
                                             text_kws={'fontweight':'bold'}),
                           Tissue=df_col.Tissue, Sex=df_col.Sex,
                           label_kws={'visible':False},axis=1,verbose=0)

my_cmap = LinearSegmentedColormap.from_list('my_cmap', [(0, 'lightgray'), (1, 'black')])
cm2 = ClusterMapPlotter(data=snp, top_annotation=col_ha2, left_annotation=row_ha,
                        show_rownames=False, show_colnames=False,
                        row_dendrogram=False, col_dendrogram=False,
                        col_cluster_method='ward',row_cluster_method='ward',
                        col_cluster_metric='jaccard',row_cluster_metric='jaccard',
                        row_split=df_row.loc[:, ['Target', 'Group']],
                        col_split=df_col['Strain'],
                        rasterized=True, row_split_gap=0.1,
                        plot=False,cmap=my_cmap,label='SNP') # or cmap='gray' or Greys,
cmlist=[cm2,cm1]

plt.figure(figsize=(6,9))
ax,legend_axes=composite(cmlist=cmlist, main=1,legend_hpad=2,col_gap=0.1)
cm1.ax.set_title("Beta")
cm2.ax.set_title("SNP")
ax.set_title("Influence of SNP on DNAm reading",y=1.03,fontdict={'fontweight':'bold','color':'red'})
# plt.savefig("beta_snp.pdf", bbox_inches='tight')
plt.show()


#https://dingwb.github.io/PyComplexHeatmap/build/html/notebooks/composite_heatmaps.html


