import matplotlib.pylab as plt
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi']=300
from PyComplexHeatmap import *

df = pd.read_csv("df.csv", index_col="mamal")
df_rows = pd.read_csv("df_rows.csv",index_col="mamal")
df_cols = pd.read_csv("df_cols.csv",index_col="col")
print(df)
print(df_rows)
print(df_cols)

col_colors_dict = {'Tissue': {'#A40043': 'Blood',
  '#00E5FF': 'Brain',
  '#00BECC': 'Cerebellum',
  '#B2F8FF': 'Striatum',
  '#6CBF00': 'Liver',
  '#FFCCEE': 'Muscle',
  '#1E9351': 'Skin'}}

#Put annotations on the top
col_ha = HeatmapAnnotation(label=anno_label(df_cols.Family, merge=True, rotation=45),
                               Family=anno_simple(df_cols.Family, legend=True),
                               Tissue=df_cols.Tissue,label_side='right', axis=1)
plt.figure(figsize=(7, 4.5))
cm = ClusterMapPlotter(data=df, top_annotation=col_ha,
                       show_rownames=True, show_colnames=False,row_names_side='left',
                       col_split=df_cols.Family, cmap='exp1', label='AUC',
                       rasterized=True, legend=True,legend_anchor='ax_heatmap',legend_width=50)
#legend_pad control the space between heatmap and legend.
#plt.savefig("clustermap.pdf", bbox_inches='tight')
plt.show()
#Put annotations on the bottom
col_ha = HeatmapAnnotation(Tissue=anno_simple(df_cols.Tissue,height=5),
                           Family=anno_simple(df_cols.Family, legend=False,height=6),
                           label=anno_label(df_cols.Family, merge=True,rotation=-45),
                           label_side='right',axis=1)
plt.figure(figsize=(6.5, 4))
cm = ClusterMapPlotter(data=df, bottom_annotation=col_ha,
                       show_rownames=True, show_colnames=False,row_names_side='right',
                       col_split=df_cols.Family, cmap='jet', label='AUC',
                       rasterized=True, legend=True)
plt.show()
#Put annotations on the left
row_ha = HeatmapAnnotation(label=anno_label(df_cols.Family, merge=True,rotation=45),
                           Family=anno_simple(df_cols.Family, legend=True,height=5),
                           Tissue=anno_simple(df_cols.Tissue,height=5),
                           label_side='top',
                           label_kws={'rotation':45,'rotation_mode':'anchor','color':'red'},
                           axis=0)
plt.figure(figsize=(4, 6))
cm = ClusterMapPlotter(data=df.T,left_annotation=row_ha,
                       show_rownames=False, show_colnames=True,col_names_side='top',
                       row_split=df_cols.Family, row_split_gap=0,
                       cmap='exp1', label='AUC',
                       rasterized=True, legend=True,
                       xticklabels_kws={'labelrotation':45,'labelcolor':'blue'})
plt.show()

#Put annotation on the right
row_ha = HeatmapAnnotation(Tissue=df_cols.Tissue,
                           Family=anno_simple(df_cols.Family, legend=False,height=5),
                           label=anno_label(df_cols.Family, merge=True,rotation=45),
                           label_side='bottom',
                           label_kws={'rotation':-45,'color':'red'},
                           axis=0)
plt.figure(figsize=(4, 6))
cm = ClusterMapPlotter(data=df.T,right_annotation=row_ha,
                       show_rownames=False, show_colnames=True,col_names_side='bottom',
                       row_split=df_cols.Family, cmap='jet', label='AUC',
                       rasterized=True, legend=True,row_split_gap=0.1,
                       xticklabels_kws={'labelrotation':-45,'labelcolor':'blue'})
#plt.savefig("annotation.pdf", bbox_inches='tight')
plt.show()

#https://dingwb.github.io/PyComplexHeatmap/build/html/notebooks/clustermap.html
