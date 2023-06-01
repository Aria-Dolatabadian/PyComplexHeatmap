import matplotlib.pylab as plt
import pickle
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi']=300
import PyComplexHeatmap
from PyComplexHeatmap import *
data=pd.read_csv("Loyfer2023.meth.csv",sep='\t',index_col=0)
df_row=pd.read_csv("Loyfer2023.meth.rows.csv",sep='\t',index_col=0)
df_col=pd.read_csv("Loyfer2023.meth.cols.csv",sep='\t',index_col=0)
print(data.head())
print(df_row.head())
print(df_col.head())
print(df_col.Group.unique())

col_colors_dict={
                'Adipocytes':'#1E93AE','Bladder-Ep':'#FF9C00','Blood-B':'#A40043','Blood-Granul':'#FF9F7F',
                'Blood-Mono+Macro':'#FF7F00','Blood-NK':'#FF2E8D','Blood-T':'#CC0043','Bone-Osteob':'#E5E5E5',
                'Breast-Basal-Ep':'#CC4407','Breast-Luminal-Ep':'#CC843D','Colon-Ep':'#663D28','Dermal-Fibro':'#1E937C',
                'Epid-Kerat':'#1E93AE','Eryth-prog':'#40705F','Fallopian-Ep':'#009351','Gallbladder':'#E7E4BF',
                'Gastric-Ep':'#CCA300','Head-Neck-Ep':'#002929','Heart-Cardio':'#FF99AA','Heart-Fibro':'#FF99FF',
                'Kidney-Ep':'#F6FF99','Liver-Hep':'#6CBF00','Lung-Ep-Alveo':'#BA99FF','Lung-Ep-Bron':'#CCCCFF',
                'Neuron':'#9e542e','Oligodend':'#2ca02c','Pancreas-Acinar':'#DF7F00','Pancreas-Beta':'#FFD866',
                'Pancreas-Delta':'#FFCC32','Pancreas-Duct':'#7F4C33','Skeletal-Musc':'#FFCCEE','Small-Int-Ep':'#CC9951',
                'Smooth-Musc':'#1E93AE','Thyroid-Ep':'#B2BFFF'}

col_ha = HeatmapAnnotation(label=anno_label(df_col['Group'],merge=True,rotation=90,extend=True,
                                            colors=col_colors_dict,adjust_color=True,luminance=0.75,
                                            relpos=(0.5,0)), #fontsize=10
                           Group=anno_simple(df_col['Group'],colors=col_colors_dict), #legend_kws={'fontsize':4}
                           verbose=0,axis=1)
row_ha = HeatmapAnnotation(
                           Group=anno_simple(df_row['Group'],legend=True,
                                             colors=col_ha.annotations[1].color_dict),
                           verbose=0,axis=0,plot_legend=False) #label_kws={'rotation':90,'rotation_mode':'anchor','color':'black'}

plt.figure(figsize=(6, 10))
cm = ClusterMapPlotter(data=data.loc[df_row.index.tolist(),df_col.index.tolist()],
                       top_annotation=col_ha, left_annotation=row_ha,
                       row_cluster=False,col_cluster=False,
                       label='beta', row_dendrogram=False,legend_gap=7,
                       cmap='parula',rasterized=True)
plt.savefig("Loyfer2023_heatmap.pdf",bbox_inches='tight')
plt.show()

#https://dingwb.github.io/PyComplexHeatmap/build/html/notebooks/single_cell_methylation.html
