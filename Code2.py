import matplotlib.pylab as plt
plt.rcParams['figure.dpi'] = 100
plt.rcParams['savefig.dpi']=300
from PyComplexHeatmap import *

df = pd.read_csv("data2.csv",index_col="sample")
print(df)

df_box = pd.DataFrame(np.random.randn(10, 4), columns=['Gene' + str(i) for i in range(1, 5)])
df_box.index = ['sample' + str(i) for i in range(1, df_box.shape[0] + 1)]
df_bar = pd.DataFrame(np.random.uniform(0, 10, (10, 2)), columns=['TMB1', 'TMB2'])
df_bar.index = ['sample' + str(i) for i in range(1, df_box.shape[0] + 1)]
df_scatter = pd.DataFrame(np.random.uniform(0, 10, 10), columns=['Scatter'])
df_scatter.index = ['sample' + str(i) for i in range(1, df_box.shape[0] + 1)]
df_heatmap = pd.DataFrame(np.random.randn(30, 10), columns=['sample' + str(i) for i in range(1, 11)])
df_heatmap.index = ["Fea" + str(i) for i in range(1, df_heatmap.shape[0] + 1)]
df_heatmap.iloc[1, 2] = np.nan

plt.figure(figsize=(5, 10))
col_ha = HeatmapAnnotation(label=anno_label(df.AB, merge=True,rotation=15),
                           AB=anno_simple(df.AB,add_text=True),axis=1,
                           CD=anno_simple(df.CD,add_text=True,colors={'C':'red','D':'green','G':'blue'},
                                            legend_kws={'frameon':False}),
                           Exp=anno_boxplot(df_box, cmap='turbo',height=15),
                           Scatter=anno_scatterplot(df_scatter,height=15),
                           TMB_bar=anno_barplot(df_bar,height=10,legend_kws={'color_text':False,'labelcolor':'blue'}))
cm = ClusterMapPlotter(data=df_heatmap, top_annotation=col_ha, col_split=2, row_split=3, col_split_gap=0.5,
                     row_split_gap=1,label='values',row_dendrogram=True,show_rownames=False,show_colnames=False,
                     tree_kws={'row_cmap': 'Dark2'},
                       legend_gap=5,legend_hpad=2,legend_vpad=5)
#legend_gap is the gap between two legends, legend_hpad is the horizonal space between legend and heatmap, legend_vpad
# is the verticall space between the first legend and the top of axes (legend_anchor).
# cm.ax_heatmap.set_axis_off()
plt.show()

plt.figure(figsize=(6, 4))
col_ha = HeatmapAnnotation(label=anno_label(df.AB, merge=True,rotation=15),
                            AB=anno_simple(df.AB,add_text=True,legend=True), axis=1,
                            CD=anno_simple(df.CD, colors={'C': 'red', 'D': 'gray', 'G': 'yellow'},
                                           add_text=True,legend=True),
                            Exp=anno_boxplot(df_box, cmap='turbo',legend=True),
                            Scatter=anno_scatterplot(df_scatter), TMB_bar=anno_barplot(df_bar,legend=True),
                           plot=True,legend=True,legend_gap=3)
col_ha.show_ticklabels(df.index.tolist())
plt.show()


#Annotate the rows with average > 0.3
df_rows = df_heatmap.apply(lambda x:x.name if x.mean() > 0.3 else None,axis=1)
df_rows.name='Selected'

row_ha = HeatmapAnnotation(Avg=anno_simple(df_heatmap.mean(axis=1).apply(lambda x:round(x,2)),
                                           cmap='jet'), #add_text=True,,text_kws={'rotation':0,'fontsize':10,'color':'black'}
                            # Avg=anno_barplot(df_heatmap.mean(axis=1).apply(lambda x:round(x,2)),
                            #                height=10,colors='orangered'),
                           selected=anno_label(df_rows,colors='red'),
                           axis=0,verbose=0,label_kws={'rotation':0,'horizontalalignment':'left'})

col_ha = HeatmapAnnotation(label=anno_label(df.AB, merge=True,rotation=15),
                           AB=anno_simple(df.AB,add_text=True),axis=1,
                           Exp=anno_boxplot(df_box, cmap='turbo'),
                           verbose=0) #verbose=0 will turn off the log.

plt.figure(figsize=(6, 8))
cm = ClusterMapPlotter(data=df_heatmap, top_annotation=col_ha,right_annotation=row_ha,
                       col_split=2,row_split=2, col_split_gap=0.5,row_split_gap=1,
                       label='values',row_dendrogram=True,show_rownames=False,show_colnames=False,
                         tree_kws={'row_cmap': 'Set1'},verbose=0,legend_gap=7,
                       annot=True,linewidths=0.05,linecolor='gold',cmap='turbo')
plt.show()


plt.figure(figsize=(5, 4))
col_ha = HeatmapAnnotation(label=anno_label(df.AB, merge=True,rotation=15),
                            AB=anno_simple(df.AB,add_text=True,legend=True,text_kws={'color':'gold'}),
                            CD=anno_simple(df.CD,add_text=True,legend=True,text_kws={'color':'purple'}),
                            Exp=anno_boxplot(df_box, cmap='turbo',legend=True),
                            Scatter=anno_scatterplot(df_scatter), TMB_bar=anno_barplot(df_bar,legend=True),
                           plot=True,legend=True,legend_gap=5,axis=1)
col_ha.show_ticklabels(df.index.tolist())
plt.show()

plt.figure(figsize=(5, 3))
col_ha = HeatmapAnnotation(df=df,plot=True,legend=True)
plt.show()

plt.figure(figsize=(6, 4))
col_ha = HeatmapAnnotation(label=anno_label(df.AB, merge=True,rotation=15),
                            AB=anno_simple(df.AB,add_text=True,legend=True), axis=1,
                            CD=anno_simple(df.CD,add_text=True,legend=True),
                            Exp=anno_boxplot(df_box, cmap='turbo',legend=True),
                            Scatter=anno_scatterplot(df_scatter), TMB_bar=anno_barplot(df_bar,legend=True),
                           plot=True,legend=True,plot_legend=False,
                           legend_gap=5,hgap=0)
col_ha.show_ticklabels(df.index.tolist())
plt.show()

plt.figure()
row_ha.plot_legends()
plt.show()

