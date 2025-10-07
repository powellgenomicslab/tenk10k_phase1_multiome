import pickle
import numpy as np
import pandas as pd
import plotly

def heatmap(
    data_array: np.ndarray,
    celltype_titles = celltype_titles,
    row_names = None,
    column_names = None,
    cluster_columns: bool = True,
    cluster_rows: bool = True,
    colorscale = "Blues",
    linkage: str = "ward",
    **kwargs,
):
    import plotly.graph_objects as go
    import plotly.figure_factory as ff
    import scipy.cluster.hierarchy as sch
    
    link_func = lambda x: sch.linkage(x, linkage)
    fig = go.Figure()
    
    if cluster_columns:
        dendro_upper = ff.create_dendrogram(data_array.T, linkagefun=link_func, orientation='bottom')
        upper_leaves = list(map(int, dendro_upper['layout']['xaxis']['ticktext']))
        data_array = data_array[:, upper_leaves]
        if column_names is not None:
            dendro_upper['layout']['xaxis']['ticktext'] = np.array(column_names)[upper_leaves]
        for i in range(len(dendro_upper['data'])):
            dendro_upper['data'][i]['yaxis'] = 'y2'
            fig.add_trace(dendro_upper['data'][i])
        fig['layout'] = dendro_upper['layout']

    if cluster_rows:
        dendro_side = ff.create_dendrogram(data_array, linkagefun=link_func, orientation='right')
        side_leaves = list(map(int, dendro_side['layout']['yaxis']['ticktext']))
        data_array = data_array[side_leaves, :]
        if row_names is not None:
            dendro_side['layout']['yaxis']['ticktext'] = np.array(row_names)[side_leaves]
        for i in range(len(dendro_side['data'])):
            dendro_side['data'][i]['xaxis'] = 'x2'
            fig.add_trace(dendro_side['data'][i])
        fig['layout']['yaxis'] = dendro_side['layout']['yaxis']

    # Create Heatmap
    heatmap = [go.Heatmap(
        z=data_array,
        colorscale=colorscale,
        colorbar=dict(
            orientation='h',
            title=dict(
                text='ln(-log<sub>10</sub>(<i>p</i>))',
                font=dict(size=18)
            ),
            tickfont=dict(size=16)
        )
    )]
    if cluster_columns:
        heatmap[0]['x'] = dendro_upper['layout']['xaxis']['tickvals']
    if cluster_rows:
        heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']
    for data in heatmap:
        fig.add_trace(data)
        
    fig.update_layout({'showlegend':False, 'hovermode': 'closest'})
    
    if cluster_rows:
        fig.update_layout(
            xaxis={
                'domain': [.15, 1], 'mirror': False, 'showgrid': False,
                'showline': False, 'zeroline': False, 'ticks':"",
            },
            xaxis2={
                'domain': [0, .15], 'mirror': False, 'showgrid': False,
                'showline': False, 'zeroline': False, 'showticklabels': False,
                'ticks':"",
            },
        )
    if cluster_columns:
        fig.update_layout(
            yaxis={
                'domain': [0, .85], 'mirror': False, 'showgrid': False,
                'showline': False, 'zeroline': False, 'showticklabels': True,
                'side': 'right', 'ticks': ""
            },
            yaxis2={
                'domain':[.825, .975], 'mirror': False, 'showgrid': False,
                'showline': False, 'zeroline': False, 'showticklabels': False,
                'ticks':""
            },
        )
    fig.update_xaxes(tickangle=270, tickfont = dict(size=15))
    ordered_celltype_titles = [celltype_titles[i] for i in side_leaves]
    fig.update_yaxes(ticktext = ordered_celltype_titles, tickfont = dict(size=15))

    fig.update_layout({
        'plot_bgcolor': 'white',
        'paper_bgcolor': 'white'
    })
    return render_plot(fig, **kwargs)

def render_plot(
    fig: 'plotly.graph_objects.Figure',
    width: int = 600,
    height: int = 400,
    interactive: bool = True,
    show: bool = True,
    out_file = None,
    scale = None,
):
    fig.update_layout({
        "width": width,
        "height": height,
    })

    # save figure to file
    if out_file is not None:
        if out_file.endswith(".html"):
            fig.write_html(out_file, include_plotlyjs="cdn")
        else:
            fig.write_image(out_file, scale=scale)

    # show figure
    if show:
        if interactive:
            fig.show()
        else:
            from IPython.display import Image
            return Image(fig.to_image(format="png"))

    # return plot object
    if not show and not out_file: return fig

# Load the data and create a heatmap
df_filtered = pd.read_csv("/g/data/ei56/od8037/Plotting/Motif_Plot/motifs.csv").set_index("Unnamed: 0")
celltype_titles = [
    'ASDC',
    'B<sub>intermediate</sub>',
    'B<sub>memory</sub>',
    'B<sub>naive</sub>',
    'CD4<sub>CTL</sub>',
    'CD4<sub>Naive</sub>',
    'CD4<sub>Proliferating</sub>',
    'CD4<sub>TCM</sub>',
    'CD4<sub>TEM</sub>',
    'CD8<sub>Naive</sub>',
    'CD8<sub>TCM</sub>',
    'CD8<sub>TEM</sub>',
    'CD14<sub>Mono</sub>',
    'CD16<sub>Mono</sub>',
    'HSPC',
    'MAIT',
    'NK',
    'NK<sub>CD56bright</sub>',
    'NK<sub>Proliferating</sub>',
    'Plasmablast',
    'T<sub>reg</sub>',
    'cDC1',
    'cDC2',
    'dnT',
    'gdT',
    'pDC'
]

heatmap(
    df_filtered.to_numpy().T,
    row_names=df_filtered.columns,
    column_names=df_filtered.index,
    colorscale='RdBu_r', scale=15,
    height=720, width=1000, interactive=False,
    out_file="/g/data/ei56/od8037/Plotting/Motif_Plot/motif_plot.png"
)