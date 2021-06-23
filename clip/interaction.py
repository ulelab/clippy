import dash
import dash_core_components as dash_cc
import dash_html_components as dash_html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dash_bs
import plotly.graph_objects as plotlygo
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import pybedtools
import clip
import re

top_x_search_results = 20
name_delimiter = "; "
gtf_delimiter = ";"
gtf_attribute_filters = ["_name", "_id"]

class DashApp:
    def __init__(self, counts_bed, annot):
        self.read_annot(annot)
        self.counts_bed = counts_bed
        self.gene_xlink_dicts = {}
        self.gene_exon_dicts = {}
        self.app = dash.Dash(__name__, external_stylesheets=[dash_bs.themes.BOOTSTRAP])

    def read_annot(self, annot):
        annot = pd.read_table(annot, header=None, names=["chrom", "source",
            "feature_type", "start", "end", "score", "strand", "frame", "attributes"])
        # Import gene information
        print("Importing gene information...")
        self.gene_annot = annot.loc[annot.feature_type=="gene", :].copy()
        self.gene_names = self.format_gene_list(self.gene_annot.attributes.tolist())
        self.gene_annot.loc[:, 'gene_id'] = self.format_gene_list(
            self.gene_annot.attributes.tolist(), ['gene_id'])
        self.gene_annot.loc[:,'gene_names'] = self.gene_names
        # Import exon information
        print("Importing exon information...")
        self.exon_annot = annot.loc[annot.feature_type=="exon", :].copy()

    def format_gene_list(self, raw_gene_list, patterns=None):
        return([self.format_gene_name(raw_gene_name, patterns)
            for raw_gene_name in raw_gene_list])

    def format_gene_name(self, raw_gene_name, patterns=None):
        p = re.compile('(\S*).+"(\S*)"')
        gtf_dict = {
            m.group(1): m.group(2)
            for m in [
                p.search(attr.strip())
                for attr in raw_gene_name.split(gtf_delimiter)
            ] if m
        }
        if patterns == None:
            patterns = gtf_attribute_filters
        return(name_delimiter.join(set([value for key, value in gtf_dict.items()
            if any([filt in key for filt in patterns])])))

    def setup_layout(self):
        self.app.layout = dash_bs.Container([
            dash_bs.Row(dash_bs.Col(dash_html.Center(
                dash_html.H1('Clippy Interactive Parameter Search')
            ))),
            dash_bs.Row([
                dash_bs.Col([
                    dash_html.Div(id='gene-graphs')
                ], lg=9),
                dash_bs.Col(dash_html.Div([
                    dash_html.Div(dash_html.Center(dash_html.H3('Controls')), className='card-header'),
                    dash_html.Div([
                        dash_bs.ButtonGroup(
                            [
                                dash_bs.Button('Server status:', disabled=True),
                                dash_bs.Button(
                                    dash_bs.Spinner(dash_html.Div(
                                        id="graph-loading-indicator",
                                        className='m-3'
                                    )),
                                    disabled=True
                                )
                            ]
                        ),
                        dash_html.Label('Gene search'),
                        dash_html.Div([ dash_cc.Dropdown(
                            options=[
                                {'label': gene, 'value': gene}
                                for gene in self.gene_names[:top_x_search_results]],
                            id="gene-select",
                            value=[],
                            multi=True,
                            optionHeight=80
                        )], style={'marginBottom': '1.5em','marginTop': '0.5em'}),
                        dash_html.Label('Rolling mean window size'),
                        dash_html.Div([ dash_cc.Dropdown(
                            options=[
                                {'label': num, 'value': num}
                                for num in [1,5,10,15,20,50,100,200]],
                            id="n-select",
                            value=50,
                            multi=False,
                            optionHeight=20
                        )], style={'marginBottom': '1.5em','marginTop': '0.5em'}),
                        dash_html.Label('Prominence adjustment'),
                        dash_html.Div([ dash_cc.Slider(
                            id='x-slider',
                            min=0,
                            max=3,
                            step=0.1,
                            value=1,
                            marks={
                                0: '0',
                                1: '1',
                                2: '2',
                                3: '3',
                            },
                            tooltip={'always_visible': True, 'placement': 'bottom'}
                        )], style={'marginBottom': '1.5em','marginTop': '0.5em'}),
                        dash_html.Label('Relative height (broad peak threshold)'),
                        dash_html.Div([ dash_cc.Slider(
                            id='rel_height',
                            min=0,
                            max=1,
                            step=0.05,
                            value=0.8,
                            marks={
                                0: '0',
                                1: '1',
                            },
                            tooltip={'always_visible': True, 'placement': 'bottom'}
                        )], style={'marginBottom': '1.5em','marginTop': '0.5em'}),
                        dash_html.Label('Minimum counts per gene to look for peaks'),
                        dash_html.Div([ dash_cc.Slider(
                            id='min-count-slider',
                            min=1,
                            max=200,
                            step=1,
                            value=5,
                            marks={
                                1: '1',
                                200: '200',
                            },
                            tooltip={'always_visible': True, 'placement': 'bottom'}
                        )], style={'marginBottom': '1.5em','marginTop': '0.5em'}),
                                                dash_html.Label('Minimum counts per broad peak'),
                        dash_html.Div([ dash_cc.Slider(
                            id='min-peak-count-slider',
                            min=1,
                            max=200,
                            step=1,
                            value=5,
                            marks={
                                1: '1',
                                200: '200',
                            },
                            tooltip={'always_visible': True, 'placement': 'bottom'}
                        )], style={'marginBottom': '1.5em','marginTop': '0.5em'})
                    ], className='card-body')
                ], className='card bg-default sticky-top'), lg=3)
            ])
        ])

    def setup_callbacks(self):
        self.app.callback(
            Output('gene-graphs', 'children'),
            Output('graph-loading-indicator', 'children'),
            Input('gene-select', 'value'),
            Input('n-select', 'value'),
            Input('x-slider', 'value'),
            Input('rel_height', 'value'),
            Input('min-count-slider', 'value'),
            Input('min-peak-count-slider', 'value'),
            State('gene-graphs', 'children'))(self.update_figures)
        self.app.callback(
            Output('gene-select', 'options'),
            Input('gene-select', 'search_value'),
            State('gene-select', 'value'))(self.update_gene_select)

    def run(self):
        self.setup_layout()
        self.setup_callbacks()
        self.app.run_server(debug=True, host = '127.0.0.1')

    def update_gene_select(self, search_value, value):
        if not search_value:
            raise(PreventUpdate)
        p = re.compile(search_value, re.IGNORECASE)
        matching_gene_names = [
            gene_name for gene_name in self.gene_names
            if p.search(gene_name)
        ]
        return_options = [
            {'label': gene_name, 'value': gene_name}
            for gene_name in set(matching_gene_names[:top_x_search_results] + value)
        ]
        return(return_options)

    def update_figures(self, gene_list, N, X, rel_height, min_gene_count, min_peak_count, current_figures):
        # Subset the xlink BED file for each gene
        if len(gene_list) > 0:
            for gene in gene_list:
                if not gene in self.gene_xlink_dicts:
                    annot_gene = self.gene_annot.loc[self.gene_annot.gene_names==gene]
                    annot_gene = pybedtools.BedTool.from_dataframe(annot_gene)
                    self.gene_xlink_dicts[gene] = self.counts_bed.intersect(
                        annot_gene, s=True, wo=True).to_dataframe(
                        names=['chrom', 'start', 'end', 'name', 'score', 'strand',
                        'chrom2', 'source', 'feature', 'gene_start', 'gene_stop',
                        'nothing', 'strand2', 'nothing2', 'attributes', 'gene_id',
                        'gene_name', 'interval'])
                    self.gene_exon_dicts[gene] = self.exon_annot.loc[
                        self.exon_annot.attributes.str.contains(
                        self.gene_xlink_dicts[gene].gene_id[0]), :].copy()
                    self.gene_exon_dicts[gene].loc[:,'transcript_names'] = self.format_gene_list(
                        self.gene_exon_dicts[gene].attributes.tolist())
                    self.gene_xlink_dicts[gene].drop(['name', 'chrom2', 'nothing',
                        'nothing2', 'interval', 'strand2', 'source', 'feature',
                        'attributes', 'gene_id'], axis=1, inplace=True)
        if len(gene_list) == 0:
            figs = [self.peak_call_and_plot(None, N, X, rel_height, min_gene_count, min_peak_count, current_figures)]
        else:
            figs = [self.peak_call_and_plot(gene, N, X, rel_height, min_gene_count, min_peak_count, current_figures)
                for gene in gene_list]
        return(figs, 'Idle')

    def peak_call_and_plot(self, gene_name, N, X, rel_height, min_gene_count, min_peak_count, current_figures):
        # Perform the peak calling if the gene is valid
        if gene_name == None or self.gene_xlink_dicts[gene_name].shape[0] == 0:
            peaks, broad_peaks, roll_mean_smoothed_scores, peak_details = [[]]*3 + [[[]]]
        else:
            peaks, broad_peaks, roll_mean_smoothed_scores, peak_details = clip.getThePeaks(
                self.gene_xlink_dicts[gene_name], N, X, rel_height, min_gene_count, min_peak_count)
        # Plot the rolling mean and thresholds
        fig = make_subplots(rows=3, row_heights=[0.90, 0.05, 0.05], shared_xaxes=True, vertical_spacing=0.12)
        # below is code for adding relative height (broad peak) trace
        fig.add_trace(plotlygo.Scatter(
            x=np.array([
                [
                    peak_details[1]['left_ips'][idx],
                    peak_details[0][idx],
                    peak_details[1]['right_ips'][idx],
                    None
                ]
                for idx in range(len(peak_details[0]))
            ]).flatten(),
            y=np.array([
                [
                    peak_details[1]['width_heights'][idx],
                    peak_details[1]['peak_heights'][idx],
                    peak_details[1]['width_heights'][idx],
                    None
                ]
                for idx in range(len(peak_details[0]))
            ]).flatten(),
            mode='lines', name="Broad peak width",
            line=dict(color='darkorange', width=1)),
            row=1, col=1
        )
        # below is code for adding smoothed signal trace
        fig.add_trace(plotlygo.Scatter(
            x=list(range(len(roll_mean_smoothed_scores))),
            y=roll_mean_smoothed_scores,
            mode='lines',
            showlegend=False,
            line=dict(color='darkslateblue', width=2)),
            row=1, col=1)
        grid_colour = 'darkgrey'
        fig.update_xaxes(title={"text":"Position", "standoff": 0.05}, zerolinecolor=grid_colour,
            gridcolor=grid_colour, row=1, col=1)
        fig.update_yaxes(title_text='Rolling Mean Crosslink Count', zerolinecolor=grid_colour,
            gridcolor=grid_colour, row=1, col=1, fixedrange=True)
        # Add gene models
        if gene_name:
            starts = self.gene_exon_dicts[gene_name]['start'].to_numpy()
            ends = self.gene_exon_dicts[gene_name]['end'].to_numpy()
            min_value = min(starts.min(), ends.min())
            starts = starts - min_value
            ends = ends - min_value
            max_value = max(starts.max(), ends.max())

            fig.add_trace(plotlygo.Scatter(
                x=np.array([
                    [starts[idx], ends[idx], ends[idx], starts[idx], None]
                    for idx in range(len(starts))
                ]).flatten(),
                y=np.array([
                    [0, 0, 1, 1, None]
                    for idx in range(len(starts))
                ]).flatten(),
                fill='toself',
                line={'color': "#cacaca"},
                fillcolor="#cacaca",
                showlegend=False),
                row=2, col=1
            )
            fig.add_shape(type="line",
                x0=0, y0=0.5, x1=max_value, y1=0.5,
                line=dict(color="#cacaca", width=2),
                row=2, col=1
            )
        # add broad peaks as boxes
        fig.add_trace(plotlygo.Scatter(
            x=np.array([
                [
                    peak_details[1]['left_ips'][idx],
                    peak_details[1]['right_ips'][idx],
                    peak_details[1]['right_ips'][idx],
                    peak_details[1]['left_ips'][idx],
                    None
                ]
                for idx in range(len(peak_details[0]))
            ]).flatten(),
            y=np.array([
                [0, 0, 1, 1, None]
                for idx in range(len(peak_details[0]))
            ]).flatten(),
            fill='toself',
            line={'color': "darkorange"},
            fillcolor="darkorange",
            showlegend=False),
            row=3, col=1
        )
        # Remove axes from gene model figure and broad peaks track
        fig.update_xaxes(showgrid=False, zeroline=False, row=2, col=1, title={"text":"Gene model", "standoff": 1})
        fig.update_yaxes(showgrid=False, zeroline=False, visible=False, row=2, col=1)
        fig.update_xaxes(showgrid=False, zeroline=False, row=3, col=1, showticklabels=False, title={"text":"Clippy broad peaks", "standoff": 0.05})
        fig.update_yaxes(showgrid=False, zeroline=False, visible=False, row=3, col=1)
        fig.update_layout(xaxis_showticklabels=True)
        fig.update_layout(
            margin=dict(l=10, r=10, t=20, b=10),
            plot_bgcolor='rgba(0,0,0,0)',
            title={
                'text': gene_name,
                'y':0.98,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top'
            },
            legend={
                "yanchor": "top",
                "y": 0.99,
                "xanchor": "left",
                "x": 0.01
            }
        )
        if len(roll_mean_smoothed_scores) > 0:
            mean_val = np.mean(roll_mean_smoothed_scores)
            fig.add_trace(plotlygo.Scatter(
                x=list(range(len(roll_mean_smoothed_scores))),
                y=[mean_val] * len(roll_mean_smoothed_scores),
                mode='lines',
                name='Mean', 
                line=dict(color='crimson', width=2)),
                row=1, col=1)
            prominence_threshold_val = mean_val + (np.std(roll_mean_smoothed_scores)*X)
            fig.add_trace(plotlygo.Scatter(
                x=list(range(len(roll_mean_smoothed_scores))),
                y=[prominence_threshold_val] * len(roll_mean_smoothed_scores),
                mode='lines',
                name='Prominence threshold',
                line=dict(color='forestgreen', width=2)),
                row=1, col=1)
        # Add in peaks, if they have been called
        if len(peak_details[0]) > 0:
            fig.add_trace(plotlygo.Scatter(
                x=peak_details[0],
                y=[roll_mean_smoothed_scores[idx] for idx in peak_details[0]],
                mode='markers',
                name='Narrow peaks'),
                row=1, col=1)
            fig.update_traces(
                marker={
                    "size": 12,
                    "color": "mediumvioletred",
                    "line": {
                        "width": 2,
                        "color": "darkslateblue"
                    }
                },
                selector={"mode": "markers"}
            )

        current_relayout_data = None
        if current_figures:
            for graph in current_figures:
                title = graph['props']['figure']['layout']['title']
                if 'text' in title and title['text'] == gene_name:
                    if 'relayoutData' in graph['props']:
                        current_relayout_data = graph['props']['relayoutData']
        # Keep the same zoom level for the graph, if the user has changed that
        if current_relayout_data and 'xaxis.range[0]' in current_relayout_data:
            fig['layout']['xaxis']['range'] = [
                current_relayout_data['xaxis.range[0]'],
                current_relayout_data['xaxis.range[1]']
            ]
        # I think fixing y axis provides more intuitive "genome browser" experience - CC
        #if current_relayout_data and 'yaxis.range[0]' in current_relayout_data:
        #    fig['layout']['yaxis']['range'] = [
        #        current_relayout_data['yaxis.range[0]'],
        #        current_relayout_data['yaxis.range[1]']
        #    ]
        return(dash_cc.Graph(
            id='gene-graph-' + str(gene_name),
            figure=fig,
            relayoutData=current_relayout_data))
