import dash
import dash_core_components as dash_cc
import dash_html_components as dash_html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dash_bs
import plotly.express as plotlyex
import plotly.graph_objects as plotlygo
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
        self.app = dash.Dash(__name__, external_stylesheets=[dash_bs.themes.BOOTSTRAP])

    def read_annot(self, annot):
        self.annot = pd.read_table(annot, header=None, names=["chrom", "source",
            "feature_type", "start", "end", "score", "strand", "frame", "attributes"])
        self.annot = self.annot[self.annot.feature_type=="gene"]
        self.gene_names = self.format_gene_list(self.annot.attributes.tolist())
        self.annot['gene_names'] = self.gene_names

    def format_gene_list(self, raw_gene_list):
        return([self.format_gene_name(raw_gene_name)
            for raw_gene_name in raw_gene_list])

    def format_gene_name(self, raw_gene_name):
        p = re.compile('(\S*).+"(\S*)"')
        gtf_dict = {
            m.group(1): m.group(2)
            for m in [
                p.search(attr.strip())
                for attr in raw_gene_name.split(gtf_delimiter)
            ] if m
        }
        return(name_delimiter.join(set([value for key, value in gtf_dict.items()
            if any([filt in key for filt in gtf_attribute_filters])])))

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
                        dash_cc.Dropdown(
                            options=[
                                {'label': gene, 'value': gene}
                                for gene in self.gene_names[:top_x_search_results]],
                            id="gene-select",
                            value=[],
                            multi=True,
                            optionHeight=80
                        ),
                        dash_html.Label('Rolling mean window size'),
                        dash_cc.Slider(
                            id='n-slider',
                            min=1,
                            max=200,
                            step=1,
                            value=50,
                            tooltip={'always_visible': True, 'placement': 'bottom'}
                        ),
                        dash_html.Label('Prominence adjustment'),
                        dash_cc.Slider(
                            id='x-slider',
                            min=0,
                            max=3,
                            step=0.1,
                            value=1,
                            tooltip={'always_visible': True, 'placement': 'bottom'}
                        ),
                        dash_html.Label('Minimum counts per gene to look for peaks'),
                        dash_cc.Slider(
                            id='min-count-slider',
                            min=1,
                            max=200,
                            step=1,
                            value=5,
                            tooltip={'always_visible': True, 'placement': 'bottom'}
                        )
                    ], className='card-body')
                ], className='card bg-default sticky-top'), lg=3)
            ])
        ])

    def setup_callbacks(self):
        self.app.callback(
            Output('gene-graphs', 'children'),
            Output('graph-loading-indicator', 'children'),
            Input('gene-select', 'value'),
            Input('n-slider', 'value'),
            Input('x-slider', 'value'),
            Input('min-count-slider', 'value'),
            State('gene-graphs', 'children'))(self.update_figures)
        self.app.callback(
            Output('gene-select', 'options'),
            Input('gene-select', 'search_value'),
            State('gene-select', 'value'))(self.update_gene_select)

    def run(self):
        self.setup_layout()
        self.setup_callbacks()
        self.app.run_server(debug=True)

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

    def update_figures(self, gene_list, N, X, min_gene_count, current_figures):
        # Subset the xlink BED file for each gene
        if len(gene_list) > 0:
            for gene in gene_list:
                if not gene in self.gene_xlink_dicts:
                    annot_gene = self.annot.loc[self.annot.gene_names==gene]
                    annot_gene = pybedtools.BedTool.from_dataframe(annot_gene)
                    self.gene_xlink_dicts[gene] = self.counts_bed.intersect(
                        annot_gene, s=True, wo=True).to_dataframe(
                        names=['chrom', 'start', 'end', 'name', 'score', 'strand',
                        'chrom2', 'source', 'feature', 'gene_start', 'gene_stop',
                        'nothing', 'strand2', 'nothing2', 'attributes', 'gene_name',
                        'interval'])
                    self.gene_xlink_dicts[gene].drop(['name', 'chrom2', 'nothing',
                        'nothing2', 'interval', 'strand2', 'source', 'feature',
                        'attributes'], axis=1, inplace=True)
        if len(gene_list) == 0:
            figs = [self.peak_call_and_plot(None, N, X, min_gene_count, current_figures)]
        else:
            figs = [self.peak_call_and_plot(gene, N, X, min_gene_count, current_figures)
                for gene in gene_list]
        return(figs, 'Idle')

    def peak_call_and_plot(self, gene_name, N, X, min_gene_count, current_figures):
        # Perform the peak calling if the gene is valid
        if gene_name == None or self.gene_xlink_dicts[gene_name].shape[0] == 0:
            peaks, roll_mean_smoothed_scores, plotting_peaks = [[]]*3
        else:
            peaks, roll_mean_smoothed_scores, plotting_peaks = clip.getThePeaks(
                self.gene_xlink_dicts[gene_name], N, X, min_gene_count, counter=1)
        # Plot the rolling mean and thresholds
        fig = plotlyex.line(
            {
                "position": list(range(len(roll_mean_smoothed_scores))),
                "roll_mean_smoothed_scores": roll_mean_smoothed_scores
            },
            x="position", y="roll_mean_smoothed_scores"
        )
        fig.update_layout(
            title={
                'text': gene_name,
                'y':0.9,
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
                name='mean'))
            prominence_threshold_val = mean_val + (np.std(roll_mean_smoothed_scores)*X)
            fig.add_trace(plotlygo.Scatter(
                x=list(range(len(roll_mean_smoothed_scores))),
                y=[prominence_threshold_val] * len(roll_mean_smoothed_scores),
                mode='lines',
                name='prominence threshold'))
        # Add in peaks, if they have been called
        if len(plotting_peaks) > 0:
            fig.add_trace(plotlygo.Scatter(
                x=plotting_peaks,
                y=[roll_mean_smoothed_scores[idx] for idx in plotting_peaks],
                mode='markers',
                name='peaks'))
            fig.update_traces(
                marker={
                    "size": 12,
                    "line": {
                        "width": 2,
                        "color": "DarkSlateGrey"
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
        if current_relayout_data and 'yaxis.range[0]' in current_relayout_data:
            fig['layout']['yaxis']['range'] = [
                current_relayout_data['yaxis.range[0]'],
                current_relayout_data['yaxis.range[1]']
            ]
        return(dash_cc.Graph(
            id='gene-graph-' + str(gene_name),
            figure=fig,
            relayoutData=current_relayout_data))
