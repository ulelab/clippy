import dash
import dash_core_components as dash_cc
import dash_html_components as dash_html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
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
        self.app = dash.Dash(__name__)

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
                p.match(attr.strip())
                for attr in raw_gene_name.split(gtf_delimiter)
            ] if m
        }
        return(name_delimiter.join(set([value for key, value in gtf_dict.items()
            if any([filt in key for filt in gtf_attribute_filters])])))

    def setup_layout(self):
        self.app.layout = dash_html.Div([
            dash_cc.Graph(id='gene-graph'),
            dash_cc.Loading(
                id="loading",
                type="default",
                children=[
                    dash_html.Div(id="graph-loading-indicator")
                ]
            ),
            dash_html.Label('Gene search'),
            dash_cc.Dropdown(
                options=[
                    {'label': gene, 'value': gene}
                    for gene in self.gene_names[:top_x_search_results]],
                id="gene-select",
                value=[self.gene_names[0]],
                multi=True
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
        ])

    def setup_callbacks(self):
        self.app.callback(
            Output('gene-graph', 'figure'),
            Output('graph-loading-indicator', 'value'),
            Input('gene-select', 'value'),
            Input('n-slider', 'value'),
            Input('x-slider', 'value'),
            Input('min-count-slider', 'value'),
            State('gene-graph', 'relayoutData'))(self.update_figure)
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
        matching_gene_names = [
            gene_name for gene_name in self.gene_names
            if search_value in gene_name
        ]
        return_options = [
            {'label': gene_name, 'value': gene_name}
            for gene_name in set(matching_gene_names[:top_x_search_results] + value)
        ]
        return(return_options)

    def update_figure(self, gene, N, X, min_gene_count, current_figure):
        # Perform the peak calling
        if len(gene) > 0:
            gene = gene[0]
            if not gene in self.gene_xlink_dicts:
                annot_gene = self.annot.loc[self.annot.gene_names==gene]
                annot_gene = pybedtools.BedTool.from_dataframe(annot_gene)
                self.gene_xlink_dicts[gene] = self.counts_bed.intersect(annot_gene, s=True, wo=True).to_dataframe(
                    names=['chrom', 'start', 'end', 'name', 'score', 'strand','chrom2','source','feature',
                    'gene_start','gene_stop','nothing','strand2','nothing2','attributes','gene_name','interval'])
                self.gene_xlink_dicts[gene].drop(['name','chrom2','nothing','nothing2','interval','strand2','source',
                    'feature','attributes'], axis=1, inplace=True)

        # Perform the peak calling if the gene is valid
        if len(gene) == 0 or self.gene_xlink_dicts[gene].shape[0] == 0:
            peaks, roll_mean_smoothed_scores, plotting_peaks = [[]]*3
        else:
            peaks, roll_mean_smoothed_scores, plotting_peaks = clip.getThePeaks(
                self.gene_xlink_dicts[gene], N, X, min_gene_count, counter=1)

        # Plot the rolling mean and thresholds
        fig = plotlyex.line(
            {
                "position": list(range(len(roll_mean_smoothed_scores))),
                "roll_mean_smoothed_scores": roll_mean_smoothed_scores
            },
            x="position", y="roll_mean_smoothed_scores"
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
        # Keep the same zoom level for the graph, if the user has changed that
        if current_figure and 'xaxis.range[0]' in current_figure:
            fig['layout']['xaxis']['range'] = [
                current_figure['xaxis.range[0]'],
                current_figure['xaxis.range[1]']
            ]
        if current_figure and 'yaxis.range[0]' in current_figure:
            fig['layout']['yaxis']['range'] = [
                current_figure['yaxis.range[0]'],
                current_figure['yaxis.range[1]']
            ]
        return(fig, '')
