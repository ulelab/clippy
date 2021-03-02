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

top_x_search_results = 20

class DashApp:
    def __init__(self, counts_bed, annot):
        self.read_annot(annot)
        self.counts_bed = counts_bed
        self.app = dash.Dash(__name__)

    def read_annot(self, annot):
        self.annot = pd.read_table(annot, header=None, names=["chrom", "source",
            "feature_type", "start", "end", "score", "strand", "frame", "attributes"])
        self.annot = self.annot[self.annot.feature_type=="gene"]
        self.gene_names = self.annot.attributes.tolist()

    def setup_layout(self):
        self.app.layout = dash_html.Div([
            dash_cc.Loading(
                id="loading-1",
                type="default",
                children=[
                    dash_cc.Graph(id='gene-graph')
                ]
            ),
            dash_html.Label('Gene search'),
            dash_cc.Dropdown(
                options=[
                    {'label': gene, 'value': gene}
                    for gene in self.gene_names[:5]],
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
            annot_gene = self.annot.loc[self.annot.attributes==gene[0]]
            annot_gene = pybedtools.BedTool.from_dataframe(annot_gene)
            gene_xlink_overlap = self.counts_bed.intersect(annot_gene, s=True, wo=True).to_dataframe(
                names=['chrom', 'start', 'end', 'name', 'score', 'strand','chrom2','source','feature',
                'gene_start', 'gene_stop','nothing','strand2','nothing2','gene_name','interval'])
            gene_xlink_overlap.drop(['name','chrom2','nothing','nothing2','interval','strand2','source',
                'feature'], axis=1, inplace=True)
        # Check if there is no gene information or there are no 
        if len(gene) == 0 or gene_xlink_overlap.shape[0] == 0:
            peaks ,roll_mean_smoothed_scores, plotting_peaks = [[]]*3
        else:
            peaks, roll_mean_smoothed_scores, plotting_peaks = clip.getThePeaks(
                gene_xlink_overlap, N, X, min_gene_count, counter=1)

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
        return(fig)
