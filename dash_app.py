import dash
import dash_core_components as dash_cc
import dash_html_components as dash_html
from dash.dependencies import Input, Output
import plotly.express as px
import pandas as pd

class DashApp:
    def __init__(self, counts_bed):
        self.counts_bed = counts_bed
        self.app = dash.Dash(__name__)
        self.df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/gapminderDataFiveYear.csv')

    def setup_layout(self):
        self.app.layout = dash_html.Div([
            dash_cc.Graph(id='graph-with-slider'),
            dash_cc.Slider(
                id='year-slider',
                min=self.df['year'].min(),
                max=self.df['year'].max(),
                value=self.df['year'].min(),
                marks={str(year): str(year) for year in self.df['year'].unique()},
                step=None
            )
        ])

    def setup_callbacks(self):
        self.app.callback(
            Output('graph-with-slider', 'figure'),
            Input('year-slider', 'value'))(self.update_figure)

    def run(self):
        self.setup_layout()
        self.setup_callbacks()
        self.app.run_server(debug=True)

    def update_figure(self, selected_year):
        filtered_df = self.df[self.df.year == selected_year]
        fig = px.scatter(filtered_df, x="gdpPercap", y="lifeExp",
                        size="pop", color="continent", hover_name="country",
                        log_x=True, size_max=55)

        fig.update_layout(transition_duration=500)
        return(fig)
