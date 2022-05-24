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
# Exposed here to allow for optimisation with regards to the xlinks saved
max_window_size = 1000


class DashApp:
    def __init__(self, counts_bed, annot, genome_file):
        self.genome_file = genome_file
        self.read_annot(annot)
        self.counts_bed = counts_bed
        self.gene_xlink_dicts = {}
        self.gene_exon_dicts = {}
        self.gene_overlap_dict = {}
        self.app = dash.Dash(__name__, external_stylesheets=[dash_bs.themes.BOOTSTRAP])
        self.base_command_list = [
            "./clip.py",
            "-i",
            counts_bed.__dict__["fn"],
            "-a",
            annot,
            "-g",
            genome_file,
            "-o",
            "OUTPUT_PREFIX",
        ]

    def read_annot(self, annot):
        annot = pd.read_table(
            annot,
            header=None,
            names=[
                "chrom",
                "source",
                "feature_type",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attributes",
            ],
            comment="#",
        )
        # Import gene information
        print("Importing gene information...")
        self.gene_annot = annot.loc[annot.feature_type == "gene", :].copy()
        self.gene_names = self.format_gene_list(self.gene_annot.attributes.tolist())
        self.gene_annot.loc[:, "gene_id"] = self.format_gene_list(
            self.gene_annot.attributes.tolist(), ["gene_id"]
        )
        self.gene_annot.loc[:, "gene_names"] = self.gene_names
        # Import exon information
        print("Importing exon information...")
        self.exon_annot = annot.loc[annot.feature_type == "exon", :].copy()
        # Create maximal flanks
        self.maximal_flanks = clip.get_flanking_regions(
            pybedtools.BedTool.from_dataframe(self.gene_annot.iloc[:, :-2]).sort(),
            self.genome_file,
            max_window_size,
            max_window_size,
        )

    def format_gene_list(self, raw_gene_list, patterns=None):
        return [
            self.format_gene_name(raw_gene_name, patterns)
            for raw_gene_name in raw_gene_list
        ]

    def format_gene_name(self, raw_gene_name, patterns=None):
        p = re.compile(r'(\S*).+"(\S*)"')
        gtf_dict = {
            m.group(1): m.group(2)
            for m in [
                p.search(attr.strip()) for attr in raw_gene_name.split(gtf_delimiter)
            ]
            if m
        }
        if patterns == None:
            patterns = gtf_attribute_filters
        return name_delimiter.join(
            set(
                [
                    value
                    for key, value in gtf_dict.items()
                    if any([filt in key for filt in patterns])
                ]
            )
        )

    def setup_layout(self):
        self.app.layout = dash_bs.Container(
            [
                dash_bs.Row(
                    dash_bs.Col(
                        dash_html.Center(
                            dash_html.H1("Clippy Interactive Parameter Search")
                        )
                    )
                ),
                dash_bs.Row(
                    [
                        dash_bs.Col(
                            dash_html.H5("Run this on the command line with:"), lg=3
                        ),
                        dash_bs.Col(
                            dash_html.Div(
                                dash_html.Code(
                                    "<clippy command>",
                                    id="clippy-command",
                                    style={"code": {"color": "#000000"}},
                                ),
                                className="alert alert-secondary",
                                role="alert",
                            ),
                            lg=9,
                        ),
                    ]
                ),
                dash_bs.Row(
                    [
                        dash_bs.Col([dash_html.Div(id="gene-graphs")], lg=9),
                        dash_bs.Col(
                            dash_html.Div(
                                [
                                    dash_html.Div(
                                        dash_html.Center(dash_html.H3("Controls")),
                                        className="card-header",
                                    ),
                                    dash_html.Div(
                                        [
                                            dash_bs.ButtonGroup(
                                                [
                                                    dash_bs.Button(
                                                        "Server status:", disabled=True
                                                    ),
                                                    dash_bs.Button(
                                                        dash_bs.Spinner(
                                                            dash_html.Div(
                                                                id="graph-loading-indicator",
                                                                className="m-3",
                                                            )
                                                        ),
                                                        disabled=True,
                                                    ),
                                                ]
                                            ),
                                            dash_html.Label("Gene search"),
                                            dash_html.Div(
                                                [
                                                    dash_cc.Dropdown(
                                                        options=[
                                                            {
                                                                "label": gene,
                                                                "value": gene,
                                                            }
                                                            for gene in self.gene_names[
                                                                :top_x_search_results
                                                            ]
                                                        ],
                                                        id="gene-select",
                                                        value=[],
                                                        multi=True,
                                                        optionHeight=80,
                                                    )
                                                ],
                                                style={
                                                    "marginBottom": "1.5em",
                                                    "marginTop": "0.5em",
                                                },
                                            ),
                                            dash_html.Label("Rolling mean window size"),
                                            dash_html.Div(
                                                [
                                                    dash_cc.Dropdown(
                                                        options=[
                                                            {"label": num, "value": num}
                                                            for num in [
                                                                1,
                                                                5,
                                                                10,
                                                                15,
                                                                20,
                                                                50,
                                                                100,
                                                                200,
                                                            ]
                                                        ],
                                                        id="n-select",
                                                        value=15,
                                                        multi=False,
                                                        optionHeight=20,
                                                    )
                                                ],
                                                style={
                                                    "marginBottom": "1.5em",
                                                    "marginTop": "0.5em",
                                                },
                                            ),
                                            dash_html.Label("Prominence adjustment"),
                                            dash_html.Div(
                                                [
                                                    dash_cc.Slider(
                                                        id="x-slider",
                                                        min=0,
                                                        max=3,
                                                        step=0.1,
                                                        value=1,
                                                        marks={
                                                            0: "0",
                                                            1: "1",
                                                            2: "2",
                                                            3: "3",
                                                        },
                                                        tooltip={
                                                            "always_visible": True,
                                                            "placement": "bottom",
                                                        },
                                                    )
                                                ],
                                                style={
                                                    "marginBottom": "1.5em",
                                                    "marginTop": "0.5em",
                                                },
                                            ),
                                            dash_html.Label(
                                                "Minimum height adjustment"
                                            ),
                                            dash_html.Div(
                                                [
                                                    dash_cc.Slider(
                                                        id="min-height-adjust-slider",
                                                        min=0,
                                                        max=10,
                                                        step=0.1,
                                                        value=1,
                                                        tooltip={
                                                            "always_visible": True,
                                                            "placement": "bottom",
                                                        },
                                                    )
                                                ],
                                                style={
                                                    "marginBottom": "1.5em",
                                                    "marginTop": "0.5em",
                                                },
                                            ),
                                            dash_html.Label(
                                                "Relative height (peak threshold)"
                                            ),
                                            dash_html.Div(
                                                [
                                                    dash_cc.Slider(
                                                        id="rel_height",
                                                        min=0,
                                                        max=1,
                                                        step=0.05,
                                                        value=0.8,
                                                        marks={0: "0", 1: "1",},
                                                        tooltip={
                                                            "always_visible": True,
                                                            "placement": "bottom",
                                                        },
                                                    )
                                                ],
                                                style={
                                                    "marginBottom": "1.5em",
                                                    "marginTop": "0.5em",
                                                },
                                            ),
                                            dash_html.Label(
                                                "Minimum counts per gene to look for peaks"
                                            ),
                                            dash_html.Div(
                                                [
                                                    dash_cc.Slider(
                                                        id="min-count-slider",
                                                        min=1,
                                                        max=200,
                                                        step=1,
                                                        value=5,
                                                        marks={1: "1", 200: "200",},
                                                        tooltip={
                                                            "always_visible": True,
                                                            "placement": "bottom",
                                                        },
                                                    )
                                                ],
                                                style={
                                                    "marginBottom": "1.5em",
                                                    "marginTop": "0.5em",
                                                },
                                            ),
                                            dash_html.Label("Minimum counts per peak"),
                                            dash_html.Div(
                                                [
                                                    dash_cc.Slider(
                                                        id="min-peak-count-slider",
                                                        min=1,
                                                        max=200,
                                                        step=1,
                                                        value=5,
                                                        marks={1: "1", 200: "200",},
                                                        tooltip={
                                                            "always_visible": True,
                                                            "placement": "bottom",
                                                        },
                                                    )
                                                ],
                                                style={
                                                    "marginBottom": "1.5em",
                                                    "marginTop": "0.5em",
                                                },
                                            ),
                                            dash_html.Label("Alternative features"),
                                            dash_html.Div(
                                                [
                                                    dash_cc.Textarea(
                                                        id="alt-features-input",
                                                        placeholder="e.g. scaRNA-gene_type-scaRNA",
                                                    )
                                                ],
                                                style={
                                                    "marginBottom": "1.5em",
                                                    "marginTop": "0.5em",
                                                },
                                            ),
                                            dash_html.Label("Upstream extension"),
                                            dash_html.Div(
                                                [
                                                    dash_cc.Slider(
                                                        id="up-ext-slider",
                                                        min=0,
                                                        max=max_window_size,
                                                        step=50,
                                                        value=0,
                                                        tooltip={
                                                            "always_visible": True,
                                                            "placement": "bottom",
                                                        },
                                                    )
                                                ],
                                                style={
                                                    "marginBottom": "1.5em",
                                                    "marginTop": "0.5em",
                                                },
                                            ),
                                            dash_html.Label("Downstream extension"),
                                            dash_html.Div(
                                                [
                                                    dash_cc.Slider(
                                                        id="down-ext-slider",
                                                        min=0,
                                                        max=max_window_size,
                                                        step=50,
                                                        value=0,
                                                        tooltip={
                                                            "always_visible": True,
                                                            "placement": "bottom",
                                                        },
                                                    )
                                                ],
                                                style={
                                                    "marginBottom": "1.5em",
                                                    "marginTop": "0.5em",
                                                },
                                            ),
                                            dash_html.Div(
                                                [
                                                    dash_cc.Checklist(
                                                        id="exon-intron-bool",
                                                        options=[
                                                            {
                                                                "label": "Separate exon thresholds",
                                                                "value": 1,
                                                            }
                                                        ],
                                                        value=[1],
                                                    )
                                                ],
                                                style={
                                                    "marginBottom": "1.5em",
                                                    "marginTop": "0.5em",
                                                },
                                            ),
                                        ],
                                        className="card-body",
                                    ),
                                ],
                                className="card bg-default sticky-top",
                            ),
                            lg=3,
                        ),
                    ]
                ),
            ]
        )

    def setup_callbacks(self):
        self.app.callback(
            Output("gene-graphs", "children"),
            Output("clippy-command", "children"),
            Output("graph-loading-indicator", "children"),
            Input("gene-select", "value"),
            Input("n-select", "value"),
            Input("x-slider", "value"),
            Input("rel_height", "value"),
            Input("min-count-slider", "value"),
            Input("min-peak-count-slider", "value"),
            Input("alt-features-input", "value"),
            Input("up-ext-slider", "value"),
            Input("down-ext-slider", "value"),
            Input("exon-intron-bool", "value"),
            Input("min-height-adjust-slider", "value"),
            State("gene-graphs", "children"),
        )(self.update_figures)
        self.app.callback(
            Output("gene-select", "options"),
            Input("gene-select", "search_value"),
            State("gene-select", "value"),
        )(self.update_gene_select)

    def run(self):
        self.setup_layout()
        self.setup_callbacks()
        self.app.run_server(debug=True, host="127.0.0.1")

    def update_gene_select(self, search_value, value):
        if not search_value:
            raise (PreventUpdate)
        p = re.compile(search_value, re.IGNORECASE)
        matching_gene_names = [
            gene_name for gene_name in self.gene_names if p.search(gene_name)
        ]
        return_options = [
            {"label": gene_name, "value": gene_name}
            for gene_name in set(matching_gene_names[:top_x_search_results] + value)
        ]
        return return_options

    def get_maximal_overlapping_genes(self, gene):
        """
        Gets the genes which, with maximal flank sizes, are overlapping the
        gene of interest. By getting all genes at this extreme, it makes
        operations later on more efficient as they don't have to be done for
        all genes, just the gene of intereract and any overlapping
        """
        # Merge all gene annotation with all maximal flanks, calculated when the
        # Dash app was initiated
        annot_bed_with_flanks = (
            pybedtools.BedTool.from_dataframe(self.gene_annot.iloc[:, :9])
            .cat(self.maximal_flanks, postmerge=False)
            .sort()
        )
        # Convert to a pandas dataframe, and add gene_id and gene_names
        annot_bed_with_flanks_df = annot_bed_with_flanks.to_dataframe()
        annot_bed_with_flanks_df.loc[:, "gene_id"] = self.format_gene_list(
            annot_bed_with_flanks_df.attributes.tolist(), ["gene_id"]
        )
        annot_bed_with_flanks_df.loc[:, "gene_names"] = self.format_gene_list(
            annot_bed_with_flanks_df.attributes.tolist()
        )
        # Select the gene and flanks of the gene of interest and convert to
        # Bedtool object
        gene_with_maximal_flanks = pybedtools.BedTool.from_dataframe(
            annot_bed_with_flanks_df.loc[
                annot_bed_with_flanks_df.gene_names == gene
            ].iloc[:, :9]
        ).sort()
        # Find the genes and flanks which overlap with the gene of interest
        maximal_flank_overlaps = annot_bed_with_flanks.intersect(
            gene_with_maximal_flanks, s=True, wa=True
        ).to_dataframe()
        # Return the gene_annot of overlapping genes
        return (
            gene_with_maximal_flanks,
            self.gene_annot.loc[
                self.gene_annot.gene_names.isin(
                    self.format_gene_list(maximal_flank_overlaps.attributes.tolist())
                )
            ].copy(True),
        )

    def update_figures(
        self,
        gene_list,
        N,
        X,
        rel_height,
        min_gene_count,
        min_peak_count,
        alt_feature_search,
        up_ext,
        down_ext,
        exon_intron_bool,
        min_height_adjust,
        current_figures,
    ):
        # Subset the xlink BED file for each gene
        if len(gene_list) > 0:
            for gene in gene_list:
                if not gene in self.gene_xlink_dicts:
                    # If the gene has been selected for the first time (e.g.
                    # the user hasn't searched for it before) then add the gene
                    # to the various data structures

                    # Get the gene with maximal flanks and set the overlap
                    # dictionary entry. The latter stores genes which can
                    # potentially overlap the gene if the flanks are large
                    # enough
                    (
                        gene_with_maximal_flanks,
                        self.gene_overlap_dict[gene],
                    ) = self.get_maximal_overlapping_genes(gene)
                    # Find the crosslinks which overlap the gene with maximal
                    # flanks
                    self.gene_xlink_dicts[gene] = self.counts_bed.intersect(
                        gene_with_maximal_flanks, s=True, wo=True
                    ).to_dataframe(
                        names=[
                            "chrom",
                            "start",
                            "end",
                            "name",
                            "score",
                            "strand",
                            "chrom2",
                            "source",
                            "feature",
                            "gene_start",
                            "gene_stop",
                            "nothing",
                            "strand2",
                            "nothing2",
                            "attributes",
                            "interval",
                        ]
                    )
                    # If there are any crosslinks in the gene
                    if len(self.gene_xlink_dicts[gene]) > 0:
                        # Find the gene's exons
                        gene_id = self.format_gene_list(
                            self.gene_xlink_dicts[gene].attributes.tolist(),
                            ["gene_id"],
                        )[0]
                        self.gene_exon_dicts[gene] = self.exon_annot.loc[
                            self.exon_annot.attributes.str.contains(
                                '"{}"'.format(gene_id)
                            ),
                            :,
                        ].copy()
                        self.gene_exon_dicts[gene].loc[
                            :, "transcript_names"
                        ] = self.format_gene_list(
                            self.gene_exon_dicts[gene].attributes.tolist()
                        )
                        # Drop unnecessary columns from the crosslink table
                        self.gene_xlink_dicts[gene].drop(
                            [
                                "name",
                                "chrom2",
                                "nothing",
                                "nothing2",
                                "interval",
                                "strand2",
                                "source",
                                "feature",
                                "attributes",
                            ],
                            axis=1,
                            inplace=True,
                        )
        if len(gene_list) == 0:
            figs = [
                self.peak_call_and_plot(
                    None,
                    N,
                    X,
                    rel_height,
                    min_gene_count,
                    min_peak_count,
                    alt_feature_search,
                    up_ext,
                    down_ext,
                    exon_intron_bool,
                    min_height_adjust,
                    current_figures,
                )
            ]
        else:
            figs = [
                self.peak_call_and_plot(
                    gene,
                    N,
                    X,
                    rel_height,
                    min_gene_count,
                    min_peak_count,
                    alt_feature_search,
                    up_ext,
                    down_ext,
                    exon_intron_bool,
                    min_height_adjust,
                    current_figures,
                )
                for gene in gene_list
            ]
        # Setup the command list
        command_list = self.base_command_list[:]
        command_list += ["-n", str(N)]
        command_list += ["-up", str(up_ext)]
        command_list += ["-down", str(down_ext)]
        command_list += ["-x", str(X)]
        command_list += ["-hc", str(rel_height)]
        command_list += ["-mg", str(min_gene_count)]
        command_list += ["-mb", str(min_peak_count)]
        command_list += ["-mx", str(min_height_adjust)]
        if alt_feature_search:
            command_list += ["-alt", alt_feature_search]
        if len(exon_intron_bool) == 0:
            command_list.append("--no_exon_info")
        return (figs, " ".join(command_list), "Idle")

    def peak_call_and_plot(
        self,
        gene_name,
        N,
        X,
        rel_height,
        min_gene_count,
        min_peak_count,
        alt_feature_search,
        up_ext,
        down_ext,
        exon_intron_bool,
        min_height_adjust,
        current_figures,
    ):
        # Perform the peak calling if the gene is valid
        if gene_name == None or self.gene_xlink_dicts[gene_name].shape[0] == 0:
            (
                peaks,
                broad_peaks,
                roll_mean_smoothed_scores,
                peak_details,
                heights,
                prominences,
            ) = ([[]] * 3 + [[[]]] + [[]] * 2)
        else:
            annot_alt_features = None
            if alt_feature_search:
                annot_alt_features = clip.return_alt_features(
                    alt_feature_search, self.gene_annot.iloc[:, :-2],
                )
            annot_gene = self.gene_annot.loc[
                self.gene_annot.gene_names == gene_name
            ].copy()
            assert len(annot_gene) == 1

            # quoting=3 is necessary here has otherwise pybedtools double quotes
            # the attribute field, causing the attribute matching to fail. I
            # think this is a bug in pybedtools, caused by it using "to_csv" to
            # create the file read by BedTools
            gene_flanks_bed = clip.get_flanking_regions(
                pybedtools.BedTool.from_dataframe(
                    self.gene_overlap_dict[gene_name].iloc[:, :9], quoting=3
                ),
                self.genome_file,
                max(N, up_ext),
                max(N, down_ext),
            )

            # Filters for the flanks which match the gene of interest
            def attributes_match(interval):
                if interval.fields[8] == annot_gene.attributes.iloc[0]:
                    return interval

            annot_gene_bed = (
                pybedtools.BedTool.from_dataframe(annot_gene.iloc[:, :9])
                .cat(gene_flanks_bed.each(attributes_match).saveas(), postmerge=False)
                .sort()
            )

            gene_with_flanks_df = annot_gene_bed.to_dataframe(
                names=[
                    "chrom",
                    "source",
                    "feature_type",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "frame",
                    "attributes",
                ],
            )

            xlinks = (
                pybedtools.BedTool.from_dataframe(
                    pd.DataFrame.from_dict(
                        {
                            "chrom": self.gene_xlink_dicts[gene_name].chrom,
                            "start": self.gene_xlink_dicts[gene_name].start,
                            "end": self.gene_xlink_dicts[gene_name].end,
                            "name": ".",
                            "score": self.gene_xlink_dicts[gene_name].score,
                            "strand": self.gene_xlink_dicts[gene_name].strand,
                        }
                    )
                )
                .intersect(annot_gene_bed, s=True, wo=True)
                .to_dataframe(
                    names=[
                        "chrom",
                        "start",
                        "end",
                        "name",
                        "score",
                        "strand",
                        "chrom2",
                        "source",
                        "feature",
                        "gene_start",
                        "gene_stop",
                        "nothing",
                        "strand2",
                        "nothing2",
                        "attributes",
                        "interval",
                    ]
                )
                .drop(
                    [
                        "name",
                        "chrom2",
                        "nothing",
                        "nothing2",
                        "interval",
                        "strand2",
                        "source",
                        "feature",
                        "attributes",
                    ],
                    axis=1,
                )
            )
            xlinks.loc[:, "gene_name"] = gene_name
            annot_exon = self.gene_exon_dicts[gene_name]
            if len(exon_intron_bool) == 0:
                annot_exon = None
            (
                peaks,
                broad_peaks,
                roll_mean_smoothed_scores,
                peak_details,
                heights,
                prominences,
            ) = clip.single_gene_get_peaks(
                xlinks,
                N,
                X,
                rel_height,
                min_gene_count,
                min_peak_count,
                annot_exon,
                annot_alt_features,
                gene_with_flanks_df,
                min_height_adjust,
            )
            if not isinstance(peaks, np.ndarray):
                (
                    peaks,
                    broad_peaks,
                    roll_mean_smoothed_scores,
                    peak_details,
                    heights,
                    prominences,
                ) = ([[]] * 3 + [[[]]] + [[]] * 2)
        # Plot the rolling mean and thresholds
        fig = make_subplots(
            rows=3,
            row_heights=[0.90, 0.05, 0.05],
            shared_xaxes=True,
            vertical_spacing=0.12,
        )
        # below is code for adding relative height (broad peak) trace
        fig.add_trace(
            plotlygo.Scatter(
                x=np.array(
                    [
                        [
                            peak_details[1]["left_ips"][idx],
                            peak_details[0][idx],
                            peak_details[1]["right_ips"][idx],
                            None,
                        ]
                        for idx in range(len(peak_details[0]))
                    ]
                ).flatten(),
                y=np.array(
                    [
                        [
                            peak_details[1]["width_heights"][idx],
                            peak_details[1]["peak_heights"][idx],
                            peak_details[1]["width_heights"][idx],
                            None,
                        ]
                        for idx in range(len(peak_details[0]))
                    ]
                ).flatten(),
                mode="lines",
                name="Peak width",
                line=dict(color="darkorange", width=1),
            ),
            row=1,
            col=1,
        )
        # below is code for adding smoothed signal trace
        fig.add_trace(
            plotlygo.Scatter(
                x=list(range(len(roll_mean_smoothed_scores))),
                y=roll_mean_smoothed_scores,
                mode="lines",
                showlegend=False,
                line=dict(color="darkslateblue", width=2),
            ),
            row=1,
            col=1,
        )
        grid_colour = "darkgrey"
        fig.update_xaxes(
            title={"text": "Position", "standoff": 0.05},
            zerolinecolor=grid_colour,
            gridcolor=grid_colour,
            row=1,
            col=1,
        )
        fig.update_yaxes(
            title_text="Rolling Mean Crosslink Count",
            zerolinecolor=grid_colour,
            gridcolor=grid_colour,
            row=1,
            col=1,
        )
        # Add gene models
        if gene_name and not self.gene_xlink_dicts[gene_name].shape[0] == 0:
            flank_start = min(gene_with_flanks_df.start)
            flank_stop = max(gene_with_flanks_df.end)

            gene_start = annot_gene.start.iloc[0]
            gene_stop = annot_gene.end.iloc[0]

            exon_starts = self.gene_exon_dicts[gene_name]["start"].to_numpy()
            exon_ends = self.gene_exon_dicts[gene_name]["end"].to_numpy()
            exon_starts = exon_starts - flank_start
            exon_ends = exon_ends - flank_start

            fig.add_trace(
                plotlygo.Scatter(
                    x=np.array(
                        [
                            [
                                exon_starts[idx],
                                exon_ends[idx],
                                exon_ends[idx],
                                exon_starts[idx],
                                None,
                            ]
                            for idx in range(len(exon_starts))
                        ]
                    ).flatten(),
                    y=np.array(
                        [[0, 0, 1, 1, None] for idx in range(len(exon_starts))]
                    ).flatten(),
                    fill="toself",
                    mode="lines",
                    line={"color": "#cacaca"},
                    fillcolor="#cacaca",
                    showlegend=False,
                ),
                row=2,
                col=1,
            )

            fig.add_shape(
                type="line",
                x0=gene_start - flank_start,
                y0=0.5,
                x1=gene_stop - flank_start,
                y1=0.5,
                line=dict(color="#cacaca", width=2),
                row=2,
                col=1,
                layer="below",
            )

            # add arrow corresponding to gene strand
            # get strand - how?????
            direction = self.gene_xlink_dicts[gene_name]["strand"][0]
            if direction == "+":
                marksymb = "triangle-right"
            else:
                marksymb = "triangle-left"
            fig.add_trace(
                plotlygo.Scatter(
                    x=np.array(
                        [
                            (gene_start - flank_start)
                            + (gene_stop - gene_start) * 0.25,
                            (gene_start - flank_start) + (gene_stop - gene_start) * 0.5,
                            (gene_start - flank_start)
                            + (gene_stop - gene_start) * 0.75,
                        ]
                    ),
                    y=np.array([0.5, 0.5, 0.5]),
                    fill="toself",
                    mode="markers",
                    marker_symbol=marksymb,
                    marker_size=20,
                    marker_color="black",
                    showlegend=False,
                ),
                row=2,
                col=1,
            )

            filtered_broad_peaks = broad_peaks

            # add broad peaks as boxes
            fig.add_trace(
                plotlygo.Scatter(
                    x=np.array(
                        [
                            [
                                float(filtered_broad_peaks[idx][1]) - flank_start,
                                float(filtered_broad_peaks[idx][2]) - flank_start,
                                float(filtered_broad_peaks[idx][2]) - flank_start,
                                float(filtered_broad_peaks[idx][1]) - flank_start,
                                None,
                            ]
                            for idx in range(len(filtered_broad_peaks))
                        ]
                    ).flatten(),
                    y=np.array(
                        [[0, 0, 1, 1, None] for idx in range(len(peak_details[0]))]
                    ).flatten(),
                    fill="toself",
                    mode="lines",
                    line={"color": "darkorange"},
                    fillcolor="darkorange",
                    showlegend=False,
                ),
                row=3,
                col=1,
            )
        # Remove axes from gene model figure and broad peaks track
        fig.update_xaxes(
            showgrid=False,
            zeroline=False,
            row=2,
            col=1,
            title={"text": "Gene model", "standoff": 1},
        )
        fig.update_yaxes(showgrid=False, zeroline=False, visible=False, row=2, col=1)
        fig.update_xaxes(
            showgrid=False,
            zeroline=False,
            row=3,
            col=1,
            showticklabels=False,
            title={"text": "Clippy peaks", "standoff": 0.05},
        )
        fig.update_yaxes(showgrid=False, zeroline=False, visible=False, row=3, col=1)
        fig.update_layout(xaxis_showticklabels=True)
        # Calculate the xlink number
        xlink_number = 0
        if (
            gene_name in self.gene_xlink_dicts
            and not self.gene_xlink_dicts[gene_name].shape[0] == 0
        ):
            xlink_number = self.gene_xlink_dicts[gene_name]["score"].sum()
        plot_title = (
            gene_name + " ; Total xlinks = {}".format(xlink_number)
            if gene_name is not None
            else gene_name
        )
        fig.update_layout(
            margin=dict(l=10, r=10, t=20, b=10),
            plot_bgcolor="rgba(0,0,0,0)",
            title={
                "text": plot_title,
                "y": 0.98,
                "x": 0.5,
                "xanchor": "center",
                "yanchor": "top",
            },
            legend={"yanchor": "top", "y": 0.99, "xanchor": "left", "x": -0.5},
        )
        if len(roll_mean_smoothed_scores) > 0:
            mean_val = np.mean(roll_mean_smoothed_scores)
            fig.add_trace(
                plotlygo.Scatter(
                    x=list(range(len(roll_mean_smoothed_scores))),
                    y=heights,
                    mode="lines",
                    name="Mean",
                    line=dict(color="crimson", width=2),
                ),
                row=1,
                col=1,
            )
            fig.add_trace(
                plotlygo.Scatter(
                    x=list(range(len(roll_mean_smoothed_scores))),
                    y=heights
                    + prominences,  # [prominence_threshold_val] * len(roll_mean_smoothed_scores),
                    mode="lines",
                    name="Prominence threshold",
                    line=dict(color="forestgreen", width=2),
                ),
                row=1,
                col=1,
            )
        # Add in peaks, if they have been called
        if len(peak_details[0]) > 0:
            fig.add_trace(
                plotlygo.Scatter(
                    x=peak_details[0],
                    y=[roll_mean_smoothed_scores[idx] for idx in peak_details[0]],
                    mode="markers",
                    marker_size=12,
                    marker_color="mediumvioletred",
                    marker_line={"width": 2, "color": "darkslateblue"},
                    name="Peak summits",
                ),
                row=1,
                col=1,
            )
            # fig.update_traces(
            #    marker={
            #        "size": 12,
            #        "color": "mediumvioletred",
            #        "line": {
            #            "width": 2,
            #            "color": "darkslateblue"
            #        }
            #    },
            #    selector={"mode": "markers"}
            # )

        current_relayout_data = None
        if current_figures:
            for graph in current_figures:
                title = graph["props"]["figure"]["layout"]["title"]
                if "text" in title and title["text"] == gene_name:
                    if "relayoutData" in graph["props"]:
                        current_relayout_data = graph["props"]["relayoutData"]
        # Keep the same zoom level for the graph, if the user has changed that
        if current_relayout_data and "xaxis.range[0]" in current_relayout_data:
            fig["layout"]["xaxis"]["range"] = [
                current_relayout_data["xaxis.range[0]"],
                current_relayout_data["xaxis.range[1]"],
            ]
        # I think fixing y axis provides more intuitive "genome browser" experience - CC
        if current_relayout_data and "yaxis.range[0]" in current_relayout_data:
            fig["layout"]["yaxis"]["range"] = [
                current_relayout_data["yaxis.range[0]"],
                current_relayout_data["yaxis.range[1]"],
            ]
        return dash_cc.Graph(
            id="gene-graph-" + str(gene_name),
            figure=fig,
            relayoutData=current_relayout_data,
        )
