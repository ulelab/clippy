import clip.interaction
import pybedtools
import os


def test_web_app(dash_duo, rootdir):
    counts_bed = pybedtools.BedTool(
        os.path.join(rootdir, "tests", "data", "crosslinkcounts.bed")
    )
    app = clip.interaction.DashApp(
        counts_bed, os.path.join(rootdir, "tests", "data", "annot.gff")
    )
    app.setup_layout()
    app.setup_callbacks()
    dash_duo.start_server(app.app)
    dash_duo.find_element("input#gene-select").send_keys("icr")
    dash_duo.find_element(".VirtualizedSelectOption").click()
    dash_duo.wait_for_text_to_equal("text.gtitle", "ICR1 ; Total xlinks = 104")
    assert dash_duo.get_logs() == []
