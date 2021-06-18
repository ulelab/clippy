import clip
import pandas as pd


def test_parse_alt_features():
    assert clip.parse_alt_features(
        "lncRNA-gene_type-lncRNA,"
        "lncRNA-gene_name-SCARNA[0-9]+,"
        "snoRNA-gene_type-snoRNA"
    ) == {
        "lncRNA": [
            {"key": "gene_type", "regex": "lncRNA"},
            {"key": "gene_name", "regex": "SCARNA[0-9]+"},
        ],
        "snoRNA": [{"key": "gene_type", "regex": "snoRNA"},],
    }


def test_test_alt_features():
    assert clip.test_alt_features(
        'gene_id "ENSG00000254911.3"; gene_type "snoRNA"; gene_name "SCARNA9"; '
        'level 2; hgnc_id "HGNC:32566"; tag "ncRNA_host"; '
        'tag "overlapping_locus"; havana_gene "OTTHUMG00000167450.2";',
        {
            "lncRNA": [
                {"key": "gene_type", "regex": "lncRNA"},
                {"key": "gene_name", "regex": "SCARNA[0-9]+"},
            ],
            "snoRNA": [{"key": "gene_type", "regex": "snoRNA"},],
        },
    ) == set(["lncRNA", "snoRNA"])
    assert clip.test_alt_features(
        'gene_id "ENSG00000254911.3"; gene_type "lncRNA"; gene_name "SCARNA9"; '
        'level 2; hgnc_id "HGNC:32566"; tag "ncRNA_host"; '
        'tag "overlapping_locus"; havana_gene "OTTHUMG00000167450.2";',
        {
            "lncRNA": [
                {"key": "gene_type", "regex": "lncRNA"},
                {"key": "gene_name", "regex": "SCARNA[0-9]+"},
            ],
            "snoRNA": [{"key": "gene_type", "regex": "snoRNA"},],
        },
    ) == set(["lncRNA"])


def test_return_alt_features():
    annot_gene = pd.DataFrame.from_dict(
        {
            "attributes": [
                'gene_id "ENSG00000166004.15"; gene_type "protein_coding"; '
                'gene_name "CEP295"; level 1; hgnc_id "HGNC:29366"; '
                'tag "ncRNA_host"; tag "overlapping_locus"; '
                'havana_gene "OTTHUMG00000167449.5";',
                'gene_id "ENSG00000254911.3"; gene_type "lncRNA"; '
                'gene_name "SCARNA9"; level 2; hgnc_id "HGNC:32566"; '
                'tag "ncRNA_host"; tag "overlapping_locus"; '
                'havana_gene "OTTHUMG00000167450.2";',
                'gene_id "ENSG00000273885.1"; gene_type "snoRNA"; '
                'gene_name "ENSG10010138116.1"; level 3;',
            ]
        }
    )
    alt_features = clip.return_alt_features(
        "lncRNA-gene_type-lncRNA,"
        "lncRNA-gene_name-SCARNA[0-9]+,"
        "snoRNA-gene_type-snoRNA",
        annot_gene,
    )
    assert alt_features["lncRNA"].equals(
        pd.DataFrame.from_dict(
            {
                "attributes": [
                    'gene_id "ENSG00000254911.3"; gene_type "lncRNA"; '
                    'gene_name "SCARNA9"; level 2; hgnc_id "HGNC:32566"; '
                    'tag "ncRNA_host"; tag "overlapping_locus"; '
                    'havana_gene "OTTHUMG00000167450.2";',
                ]
            }
        )
    )
    assert alt_features["snoRNA"].equals(
        pd.DataFrame.from_dict(
            {
                "attributes": [
                    'gene_id "ENSG00000273885.1"; gene_type "snoRNA"; '
                    'gene_name "ENSG10010138116.1"; level 3;',
                ]
            }
        ),
    )
