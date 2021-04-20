import clip
import pytest

def test_no_argument():
    with pytest.raises(SystemExit) as e:
        clip.parse_arguments([])
    assert(e.type == SystemExit)
    assert(e.value.code == 2)

def test_minimal_argument():
    args = clip.parse_arguments(['-i', 'test.bed', '-a', 'test.gtf', '-o', 'test'])
    assert(isinstance(args, tuple))
    assert(args == ('test.bed', 'test.gtf', 50, 1, 5,
        'test_rollmean50_stdev1_minGeneCount5.bed', None, 10, 5, False))

def test_optional_argument():
    args = clip.parse_arguments(['-i', 'test.bed', '-a', 'test.gtf', '-o', 'test', '-n', '100', '-int'])
    assert(isinstance(args, tuple))
    assert(args == ('test.bed', 'test.gtf', 100, 1, 5,
        'test_rollmean100_stdev1_minGeneCount5.bed', None, 10, 5, True))
