import clip
import os

def test_full_run(rootdir, monkeypatch, tmp_path):
    monkeypatch.setattr("sys.argv", ['clip.py',
        '-i', os.path.join(rootdir, 'tests', 'data', 'crosslinkcounts.bed'),
        '-a', os.path.join(rootdir, 'tests', 'data', 'annot.gff'),
        '-o', os.path.join(tmp_path, 'clippy'),
        '-t', str(os.cpu_count())])
    clip.main()
    for file_name in os.listdir(tmp_path):
        test_lines = open(os.path.join(tmp_path, file_name)).readlines()
        ref_lines = open(os.path.join(rootdir, 'tests', 'data', 'output', file_name)).readlines()
        assert(test_lines == ref_lines)
