def find_diff_entries(file_1_path, file_2_path, output_1_path, output_2_path):
    file_1_lines = set(open(file_1_path).readlines())
    file_2_lines = set(open(file_2_path).readlines())
    with open(output_1_path, "w") as out_f:
        for line in file_1_lines.difference(file_2_lines):
            out_f.write(line)
    with open(output_2_path, "w") as out_f:
        for line in file_2_lines.difference(file_1_lines):
            out_f.write(line)


find_diff_entries(
    "tests/data/output/clippy_rollmean50_stdev1_minGeneCount5.bed",
    "TEST_rollmean50_stdev1_minGeneCount5.bed",
    "original_single.bed",
    "new_single.bed",
)

find_diff_entries(
    "tests/data/output/clippy_rollmean50_stdev1_minGeneCount5_broadPeaks.bed",
    "TEST_rollmean50_stdev1_minGeneCount5_broadPeaks.bed",
    "original_broad.bed",
    "new_broad.bed",
)
