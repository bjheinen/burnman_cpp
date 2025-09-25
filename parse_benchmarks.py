#!/usr/bin/env python3
# ----------------- START OF LICENSE SECTION -----------------
# Copyright (c) 2025 Benedict Heinen
#
# This file is part of burnman_cpp and is licensed under the
# GNU General Public License v3.0 or later. See the LICENSE file
# or <https://www.gnu.org/licenses/> for details.
#
# burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
#
# ------------------ END OF LICENSE SECTION ------------------
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import pandas as pd
import re
import os

def unique_filename(path: str) -> str:
    base, ext = os.path.splitext(path)
    counter = 2
    new_path = path
    while os.path.exists(new_path):
        new_path = f"{base}_{counter}{ext}"
        counter += 1
    return new_path

def bar_plot(df, show=0, save_path=None, xlabel='Time (ns)'):
    fig, ax = plt.subplots(figsize=(16, 8))
    colors = plt.cm.tab20(np.linspace(0, 1, df.shape[0]))
    df.plot(kind='barh', ax=ax, color=colors)
    ax.yaxis.label.set_visible(False)
    ax.set_xlabel(xlabel)
    ax.set_title("Mineral Property Benchmarks")
    ax.legend(title="Property", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    if save_path:
        save_path = unique_filename(save_path)
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
    if show:
        plt.show(block=False)

def format_line(label, value, box_width=56, label_width=18, val_width=4):
    # | + label + " : " + value + fill + |
    value_str = str(value)
    label_width = max(len(label), label_width)
    val_width = max(len(value_str), val_width)
    spaces = box_width - label_width - val_width - 7
    return f"| {label:<{label_width}} : {value_str:>{val_width}}{' ' * spaces} |\n"

def format_title_line(label=None, box_width=56):
    if not label:
        return f"+{'-' * (box_width - 2)}+"
    else:
        fill = box_width - len(label) - 4
        left = fill // 2
        right = fill - left
        return f"\n +{'-' * left} {label} {'-' * right}+\n"

def parse_mineral_benchmarks(case, verbose=1):
    data = []
    benchmarks = case.findall("BenchmarkResults")
    print(
        format_title_line("Processing Test Case"),
        format_line("Benchmark Suite", case.get('name')),
        format_line("Benchmarks", len(benchmarks)),
        format_title_line()
    )
    # Loop through each benchmark
    for bench in case.findall("BenchmarkResults"):
        # Get name - e.g. 'get_density [MGD3]'
        name_field = bench.get("name")
        # Strip benchmark name and EOS tag for grouping
        match = re.match(r'(.+?)\s*\[([^\]]+)\](.*)', name_field)
        if match:
            bench_name = match.group(1).strip()
            eos_tag = match.group(2).strip()
            extra = match.group(3).strip()
            if verbose:
                print(f"Processing benchmark '{bench_name}' for EOS '{eos_tag}'...")
            if extra:
                print(f"    Warning: Stripping extra text after EOS tag -> '{extra}'")
        else:
            print(f"Warning: Skipping benchmark -> '{name_field}'")
            continue
        # Get xml elements for mean / stddev
        mean_elem = bench.find("mean")
        std_elem = bench.find("standardDeviation")
        # Extract times (ns)
        mean_ns = float(mean_elem.get("value")) if mean_elem is not None else None
        std_ns = float(std_elem.get("value")) if std_elem is not None else None
        data.append(
            {
            "Benchmark" : bench_name,
            "EOS": eos_tag,
            "Mean": mean_ns,
            "Std_dev": std_ns
            }
        )
    # Convert to pandas
    df = pd.DataFrame(data)
    # Group and average multi entries
    df_grouped = df.groupby(["Benchmark", "EOS"], as_index=False).mean()
    # Pivot: Benchmarks as rows, EOS as columns
    df_mean = df_grouped.pivot(index="Benchmark", columns="EOS", values="Mean").reset_index()
    df_std = df_grouped.pivot(index="Benchmark", columns="EOS", values="Std_dev").reset_index()
    print(
        format_title_line("Finished Processing Suite"),
        format_line("Benchmark Suite", case.get('name')),
        format_line("Benchmarks", df_mean.shape[0]),
        format_line("EOS Variants", df_mean.shape[1]-1),
        format_title_line()
    )
    # Merge mean and std_dev data frames
    df_flat = pd.concat([df_mean, df_std.iloc[:, 1:].add_suffix('_std')], axis=1)
    return df_flat

def plot_mineral_benchmarks(data, show_plots=1, save_plots=0, out_file_ext=''):
    # Ignore std_devs for now
    df = data[data.columns[~data.columns.str.endswith('_std')]]
    # Set 'Benchmark' column to the index
    df.set_index('Benchmark', inplace=True)
    fn_a = None
    fn_b = None
    if save_plots:
        fn_a = 'mineral_benchmarks_per_eos_' + out_file_ext + '.pdf'
        fn_b = 'mineral_benchmarks_per_property_' + out_file_ext + '.pdf'
    bar_plot(df.T, show=show_plots, save_path=fn_a)
    bar_plot(df, show=show_plots, save_path=fn_b)

def plot_mineral_baseline_comparison(data, show_plots=1, save_plots=0, out_file_ext=''):
    fn_a = None
    fn_b = None
    fn_c = None
    fn_d = None
    if save_plots:
        fn_a = 'mineral_benchmarks_pct_change_per_eos_' + out_file_ext + '.pdf'
        fn_b = 'mineral_benchmarks_pct_change_per_property_' + out_file_ext + '.pdf'
        fn_c = 'mineral_benchmarks_pct_change_per_property_averages' + out_file_ext + '.pdf'
        fn_d = 'mineral_benchmarks_pct_change_per_eos_averages' + out_file_ext + '.pdf'
    # EOS columns (exclude Benchmark + averages)
    eos_cols = [c for c in data.columns if c not in ["Benchmark", "average_pct_change_property"]]
    raw_data = data[data["Benchmark"] != "average_pct_change_EOS"]
    col_avg = data[data["Benchmark"] == "average_pct_change_EOS"]
    col_avg = col_avg[eos_cols].T
    col_avg.index.name = "EOS"
    col_avg.columns = ["Average % Change"]
    row_avg = raw_data.copy().set_index("Benchmark")
    row_avg = row_avg[["average_pct_change_property"]]
    raw_data.set_index("Benchmark", inplace=True)
    raw_data = raw_data[eos_cols]
    bar_plot(raw_data, show=show_plots, save_path=fn_a, xlabel='% Change')
    bar_plot(raw_data.T, show=show_plots, save_path=fn_b, xlabel='% Change')
    bar_plot(row_avg, show=show_plots, save_path=fn_c, xlabel='Average % Change')
    bar_plot(col_avg, show=show_plots, save_path=fn_d, xlabel='Average % Change')

def compare_mineral_benchmarks(data, baseline_fn='mineral_benchmarks_baseline.csv'):
    # Load baseline benchmark
    baseline = pd.read_csv(baseline_fn)
    # Get eos cols
    eos_cols = [c for c in baseline.columns if c != "Benchmark" and not c.endswith("_std")]
    # Merge on Benchmark name, ignore std_dev
    merged = baseline[["Benchmark"] + eos_cols].merge(
        data[["Benchmark"] + eos_cols],
        on="Benchmark",
        suffixes=("_baseline", "_new")
    )
    # Calculate percentage change
    pct_changes = {"Benchmark": merged["Benchmark"]}
    for col in eos_cols:
        pct_changes[col] = 100 * (merged[col + "_new"] - merged[col + "_baseline"]) / merged[col + "_baseline"]
    # Convert back to data frame
    df_pct = pd.DataFrame(pct_changes)
    df_pct["average_pct_change_property"] = df_pct[eos_cols].mean(axis=1)
    avg_row = ["average_pct_change_EOS"] + list(df_pct[eos_cols].mean(axis=0)) + [df_pct["average_pct_change_property"].mean()]
    df_pct.loc[len(df_pct)] = avg_row
    max_change = df_pct[eos_cols].values.max()
    min_change = df_pct[eos_cols].values.min()
    avg_change = df_pct.loc[df_pct["Benchmark"] == "average_pct_change_EOS", "average_pct_change_property"].values[0]
    print(
        format_title_line("Baseline Comparison"),
        format_line("Baseline data", baseline_fn),
        format_line("Max change", (str(max_change)+'%')),
        format_line("Min change", (str(min_change)+'%')),
        format_line("Avg change", (str(avg_change)+'%')),
        format_title_line()
    )
    return df_pct

def parse_catch2_benchmark_xml(filename, out_file_ext, save_data=1, plot_data=1, show_plots=1, save_plots=1):
    tree = ET.parse(filename)
    root = tree.getroot()
    # Get Test Cases
    test_cases = root.findall("TestCase")
    # Parse OverallResultsCases and print info
    success_info = root.find("OverallResultsCases")
    print(
        format_title_line("Benchmark Test Cases"),
        format_line("Successes", success_info.get('successes')),
        format_line("Failures", success_info.get('failures')),
        format_line("Expected Failures", success_info.get('expectedFailures')),
        format_line("Skips", success_info.get('skips')),
        format_line("Test Cases", len(test_cases)),
        format_line("Total Benchmarks", sum([len(case.findall('BenchmarkResults')) for case in test_cases])),
        format_title_line()
    )
    # Loop through TestCases and process accordingly
    for case in test_cases:
        # Parse mineral benchmarks (property for each EOS)
        if case.get("name") == 'Mineral per-property benchmarks':
            eos_data = parse_mineral_benchmarks(case)
            compare_data = compare_mineral_benchmarks(eos_data)
            if save_data:
                min_fname = unique_filename("mineral_benchmarks_" + out_file_ext + ".csv")
                comp_fname = unique_filename('mineral_benchmarks_compare' + out_file_ext + '.csv')
                eos_data.to_csv(min_fname, index=False)
                compare_data.to_csv(comp_fname, index=False)
                print(f"Data saved: {min_fname} | {comp_fname}")
            if plot_data:
                plot_mineral_benchmarks(eos_data, show_plots=show_plots, save_plots=save_plots, out_file_ext=out_file_ext)
                plot_mineral_baseline_comparison(compare_data, show_plots=show_plots, save_plots=save_plots, out_file_ext=out_file_ext)
        else:
            print("Parsing for benchmark suite type not implemented!")

if __name__ == "__main__":
    import sys
    # Get XML file and output fname to append to pass as args
    n_args = len(sys.argv)
    xml_file = sys.argv[1]
    if n_args < 3:
        save = 0
        out_ext = None
    else:
        save = 1
        out_ext = sys.argv[2]
    parse_catch2_benchmark_xml(xml_file, out_ext, save_data=save, save_plots=save)
    input("Press <Enter> to close plots and exit...")
