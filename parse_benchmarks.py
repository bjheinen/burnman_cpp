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

def bar_plot(df, show=0, save_path=None, xlabel='Time (ns)', title='Mineral Property Benchmarks', plot_legend=1):
    fig, ax = plt.subplots(figsize=(16, 8))
    colors = plt.cm.tab20(np.linspace(0, 1, df.shape[0]))
    df.plot(kind='barh', ax=ax, color=colors)
    ax.yaxis.label.set_visible(False)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    if plot_legend:
        ax.legend(title="Property", bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        ax.get_legend().remove()
    plt.tight_layout()
    if save_path:
        save_path = unique_filename(save_path)
        plt.savefig(save_path, bbox_inches='tight', dpi=300)
    if show:
        plt.show(block=False)

def format_line(label, value, note='', box_width=56, label_width=18, val_width=4):
    # | + label + " : " + value + fill + |
    value_str = str(value)
    label_width = max(len(label), label_width)
    val_width = max(len(value_str), val_width)
    note = (" " + note) if note else note
    note_width = len(note)
    spaces = box_width - label_width - val_width - note_width - 7
    return f"| {label:<{label_width}} : {value_str:>{val_width}}{note}{' ' * spaces} |\n"

def format_title_line(label=None, box_width=56):
    if not label:
        return f"+{'-' * (box_width - 2)}+"
    else:
        fill = box_width - len(label) - 4
        left = fill // 2
        right = fill - left
        return f"\n +{'-' * left} {label} {'-' * right}+\n"

def parse_mineral_benchmarks(case, verbose=0):
    data = []
    high_std_count = 0
    benchmarks = case.findall("BenchmarkResults")
    n_benchmarks = len(benchmarks)
    print(
        format_title_line("Processing Test Case"),
        format_line("Benchmark Suite", case.get('name')),
        format_line("Benchmarks", n_benchmarks),
        format_title_line()
    )
    # Loop through each benchmark
    for bench in benchmarks:
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

        std_pct = std_ns / mean_ns * 100.0
        threshold = 5.0
        if std_pct > threshold:
            high_std_count += 1
            print(
                f"    Warning: High variability in benchmark '{bench_name} [{eos_tag}]'\n"
                f"             (std_dev = {std_pct:.1f}%)"
            )
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
    uncertain_pct = high_std_count / n_benchmarks * 100.0
    print(
        format_title_line("Finished Processing Suite"),
        format_line("Benchmark Suite", case.get('name')),
        format_line("Benchmarks", df_mean.shape[0]),
        format_line("EOS Variants", df_mean.shape[1]-1),
        format_line("Total benchmarks", len(benchmarks)),
        format_line("Uncertain", high_std_count, note=f'({uncertain_pct:.1f}%)'),
        format_title_line()
    )
    # Merge mean and std_dev data frames
    df_flat = pd.concat([df_mean, df_std.iloc[:, 1:].add_suffix('_std')], axis=1)
    return df_flat

def parse_prop_mod_benchmarks(case, verbose=0):
    data = []
    high_std_count = 0
    benchmarks = case.findall("BenchmarkResults")
    n_benchmarks = len(benchmarks)
    print(
        format_title_line("Processing Test Case"),
        format_line("Benchmark Suite", case.get('name')),
        format_line("Benchmarks", n_benchmarks),
        format_title_line()
    )
    # Loop through each benchmark
    for bench in benchmarks:
        # Get name
        name_field = bench.get("name")
        bench_name = name_field.strip()
        if verbose:
            print(f"Processing benchmark '{bench_name}'...")
        # Get xml elements for mean / stddev
        mean_elem = bench.find("mean")
        std_elem = bench.find("standardDeviation")
        # Extract times (ns)
        mean_ns = float(mean_elem.get("value")) if mean_elem is not None else None
        std_ns = float(std_elem.get("value")) if std_elem is not None else None
        std_pct = std_ns / mean_ns * 100.0
        threshold = 5.0
        if std_pct > threshold:
            high_std_count += 1
            print(
                f"    Warning: High variability in benchmark '{bench_name}'\n"
                f"             (std_dev = {std_pct:.1f}%)"
            )
        data.append(
            {
            "Benchmark" : bench_name,
            "Mean": mean_ns,
            "Std_dev": std_ns
            }
        )
    # Convert to pandas
    df = pd.DataFrame(data)
    uncertain_pct = high_std_count / n_benchmarks * 100.0
    print(
        format_title_line("Finished Processing Suite"),
        format_line("Benchmark Suite", case.get("name")),
        format_line("Total Benchmarks", df.shape[0]),
        format_line("Uncertain", high_std_count, note=f'({uncertain_pct:.1f}%)'),
        format_title_line()
    )
    return df

def compare_prop_mod_benchmarks(data, baseline_fn='property_modifier_benchmarks_baseline.csv'):
    if not os.path.exists(baseline_fn):
        print(
            format_title_line("Baseline Comparison"),
            format_line("Error!", "File not found"),
            format_line("Baseline file", baseline_fn),
            format_title_line() + "\n"
        )
        return
    # Load baseline benchmark
    baseline = pd.read_csv(baseline_fn)
    merged = baseline.merge(data, on="Benchmark", suffixes=("_baseline", "_new"))
    merged["Pct Change from Baseline"] = 100 * (merged["Mean_new"] - merged["Mean_baseline"]) / merged["Mean_baseline"]
    comparison = merged[["Benchmark", "Pct Change from Baseline"]].copy()
    comparison.loc[len(comparison)] = ["Average", comparison["Pct Change from Baseline"].mean()]
    mask = comparison["Benchmark"] != "Average"
    changes = comparison.loc[mask, "Pct Change from Baseline"]
    min_change = changes.min()
    max_change = changes.max()
    avg_change = changes.mean()
    print(
        format_title_line("Baseline Comparison"),
        format_line("Baseline data", baseline_fn),
        format_line("Max change", f"{max_change:.1f}%", val_width=6),
        format_line("Min change", f"{min_change:.1f}%", val_width=6),
        format_line("Avg change", f"{avg_change:.1f}%", val_width=6),
        format_title_line()
    )
    return comparison

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
    bar_plot(row_avg, show=show_plots, save_path=fn_c, xlabel='Average % Change', plot_legend=0)
    bar_plot(col_avg, show=show_plots, save_path=fn_d, xlabel='Average % Change', plot_legend=0)

def plot_prop_mod_benchmarks(data, show_plots=1, save_plots=0, out_file_ext=''):
    df = data[["Benchmark", "Mean"]].copy()
    df.set_index("Benchmark", inplace=True)
    fn = None
    if save_plots:
        fn = 'property_modifier_benchmarks_' + out_file_ext + '.pdf'
    bar_plot(df, show=show_plots, save_path=fn, title='Property Modifier Benchmarks', plot_legend=0)

def plot_prop_mod_baseline_comparison(data, show_plots=1, save_plots=0, out_file_ext=''):
    fn = None
    if save_plots:
        fn = 'property_modifier_benchmarks_pct_change' + out_file_ext + '.pdf'
    mask = data["Benchmark"] != "Average"
    plot_data = data[mask].copy()
    plot_data.set_index("Benchmark", inplace=True)
    bar_plot(plot_data, show=show_plots, save_path=fn, xlabel='% Change', title='Property Modifier Benchmarks', plot_legend=0)

def compare_mineral_benchmarks(data, baseline_fn='mineral_benchmarks_baseline.csv'):
    if not os.path.exists(baseline_fn):
        print(
            format_title_line("Baseline Comparison"),
            format_line("Error!", "File not found"),
            format_line("Baseline file", baseline_fn),
            format_title_line() + "\n"
        )
        return
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
        format_line("Max change", f"{max_change:.1f}%", val_width=6),
        format_line("Min change", f"{min_change:.1f}%", val_width=6),
        format_line("Avg change", f"{avg_change:.1f}%", val_width=6),
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
                comp_fname = unique_filename('mineral_benchmarks_compare_' + out_file_ext + '.csv')
                eos_data.to_csv(min_fname, index=False)
                msg = (
                    format_title_line("Data Saved!")
                    + ' ' + format_line("Data file", min_fname)
                )
                if compare_data is not None:
                    compare_data.to_csv(comp_fname, index=False)
                    msg += ' ' + format_line("Comparison", comp_fname)
                msg += ' ' + format_title_line()
                print(msg)
            if plot_data:
                plot_mineral_benchmarks(eos_data, show_plots=show_plots, save_plots=save_plots, out_file_ext=out_file_ext)
                if compare_data is not None:
                    plot_mineral_baseline_comparison(compare_data, show_plots=show_plots, save_plots=save_plots, out_file_ext=out_file_ext)
        elif case.get("name") == "Property modifier benchmarks":
            prop_mod_data = parse_prop_mod_benchmarks(case)
            prop_mod_compare_data = compare_prop_mod_benchmarks(prop_mod_data)
            # save check and save
            if save_data:
                prop_mod_fname = unique_filename("property_modifier_benchmarks_" + out_file_ext + ".csv")
                comp_fname = unique_filename("property_modifier_benchmarks_compare_" + out_file_ext + ".csv")
                prop_mod_data.to_csv(prop_mod_fname, index=False)
                msg = (
                    format_title_line("Data Saved!")
                    + ' ' + format_line("Data file", prop_mod_fname)
                )
                if prop_mod_compare_data is not None:
                    prop_mod_compare_data.to_csv(comp_fname, index=False)
                    msg += ' ' + format_line("Comparison", comp_fname)
                msg += ' ' + format_title_line()
                print(msg)
            if plot_data:
                plot_prop_mod_benchmarks(prop_mod_data, show_plots=show_plots, save_plots=save_plots, out_file_ext=out_file_ext)
                if prop_mod_compare_data is not None:
                    plot_prop_mod_baseline_comparison(prop_mod_compare_data, show_plots=show_plots, save_plots=save_plots, out_file_ext=out_file_ext)

        else:
            print(
                format_title_line("Warning!"),
                "    Parsing for benchmark suite type not implemented!\n",
                format_title_line()
            )

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
