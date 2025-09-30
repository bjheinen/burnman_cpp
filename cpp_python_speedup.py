import pandas as pd

def compare_simple_benchmarks(csv_cpp, csv_py):
    df_cpp = pd.read_csv(csv_cpp)
    df_py = pd.read_csv(csv_py)
    # Merge on Benchmark name
    df = pd.merge(df_cpp, df_py, on="Benchmark", suffixes=("_cpp", "_py"))
    # Speedup = how many times C++ is faster than Python
    df["Speedup"] = df["Mean_py"] / df["Mean_cpp"]
    # Optional: sort by speedup descending
    df = df.sort_values("Speedup", ascending=False)
    return df

def compare_multi_benchmarks(csv_cpp, csv_py):
    df_cpp = pd.read_csv(csv_cpp)
    df_py = pd.read_csv(csv_py)
    eos_cols = [c for c in df_cpp.columns if not c.endswith("_std") and c != "Benchmark"]
    df_multi = pd.merge(df_cpp, df_py, on="Benchmark", suffixes=("_cpp", "_py"))
    for eos in eos_cols:
        df_multi[f"{eos}_speedup"] = df_multi[f"{eos}_py"] / df_multi[f"{eos}_cpp"]
    speedup_cols = [f"{eos}_speedup" for eos in eos_cols]
    avg_speedup_per_eos = df_multi[speedup_cols].mean(axis=0)
    df_multi["avg_speedup"] = df_multi[speedup_cols].mean(axis=1)
    df_multi = df_multi.sort_values("avg_speedup", ascending=False)
    return df_multi

def compute_cpp_speedup(save_data=0):
    mineral_benchmarks = compare_multi_benchmarks('benchmarks/mineral_benchmarks_baseline.csv', 'python_mineral_benchmarks.csv')
    assemblage_benchmarks = compare_simple_benchmarks('benchmarks/assemblage_benchmarks_baseline.csv', 'python_assemblage_benchmarks.csv')
    if save_data:
        mineral_benchmarks.to_csv('mineral_benchmark_speedups.csv')
        assemblage_benchmarks.to_csv('assemblage_benchmark_speedups.csv')
    avg_mineral_speedup = mineral_benchmarks["avg_speedup"].mean()
    avg_assemblage_speedup = assemblage_benchmarks["Speedup"].mean()
    avg_speedup = (avg_mineral_speedup + avg_assemblage_speedup) / 2.0
    print(
        f" +-------------------------------+\n",
        f"  Mineral Speed-up    : {avg_mineral_speedup:.1f}x\n",
        f"  Assemblage Speed-up : {avg_assemblage_speedup:.1f}x\n",
        f"  Average Speed-up    : {avg_speedup:.1f}x\n",
        f"+-------------------------------+\n"
    )

if __name__ == "__main__":
    compute_cpp_speedup()
