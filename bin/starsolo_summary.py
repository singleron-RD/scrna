#!/usr/bin/env python

import argparse
import csv
import json
from collections import defaultdict

import pandas as pd
import utils


def parse_summary(f):
    parsed_data = {}
    reader = csv.reader(f)
    for row in reader:
        parsed_data[row[0]] = row[1]
    return parsed_data


MAX_CELL = 2 * 10**5

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Starsolo summary")
    parser.add_argument("--read_stats", help="cellReadsStats file")
    parser.add_argument("--barcodes", help="barcode file")
    parser.add_argument("--summary", help="summary file")
    parser.add_argument("--sample", help="sample name")
    args = parser.parse_args()
    dtypes = defaultdict(lambda: "int")
    dtypes["CB"] = "object"
    df = pd.read_csv(
        args.read_stats, sep="\t", header=0, index_col=0, skiprows=[1], dtype=dtypes
    )  # skip first line cb not pass whitelist
    umi_count = df["nUMIunique"]  # keep dataframe format
    df = df.loc[
        :, ["cbMatch", "cbPerfect", "genomeU", "genomeM", "exonic", "intronic", "exonicAS", "intronicAS", "countedU"]
    ]
    s = df.sum()
    # json does not recognize NumPy data types. TypeError: Object of type int64 is not JSON serializable
    valid = int(s["cbMatch"])
    perfect = int(s["cbPerfect"])
    corrected = valid - perfect
    genome_uniq = int(s["genomeU"])
    genome_multi = int(s["genomeM"])
    mapped = genome_uniq + genome_multi
    exonic = int(s["exonic"])
    intronic = int(s["intronic"])
    antisense = int(s["exonicAS"] + s["intronicAS"])
    intergenic = mapped - exonic - intronic - antisense
    counted_uniq = int(s["countedU"])
    data_dict = {
        "sample": args.sample,
        "valid": valid,
        "perfect": perfect,
        "corrected": corrected,
        "genomeU": genome_uniq,
        "genomeM": genome_multi,
        "mapped": mapped,
        "exonic": exonic,
        "intronic": intronic,
        "antisense": antisense,
        "intergenic": intergenic,
        "countedU": counted_uniq,
    }
    read_stats_file = args.sample + ".read_stats.json"
    with open(read_stats_file, "w") as f:
        json.dump(data_dict, f)

    # summary
    parsed_data = parse_summary(open(args.summary))
    summary_file = args.sample + ".summary.json"
    with open(summary_file, "w") as f:
        json.dump(parsed_data, f)

    # UMI count
    umi_count.loc[lambda x: x > 0]
    umi_count = umi_count.sort_values(ascending=False)
    cbs = set(utils.read_one_col(args.barcodes))
    plot_data = {}
    n = len(umi_count)
    first_noncell = n - 1
    for i, bc in enumerate(umi_count.index):
        if bc not in cbs:
            first_noncell = i
            break
    last_cell = 0
    for i in range(min(n - 1, MAX_CELL), -1, -1):
        bc = umi_count.index[i]
        if bc in cbs:
            last_cell = i
            break
    pure = args.sample + ".cells.pure" + f"({first_noncell}/{first_noncell}, 100%)"
    bg = args.sample + ".cells.background"
    plot_data[pure] = {}
    plot_data[bg] = {}
    for i in range(first_noncell):
        plot_data[pure][i + 1] = int(umi_count.iloc[i])

    n_mix = last_cell - first_noncell + 1
    if n_mix != 0:
        n_total = len(cbs)
        n_mix_cell = n_total - first_noncell
        mix_rate = round(n_mix_cell / n_mix * 100, 2)
        mix = args.sample + ".cells.mix" + f"({n_mix_cell}/{n_mix}, {mix_rate}%)"
        plot_data[mix] = {}
        for i in range(first_noncell, last_cell + 1):
            plot_data[mix][i + 1] = int(umi_count.iloc[i])

    for i in range(last_cell + 1, min(MAX_CELL, n), 10):
        plot_data[bg][i + 1] = int(umi_count.iloc[i])
    # do not record every umi count
    for i in range(MAX_CELL, n, 1000):
        plot_data[bg][i + 1] = int(umi_count.iloc[i])

    umi_file = args.sample + ".umi_count.json"
    with open(umi_file, "w") as f:
        json.dump(plot_data, f)
