#!/usr/bin/env python

import argparse
from collections import defaultdict

import pandas as pd
import utils
from __init__ import ASSAY

MAX_CELL = 2 * 10**5
FEATURE_FILE_NAME = "features.tsv.gz"
BARCODE_FILE_NAME = "barcodes.tsv.gz"


class StarsoloSummary:
    def __init__(self, args):
        self.args = args
        self.cbs = utils.read_one_col(self.args.barcodes)

        self.stats = {}

    def parse_read_stats(self):
        dtypes = defaultdict(lambda: "int")
        dtypes["CB"] = "object"
        df = pd.read_csv(
            self.args.read_stats, sep="\t", header=0, index_col=0, skiprows=[1], dtype=dtypes
        )  # skip first line cb not pass whitelist
        rbs = list(df.index)
        umi_count = list(df["nUMIunique"])
        df = df.loc[
            :,
            [
                "cbMatch",
                "cbPerfect",
                "genomeU",
                "genomeM",
                "exonic",
                "intronic",
                "exonicAS",
                "intronicAS",
                "countedU",
                "nUMIunique",
                "nGenesUnique",
            ],
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
            "Corrected Barcodes": corrected / valid,
            "Reads Mapped To Unique Loci": genome_uniq / valid,
            "Reads Mapped To Multiple Loci": genome_multi / valid,
            "Reads Mapped Uniquely To Transcriptome": counted_uniq / valid,
            "Mapped Reads Assigned To Exonic Regions": exonic / mapped,
            "Mapped Reads Assigned To Intronic Regions": intronic / mapped,
            "Mapped Reads Assigned To Intergenic Regions": intergenic / mapped,
            "Mapped Reads Assigned Antisense To Gene": antisense / mapped,
        }
        for k in data_dict:
            data_dict[k] = utils.get_frac(data_dict[k])
        self.stats.update(data_dict)

        n_cells = len(self.cbs)
        reads_cell = df.loc[self.cbs, "countedU"].sum()
        fraction_reads_in_cells = utils.get_frac(float(reads_cell / counted_uniq))
        mean_used_reads_per_cell = int(reads_cell / n_cells)
        median_umi_per_cell = int(df.loc[self.cbs, "nUMIunique"].median())
        median_genes_per_cell = int(df.loc[self.cbs, "nGenesUnique"].median())
        data_dict = {
            "Estimated Number of Cells": n_cells,
            "Fraction of Reads in Cells": fraction_reads_in_cells,
            "Mean Used Reads per Cell": mean_used_reads_per_cell,
            "Median UMI per Cell": median_umi_per_cell,
            "Median Genes per Cell": median_genes_per_cell,
        }
        self.stats.update(data_dict)
        return rbs, umi_count

    def parse_summary(self):
        data = utils.csv2dict(self.args.summary)
        origin_new = {
            "Number of Reads": "Raw Reads",
            "Reads With Valid Barcodes": "Valid Reads",
            "Sequencing Saturation": "Saturation",
        }
        parsed_data = {}
        for origin, new in origin_new.items():
            parsed_data[new] = data[origin]
        frac_names = {"Valid Reads", "Saturation"}
        for k in frac_names:
            parsed_data[k] = utils.get_frac(parsed_data[k])
        for k in set(origin_new.values()) - frac_names:
            parsed_data[k] = int(parsed_data[k])
        self.stats.update(parsed_data)

    def run(self):
        self.cbs = utils.read_one_col(self.args.barcodes)
        rbs, umis = self.parse_read_stats()
        self.parse_summary()
        plot_data = utils.get_umi_count(rbs, umis, self.cbs, self.args.sample)
        utils.write_multiqc(plot_data, args.sample, ASSAY, "umi_count")
        utils.write_multiqc(self.stats, args.sample, ASSAY, "starsolo_summary.stats")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Starsolo summary")
    parser.add_argument("--read_stats", help="cellReadsStats file")
    parser.add_argument("--barcodes", help="barcode file")
    parser.add_argument("--summary", help="summary file")
    parser.add_argument("--sample", help="sample name")
    args = parser.parse_args()

    StarsoloSummary(args).run()
