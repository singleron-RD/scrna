#!/usr/bin/env python

import argparse
import gzip
import json
import random
import statistics
from collections import defaultdict

import pysam


def openfile(file_name, mode="rt", **kwargs):
    """open gzip or plain file"""
    if file_name.endswith(".gz"):
        file_obj = gzip.open(file_name, mode=mode, **kwargs)
    else:
        file_obj = open(file_name, mode=mode, **kwargs)
    return file_obj


def read_one_col(fn):
    """read one column file into list"""
    with openfile(fn) as f:
        return [x.strip() for x in f]


def get_records(bam_file):
    a = []
    cb_int = {}
    ub_int = {}
    gx_int = {}
    name_int = {}
    n_cb = n_ub = n_gx = n_read = 0
    with pysam.AlignmentFile(bam_file) as bam:
        for record in bam:
            if record.query_name in name_int:  # already added
                continue
            cb = record.get_tag("CB")
            ub = record.get_tag("UB")
            gx = record.get_tag("GX")

            if all(x != "-" for x in (cb, ub, gx)):
                # use int instead of str to avoid memory hog
                if cb not in cb_int:
                    n_cb += 1
                    cb_int[cb] = n_cb
                if ub not in ub_int:
                    n_ub += 1
                    ub_int[ub] = n_ub
                if gx not in gx_int:
                    n_gx += 1
                    gx_int[gx] = n_gx
                a.append((cb_int[cb], ub_int[ub], gx_int[gx]))
                n_read += 1
                name_int[record.query_name] = n_read
    return a, cb_int


def sub_saturation(a):
    """get saturation and median gene"""
    n = len(a)
    fraction_saturation = {0.0: 0.0}
    for fraction in range(1, 11):
        fraction /= 10.0
        nread = int(n * fraction)
        uniq = len(set(a[:nread]))
        saturation = 1 - float(uniq) / nread
        saturation = round(saturation * 100, 2)
        fraction_saturation[fraction] = saturation
    return fraction_saturation


def sub_gene(a, barcodes):
    """get median gene for each fraction"""
    nread_fraction = {}
    n = len(a)
    for fraction in range(11):
        fraction /= 10.0
        nread = int(n * fraction)
        nread_fraction[nread] = fraction

    fraction_mg = {0.0: 0}
    cb_gx = defaultdict(set)
    for i, (cb, _, gx) in enumerate(a, start=1):
        if cb in barcodes:
            cb_gx[cb].add(gx)
        if i in nread_fraction:
            fraction = nread_fraction[i]
            fraction_mg[fraction] = int(statistics.median([len(x) for x in cb_gx.values()]))
    return fraction_mg


def main(args):
    """main function"""
    a, cb_dict = get_records(args.bam)
    barcodes = read_one_col(args.cell_barcode)
    barcodes = set(cb_dict[x] for x in barcodes)
    random.seed(0)
    random.shuffle(a)

    fraction_saturation = sub_saturation(a)
    fraction_mg = sub_gene(a, barcodes)
    saturation_file = f"{args.sample}.scrna.saturation.json"
    median_gene_file = f"{args.sample}.scrna.median_gene.json"
    # write json
    with open(saturation_file, "w") as f:
        f.write(json.dumps(fraction_saturation))
    with open(median_gene_file, "w") as f:
        f.write(json.dumps(fraction_mg))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="saturation")
    parser.add_argument("-b", "--bam", help="bam file", required=True)
    parser.add_argument("-c", "--cell_barcode", help="barcode file", required=True)
    parser.add_argument("-s", "--sample", help="sample name", required=True)
    args = parser.parse_args()
    main(args)
