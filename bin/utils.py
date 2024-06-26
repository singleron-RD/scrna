import csv
import gzip
import json
import logging
import sys
import time
from datetime import timedelta
from functools import wraps


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


def fastq_str(name, seq, qual):
    """return fastq read string"""
    return f"@{name}\n{seq}\n+\n{qual}\n"


def get_logger(name, level=logging.INFO):
    """out to stderr"""
    logger = logging.getLogger(name)
    logger.setLevel(level)
    log_formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setFormatter(log_formatter)
    logger.addHandler(console_handler)
    return logger


def write_json(data, fn):
    with open(fn, "w") as f:
        json.dump(data, f, indent=4)


def get_frac(raw_frac: float):
    return round(float(raw_frac) * 100, 2)


def csv2dict(csv_file):
    data = {}
    reader = csv.reader(openfile(csv_file))
    for row in reader:
        data[row[0]] = row[1]
    return data


def add_log(func):
    """
    logging start and done.
    """
    log_formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

    module = func.__module__
    name = func.__name__
    logger_name = f"{module}.{name}"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)

    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setFormatter(log_formatter)
    logger.addHandler(console_handler)

    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.info("start...")
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info("done. time used: %s", used)
        return result

    wrapper.logger = logger
    return wrapper


def write_multiqc(data, sample, assay, step):
    fn = f"{sample}.{assay}.{step}.json"
    write_json(data, fn)


MAX_CELL = 10**5


def get_umi_count(rbs, umis, cbs, sample):
    """
    Args:
        rbs: raw barcodes
        umis: umi count
        cbs: cell barcodes
    """
    a = [(umi, bc) for umi, bc in zip(umis, rbs) if umi > 0]
    a.sort(reverse=True)
    cbs = set(cbs)
    plot_data = {}
    n = len(a)
    first_noncell = n - 1
    for i, (umi, bc) in enumerate(a):
        if bc not in cbs:
            first_noncell = i
            break
    print(f"first non-cell barcode rank: {first_noncell}")
    last_cell = 0
    for i in range(min(n - 1, MAX_CELL), -1, -1):
        bc = a[i][1]
        if bc in cbs:
            last_cell = i
            break
    pure = sample + ".cells.pure" + f"({first_noncell}/{first_noncell}, 100%)"
    bg_cells = n - first_noncell
    bg = sample + ".cells.background" + f"(0/{bg_cells}, 0%)"
    plot_data[pure] = {}
    plot_data[bg] = {}
    for i in range(first_noncell):
        plot_data[pure][i + 1] = int(a[i][0])

    n_mix = last_cell - first_noncell + 1
    if n_mix != 0:
        n_total = len(cbs)
        n_mix_cell = n_total - first_noncell
        mix_rate = round(n_mix_cell / n_mix * 100, 2)
        mix = sample + ".cells.mix" + f"({n_mix_cell}/{n_mix}, {mix_rate}%)"
        plot_data[mix] = {}
        for i in range(first_noncell, last_cell + 1):
            plot_data[mix][i + 1] = int(a[i][0])

    for i in range(last_cell + 1, min(MAX_CELL, n), 10):
        plot_data[bg][i + 1] = int(a[i][0])
    # do not record every umi count
    for i in range(MAX_CELL, n, 1000):
        plot_data[bg][i + 1] = int(a[i][0])
    return plot_data
