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
