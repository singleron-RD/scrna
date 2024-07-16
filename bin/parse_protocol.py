import itertools
import json
import os
import re
import sys
from collections import defaultdict

import pyfastx
import utils

logger = utils.get_logger(__name__)


def get_seq_str(seq, sub_pattern):
    """
    join seq slices.

    Args:
        seq: usually R1 read
        sub_pattern: [slice(0,8),slice(16,24)]

    Returns:
        joined intervals seq

    >>> sub_pattern_dict = [slice(0,2)]
    >>> seq = "A" * 2 + "T" * 2
    >>> get_seq_str(seq, sub_pattern_dict)
    'AA'
    """
    return "".join([seq[x] for x in sub_pattern])


def findall_mismatch(seq, n_mismatch=1, bases="ACGTN"):
    """
    choose locations where there's going to be a mismatch using combinations
    and then construct all satisfying lists using product

    Return:
    all mismatch <= n_mismatch set.

    >>> answer = set(["TCG", "AAG", "ACC", "ATG", "ACT", "ACN", "GCG", "ANG", "ACA", "ACG", "CCG", "AGG", "NCG"])
    >>> seq_set = findall_mismatch("ACG")
    >>> seq_set == answer
    True
    """
    seq_set = set()
    seq_len = len(seq)
    if n_mismatch > seq_len:
        n_mismatch = seq_len
    for locs in itertools.combinations(range(seq_len), n_mismatch):
        seq_locs = [[base] for base in seq]
        for loc in locs:
            seq_locs[loc] = list(bases)
        for poss in itertools.product(*seq_locs):
            seq_set.add("".join(poss))
    return seq_set


def get_mismatch_dict(seq_list, n_mismatch=1):
    """
    Return:
    mismatch dict. Key: mismatch seq, value: seq in seq_list

    >>> seq_list = ["AACGTGAT", "AAACATCG"]
    >>> mismatch_dict = get_mismatch_dict(seq_list)
    >>> mismatch_dict["AACGTGAA"] == "AACGTGAT"
    True
    """
    mismatch_dict = {}
    for seq in seq_list:
        seq = seq.strip()
        if seq == "":
            continue
        for mismatch_seq in findall_mismatch(seq, n_mismatch):
            mismatch_dict[mismatch_seq] = seq
    return mismatch_dict


def parse_pattern(pattern, allowed="CLUNT"):
    """
    >>> pattern_dict = parse_pattern("C8L16C8L16C8L1U12T18")
    >>> pattern_dict['C']
    [slice(0, 8, None), slice(24, 32, None), slice(48, 56, None)]
    >>> pattern_dict['L']
    [slice(8, 24, None), slice(32, 48, None), slice(56, 57, None)]
    """
    pattern_dict = {}
    p = re.compile(r"([A-Z])(\d+)")
    tmp = p.findall(pattern)
    if not tmp:
        sys.exit(f"Invalid pattern: {pattern}")
    start = 0
    for x, length in tmp:
        if x not in allowed:
            sys.exit(f"Invalid pattern: {pattern}")
        if x not in pattern_dict:
            pattern_dict[x] = []
        end = start + int(length)
        pattern_dict[x].append(slice(start, end))
        start = end
    return pattern_dict


def get_raw_mismatch(files: list, n_mismatch: int):
    """
    Args:
        files: whitelist file paths
        n_mismatch: allowed number of mismatch bases
    Returns:
        raw_list
        mismatch_list
    """
    raw_list, mismatch_list = [], []
    for f in files:
        barcodes = utils.read_one_col(f)
        raw_list.append(set(barcodes))
        barcode_mismatch_dict = get_mismatch_dict(barcodes, n_mismatch)
        mismatch_list.append(barcode_mismatch_dict)

    return raw_list, mismatch_list


def check_seq_mismatch(seq_list, raw_list, mismatch_list):
    """
    Returns
        valid: True if seq in mismatch_list
        corrected: True if seq in mismatch_list but not in raw_list
        res: joined seq

    >>> seq_list = ['ATA', 'AAT', 'ATA']
    >>> correct_set_list = [{'AAA'},{'AAA'},{'AAA'}]
    >>> mismatch_dict_list = [get_mismatch_dict(['AAA'])] * 3

    >>> check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list)
    (True, True, 'AAA_AAA_AAA')

    >>> seq_list = ['AAA', 'AAA', 'AAA']
    >>> check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list)
    (True, False, 'AAA_AAA_AAA')
    """
    valid = True
    corrected = False
    res = []
    for index, seq in enumerate(seq_list):
        if seq not in raw_list[index]:
            if seq not in mismatch_list[index]:
                valid = False
                res = []
            else:
                corrected = True
                res.append(mismatch_list[index][seq])
        else:
            res.append(seq)

    return valid, corrected, "_".join(res)


def get_protocol_dict(assets_dir):
    """
    Return:
    protocol_dict. Key: protocol name, value: protocol dict

    >>> protocol_dict = get_protocol_dict("./assets/")
    >>> protocol_dict["GEXSCOPE-MicroBead"]["pattern_dict"]
    {'C': [slice(0, 12, None)], 'U': [slice(12, 20, None)]}
    """
    json_file = os.path.join(assets_dir, "protocols.json")
    protocol_dict = json.load(open(json_file))
    whitelist_dir = os.path.join(assets_dir, "whitelist")
    # add folder prefix
    for protocol in protocol_dict:
        cur = protocol_dict[protocol]
        bc = cur.get("bc", [])
        linker = cur.get("linker", [])
        if bc:
            cur["bc"] = [os.path.join(whitelist_dir, protocol, x) for x in bc]
        if linker:
            cur["linker"] = [os.path.join(whitelist_dir, protocol, x) for x in linker]
        cur["pattern_dict"] = parse_pattern(cur["pattern"])
    return protocol_dict


class Auto:
    """
    Auto detect singleron protocols from R1-read
    GEXSCOPE-MicroBead
    GEXSCOPE-V1
    GEXSCOPE-V2
    GEXSCOPE-V3
    """

    def __init__(self, fq1_list, sample, assets_dir="assets/", max_read=10000):
        """
        Args:
            assets_dir: Expects file 'protocols.json' and 'whitelist/{protocol}' folder under assets_dir

        Returns:
            protocol, protocol_dict[protocol]
        """
        self.fq1_list = fq1_list
        self.max_read = max_read
        self.sample = sample
        self.protocol_dict = get_protocol_dict(assets_dir)
        self.mismatch_dict = {}
        for protocol in self.protocol_dict:
            if "bc" in self.protocol_dict[protocol]:
                self.mismatch_dict[protocol] = get_raw_mismatch(self.protocol_dict[protocol]["bc"], 1)

    def run(self):
        """
        Returns:
            protocol, protocol_dict[protocol]
        """
        protocol = self.get_protocol()
        return protocol, self.protocol_dict[protocol]

    def get_protocol(self):
        """check protocol in the fq1_list"""
        fq_protocol = {}
        for fastq1 in self.fq1_list:
            protocol = self.get_fq_protocol(fastq1)
            fq_protocol[fastq1] = protocol
        if len(set(fq_protocol.values())) != 1:
            sys.exit(f"Error: multiple protocols are not allowed for one sample: {self.sample}! \n" + str(fq_protocol))
        protocol = list(fq_protocol.values())[0]
        return protocol

    def is_protocol(self, seq, protocol):
        """check if seq matches the barcode of protocol"""
        raw_list, mismatch_list = self.mismatch_dict[protocol]
        bc_list = [seq[x] for x in self.protocol_dict[protocol]["pattern_dict"]["C"]]
        valid, _corrected, _res = check_seq_mismatch(bc_list, raw_list, mismatch_list)
        return valid

    def seq_protocol(self, seq):
        """
        Returns: protocol or None

        >>> import tempfile
        >>> runner = Auto([], "fake_sample")
        >>> seq = "TCGACTGTC" + "ATCCACGTGCTTGAGA" + "TTCTAGGAT" + "TCAGCATGCGGCTACG" + "TGCACGAGA" + "C" + "CATATCAATGGG" + "TTTTTTTTTT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-V2'
        >>> seq = "NCAGATTC" + "TCGGTGACAGCCATAT" + "GTACGCAA" + "CGTAGTCAGAAGCTGA" + "CTGAGCCA" + "C" + "TCCGAAGCCCAT" + "TTTTTTTTTT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-V1'
        >>> seq = "NCAGATTC" + "TCGGTGACAGCCATAT" + "GTACGCAA" + "CGTAGTCAGAAGCTGA" + "CTGAGCCA"  + "TCCGAAGCC" + "CTGTCT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-V1'
        >>> seq = "NCAGATTC" + "TCGGTGACAGCCATAT" + "GTACGCAA" + "CGTAGTCAGAAGCTGA" + "CTGAGCCA"  + "TCCGAAGCC"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-V1'
        >>> seq = "GCTAGCTAGCTA" + "GCTAGCTA" + "C" + "TTTTTTTTTT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-MicroBead'
        """
        for start in range(9, 12):
            if seq[start : start + 8] == "ATCGATCG":
                return "GEXSCOPE-V3"

        for protocol in ["GEXSCOPE-V2", "GEXSCOPE-V1"]:
            if self.is_protocol(seq, protocol):
                return protocol

        # check if it is MicroBead
        if seq[16:20] != "TTTT" and seq[22:26] == "TTTT":
            return "GEXSCOPE-MicroBead"

    def get_fq_protocol(self, fq1):
        results = defaultdict(int)

        fq = pyfastx.Fastx(fq1)
        n = 0
        for name, seq, qual in fq:
            n += 1
            protocol = self.seq_protocol(seq)
            if protocol:
                results[protocol] += 1
            if n == self.max_read:
                break
        sorted_counts = sorted(results.items(), key=lambda x: x[1], reverse=True)
        logger.info(sorted_counts)

        protocol, read_counts = sorted_counts[0]
        percent = float(read_counts) / n
        if percent < 0.5:
            logger.warning("Valid protocol read counts percent < 0.5")
        if percent < 0.1:
            logger.error("Valid protocol read counts percent < 0.1")
            raise Exception("Auto protocol detection failed! ")
        logger.info(f"{fq1}: {protocol}")

        return protocol
