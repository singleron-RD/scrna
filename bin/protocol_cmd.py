#!/usr/bin/env python

import argparse
import sys

import parse_protocol
import utils
from __init__ import ASSAY

logger = utils.get_logger(__name__)
SOLOFEATURE = "GeneFull_Ex50pAS"


class Starsolo:
    def __init__(self, args):
        self.args = args
        fq1_list = args.fq1.split(",")
        fq2_list = args.fq2.split(",")
        fq1_number = len(fq1_list)
        fq2_number = len(fq2_list)
        if fq1_number != fq2_number:
            sys.exit("fastq1 and fastq2 do not have same file number!")

        self.read_command = "cat"
        if str(fq1_list[0]).endswith(".gz"):
            self.read_command = "zcat"

        if args.protocol == "new":
            protocol = "new"
            pattern = args.pattern
            whitelist_str = args.whitelist
        else:
            if args.protocol == "auto":
                runner = parse_protocol.Auto(fq1_list, args.sample)
                protocol, protocol_meta = runner.run()
            else:
                protocol = args.protocol
                protocol_meta = parse_protocol.get_protocol_dict(args.assets_dir)[protocol]
            pattern = protocol_meta["pattern"]
            whitelist_str = " ".join(protocol_meta.get("bc", []))
        self.protocol = protocol

        if protocol == "GEXSCOPE-V3":
            v3_linker = "ATCG" * 2
            bc = "N" * 9
            pattern_args = (
                "--soloType CB_UMI_Complex "
                "--soloCBposition 2_0_2_8 2_13_2_21 3_1_3_9 "
                "--soloUMIposition 3_10_3_21 "
                f"--soloAdapterSequence {bc}{v3_linker}{bc}{v3_linker} "
                "--soloAdapterMismatchesNmax 1 "
                "--soloCBmatchWLtype EditDist_2 "
            )
        else:
            pattern_args = Starsolo.get_solo_pattern(pattern)
        if not whitelist_str:
            whitelist_str = args.whitelist if args.whitelist else "None"
        whitelist_str = whitelist_str.strip()
        if whitelist_str.startswith("http"):
            whitelist_str = whitelist_str.split("/")[-1]
        if whitelist_str.endswith(".gz"):
            whitelist_str = f"<(gzip -cdf {whitelist_str})"
        self.cb_umi_args = pattern_args + f" --soloCBwhitelist {whitelist_str} "

        # out cmd
        self.cmd_fn = args.sample + ".starsolo_cmd.txt"

    @staticmethod
    def get_solo_pattern(pattern) -> str:
        """
        Returns:
            starsolo_cb_umi_args
        """
        pattern_dict = parse_protocol.parse_pattern(pattern)
        if len(pattern_dict["U"]) != 1:
            sys.exit(f"Error: Wrong pattern:{pattern}. \n Solution: fix pattern so that UMI only have 1 position.\n")
        ul = pattern_dict["U"][0].start
        ur = pattern_dict["U"][0].stop
        umi_len = ur - ul

        if len(pattern_dict["C"]) == 1:
            solo_type = "CB_UMI_Simple"
            start, stop = pattern_dict["C"][0].start, pattern_dict["C"][0].stop
            cb_start = start + 1
            cb_len = stop - start
            umi_start = ul + 1
            cb_str = f"--soloCBstart {cb_start} --soloCBlen {cb_len} --soloCBmatchWLtype 1MM "
            umi_str = f"--soloUMIstart {umi_start} --soloUMIlen {umi_len} "
        else:
            solo_type = "CB_UMI_Complex"
            cb_pos = " ".join([f"0_{x.start}_0_{x.stop-1}" for x in pattern_dict["C"]])
            umi_pos = f"0_{ul}_0_{ur-1}"
            cb_str = f"--soloCBposition {cb_pos} --soloCBmatchWLtype EditDist_2 "
            umi_str = f"--soloUMIposition {umi_pos} --soloUMIlen {umi_len} "

        starsolo_cb_umi_args = " ".join([f"--soloType {solo_type} ", cb_str, umi_str])
        return starsolo_cb_umi_args

    def run(self):
        """
        If UMI+CB length is not equal to the barcode read length, specify barcode read length with --soloBarcodeReadLength.
        To avoid checking of barcode read length, specify soloBarcodeReadLength 0
        """
        cmd = " ".join([self.cb_umi_args, f"--readFilesCommand {self.read_command}"])
        logger.info(cmd)
        with open(self.cmd_fn, "w") as f:
            f.write(cmd)

        utils.write_multiqc({"Protocol": self.protocol}, self.args.sample, ASSAY, "protocol.stats")


if __name__ == "__main__":
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--fq1", required=True)
    parser.add_argument("--fq2", required=True)
    parser.add_argument("--assets_dir", required=True)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--whitelist")
    parser.add_argument("--pattern")
    # add version
    parser.add_argument("--version", action="version", version="1.0")

    args = parser.parse_args()

    runner = Starsolo(args)
    runner.run()
