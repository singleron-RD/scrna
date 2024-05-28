#!/usr/bin/env python

import argparse
import sys

import parse_protocol
import utils

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

        if args.protocol == "auto":
            runner = parse_protocol.Auto(fq1_list, args.sample)
            protocol, protocol_meta = runner.run()
        else:
            protocol = args.protocol
            protocol_meta = parse_protocol.get_protocol_dict(args.assets_dir)[protocol]
        self.protocol = protocol

        if protocol == "new":
            pattern = args.pattern
            whitelist_str = args.whitelist
        else:
            pattern = protocol_meta["pattern"]
            whitelist_str = " ".join(protocol_meta.get("bc", []))

        pattern_args = Starsolo.get_solo_pattern(pattern)
        if not whitelist_str:
            whitelist_str = args.whitelist if args.whitelist else "None"
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

    def write_cmd(self):
        """
        If UMI+CB length is not equal to the barcode read length, specify barcode read length with --soloBarcodeReadLength.
        To avoid checking of barcode read length, specify soloBarcodeReadLength 0
        """
        prefix = self.args.sample + "."
        cmd = (
            "STAR \\\n"
            f"{self.cb_umi_args} \\\n"
            f"--genomeDir {self.args.genomeDir} \\\n"
            f"--readFilesIn {self.args.fq2} {self.args.fq1} \\\n"
            f"--readFilesCommand {self.read_command} \\\n"
            f"--outFileNamePrefix {prefix} \\\n"
            f"{self.args.ext_args} "
        )
        logger.info(cmd)
        with open(self.cmd_fn, "w") as f:
            f.write(cmd)

    def write_stats(self, assay):
        fn = f"{self.args.sample}.{assay}.protocol.stats.json"
        utils.write_json({"Protocol": self.protocol}, fn)


if __name__ == "__main__":
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--genomeDir", required=True)
    parser.add_argument("--fq1", required=True)
    parser.add_argument("--fq2", required=True)
    parser.add_argument("--assets_dir", required=True)
    parser.add_argument("--protocol", required=True)
    parser.add_argument("--ext_args", required=True)
    parser.add_argument("--whitelist")
    parser.add_argument("--pattern")
    # add version
    parser.add_argument("--version", action="version", version="1.0")

    args = parser.parse_args()

    runner = Starsolo(args)
    runner.write_cmd()
    runner.write_stats("scrna")
