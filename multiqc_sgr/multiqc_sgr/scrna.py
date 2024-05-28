import json
import logging
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger("multiqc")


def get_int(x):
    return str(int(x))


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name="scrna",
            anchor="scrna",
            info="mapping, demultiplexing and quantification for single cell RNA-seq",
        )
        log.info(f"Running module: {self.name}")

        stat_data = self.parse_json(self.name, "stats")
        umi_count_data = self.parse_json(self.name, "umi_count")
        saturation_data = self.parse_json(self.name, "saturation")
        median_gene_data = self.parse_json(self.name, "median_gene")
        if all(len(x) == 0 for x in [stat_data, umi_count_data, saturation_data, median_gene_data]):
            raise ModuleNoSamplesFound

        # Basic Stats Table
        self.general_stats_table(stat_data)

        # barcode rank plot
        self.add_section(name="Barcode Rank", anchor="scrna_barcode_rank", plot=self.barcode_rank_plot(umi_count_data))

        # subsample
        if saturation_data:
            self.add_section(name="Saturation", anchor="scrna_subsample", plot=self.saturation_plot(saturation_data))
        if median_gene_data:
            self.add_section(
                name="Median Gene", anchor="scrna_median_gene", plot=self.median_gene_plot(median_gene_data)
            )

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

    def parse_json(self, assay, seg):
        data_dict = defaultdict(dict)
        n = 0
        for f in self.find_log_files(f"{assay}/{seg}"):
            log.debug(f"Found file: {f['fn']}")
            n += 1
            parsed_data = json.loads(f["f"])
            if parsed_data is not None:
                x = f["s_name"]
                s_name = x[: x.find(f".{assay}")]
                if s_name in data_dict:
                    log.info(f"Duplicate sample name found! Update: {s_name}")
                self.add_data_source(f, s_name=s_name, section=seg)
                data_dict[s_name].update(parsed_data)

        data_dict = self.ignore_samples(data_dict)

        log.info(f"Found {n} {assay} {seg} reports")
        # Write parsed report data to a file
        self.write_data_file(data_dict, f"multiqc_{assay}_{seg}")
        return data_dict

    def general_stats_table(self, summary_data):
        headers = {
            "Protocol": {
                "title": "Protocol",
                "description": "Predefined pattern of barcode and UMI",
                "scale": "purple",
                "hidden": True,
            },
            "Raw Reads": {
                "title": "Raw Reads",
                "description": "Number of reads in the input file",
                "scale": "blue",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "Valid Reads": {
                "title": "Valid Reads",
                "description": "Percent of reads with valid barcodes",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Corrected Barcodes": {
                "title": "Corrected Barcodes",
                "description": "Percent of corrected barcodes",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Reads Mapped To Unique Loci": {
                "title": "Unique Reads",
                "description": "Percent of reads mapped uniquely to genome",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Reads Mapped To Multiple Loci": {
                "title": "Multi-Loci Reads",
                "description": "Percent of reads mapped to multiple loci",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Reads Mapped Uniquely To Transcriptome": {
                "title": "Counted Unique Reads",
                "description": "Percent of reads mapped uniquely to transcriptome; These reads are used for UMI counting",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Mapped Reads Assigned To Exonic Regions": {
                "title": "Exonic Reads",
                "description": "Percent of mapped reads assigned to exonic regions",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Mapped Reads Assigned To Intronic Regions": {
                "title": "Intronic Reads",
                "description": "Percent of mapped reads assigned to intronic regions",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Mapped Reads Assigned To Intergenic Regions": {
                "title": "Intergenic Reads",
                "description": "Percent of mapped reads assigned to intergenic regions",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Mapped Reads Assigned Antisense To Gene": {
                "title": "Antisense Reads",
                "description": "Percent of mapped reads assigned antisense to gene",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "hidden": True,
            },
            "Estimated Number of Cells": {
                "title": "Number of Cells",
                "description": "Estimated number of cells",
                "scale": "blue",
                "format": "{:,.0f}",
            },
            "Fraction Reads in Cells": {
                "title": "Reads in Cells",
                "description": "Percent of unique reads in cells",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
            "Mean Used Reads per Cell": {
                "title": "Mean Reads",
                "description": "Mean number of reads used per cell",
                "scale": "blue",
                "format": "{:,.0f}",
            },
            "Median UMI per Cell": {
                "title": "Median UMI",
                "description": "Median number of UMIs per cell",
                "scale": "blue",
                "format": "{:,.0f}",
            },
            "Median Genes per Cell": {
                "title": "Median Genes",
                "description": "Median number of genes per cell",
                "scale": "blue",
                "format": "{:,.0f}",
            },
            "Total Genes": {
                "title": "Total Genes",
                "description": "Total number of genes detected",
                "scale": "blue",
                "format": "{:,.0f}",
            },
            "Saturation": {
                "title": "Saturation",
                "description": "Percent of reads originating from an already-observed UMI",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
            },
        }
        self.general_stats_addcols(summary_data, headers=headers)

    def barcode_rank_plot(self, umi_count_data):
        plot_data = {}
        colors = {}
        for sample in umi_count_data:
            for sub in umi_count_data[sample]:
                cur = umi_count_data[sample][sub]
                if not cur:
                    continue
                new = {}
                for k, v in cur.items():
                    new[int(k)] = v
                plot_data[sub] = new
                if "pure" in sub:
                    colors[sub] = "darkblue"
                elif "mix" in sub:
                    colors[sub] = "lightblue"
                elif "background" in sub:
                    colors[sub] = "lightgray"

        # Config for the plot
        pconfig = {
            "id": "scrna_barcode_rank_plot",
            "title": "scrna: Barcode Rank",
            "ylab": "UMI counts",
            "xlab": "Barcode Rank",
            "yLog": True,
            "xLog": True,
            "colors": colors,
            "ymin": 0,
            "height": 750,
        }

        return linegraph.plot(plot_data, pconfig)

    def saturation_plot(self, saturation_data):
        # Config for the plot
        pconfig = {
            "id": "scrna_saturation_plot",
            "title": "scrna: Saturation",
            "ylab": "Saturation",
            "xlab": "Percent of Reads",
            "height": 750,
        }

        return linegraph.plot(saturation_data, pconfig)

    def median_gene_plot(self, median_gene_data):
        # Config for the plot
        pconfig = {
            "id": "scrna_median_gene_plot",
            "title": "scrna: Median Gene",
            "ylab": "Median Gene",
            "xlab": "Percent of Reads",
            "height": 750,
        }
        return linegraph.plot(median_gene_data, pconfig)
