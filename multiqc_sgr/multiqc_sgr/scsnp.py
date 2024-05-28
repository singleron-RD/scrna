import json
import logging
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

# Initialise the logger
log = logging.getLogger("multiqc")


def get_frac(x):
    return round(float(x) * 100, 2)


def get_int(x):
    return str(int(x))


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name="scsnp",
            anchor="scsnp",
            info="Single cell amplicon variant calling.",
        )

        # Find and load any STAR reports
        stats_data = self.parse_json("stats")
        gene_data = self.parse_json("gene")
        count_data = self.parse_json("count")
        meta_data = self.parse_json("meta")
        if all(len(x) == 0 for x in [stats_data, count_data, meta_data]):
            raise ModuleNoSamplesFound

        self.general_stats(stats_data)

        self.add_section(
            name="Gene",
            anchor="scsnp_gene",
            helptext="Cell Reads assigned to each gene",
            plot=self.gene_bar(gene_data),
        )

        # assgin plot
        # https://pcingola.github.io/SnpEff/snpsift/extractfields/
        helptext = """
                Effect: Effect in Sequence ontology terms (e.g. 'missense_variant', 'synonymous_variant', 'stop_gained', etc.)

                Impact: There are four levels of impacts predicted by SnpEff:

                * **High**: High impact (like stop codon)
                * **Moderate**: Middle impact (like same type of amino acid substitution)
                * **Low**: Low impact (ie silence mutation)
                * **Modifier**: No impact

                Rank: Exon or Intron rank (i.e. exon number in a transcript)

                HGVS_C: Variant in [HGVS DNA notation](https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/)

                HGVS_P: Variant in [HGVS protein notation](https://hgvs-nomenclature.org/stable/recommendations/protein/deletion/)

                """
        self.add_section(
            name="Variant after filtering", anchor="scsnp_meta", helptext=helptext, plot=self.meta_table(meta_data)
        )

        helptext = """
            The number of cells with each genotype.
            
            Genotypes: From the [VCF version 4.1](https://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41): 
            
            GT : genotype, encoded as allele values separated by either of / or |. The allele values are 0 for the reference allele (what is in the REF field), 1 for the first allele listed in ALT, 2 for the second allele list in ALT and so on. For diploid calls examples could be 0/1, 1 | 0, or 1/2, etc.
            
            'NA' means not available(no reads at this position were found).
            """
        self.add_section(name="Variant count", anchor="scsnp_count", helptext=helptext, plot=self.count_bar(count_data))

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

    def parse_json(self, seg):
        data_dict = defaultdict(dict)
        for f in self.find_log_files(f"scsnp/{seg}"):
            parsed_data = json.loads(f["f"])
            if parsed_data is not None:
                x = f["s_name"]
                s_name = x[: x.find(".scsnp")]
                if s_name in data_dict:
                    log.info(f"Duplicate sample name found! Update: {s_name}")
                self.add_data_source(f, s_name=s_name, section=seg)
                data_dict[s_name].update(parsed_data)

        data_dict = self.ignore_samples(data_dict)

        log.info(f"Found {len(data_dict)} scsnp {seg} reports")
        # Write parsed report data to a file
        self.write_data_file(data_dict, f"multiqc_scsnp_{seg}")
        return data_dict

    def general_stats(self, stats_data):
        headers = {
            "Protocol": {
                "title": "Protocol",
                "description": "barcode and UMI protocol",
            },
            "Raw Reads": {"title": "Raw Reads", "description": "Number of reads from fastq files", "format": "{:,.0f}"},
            "Valid Reads": {
                "title": "Valid Reads",
                "description": "fraction of reads with valid barcode and UMI",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "modify": get_frac,
                "format": "{:,.2f}",
            },
            "Mean Used Reads per Cell": {
                "title": "Mean Used Reads",
                "description": "Mean Used Reads per Cell: Reads mapped to target genes and remove PCR duplicates",
                "format": "{:,.0f}",
            },
            "Saturation": {
                "title": "Saturation",
                "description": "fraction of reads originating from an already-observed reads. saturation = 1 - distinct_reads / total_reads",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "green",
                "modify": get_frac,
                "format": "{:,.2f}",
            },
            "Number of variant after filtering": {
                "title": "N variant",
                "description": "Number of variant after filtering",
                "format": "{:,.0f}",
            },
        }
        self.general_stats_addcols(stats_data, headers=headers)

    def gene_bar(self, gene_data):
        pconfig = {
            "id": "gene_read",
            "title": "Read count per gene",
        }
        return bargraph.plot(gene_data, pconfig=pconfig)

    def meta_table(self, meta_data):
        name_meta = {}
        for sample in meta_data:
            for name in meta_data[sample]:
                name_meta[name] = meta_data[sample][name]
        table_config = {
            "id": "scsnp_meta",
            "title": "Variant after filtering",
            "col1_header": "Variant",
        }
        return table.plot(name_meta, pconfig=table_config)

    def count_bar(self, count_data):
        v_sample = defaultdict(dict)
        for sample in count_data:
            for v in count_data[sample]:
                v_sample[v].update({sample: count_data[sample][v]})
        cats = sorted(list(v_sample.keys()))
        plot_data = [v_sample[v] for v in cats]
        pconfig = {
            "id": "variant_count",
            "title": "Variant count",
            "data_labels": cats,
        }
        return bargraph.plot(plot_data, pconfig=pconfig)
