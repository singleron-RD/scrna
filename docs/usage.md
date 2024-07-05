# singleron-RD/scrna: Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below. An example `samplesheet.csv` can be found in the [test data repository](https://github.com/singleron-RD/scrna_test_data/tree/master/GEXSCOPE-V2-human).

```bash
--input '[path to samplesheet file]'
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                                                |
| `fastq_2` | Full path to FastQ file reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                                                |

> [!NOTE]
> fastq_1 and fastq_2 must be full path. Relative path are not allowed.

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. 

### Create `samplesheet.csv` using helper script

When you have many samples, manually creating `samplesheet.csv` can be tedious and error-prone. There is a python script [manifest.py](https://github.com/singleron-RD/sccore/blob/main/sccore/cli/manifest.py) that can help you create a `samplesheet.csv` file.

```
pip install sccore
manifest -m manifest.csv -f /workspaces/scrna_test_data/GEXSCOPE-V2
```

Recursively search the specified folders for fastq files and (optional) matched barcode files.

`-m --manifest` Path to the manifest CSV file containing mappings between fastq file prefixes and sample names. An example `manifest.csv` can be found in the [test data repository](https://github.com/singleron-RD/scrna_test_data/tree/master/GEXSCOPE-V2-human).

`-f --folders` Comma-separated paths to folders to search for fastq files. If `--match` is used, all `barcode.tsv.gz` files with sample name in the full path will also be searched.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run singleron-RD/scrna \
 --input ./samplesheet.csv \
 --outdir ./results \
 --star_genome path_to_star_genome_index \
 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run singleron-RD/scrna -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
star_genome: 'path_to_star_genome_index'
<...>
```

If you prefer a web-based graphical interface or an interactive command-line wizard tool to generate the pipeline parameters, you can use [nf-core launch](https://oldsite.nf-co.re/tools/#launch-a-pipeline):

```
pip install nf-core
nf-core launch singleron-RD/scrna
```

### Create genome index

Since indexing is an expensive process in time and resources you should ensure that it is only done once, by retaining the indices generated from each batch of reference files.

When running the data of a certain species for the first time, you can provide `fasta`, `gtf` and `genome_name` instead of `star_genome`. For example,

```yaml
fasta: "https://raw.githubusercontent.com/singleron-RD/test_genome/master/human.GRCh38.99.MT/human.GRCh38.99.MT.fasta"
gtf: "https://raw.githubusercontent.com/singleron-RD/test_genome/master/human.GRCh38.99.MT/human.GRCh38.99.MT.gtf"
genome_name: "human.GRCh38.99.MT"
```

The STAR index files will be saved in `{outdir}/star_genome/{genome_name}/`.
When running data from the same genome later, you can provide `star_genome` to skip the indexing:

```yaml
star_genome: "/workspaces/test/outs/star_genome/human.GRCh38.99.MT/"
```

### Cell-calling algorithm

STARsolo implements two [cell-calling algorithms](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#cell-filtering-calling): Knee filtering(`cellranger2.2`) and EmptyDrop-like filtering(`EmptyDrops_CR`). EmptyDrop-like filtering considers more barcodes with low UMI as real cells, which helps to recover immune cells with low RNA content, but there is also a risk of including more background barcodes.

The cell-calling algorithm is controlled by the `soloCellFilter` parameter, for example
```yaml
soloCellFilter: EmptyDrops_CR
```

### Running the pipeline with test data

This pipeline contains a small test data. The test config file can be found [here](../conf/test.config).
Run the following command to test

```
nextflow run singleron-RD/scrna -profile test,docker --outdir results
```

> [!NOTE]
> This command might fail if you have [trouble connecting to raw.githubusercontent.com](https://github.com/openvinotoolkit/openvino/issues/8492).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull singleron-RD/scrna
```

> [!NOTE]
> This command might fail if you have trouble connecting to github. In this case, you can manually git clone the master branch and run with the path to the folder.
> ```
> git clone https://github.com/singleron-RD/scrna.git
> nextflow run /workspace/pipeline/scrna ...
> ```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [singleron-RD/scrna releases page](https://github.com/singleron-RD/scrna/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!NOTE]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-bg`

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](../conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original). If it still fails after the attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
