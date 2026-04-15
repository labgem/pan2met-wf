# labgem/pan2met-wf

[![GitHub Actions CI Status](https://github.com/labgem/pan2met-wf/actions/workflows/nf-test.yml/badge.svg)](https://github.com/labgem/pan2met-wf/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/labgem/pan2met-wf/actions/workflows/linting.yml/badge.svg)](https://github.com/labgem/pan2met-wf/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**labgem/pan2met-wf** is a bioinformatics pipeline that predicts metabolic network of procaryotes at pangenome scale.
It takes, either a set of genomes, a [PPanGGOLiN](https://github.com/labgem/PPanGGOLiN) pangenome, or directly a proteome in FASTA format and predict what metabolic pathways forms the metabolism of the organism.

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

A typical workflow boils down into these main steps:

- Prepare a PPanGGOLiN pangenome ([PPanGGOLiN](https://github.com/labgem/PPanGGOLiN))
- Identify reference proteins for gene families ([PPanGGOLiN](https://github.com/labgem/PPanGGOLiN))
- Align MetaCyc proteins with the gene families reference proteins ([diamond blastp](https://github.com/bbuchfink/diamond)); associate the found homolog proteins with MetaCyc reaction identifiers.
- Align EcoCyc proteins with the gene families reference proteins ([diamond blastp](https://github.com/bbuchfink/diamond)); associate the found homolog proteins with MetaCyc reaction identifiers.
- Scan the gene families reference proteins for KEGG KOFam ([KOFamScan](https://github.com/takaram/kofam_scan)); associate the KEGG KO identified, with MetaCyc reaction identifiers, or, if none found, EC-numbers.
- Merge EcoCyc-, MetaCyc- and KOFam-based MetaCyc reaction identifiers annotations. Keep at first, the EcoCyc-based annotation. If none found, keep the MetaCyc-based annotation. Again, if none found, keep the MetaCyc reaction identifiers obtained through KOFamScan. In last resort, if no MetaCyc reaction identifier could be found, and there exist a EC-number associated with the protein, keep this EC-number.7

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

<!-- TODO nf-core: Describe the minimum required steps to execute the pipeline, e.g. how to prepare samplesheets.
     Explain what rows and columns represent. For instance (please edit as appropriate):

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

-->


Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run labgem/pan2met-wf \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

labgem/pan2met-wf was originally written by Samuel Ortion.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use labgem/pan2met-wf for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
