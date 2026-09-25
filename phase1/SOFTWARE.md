# Methods tools and where they live in this repository

Versions in the Methods are listed next to the Docker image or command actually vendored here. When the image tag does not state the binary version, the row says so. Do not treat an unpinned image as confirmation of the Methods version.

Reference builds cited in the Methods:

- GRCh38 no-alt: `ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`
- CHM13v2.0: `https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz`

## Version matches the Methods

| Tool | Methods version | In this repo |
|------|-----------------|--------------|
| hifiasm | 0.16.1 | [PBAssembleWithHifiasm](terra/wdl/long-read-pipelines/PBAssembleWithHifiasm/PBAssembleWithHifiasm.wdl), image `us.gcr.io/broad-dsp-lrma/lr-hifiasm:0.16.1` |
| Sniffles2 | 2.0.6 | Task inside [PBCCSWholeGenome](terra/wdl/long-read-pipelines/PBCCSWholeGenome/PBCCSWholeGenome.wdl), image `us.gcr.io/broad-dsp-lrma/lr-sniffles2:2.0.6` |
| PAV | 1.2.1 | [pav](terra/wdl/pav-wdl/pav/pav.wdl), image `us.gcr.io/broad-dsp-lrma/lr-pav:1.2.1.aou` |
| GLnexus | 1.4.1 | [LRJointCallGVCFs](terra/wdl/long-read-pipelines/LRJointCallGVCFs/LRJointCallGVCFs.wdl), image `ghcr.io/dnanexus-rnd/glnexus:v1.4.1` |
| SHAPEIT4 | 4.2.2 | Command `shapeit4.2` in [ChromosomePhasedPanelCreationFromHiPhase](terra/wdl/lrma-aou1-panel-creation/ChromosomePhasedPanelCreationFromHiPhase/ChromosomePhasedPanelCreationFromHiPhase.wdl). Image tag is `hangsuunc/shapeit4:v1`, which does not encode the 4.2.2 patch. |

## Image tag does not confirm the Methods version

These workflows run the tool, but the Docker tag is an internal image name. Confirm the binary inside the image before treating it as the Methods version.

| Tool | Methods version | Image in this repo | Workflow |
|------|-----------------|--------------------|----------|
| pbmm2 | 1.12.0 | `us.gcr.io/broad-dsp-lrma/lr-pb:0.1.38` | Alignment task consumed by downstream calling. [PBCCSWholeGenome](terra/wdl/long-read-pipelines/PBCCSWholeGenome/PBCCSWholeGenome.wdl) starts from aligned BAMs. The Methods command is `pbmm2 align <input.bam> <ref.fa> <output.bam> --preset CCS --sample <sample_name> --strip --sort --unmapped`. |
| DeepVariant | 1.4.0 | `us.gcr.io/broad-dsp-lrma/lr-deepvariant:1.3.0` | [PBCCSWholeGenome](terra/wdl/long-read-pipelines/PBCCSWholeGenome/PBCCSWholeGenome.wdl) (`repo/wdl/tasks/VariantCalling/CCSPepper.wdl`). Pepper image is `kishwars/pepper_deepvariant:r0.4.1`. |
| PBSV | 2.6.0 | `us.gcr.io/broad-dsp-lrma/lr-sv:0.1.8` | PBSV task inside PBCCSWholeGenome. Methods commands: `pbsv discover --tandem-repeats hg38.trf.bed.gz` and `pbsv call`. |
| Kanpig | 0.3.1 | `fcunial/kanpig_experiments` | [GetRegenotypedVcfKanpigV11](terra/wdl/callset-integration/GetRegenotypedVcfKanpigV11/GetRegenotypedVcfKanpigV11.wdl). Methods parameters for the ONT check: `--hapsim 0.97 --chunksize 500 --maxpaths 1000 --gpenalty 0.04`. |
| Truvari | 4.3.1 | `us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample` | [TruvariIntrasampleV11](terra/wdl/callset-integration/TruvariIntrasampleV11/TruvariIntrasampleV11.wdl). Methods cite Truvari v4.3.1 for ONT concordance; the image tag does not state that version. |
| HiPhase | cited, version not given in the main Methods | `hangsuunc/hiphase:1.4.5` | [PhysicalPhasing](terra/wdl/lrma-aou1-panel-creation/PhysicalPhasing/PhysicalPhasing.wdl). The workflow name inside the WDL is `PhasicalPhasing`. |

## Described in the Methods, not a workflow in this checkout

| Tool | Methods version | What this repo has |
|------|-----------------|--------------------|
| TRGT | 0.3.3 | Optional VCF input to PhysicalPhasing. No caller workflow. Command in the Methods: `trgt --genome <ref.fa> --reads <bam> --repeats <trgtCatalog> --output-prefix <sample_name>_TRGT` with `--karyotype XY` or `XX`, Adotto v1.0 padded by 25 bp. |
| LongTR | 1.2 | No workflow. Methods: `--bams <bam> --regions <locus.bed> --fasta <ref.fa> --min-reads 0 --tr-vcf <sample_name>.vcf.gz` on 50 ONT samples. |
| ExpansionHunter | 5.0.0 | No workflow. Run on CDR v8 srWGS at *ATXN3*, *FMR1*, *FGF14*, *HTT*, and *RFC1*. |
| minimap2 (ONT) | 2.24 | No standalone Phase 1 workflow. Methods: `minimap2 -ayYL --MD --eqx -x map-ont`. |
| MosDepth, NanoPlot | not versioned in the main Methods | Tasks exist under vendored `terra/wdl/` imports. They are not a top-level stage in [PIPELINE.md](PIPELINE.md). |
| CYP2D6 callers | Aldy 4.4, Cyrius 1.1.1, PyPGx 0.20.0, StellarPGx 1.2.7, Augustus 3.5.0 | [rw/notebooks/12_cyp2d6.ipynb](rw/notebooks/12_cyp2d6.ipynb) is a public Workbench notebook. Several cells are `sys.argv` scripts pasted into the notebook, so it is not a step-by-step run. Short-read consensus command in the Supplementary Methods: `pypgx run-ngs-pipeline --variants variants.vcf.gz --depth-of-coverage depth.zip --control-statistics control.zip --assembly GRCh38 CYP2D6 results`. |
| AlphaFold3 | not versioned | Not in this repository. |
| Locityper (UK Biobank replication) | not versioned | [github.com/tprodanov/ukbb-scripts](https://github.com/tprodanov/ukbb-scripts). |
| Repeat-expansion figures | — | [rw/notebooks/13_repeat_expansion_figures.ipynb](rw/notebooks/13_repeat_expansion_figures.ipynb) clones `EichlerLab/AoU_WDL` at runtime. |

Upstream repositories that are not fully vendored are listed in [README.md](README.md#upstream-repositories). Workflows named in older docs but absent from the checkout (EvaluateAssemblyHap, SummarizePAVSVs, SummarizeSnifflesSVs) are listed in [terra/wdl/README.md](terra/wdl/README.md#not-vendored).
