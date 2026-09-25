version 1.0

struct VcfdistOutputs {
    File summary_vcf
    File precision_recall_summary_tsv
    File precision_recall_tsv
    File query_tsv
    File truth_tsv
    File phasing_summary_tsv
    File switchflips_tsv
    File superclusters_tsv
    File phase_blocks_tsv
}

# TODO trivial wrapper so we can reuse code, refactor later
struct VcfdistPrecisionRecallSummary {
    File precision_recall_summary_tsv
}

struct OverlapMetricsOutputs {
    File inconsistent_tsv
    File metrics_tsv
}

workflow VcfdistAndOverlapMetricsEvaluation {
    input {
        Array[String] samples
        Array[File] confident_regions_bed_files
        File truth_vcf
        File truth_vcf_idx
        File eval_vcf
        File eval_vcf_idx
        String region
        File reference_fasta
        File reference_fasta_fai

        File vcfdist_bed_file
        String? vcfdist_extra_args
        Int? vcfdist_mem_gb

        Boolean do_overlap_metrics = true
        String overlap_phase_tag
        String overlap_metrics_docker
    }

    scatter (i in range(length(samples))) {
        String sample = samples[i]

        call SubsetSampleFromVcf as SubsetSampleFromVcfEval { input:
            vcf = eval_vcf,
            vcf_idx = eval_vcf_idx,
            sample = sample,
            region = region,
            bed_file = confident_regions_bed_files[i],
            reference_fasta_fai = reference_fasta_fai
        }

        call SubsetSampleFromVcf as SubsetSampleFromVcfTruth { input:
            vcf = truth_vcf,
            vcf_idx = truth_vcf_idx,
            sample = sample,
            region = region,
            bed_file = confident_regions_bed_files[i],
            reference_fasta_fai = reference_fasta_fai
        }

        call Vcfdist { input:
            sample = sample,
            eval_vcf = SubsetSampleFromVcfEval.single_sample_vcf,
            truth_vcf = SubsetSampleFromVcfTruth.single_sample_vcf,
            bed_file = vcfdist_bed_file,
            reference_fasta = reference_fasta,
            extra_args = vcfdist_extra_args,
            mem_gb = vcfdist_mem_gb
        }
    }

    if (do_overlap_metrics) {
        call CalculateOverlapMetrics { input:
            vcf = eval_vcf,
            vcf_idx = eval_vcf_idx,
            region = region,
            phase_tag = overlap_phase_tag,
            docker = overlap_metrics_docker
        }
    }

    output {
        # per-sample
        Array[VcfdistOutputs] vcfdist_outputs_per_sample = Vcfdist.outputs
        Array[VcfdistPrecisionRecallSummary] precision_recall_summary_per_sample = Vcfdist.precision_recall_summary

        # per-cohort
        OverlapMetricsOutputs? overlap_metrics_outputs = CalculateOverlapMetrics.outputs
    }
}

task SubsetSampleFromVcf {
    input {
        File vcf
        File vcf_idx
        String sample
        String region
        File? bed_file
        File reference_fasta_fai
    }

    Int disk_size = 1 + ceil(2 * (size(vcf, "GiB")))

    command <<<
        set -euxo pipefail

        # must use -T bed_file to intersect with -r region properly
        bcftools view ~{vcf} \
            -s ~{sample} \
            -r ~{region} \
            ~{"-T " + bed_file} \
            -Oz -o ~{sample}.subset.g.vcf.gz
        bcftools reheader ~{sample}.subset.g.vcf.gz \
            --fai ~{reference_fasta_fai} \
            -o ~{sample}.subset.reheadered.g.vcf.gz
        bcftools index -t ~{sample}.subset.reheadered.g.vcf.gz
    >>>
    
    output {
        File single_sample_vcf = "~{sample}.subset.reheadered.g.vcf.gz"
        File single_sample_vcf_tbi = "~{sample}.subset.reheadered.g.vcf.gz.tbi"
    }

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task Vcfdist {
    input {
        String sample
        File eval_vcf
        File truth_vcf
        File bed_file
        File reference_fasta
        String? extra_args
        Int verbosity = 1

        Int disk_size_gb = ceil(size(truth_vcf, "GiB") + 10)
        Int mem_gb = 32
        Int cpu = 4
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        vcfdist \
            ~{eval_vcf} \
            ~{truth_vcf} \
            ~{reference_fasta} \
            -b ~{bed_file} \
            -v ~{verbosity} \
            ~{extra_args}

        for tsv in $(ls *.tsv); do mv $tsv ~{sample}.$tsv; done
        mv summary.vcf ~{sample}.summary.vcf
    >>>

    output {
        VcfdistOutputs outputs = {
            "summary_vcf": "~{sample}.summary.vcf",
            "precision_recall_summary_tsv": "~{sample}.precision-recall-summary.tsv",
            "precision_recall_tsv": "~{sample}.precision-recall.tsv",
            "query_tsv": "~{sample}.query.tsv",
            "truth_tsv": "~{sample}.truth.tsv",
            "phasing_summary_tsv": "~{sample}.phasing-summary.tsv",
            "switchflips_tsv": "~{sample}.switchflips.tsv",
            "superclusters_tsv": "~{sample}.superclusters.tsv",
            "phase_blocks_tsv": "~{sample}.phase-blocks.tsv"
        }

        VcfdistPrecisionRecallSummary precision_recall_summary = {
            "precision_recall_summary_tsv": "~{sample}.precision-recall-summary.tsv"
        }
    }

    runtime {
        docker: "timd1/vcfdist:v2.5.3"
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}

# see description at https://github.com/broadinstitute/lrma-aou1-panel-creation/pull/22#issue-2456295093
task CalculateOverlapMetrics {
    input {
        File vcf
        File vcf_idx
        String region
        String? phase_tag    # e.g., "PG" for kanpig integrated, "PS" for HiPhase, "NONE" for Shapeit4

        # docker needs bcftools, scikit-allel, and pandas
        String docker
        Int disk_size_gb = ceil(size(vcf, "GiB") + 10)
        Int mem_gb = 32
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        set -euxo pipefail

        # split to biallelic
        bcftools norm \
            -r ~{region} \
            -m- --do-not-normalize ~{vcf} \
            -Oz -o split.vcf.gz
        bcftools index -t split.vcf.gz

        # join and subset to multiallelic
        bcftools norm \
            -r ~{region} \
            -m+ --do-not-normalize ~{vcf} | \
            bcftools view --min-alleles 3 \
                -Oz -o multi.vcf.gz
        bcftools index -t multi.vcf.gz

        python - --split_vcf split.vcf.gz \
                 --joined_multi_vcf multi.vcf.gz \
                 --phase_tag ~{phase_tag} \
                 <<-'EOF'
        import argparse
        import allel
        import numpy as np
        import pandas as pd
        from collections import Counter

        def calculate_metrics(split_vcf, joined_multi_vcf, phase_tag):
            joined_multi_callset = allel.read_vcf(joined_multi_vcf)
            split_callset = allel.read_vcf(split_vcf,
                                           fields=['samples', 'calldata/GT', 'CHROM', 'POS'] + ([] if phase_tag == 'NONE' else [f'calldata/{phase_tag}']))

            # V indexes all split records, v indexes only those that are in multiallelic sites
            is_multi_V = np.isin(split_callset['variants/POS'], joined_multi_callset['variants/POS'])
            gt_vsp = split_callset['calldata/GT'][is_multi_V]
            is_phased_vs = np.ones(gt_vsp.shape[:2]) if phase_tag == 'NONE' else split_callset[f'calldata/{phase_tag}'][is_multi_V] > 0
            chrom_v = split_callset['variants/CHROM'][is_multi_V]
            pos_v = split_callset['variants/POS'][is_multi_V]
            samples_s = split_callset['samples']

            # group split records by multiallelic site (indexed by m)
            m_splits = np.unique(pos_v, return_index=True)[1][1:]
            gt_mvsp = np.split(gt_vsp, m_splits)
            is_phased_mvs = np.split(is_phased_vs, m_splits)
            chrom_mv = np.split(chrom_v, m_splits)
            pos_mv = np.split(pos_v, m_splits)
            
            # calculate metrics and record inconsistent positions/samples
            metrics = Counter({'NUM_CONSISTENT_ALLELES': 0,
                               'NUM_INCONSISTENT_ALLELES': 0,
                               'NUM_CONSISTENT_SITES': 0,
                               'NUM_INCONSISTENT_SITES': 0})
            inconsistent_records = []

            for m, site_gt_vsp in enumerate(gt_mvsp):
                # site is considered inconsistent if at least one sample has phased genotype calls
                # that yield multiple alternate alleles on the same haplotype
                is_phased_alt_vsp = ((site_gt_vsp != 0) & (site_gt_vsp != -1)) * is_phased_mvs[m][:, :, np.newaxis]
                is_inconsistent_sp = (is_phased_alt_vsp).sum(axis=0) > 1
                is_inconsistent = is_inconsistent_sp.any()
                
                if is_inconsistent:
                    # find inconsistent alleles (i.e., split records indexed by v)
                    is_inconsistent_v = is_phased_alt_vsp[:, is_inconsistent_sp].any(axis=1)

                    # record genomic positions, samples, allele indices
                    chrom = chrom_mv[m][0]
                    start = pos_mv[m][0]
                    end = pos_mv[m][-1]
                    is_inconsistent_s = is_inconsistent_sp.any(axis=1)
                    samples = ','.join(samples_s[is_inconsistent_s])
                    allele_indices = ','.join(np.where(is_inconsistent_v)[0].astype(str))
                    inconsistent_records.append([chrom, start, end, samples, allele_indices])

                    # update counts
                    num_inconsistent_alleles = is_inconsistent_v.sum()
                    metrics['NUM_INCONSISTENT_ALLELES'] += num_inconsistent_alleles
                    metrics['NUM_CONSISTENT_ALLELES'] += len(site_gt_vsp) - num_inconsistent_alleles
                    metrics['NUM_INCONSISTENT_SITES'] += 1
                else:
                    metrics['NUM_CONSISTENT_ALLELES'] += len(site_gt_vsp)
                    metrics['NUM_CONSISTENT_SITES'] += 1

            # write inconsistent records
            inconsistent_df = pd.DataFrame(inconsistent_records,
                                           columns=['CHROM', 'START', 'END', 'SAMPLES', 'ALLELE_INDICES'])
            inconsistent_df.to_csv('inconsistent.tsv', sep='\t', index=False)

            # write in/consistent counts
            metrics_df = pd.DataFrame.from_records([dict(metrics)])
            metrics_df.to_csv('metrics.tsv', sep='\t', index=False,
                              columns=['NUM_INCONSISTENT_ALLELES', 'NUM_CONSISTENT_ALLELES',
                                       'NUM_INCONSISTENT_SITES', 'NUM_CONSISTENT_SITES'])


        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--split_vcf',
                                type=str)

            parser.add_argument('--joined_multi_vcf',
                                type=str)

            parser.add_argument('--phase_tag',
                                type=str)

            args = parser.parse_args()

            calculate_metrics(args.split_vcf,
                              args.joined_multi_vcf,
                              args.phase_tag)

        if __name__ == '__main__':
            main()
        EOF
    >>>

    output {
        OverlapMetricsOutputs outputs = {
            "inconsistent_tsv": "inconsistent.tsv",
            "metrics_tsv": "metrics.tsv"
        }
    }

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }
}
