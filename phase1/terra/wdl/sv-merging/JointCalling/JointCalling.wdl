########################################################################################
# JointCalling
# Helper from the SV-merging comparisons, not on the Phase 1 launch path.
#
# Joint-calls SVs across samples with Sniffles2 multi-sample mode, from each sample's
# Sniffles2 signature (.snf) file.
#
# Inputs:  snfs  a file listing the .snf paths; n_files
# Outputs: joint_vcf, counts
#
# Upstream repository, branch, and any local patches: SOURCE.txt in this directory.
# Example inputs from a Phase 1 run: example_inputs.json in this directory.
########################################################################################

version 1.0


workflow JointCalling {
    input {
        File snfs
        Int n_files
    }

    call JointCallingImpl { 
        input: 
            snfs = snfs,
            n_files = n_files
    }
    
    output {
        File joint_vcf = JointCallingImpl.joint_vcf
        File joint_tbi = JointCallingImpl.joint_tbi
        File counts = JointCallingImpl.counts
    }
}


# Remark: sniffles2 stores all signature files in RAM, so we should give the
# procedure an upper bound on the size of a signature file.
#
task JointCallingImpl {
    input {
        File snfs
        Int n_files
    }

    command <<<
        set -euxo pipefail
        
        GSUTIL_DELAY_S="600"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        
        rm -f list.txt; touch list.txt
        rm -f counts.txt; touch counts.txt
        i=0
        while read SNF_FILE; do
            LOCAL_FILE=$(basename ${SNF_FILE})
            while : ; do
                TEST=$(gsutil -m cp ${SNF_FILE} ${LOCAL_FILE} && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${SNF_FILE}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            echo -e "${LOCAL_FILE}\t${LOCAL_FILE%.sniffles2.snf}" >> list.tsv
            i=$(($i + 1))
            if [ $i -eq ~{n_files} ]; then
                break
            fi
        done < ~{snfs}
        head list.tsv
        ${TIME_COMMAND} sniffles --threads ${N_THREADS} --combine-separate-intra --input list.tsv --vcf joint.vcf
        N_INS=$(grep "SVTYPE=INS" joint.vcf | awk '{ if ($7=="PASS") print $0; }' | wc -l)
        N_DEL=$(grep "SVTYPE=DEL" joint.vcf | awk '{ if ($7=="PASS") print $0; }' | wc -l)
        echo "${N_INS},${N_DEL}" >> counts.txt
        bgzip --threads ${N_THREADS} joint.vcf
        tabix joint.vcf.gz
    >>>

    output {
        File joint_vcf = "joint.vcf.gz"
        File joint_tbi = "joint.vcf.gz.tbi"
        File counts = "counts.txt"
    }
    runtime {
        cpu: 32
        memory: "128 GiB"
        disks: "local-disk 100 HDD"
        preemptible: 0
        docker: "fcunial/simulation"
    }
}