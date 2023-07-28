version 1.0

task Print {
    input {
        String message
    }

    command <<<
        echo ~{message}
    >>>

    output {
        String text = read_string(stdout())
    }

    runtime {
        cpu:            1
        memory:         10 + " GiB"
        disks:          "local-disk " + 10 + " HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker:         "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
    }
}