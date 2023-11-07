#!/bin/bash
#
docker build --progress=plain -t us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample .
docker push us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample
