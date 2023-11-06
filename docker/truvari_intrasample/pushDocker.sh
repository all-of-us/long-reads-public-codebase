#!/bin/bash
#
docker build --progress=plain -t fcunial/truvari_intrasample .
docker push fcunial/truvari_intrasample
