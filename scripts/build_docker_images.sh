#!/bin/bash

set -euxo pipefail

DOCKER_REPO="aou-lr"
LABEL=$1

for DOCKER_FILE in $(find docker -name Dockerfile)
do
    DIR_NAME=$(dirname $DOCKER_FILE)
    DOCKER_NAME=$(basename $DIR_NAME)

    TAG="us-central1-docker.pkg.dev/broad-dsp-lrma/$DOCKER_REPO/$DOCKER_NAME:$LABEL"

    if docker manifest inspect $TAG > /dev/null; then
        docker pull $TAG
    fi

    docker build -t $TAG $DIR_NAME && docker push $TAG
done
