#!/bin/bash

set -euxo pipefail

DOCKER_REGISTRY="us.gcr.io" # us-central1-docker.pkg.dev is not supported on RWB
GOOGLE_PROJECT="broad-dsp-lrma"
DOCKER_REPO="aou-lr"
LABEL=$1

for DOCKER_FILE in $(find docker -name Dockerfile)
do
    DIR_NAME=$(dirname $DOCKER_FILE)
    DOCKER_NAME=$(basename $DIR_NAME)

    VERSION=$(grep '^current_version' .bumpversion.cfg | cut -d' ' -f3)

    TAG_VERSION="$DOCKER_REGISTRY/$GOOGLE_PROJECT/$DOCKER_REPO/$DOCKER_NAME:$VERSION"
    TAG_LABEL="$DOCKER_REGISTRY/$GOOGLE_PROJECT/$DOCKER_REPO/$DOCKER_NAME:$LABEL"

    if docker manifest inspect $TAG_VERSION > /dev/null; then
        # Make sure the most recent version of the Docker image gets used as the cache
        docker pull $TAG_VERSION
        docker tag $TAG_VERSION $TAG_LABEL
    fi

    if docker manifest inspect $TAG_LABEL > /dev/null; then
        docker pull $TAG_LABEL
    fi

    docker build -t $TAG_LABEL $DIR_NAME
    docker push $TAG_LABEL
done
