#!/bin/bash

set -euxo pipefail

RESOURCES="resources/workflow_options/default.json"
CA="cromshell-alpha -q --hide_logo -t 30 submit -d ./wdl"

$CA wdl/pipelines/HelloWorkflow.wdl inputs/HelloWorkflow.json
