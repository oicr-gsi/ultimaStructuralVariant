#!/bin/bash
set -o nounset
set -o pipefail

diff -s "$1" "$2"
