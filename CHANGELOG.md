# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-11-19
### Added
- Germline and somatic structural variant calling for Ultima data, adapted from the Ultima
  Genomics structural_variant_pipeline workflow (healthomics-workflows v1.37.1) to run on
  Cromwell + UGE with the upstream container images under Apptainer.
- Resources come from the modules `hg38-ultima/v0`, `hg38-ugsv-resources/1.37.1`,
  `hg38-ua-index/v45-79372c0` and `ugsv/1.37.1`.
