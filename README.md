# ultimaStructuralVariant

## Overview

Structural variant calling for Ultima Genomics data, germline or somatic. The genome is split into shards; per shard, reads with soft clips, indels or mismatches are locally assembled into haplotypes. The haplotypes are realigned with the UA aligner and the better of the original and UA alignment is kept per haplotype, homopolymer-rich haplotypes are realigned, and every read is matched to the haplotype it supports. GRIDSS calls breakpoints from the supported assembly and annotates them shard by shard. Germline calls are then linked into events; somatic calls are filtered by GRIPSS against panels of normals. Adapted from the Ultima Genomics [structural_variant_pipeline](https://github.com/Ultimagen/healthomics-workflows/tree/main/workflows/structural_variant_pipeline) workflow (v1.37.1).

Somatic calling runs when `tumorCrams` is given, with `germlineCrams` as the matched normal (optional, for tumor-only); otherwise `germlineCrams` is called as a germline sample. Read-selection thresholds and the default shard count follow the upstream germline and somatic templates. The optional pangenome (Giraffe) realignment of the upstream workflow is not included.

Each task runs its tools from the upstream container images through Apptainer, so the Cromwell backend needs no container support.

## Resources

Everything the workflow reads besides the input alignments comes from four modules, which tasks load through their `modules` runtime attribute:

* `hg38-ultima/v0` sets `HG38_ULTIMA_ROOT`: the hg38 FASTA with its `.fai` and `.dict`.
* `hg38-ugsv-resources/1.37.1` sets `HG38_UGSV_RESOURCES_ROOT`: the calling regions and GRIDSS blacklist under `sv/`, the GRIPSS panels of normals, known fusions and repeat mask under `sv/gripss/`, laid out as in the Ultima public bucket.
* `hg38-ua-index/v45-79372c0` sets `HG38_UA_INDEX_ROOT`: the UA aligner index and the alt-contig list.
* `ugsv/1.37.1` sets `UGSV_ROOT` and loads apptainer: the container images, pulled by digest, under `images/`.

## Dependencies

* [apptainer 1.4.5](https://apptainer.org)
* [hg38-ultima v0](https://github.com/Ultimagen/healthomics-workflows)
* [hg38-ugsv-resources 1.37.1](https://github.com/Ultimagen/healthomics-workflows)
* [hg38-ua-index v45-79372c0](https://github.com/Ultimagen/healthomics-workflows)
* [ultimagenomics make_examples:3.5.0](https://hub.docker.com/r/ultimagenomics/make_examples)
* [ultimagenomics alignment:4.1.4](https://hub.docker.com/r/ultimagenomics/alignment)
* [ultimagenomics rematcher:1.2.0_d0f0b0a](https://hub.docker.com/r/ultimagenomics/rematcher)
* [ultimagenomics gridss:1.0.2](https://hub.docker.com/r/ultimagenomics/gridss)
* [ultimagenomics gripss:ug_2.4.1_165b492](https://hub.docker.com/r/ultimagenomics/gripss)
* [broadinstitute gatk:4.6.0.0](https://hub.docker.com/r/broadinstitute/gatk)


## Usage

### Cromwell
```
java -jar cromwell.jar run ultimaStructuralVariant.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`outputFileNamePrefix`|String|Prefix for output file names


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`tumorCrams`|Array[File]|[]|Tumor CRAM files; all must carry the same read-group sample name. Non-empty selects somatic calling
`tumorCramIndexes`|Array[File]|[]|Indexes of tumorCrams, in the same order
`germlineCrams`|Array[File]|[]|CRAM files of the germline sample, or of the matched normal in somatic calling (optional there, for tumor-only); all must carry the same read-group sample name
`germlineCramIndexes`|Array[File]|[]|Indexes of germlineCrams, in the same order
`reference`|String|"hg38"|Reference genome build, selects the genome resources: hg38
`intervalList`|String?|None|Picard interval list of the regions to call, to restrict a run such as a test; defaults to the whole-genome calling regions of the reference
`numShards`|Int?|None|Number of genomic shards to process in parallel; defaults to 50 for germline and 200 for somatic calling
`scatterIntervalsBreak`|Int|10000000|Break scatter intervals at multiples of this many bases
`symbolicVcfFormat`|Boolean|true|Also write the calls with symbolic alleles (<DEL>, <INS>, <DUP>, <INV>) in place of breakend notation


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`getSampleNames.jobMemory`|Int|2|Memory allocated for this job (GB)
`getSampleNames.cores`|Int|1|Number of cores requested for this job
`getSampleNames.timeout`|Int|1|Hours before task timeout
`getSampleNames.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`scatterIntervals.jobMemory`|Int|8|Memory allocated for this job (GB); the Java heap is set 3 GB below it
`scatterIntervals.cores`|Int|1|Number of cores requested for this job
`scatterIntervals.timeout`|Int|1|Hours before task timeout
`scatterIntervals.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`createAssembly.minBase`|Int|0|Minimum base quality for a base to count toward assembly
`createAssembly.minMapq`|Int|5|Minimum mapping quality for a read to enter assembly
`createAssembly.maxNumHaps`|Int|40|Maximum number of haplotypes assembled per region
`createAssembly.maxReadsPerPartition`|Int|1500|Maximum number of reads assembled per region
`createAssembly.jobMemory`|Int|4|Memory allocated for this job (GB)
`createAssembly.cores`|Int|1|Number of cores requested for this job
`createAssembly.timeout`|Int|12|Hours before task timeout
`createAssembly.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`mergeBams.jobMemory`|Int|16|Memory allocated for this job (GB)
`mergeBams.cores`|Int|8|Number of cores requested for this job, also the samtools thread count
`mergeBams.timeout`|Int|4|Hours before task timeout
`mergeBams.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`alignWithUa.uaExtraArgs`|String|"--seed-score-ratio 0.5 --vector --soft-clipping --realignment-tag re --mismatch-cost -6"|Further ua arguments
`alignWithUa.jobMemory`|Int|48|Memory allocated for this job (GB); ua holds the whole index in memory
`alignWithUa.cores`|Int|16|Number of cores requested for this job, also the ua and samtools thread count
`alignWithUa.timeout`|Int|12|Hours before task timeout
`alignWithUa.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`chooseBestAlignment.realignMapq`|Int|60|Minimum mapping quality for a realignment to replace the original
`chooseBestAlignment.jobMemory`|Int|16|Memory allocated for this job (GB)
`chooseBestAlignment.cores`|Int|8|Number of cores requested for this job, also the samtools thread count
`chooseBestAlignment.timeout`|Int|6|Hours before task timeout
`chooseBestAlignment.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`alignLongHomopolymers.homopolymerLength`|Int|10|Minimum homopolymer length for a haplotype to be realigned
`alignLongHomopolymers.jobMemory`|Int|16|Memory allocated for this job (GB)
`alignLongHomopolymers.cores`|Int|8|Number of cores requested for this job
`alignLongHomopolymers.timeout`|Int|6|Hours before task timeout
`alignLongHomopolymers.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`matchReadsToHaplotypes.minMapq`|Int|5|Minimum mapping quality for a read to be matched
`matchReadsToHaplotypes.maxReadsPerWorkingArea`|Int|50000|Maximum number of reads matched per working area
`matchReadsToHaplotypes.jobMemory`|Int|12|Memory allocated for this job (GB)
`matchReadsToHaplotypes.cores`|Int|2|Number of cores requested for this job, also the sv_rematch thread count
`matchReadsToHaplotypes.timeout`|Int|12|Hours before task timeout
`matchReadsToHaplotypes.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`mergeMatchedReads.jobMemory`|Int|4|Memory allocated for this job (GB)
`mergeMatchedReads.cores`|Int|1|Number of cores requested for this job
`mergeMatchedReads.timeout`|Int|4|Hours before task timeout
`mergeMatchedReads.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`addReadSupportToAssembly.jobMemory`|Int|16|Memory allocated for this job (GB)
`addReadSupportToAssembly.cores`|Int|1|Number of cores requested for this job
`addReadSupportToAssembly.timeout`|Int|12|Hours before task timeout
`addReadSupportToAssembly.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`identifyVariants.gridssConfig`|String|"variantcalling.breakpointLowQuality = 120.0 variantcalling.breakendLowQuality = 120.0 variantcalling.writeFiltered = true variantcalling.requireReadPair = false variantcalling.requireSplitRead = false variantcalling.callUnassembledBreakends = true variantcalling.callUnassembledBreakpoints = true variantcalling.breakendMaxAssemblySupportBias = 100 variantcalling.breakpointMaxAssemblySupportBias = 100 variantcalling.includeSingleAssemblyFilter = false variantcalling.requiredReadAndAssemblyBreakpointOverlap = 5"|GRIDSS configuration as space-separated key = value settings
`identifyVariants.jobMemory`|Int|16|Memory allocated for this job (GB); the Java heap is set 2 GB below it
`identifyVariants.cores`|Int|2|Number of cores requested for this job
`identifyVariants.timeout`|Int|12|Hours before task timeout
`identifyVariants.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`prefilterCandidates.filterExpression`|String|"QUAL>120"|bcftools include expression a candidate must pass to be annotated
`prefilterCandidates.jobMemory`|Int|8|Memory allocated for this job (GB)
`prefilterCandidates.cores`|Int|1|Number of cores requested for this job
`prefilterCandidates.timeout`|Int|2|Hours before task timeout
`prefilterCandidates.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`annotateVariants.gridssConfig`|String|"variantcalling.breakpointLowQuality = 120.0 variantcalling.breakendLowQuality = 120.0 variantcalling.writeFiltered = true variantcalling.requireReadPair = false variantcalling.requireSplitRead = false variantcalling.callUnassembledBreakends = true variantcalling.callUnassembledBreakpoints = true variantcalling.breakendMaxAssemblySupportBias = 100 variantcalling.breakpointMaxAssemblySupportBias = 100 variantcalling.includeSingleAssemblyFilter = false variantcalling.requiredReadAndAssemblyBreakpointOverlap = 5"|GRIDSS configuration as space-separated key = value settings
`annotateVariants.germlineMemory`|Int|8|Memory allocated for this job in germline calling (GB); the Java heap is set 4 GB below it
`annotateVariants.somaticMemory`|Int|128|Memory allocated for this job in somatic calling (GB); the Java heap is set 4 GB below it
`annotateVariants.cores`|Int|1|Number of cores requested for this job
`annotateVariants.timeout`|Int|12|Hours before task timeout
`annotateVariants.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`mergeVcfs.jobMemory`|Int|10|Memory allocated for this job (GB); the Java heap is set 2 GB below it
`mergeVcfs.cores`|Int|1|Number of cores requested for this job
`mergeVcfs.timeout`|Int|2|Hours before task timeout
`mergeVcfs.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`linkGermlineBreakpoints.jobMemory`|Int|64|Memory allocated for this job (GB)
`linkGermlineBreakpoints.cores`|Int|1|Number of cores requested for this job
`linkGermlineBreakpoints.timeout`|Int|12|Hours before task timeout
`linkGermlineBreakpoints.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`filterSomatic.excludeFilters`|String|"SHORT_SR_NORMAL"|GRIPSS filters not to apply, separated by ';'
`filterSomatic.minNormalCoverage`|Int?|None|GRIPSS minimum normal coverage; GRIPSS default when unset
`filterSomatic.jobMemory`|Int|64|Memory allocated for this job (GB); the Java heap is set 4 GB below it
`filterSomatic.cores`|Int|1|Number of cores requested for this job
`filterSomatic.timeout`|Int|6|Hours before task timeout
`filterSomatic.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images
`convertVcfFormat.jobMemory`|Int|64|Memory allocated for this job (GB)
`convertVcfFormat.cores`|Int|8|Number of cores requested for this job, also the number of R worker processes
`convertVcfFormat.timeout`|Int|12|Hours before task timeout
`convertVcfFormat.modules`|String|"ugsv/1.37.1"|Environment modules to load; ugsv provides apptainer and the container images


### Outputs

Output | Type | Description | Labels
---|---|---|---
`svVcf`|File|Structural variant calls in breakend notation: linked germline calls, or GRIPSS-filtered somatic calls|vidarr_label: svVcf
`svVcfIndex`|File|Tabix index of svVcf|vidarr_label: svVcfIndex
`symbolicSvVcf`|File?|svVcf with symbolic alleles (<DEL>, <INS>, <DUP>, <INV>), when symbolicVcfFormat is true|vidarr_label: symbolicSvVcf
`symbolicSvVcfIndex`|File?|Tabix index of symbolicSvVcf|vidarr_label: symbolicSvVcfIndex
`annotatedVcf`|File|All GRIDSS candidate breakpoints with annotations, before linking or somatic filtering|vidarr_label: annotatedVcf
`annotatedVcfIndex`|File|Tabix index of annotatedVcf|vidarr_label: annotatedVcfIndex
`supportedAssemblyBam`|File|Realigned haplotypes annotated with the reads that support them, the assembly GRIDSS calls from|vidarr_label: supportedAssemblyBam
`supportedAssemblyBamIndex`|File|Index of supportedAssemblyBam|vidarr_label: supportedAssemblyBamIndex


## Commands
This section lists command(s) run by ultimaStructuralVariant workflow

* Running ultimaStructuralVariant

```
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    tumor=(~{sep=" " tumorCrams})
    germline=(~{sep=" " germlineCrams})
    for f in "${tumor[@]}" "${germline[@]}"; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # Distinct SM values over the read groups of all the given alignment files
    sample_names() {
      for f in "$@"; do
        run samtools view -H "$(readlink -f "$f")" \
          | awk -F'\t' '$1 == "@RG" { for (i = 2; i <= NF; i++) if ($i ~ /^SM:/) print substr($i, 4) }'
      done | sort -u
    }

    # Each group given must name exactly one sample; an empty group writes an empty name
    one_sample_name() {
      local group="$1"
      shift
      : > "${group}_sample_name.txt"
      [[ "$#" -gt 0 ]] || return 0
      sample_names "$@" > "${group}_names.txt"
      if [[ $(wc -l < "${group}_names.txt") -ne 1 ]]; then
        echo "Expected exactly one sample name across the ${group} files, found:" >&2
        cat "${group}_names.txt" >&2
        exit 1
      fi
      tr -d '\n' < "${group}_names.txt" > "${group}_sample_name.txt"
    }

    if [[ "${#tumor[@]}" -eq 0 && "${#germline[@]}" -eq 0 ]]; then
      echo "No input alignments: give tumorCrams, germlineCrams or both" >&2
      exit 1
    fi
    one_sample_name tumor "${tumor[@]}"
    one_sample_name germline "${germline[@]}"
    if [[ -s tumor_sample_name.txt && "$(cat tumor_sample_name.txt)" == "$(cat germline_sample_name.txt)" ]]; then
      echo "Tumor and germline files share the sample name $(cat tumor_sample_name.txt)" >&2
      exit 1
    fi

    # GRIDSS sample order: tumor first, then germline
    : > sample_names.txt
    for f in tumor_sample_name.txt germline_sample_name.txt; do
      if [[ -s "$f" ]]; then
        cat "$f" >> sample_names.txt
        echo >> sample_names.txt
      fi
    done
```
```
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    add_bind_root ~{intervalList}
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # With SCATTER_COUNT=1 IntervalListTools writes OUTPUT as a single file rather than
    # a directory of temp_NNNN_of_N subdirectories, so point it at that layout directly.
    mkdir out beds
    output_arg=out
    if [[ ~{numShards} -le 1 ]]; then
      mkdir out/temp_0001_of_1
      output_arg=out/temp_0001_of_1/scattered.interval_list
    fi

    run gatk --java-options "-Xmx~{jobMemory - 3}g" IntervalListTools \
      SCATTER_COUNT=~{numShards} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF=~{breakBandsAtMultiplesOf} \
      INPUT="$(readlink -f ~{intervalList})" \
      OUTPUT="${output_arg}"

    # Number the shards in scatter order so every BED has a unique name
    i=0
    for f in out/*/*.interval_list; do
      i=$((i + 1))
      grep -v '^@' "$f" | awk 'BEGIN { OFS = "\t" } { print $1, $2 - 1, $3 }' > "$(printf 'beds/shard_%04d.bed' "$i")"
    done
```
```
    set -euo pipefail

    tumor=(~{sep=" " tumorCrams})
    tumor_idx=(~{sep=" " tumorCramIndexes})
    germline=(~{sep=" " germlineCrams})
    germline_idx=(~{sep=" " germlineCramIndexes})

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    for f in "${tumor[@]}" "${tumor_idx[@]}" "${germline[@]}" "${germline_idx[@]}" ~{shardBed} ~{refFasta}; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # Cromwell may localize an alignment and its index to different directories; link
    # each pair side by side. The counter prefix keeps files with the same basename apart.
    mkdir alignments
    n=0
    link_group() {
      local -n files="$1" indexes="$2"
      local i name links=() index_links=()
      for i in "${!files[@]}"; do
        name="${n}_$(basename "${files[$i]}")"
        n=$((n + 1))
        ln -s "${files[$i]}" "alignments/${name}"
        ln -s "${indexes[$i]}" "alignments/${name}.${indexes[$i]##*.}"
        links+=("${PWD}/alignments/${name}")
        index_links+=("${PWD}/alignments/${name}.${indexes[$i]##*.}")
      done
      (IFS=,; echo "${links[*]}") > "alignments/$1.inputs"
      (IFS=,; echo "${index_links[*]}") > "alignments/$1.indexes"
    }

    # The assembler takes one comma-separated list per sample, samples separated by ';',
    # tumor first; tumor reads also switch on somatic active areas
    inputs=()
    indexes=()
    somatic_args=""
    if [[ "${#tumor[@]}" -gt 0 ]]; then
      link_group tumor tumor_idx
      inputs+=("$(cat alignments/tumor.inputs)")
      indexes+=("$(cat alignments/tumor.indexes)")
      somatic_args="--somatic-active-areas"
    fi
    if [[ "${#germline[@]}" -gt 0 ]]; then
      link_group germline germline_idx
      inputs+=("$(cat alignments/germline.inputs)")
      indexes+=("$(cat alignments/germline.indexes)")
    fi
    INPUT_STRING="$(IFS=';'; echo "${inputs[*]}")"
    INDEX_STRING="$(IFS=';'; echo "${indexes[*]}")"
    REF="$(readlink -f ~{refFasta})"
    BED="$(readlink -f ~{shardBed})"
    export INPUT_STRING INDEX_STRING REF BED somatic_args

    run bash -euo pipefail <<'CONTAINER'
    # shellcheck disable=SC2086
    tool \
      --input "${INPUT_STRING}" \
      --cram-index "${INDEX_STRING}" \
      ${somatic_args} \
      --output ~{outputPrefix} \
      --ref "${REF}" \
      --bed "${BED}" \
      --min-base ~{minBase} \
      --min-mapq ~{minMapq} \
      --median-coverage 1 \
      --background-median-coverage 1 \
      --max-reads-per-region ~{maxReadsPerPartition} \
      --min-feature-length '~{minIndelScSize}' \
      --min-mismatch-count '~{minMismatchCount}' \
      --max-num-haps ~{maxNumHaps} \
      --haps-max-sv-events 1 \
      --haps-hp-dedup true \
      --haps-read-witness true \
      --haps-weight-boost false \
      --haps-low-haps-rescue false \
      --haps-chain-graph-add false \
      --haps-allow-nonunique-ref false \
      --haps-min-path-cover false \
      --haps-merge-dangling false \
      --haps-contig-consensus false \
      --prog \
      --sv \
      --htslib \
      --realigned-sam
    samtools sort -o ~{outputPrefix}_hap_out_sorted.bam ~{outputPrefix}_hap_out.sam
    samtools index ~{outputPrefix}_hap_out_sorted.bam
    rm ~{outputPrefix}_hap_out.sam
    CONTAINER
```
```
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    bams=()
    for f in ~{sep=" " bams}; do
      bams+=("$(readlink -f "$f")")
      add_bind_root "${bams[-1]}"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    printf '%s\n' "${bams[@]}" > bams.txt
    run bash -euo pipefail <<'CONTAINER'
    samtools merge -c -@ ~{cores} -b bams.txt ~{outputPrefix}.bam
    samtools index ~{outputPrefix}.bam
    CONTAINER
```
```
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    bams=()
    for f in ~{sep=" " bams}; do
      bams+=("$(readlink -f "$f")")
      add_bind_root "${bams[-1]}"
    done
    add_bind_root ~{uaIndex}
    add_bind_root ~{refAlt}
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    printf '%s\n' "${bams[@]}" > bams.txt
    UA_INDEX="$(readlink -f ~{uaIndex})"
    REF_ALT="$(readlink -f ~{refAlt})"
    export UA_INDEX REF_ALT

    run bash -euo pipefail <<'CONTAINER'
    # shellcheck disable=SC2086
    samtools merge -@ ~{cores} -c -O SAM -b bams.txt - \
      | ua \
          --index "${UA_INDEX}" \
          --align true \
          --progress \
          --tp reference \
          --alt="${REF_ALT}" \
          --stat=~{outputPrefix}.%s.json \
          --nthread ~{cores} \
          --sam-input - \
          --sam-output - \
          ~{uaExtraArgs} \
      | samtools view -@ ~{cores} -o ~{outputPrefix}.bam -
    CONTAINER
```
```
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    bams=()
    for f in ~{sep=" " bams}; do
      bams+=("$(readlink -f "$f")")
      add_bind_root "${bams[-1]}"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    printf '%s\n' "${bams[@]}" > bams.txt

    # The chooser needs every source sorted by read name. The tie-breaker index points past
    # the two sources, as upstream does when the optional third (pangenome) source is off.
    run bash -euo pipefail <<'CONTAINER'
    sources=()
    i=0
    while read -r bam; do
      i=$((i + 1))
      samtools sort -@ ~{cores} -n -o "source_${i}.bam" "${bam}"
      sources+=(--alignment_source "source_${i}.bam")
    done < bams.txt
    python3 /opt/gridss/choose_best_haplotype_realignment.py \
      "${sources[@]}" \
      --output ~{outputPrefix}_unsorted.bam \
      --min_mapping_quality ~{realignMapq} \
      --tie_breaker_idx 2 \
      --overwrite_mapq 2,3,30
    samtools sort -@ ~{cores} -o ~{outputPrefix}.bam ~{outputPrefix}_unsorted.bam
    samtools index ~{outputPrefix}.bam
    rm source_*.bam ~{outputPrefix}_unsorted.bam
    CONTAINER
```
```
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    for f in ~{inputBam} ~{inputBamIndex} ~{refFasta} ~{refFai} ~{refDict}; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # The index is looked up next to the BAM, so pass the BAM by a path that has both
    ln -s ~{inputBam} input.bam
    ln -s ~{inputBamIndex} input.bam.bai

    run python3 /opt/gridss/align_long_homopolymers.py \
      --input "${PWD}/input.bam" \
      --output ~{outputPrefix}.bam \
      --reference "$(readlink -f ~{refFasta})" \
      --homopolymer_length ~{homopolymerLength}
```
```
    set -euo pipefail

    tumor=(~{sep=" " tumorCrams})
    tumor_idx=(~{sep=" " tumorCramIndexes})
    germline=(~{sep=" " germlineCrams})
    germline_idx=(~{sep=" " germlineCramIndexes})

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    for f in "${tumor[@]}" "${tumor_idx[@]}" "${germline[@]}" "${germline_idx[@]}" \
        ~{assemblyBam} ~{assemblyBamIndex} ~{shardBed} ~{refFasta}; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # Cromwell may localize an alignment and its index to different directories; link
    # each pair side by side. The counter prefix keeps files with the same basename apart.
    mkdir alignments
    n=0
    link_group() {
      local -n files="$1" indexes="$2"
      local i name links=()
      for i in "${!files[@]}"; do
        name="${n}_$(basename "${files[$i]}")"
        n=$((n + 1))
        ln -s "${files[$i]}" "alignments/${name}"
        ln -s "${indexes[$i]}" "alignments/${name}.${indexes[$i]##*.}"
        links+=("${PWD}/alignments/${name}")
      done
      (IFS=,; echo "${links[*]}") > "alignments/$1.inputs"
    }
    ln -s ~{assemblyBam} alignments/assembly.bam
    ln -s ~{assemblyBamIndex} alignments/assembly.bam.bai

    # One comma-separated list per sample, samples separated by ';', tumor first
    inputs=()
    if [[ "${#tumor[@]}" -gt 0 ]]; then
      link_group tumor tumor_idx
      inputs+=("$(cat alignments/tumor.inputs)")
    fi
    if [[ "${#germline[@]}" -gt 0 ]]; then
      link_group germline germline_idx
      inputs+=("$(cat alignments/germline.inputs)")
    fi

    run sv_rematch \
      -b "$(readlink -f ~{shardBed})" \
      -j ~{cores} -t -a -v \
      -l '~{minIndelScSize}' \
      -s '~{minMismatchCount}' \
      -m ~{minMapq} \
      -L ~{maxReadsPerWorkingArea} \
      "$(IFS=';'; echo "${inputs[*]}")" \
      "${PWD}/alignments/assembly.bam" \
      "$(readlink -f ~{refFasta})" \
      ~{outputPrefix}_rematched_hap.txt
```
```
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    : > shards.txt
    for f in ~{sep=" " matchedReads}; do
      readlink -f "$f" >> shards.txt
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # A read can be matched in more than one shard: keep its best-scoring haplotype, then
    # order the assignments by haplotype
    run bash -euo pipefail <<'CONTAINER'
    export LC_ALL=C
    mapfile -t shards < shards.txt
    sort -T . -m -k1,1 -k3,3rn "${shards[@]}" \
      | awk '$1 != prev { print; prev = $1 }' \
      | sort -T . -k2,2 > ~{outputPrefix}.all_matches_by_hap.txt
    CONTAINER
```
```
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    for f in ~{assemblyBam} ~{matchedReads} ~{refFasta}; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    ASSEMBLY="$(readlink -f ~{assemblyBam})"
    MATCHED="$(readlink -f ~{matchedReads})"
    REF="$(readlink -f ~{refFasta})"
    export ASSEMBLY MATCHED REF

    run bash -euo pipefail <<'CONTAINER'
    samtools sort -N -o assembly_name_sorted.bam "${ASSEMBLY}"
    sv_rematch -v -S -M "${MATCHED}" \
      assembly_name_sorted.bam \
      "${REF}" \
      ~{outputPrefix}_assembly.support.name_sorted.bam
    samtools sort -o ~{outputPrefix}_assembly.support.bam ~{outputPrefix}_assembly.support.name_sorted.bam
    samtools index ~{outputPrefix}_assembly.support.bam
    rm assembly_name_sorted.bam ~{outputPrefix}_assembly.support.name_sorted.bam
    CONTAINER
```
```
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    for f in ~{assemblyBam} ~{assemblyBamIndex} ~{refFasta} ~{refFai} ~{refDict} ~{blacklistBed}; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # GRIDSS writes working files next to its inputs and a cache next to the reference, so
    # give it links in the task directory rather than paths in read-only locations
    mkdir ref tmp
    ln -s "$(readlink -f ~{refFasta})" ref/
    ln -s "$(readlink -f ~{refFai})" ref/
    ln -s "$(readlink -f ~{refDict})" ref/
    ln -s ~{assemblyBam} assembly.bam
    ln -s ~{assemblyBamIndex} assembly.bam.bai
    REF="${PWD}/ref/$(basename "$(readlink -f ~{refFasta})")"
    BLACKLIST="$(readlink -f ~{blacklistBed})"
    export REF BLACKLIST

    run bash -euo pipefail <<'CONTAINER'
    python3 - <<'PYTHON'
    import re
    text = "~{gridssConfig}"
    pattern = r'(\w+(?:\.\w+)+\s*=\s*(?:\d+(?:\.\d+)?|true|false|\w+))'
    with open("gridss.config", "w") as fh:
        fh.write("\n".join(re.findall(pattern, text)))
    PYTHON
    java -Xmx~{jobMemory - 2}g -Djava.io.tmpdir="${PWD}/tmp" \
      -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.IdentifyVariants \
      ~{sep=" " prefix("SAMPLE_NAMES=", sampleNames)} \
      R="${REF}" \
      O=~{outputPrefix}.vcf \
      BLACKLIST="${BLACKLIST}" \
      ASSEMBLY=assembly.bam \
      C=gridss.config
    bcftools view ~{outputPrefix}.vcf -Oz -o ~{outputPrefix}.vcf.gz
    bcftools index -t ~{outputPrefix}.vcf.gz
    CONTAINER
```
```
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    add_bind_root ~{inputVcf}
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    VCF="$(readlink -f ~{inputVcf})"
    export VCF

    run bash -euo pipefail <<'CONTAINER'
    bcftools view -Oz -i "~{filterExpression}" -o ~{outputPrefix}.vcf.gz "${VCF}"
    bcftools index -t ~{outputPrefix}.vcf.gz
    CONTAINER
```
```
    set -euo pipefail

    tumor=(~{sep=" " tumorCrams})
    tumor_idx=(~{sep=" " tumorCramIndexes})
    germline=(~{sep=" " germlineCrams})
    germline_idx=(~{sep=" " germlineCramIndexes})

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    for f in "${tumor[@]}" "${tumor_idx[@]}" "${germline[@]}" "${germline_idx[@]}" \
        ~{inputVcf} ~{inputVcfIndex} ~{assemblyBam} ~{assemblyBamIndex} ~{shardBed} \
        ~{refFasta} ~{refFai} ~{refDict} ~{blacklistBed}; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # Cromwell may localize a file and its index to different directories; link each pair
    # side by side. The counter prefix keeps files with the same basename apart.
    mkdir alignments ref tmp
    n=0
    link_group() {
      local -n files="$1" indexes="$2"
      local i name args=""
      for i in "${!files[@]}"; do
        name="${n}_$(basename "${files[$i]}")"
        n=$((n + 1))
        ln -s "${files[$i]}" "alignments/${name}"
        ln -s "${indexes[$i]}" "alignments/${name}.${indexes[$i]##*.}"
        args+=" -I ${PWD}/alignments/${name}"
      done
      echo "${args}"
    }
    TUMOR_ARGS="$(link_group tumor tumor_idx)"
    GERMLINE_ARGS="$(link_group germline germline_idx)"
    ln -s ~{assemblyBam} alignments/assembly.bam
    ln -s ~{assemblyBamIndex} alignments/assembly.bam.bai
    ln -s ~{inputVcf} candidates.vcf.gz
    ln -s ~{inputVcfIndex} candidates.vcf.gz.tbi

    # GRIDSS writes working files next to its inputs and a cache next to the reference, so
    # give it links in the task directory rather than paths in read-only locations
    ln -s "$(readlink -f ~{refFasta})" ref/
    ln -s "$(readlink -f ~{refFai})" ref/
    ln -s "$(readlink -f ~{refDict})" ref/
    REF="${PWD}/ref/$(basename "$(readlink -f ~{refFasta})")"
    BED="$(readlink -f ~{shardBed})"
    BLACKLIST="$(readlink -f ~{blacklistBed})"
    export TUMOR_ARGS GERMLINE_ARGS REF BED BLACKLIST

    # Annotate this shard's candidates from the reads and assembly in the shard only: GRIDSS
    # reads every alignment it is given in full to collect its metrics
    run bash -euo pipefail <<'CONTAINER'
    python3 - <<'PYTHON'
    import re
    text = "~{gridssConfig}"
    pattern = r'(\w+(?:\.\w+)+\s*=\s*(?:\d+(?:\.\d+)?|true|false|\w+))'
    with open("gridss.config", "w") as fh:
        fh.write("\n".join(re.findall(pattern, text)))
    PYTHON

    bcftools view -R "${BED}" -o shard.vcf candidates.vcf.gz

    subset() {
      local out="$1"
      shift
      # shellcheck disable=SC2068
      gatk --java-options "-Xmx4g -Djava.io.tmpdir=${PWD}/tmp" PrintReads $@ \
        -O /dev/stdout -L "${BED}" -R "${REF}" \
        | samtools view -C -T "${REF}" -o "${out}" --output-fmt-option embed_ref=1 -
      samtools index "${out}"
    }
    inputs=()
    if [[ -n "${TUMOR_ARGS}" ]]; then
      subset input_tumor.cram ${TUMOR_ARGS}
      inputs+=(I=input_tumor.cram)
    fi
    if [[ -n "${GERMLINE_ARGS}" ]]; then
      subset input_germline.cram ${GERMLINE_ARGS}
      inputs+=(I=input_germline.cram)
    fi
    subset assembly_shard.cram -I alignments/assembly.bam

    java -Xmx~{jobMemory - 4}g -Xms~{floor(jobMemory / 2)}g -Djava.io.tmpdir="${PWD}/tmp" \
      -cp /opt/gridss/gridss--gridss-jar-with-dependencies.jar gridss.AnnotateVariants \
      INPUT_VCF=shard.vcf \
      R="${REF}" \
      "${inputs[@]}" \
      BLACKLIST="${BLACKLIST}" \
      C=gridss.config \
      ASSEMBLY=assembly_shard.cram \
      OUTPUT_VCF=~{outputPrefix}.ann.vcf \
      THREADS=1
    rm -f input_tumor.cram input_tumor.cram.crai input_germline.cram input_germline.cram.crai \
      assembly_shard.cram assembly_shard.cram.crai
    CONTAINER
```
```
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    inputs=()
    for f in ~{sep=" " vcfs}; do
      inputs+=("INPUT=$(readlink -f "$f")")
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    run gatk --java-options "-Xmx~{jobMemory - 2}g" MergeVcfs \
      "${inputs[@]}" \
      OUTPUT=~{outputFileName}
```
```
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    for f in ~{inputVcf} ~{inputVcfIndex} ~{refFasta} ~{refFai} ~{refDict}; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # Writes <fulloutput>.bgz and its .tbi
    run Rscript /opt/gridss/link_breakpoints.R \
      --input "$(readlink -f ~{inputVcf})" \
      --fulloutput ~{outputPrefix}_linked.vcf \
      --ref "$(readlink -f ~{refFasta})" \
      --scriptdir /opt/gridss/
```
```
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    for f in ~{inputVcf} ~{inputVcfIndex} ~{refFasta} ~{refFai} ~{refDict} \
        ~{ponSglFile} ~{ponSvFile} ~{knownHotspotFile} ~{repeatMaskFile}; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    # The index is looked up next to the VCF, so pass the VCF by a path that has both
    ln -s ~{inputVcf} input.vcf.gz
    ln -s ~{inputVcfIndex} input.vcf.gz.tbi

    normal_args=()
    if [[ -n "~{germlineSampleName}" ]]; then
      normal_args=(-reference "~{germlineSampleName}")
    fi

    mkdir gripss_output
    run java -Xmx~{jobMemory - 4}g -jar /opt/wtsi-cgp/java/gripss.jar \
      -vcf "${PWD}/input.vcf.gz" \
      -sample "~{tumorSampleName}" \
      "${normal_args[@]}" \
      -ref_genome_version ~{refGenomeVersion} \
      -ref_genome "$(readlink -f ~{refFasta})" \
      -pon_sgl_file "$(readlink -f ~{ponSglFile})" \
      -pon_sv_file "$(readlink -f ~{ponSvFile})" \
      -known_hotspot_file "$(readlink -f ~{knownHotspotFile})" \
      -repeat_mask_file "$(readlink -f ~{repeatMaskFile})" \
      ~{"-min_normal_coverage " + minNormalCoverage} \
      -exclude_filters "~{excludeFilters}" \
      -output_dir gripss_output

    mv "gripss_output/~{tumorSampleName}.gripss.vcf.gz" ~{outputPrefix}.gripss.vcf.gz
    mv "gripss_output/~{tumorSampleName}.gripss.vcf.gz.tbi" ~{outputPrefix}.gripss.vcf.gz.tbi
```
```
    set -euo pipefail

    # Bind the filesystem root of every path the tools open: the task directory and
    # wherever each input resolves to. Roots rather than directories, because the
    # container runtime skips a bind nested inside another bind.
    bind_args=(--pwd "$PWD")
    declare -A roots=()
    add_bind_root() {
      [[ -e "$1" ]] || { echo "not found: $1" >&2; exit 1; }
      local root
      root="/$(readlink -f "$1" | cut -d/ -f2)"
      if [[ -z "${roots[$root]:-}" ]]; then
        roots[$root]=1
        bind_args+=(-B "$root")
      fi
    }
    add_bind_root "$PWD"
    for f in ~{inputVcf} ~{inputVcfIndex} ~{refFasta} ~{refFai} ~{refDict}; do
      add_bind_root "$f"
    done
    run() { apptainer exec "${bind_args[@]}" ~{image} "$@"; }

    # Without squashfuse, apptainer unpacks a SIF into a temporary sandbox on every exec and
    # can fail to delete it. Keep those sandboxes in the task directory and remove them on exit,
    # without letting a failed cleanup change the task's exit status.
    export APPTAINER_TMPDIR="${PWD}/apptainer_tmp"
    mkdir -p "${APPTAINER_TMPDIR}"
    trap 'chmod -R u+w "${APPTAINER_TMPDIR}" 2> /dev/null || true; rm -rf "${APPTAINER_TMPDIR}" || true' EXIT

    VCF="$(readlink -f ~{inputVcf})"
    REF="$(readlink -f ~{refFasta})"
    export VCF REF

    # The R script writes <output_vcf>.bgz
    run bash -euo pipefail <<'CONTAINER'
    Rscript /opt/gridss/convert_vcf_format.R \
      --input_vcf "${VCF}" \
      --output_vcf ~{outputPrefix}.vcf \
      --reference "${REF}" \
      --n_jobs ~{cores}
    bcftools view -Oz -o ~{outputPrefix}.vcf.gz ~{outputPrefix}.vcf.bgz
    bcftools index -t ~{outputPrefix}.vcf.gz
    rm -f ~{outputPrefix}.vcf.bgz ~{outputPrefix}.vcf.bgz.tbi
    CONTAINER
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
