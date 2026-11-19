version 1.0

# Adapted from the Ultima Genomics structural_variant_pipeline workflow (healthomics-workflows
# v1.37.1, Apache License 2.0). Tools run from the upstream container images through Apptainer.
# Paths are rooted at environment variables set by the modules each task loads:
# HG38_ULTIMA_ROOT (hg38-ultima) for the reference, HG38_UGSV_RESOURCES_ROOT
# (hg38-ugsv-resources) for the calling regions, blacklist and GRIPSS panels,
# HG38_UA_INDEX_ROOT (hg38-ua-index) for the UA aligner index, UGSV_ROOT (ugsv) for the
# container images.

struct GenomeResources {
  String genomeModule
  String refFasta
  String refFai
  String refDict
  String refAlt
  String refGenomeVersion
  String uaIndex
  String intervalList
  String blacklistBed
  String ponSglFile
  String ponSvFile
  String knownHotspotFile
  String repeatMaskFile
}

struct ModeParameters {
  String minIndelScSize
  String minMismatchCount
  Int numShards
}

workflow ultimaStructuralVariant {
  input {
    Array[File] tumorCrams = []
    Array[File] tumorCramIndexes = []
    Array[File] germlineCrams = []
    Array[File] germlineCramIndexes = []
    String outputFileNamePrefix
    String reference = "hg38"
    String? intervalList
    Int? numShards
    Int scatterIntervalsBreak = 10000000
    Boolean symbolicVcfFormat = true
  }

  parameter_meta {
    tumorCrams: "Tumor CRAM files; all must carry the same read-group sample name. Non-empty selects somatic calling"
    tumorCramIndexes: "Indexes of tumorCrams, in the same order"
    germlineCrams: "CRAM files of the germline sample, or of the matched normal in somatic calling (optional there, for tumor-only); all must carry the same read-group sample name"
    germlineCramIndexes: "Indexes of germlineCrams, in the same order"
    outputFileNamePrefix: "Prefix for output file names"
    reference: "Reference genome build, selects the genome resources: hg38"
    intervalList: "Picard interval list of the regions to call, to restrict a run such as a test; defaults to the whole-genome calling regions of the reference"
    numShards: "Number of genomic shards to process in parallel; defaults to 50 for germline and 200 for somatic calling"
    scatterIntervalsBreak: "Break scatter intervals at multiples of this many bases"
    symbolicVcfFormat: "Also write the calls with symbolic alleles (<DEL>, <INS>, <DUP>, <INV>) in place of breakend notation"
  }

  Boolean isSomatic = length(tumorCrams) > 0

  String gatkImage = "$UGSV_ROOT/images/broadinstitute_gatk_4.6.0.0.sif"
  String makeExamplesImage = "$UGSV_ROOT/images/ultimagenomics_make_examples_3.5.0.sif"
  String alignmentImage = "$UGSV_ROOT/images/ultimagenomics_alignment_4.1.4.sif"
  String rematcherImage = "$UGSV_ROOT/images/ultimagenomics_rematcher_1.2.0_d0f0b0a.sif"
  String gridssImage = "$UGSV_ROOT/images/ultimagenomics_gridss_1.0.2.sif"
  String gripssImage = "$UGSV_ROOT/images/ultimagenomics_gripss_ug_2.4.1_165b492.sif"

  Map[String, GenomeResources] genomeResources = {
    "hg38": {
      "genomeModule": "hg38-ultima/v0 hg38-ugsv-resources/1.37.1 hg38-ua-index/v45-79372c0",
      "refFasta": "$HG38_ULTIMA_ROOT/Homo_sapiens_assembly38.fasta",
      "refFai": "$HG38_ULTIMA_ROOT/Homo_sapiens_assembly38.fasta.fai",
      "refDict": "$HG38_ULTIMA_ROOT/Homo_sapiens_assembly38.dict",
      "refAlt": "$HG38_UA_INDEX_ROOT/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt",
      "refGenomeVersion": "38",
      "uaIndex": "$HG38_UA_INDEX_ROOT/hg38/UA/b38-v45-79372c0.uai",
      "intervalList": "$HG38_UGSV_RESOURCES_ROOT/sv/wgs_calling_regions.hg38_no_cytoBandIdeo_acen.interval_list",
      "blacklistBed": "$HG38_UGSV_RESOURCES_ROOT/sv/gridss/ENCFF356LFX.bed",
      "ponSglFile": "$HG38_UGSV_RESOURCES_ROOT/sv/gripss/sgl_pon.38.bed.gz",
      "ponSvFile": "$HG38_UGSV_RESOURCES_ROOT/sv/gripss/sv_pon.38.bedpe.gz",
      "knownHotspotFile": "$HG38_UGSV_RESOURCES_ROOT/sv/gripss/known_fusions.38.sorted.bedpe",
      "repeatMaskFile": "$HG38_UGSV_RESOURCES_ROOT/sv/gripss/repeat_mask_data.38.fa.gz"
    }
  }

  # Read-selection thresholds and shard counts of the upstream germline and somatic templates.
  # In somatic mode a value "x;y" applies x to the tumor and y to the normal.
  Map[String, ModeParameters] modeParameters = {
    "germline": {
      "minIndelScSize": "10",
      "minMismatchCount": "5",
      "numShards": 50
    },
    "somatic": {
      "minIndelScSize": "20;0",
      "minMismatchCount": "5;0",
      "numShards": 200
    }
  }

  GenomeResources genome = genomeResources[reference]
  ModeParameters mode = modeParameters[if isSomatic then "somatic" else "germline"]

  call getSampleNames {
    input:
      tumorCrams = tumorCrams,
      germlineCrams = germlineCrams,
      image = alignmentImage
  }

  call scatterIntervals {
    input:
      intervalList = select_first([intervalList, genome.intervalList]),
      numShards = select_first([numShards, mode.numShards]),
      breakBandsAtMultiplesOf = scatterIntervalsBreak,
      genomeModule = genome.genomeModule,
      image = gatkImage
  }

  scatter (shardBed in scatterIntervals.shardBeds) {
    call createAssembly {
      input:
        tumorCrams = tumorCrams,
        tumorCramIndexes = tumorCramIndexes,
        germlineCrams = germlineCrams,
        germlineCramIndexes = germlineCramIndexes,
        shardBed = shardBed,
        outputPrefix = outputFileNamePrefix + "_assembly." + basename(shardBed, ".bed"),
        minIndelScSize = mode.minIndelScSize,
        minMismatchCount = mode.minMismatchCount,
        genomeModule = genome.genomeModule,
        refFasta = genome.refFasta,
        image = makeExamplesImage
    }
  }

  call mergeBams {
    input:
      bams = createAssembly.assembly,
      outputPrefix = outputFileNamePrefix + "_merged_assembly",
      image = gridssImage
  }

  call alignWithUa {
    input:
      bams = createAssembly.assembly,
      outputPrefix = outputFileNamePrefix + "_assembly_file_ua_aligned",
      genomeModule = genome.genomeModule,
      uaIndex = genome.uaIndex,
      refAlt = genome.refAlt,
      image = alignmentImage
  }

  call chooseBestAlignment {
    input:
      bams = [mergeBams.bam, alignWithUa.bam],
      outputPrefix = outputFileNamePrefix + "_assembly_realigned",
      image = gridssImage
  }

  call alignLongHomopolymers {
    input:
      inputBam = chooseBestAlignment.bam,
      inputBamIndex = chooseBestAlignment.bamIndex,
      outputPrefix = outputFileNamePrefix + "_assembly_realigned_long_homopolymers_aligned",
      genomeModule = genome.genomeModule,
      refFasta = genome.refFasta,
      refFai = genome.refFai,
      refDict = genome.refDict,
      image = gridssImage
  }

  scatter (shardBed in scatterIntervals.shardBeds) {
    call matchReadsToHaplotypes {
      input:
        tumorCrams = tumorCrams,
        tumorCramIndexes = tumorCramIndexes,
        germlineCrams = germlineCrams,
        germlineCramIndexes = germlineCramIndexes,
        assemblyBam = alignLongHomopolymers.bam,
        assemblyBamIndex = alignLongHomopolymers.bamIndex,
        shardBed = shardBed,
        outputPrefix = outputFileNamePrefix + "." + basename(shardBed, ".bed"),
        minIndelScSize = mode.minIndelScSize,
        minMismatchCount = mode.minMismatchCount,
        genomeModule = genome.genomeModule,
        refFasta = genome.refFasta,
        image = rematcherImage
    }
  }

  call mergeMatchedReads {
    input:
      matchedReads = matchReadsToHaplotypes.readHaplotypeAssignment,
      outputPrefix = outputFileNamePrefix,
      image = alignmentImage
  }

  call addReadSupportToAssembly {
    input:
      assemblyBam = alignLongHomopolymers.bam,
      matchedReads = mergeMatchedReads.readHaplotypeAssignment,
      outputPrefix = outputFileNamePrefix,
      genomeModule = genome.genomeModule,
      refFasta = genome.refFasta,
      image = rematcherImage
  }

  call identifyVariants {
    input:
      sampleNames = getSampleNames.sampleNames,
      assemblyBam = addReadSupportToAssembly.bam,
      assemblyBamIndex = addReadSupportToAssembly.bamIndex,
      outputPrefix = outputFileNamePrefix,
      genomeModule = genome.genomeModule,
      refFasta = genome.refFasta,
      refFai = genome.refFai,
      refDict = genome.refDict,
      blacklistBed = genome.blacklistBed,
      image = gridssImage
  }

  call prefilterCandidates {
    input:
      inputVcf = identifyVariants.vcf,
      inputVcfIndex = identifyVariants.vcfIndex,
      outputPrefix = outputFileNamePrefix + ".prefiltered",
      image = gridssImage
  }

  scatter (shardBed in scatterIntervals.shardBeds) {
    call annotateVariants {
      input:
        inputVcf = prefilterCandidates.vcf,
        inputVcfIndex = prefilterCandidates.vcfIndex,
        tumorCrams = tumorCrams,
        tumorCramIndexes = tumorCramIndexes,
        germlineCrams = germlineCrams,
        germlineCramIndexes = germlineCramIndexes,
        assemblyBam = addReadSupportToAssembly.bam,
        assemblyBamIndex = addReadSupportToAssembly.bamIndex,
        shardBed = shardBed,
        outputPrefix = outputFileNamePrefix + "." + basename(shardBed, ".bed"),
        isSomatic = isSomatic,
        genomeModule = genome.genomeModule,
        refFasta = genome.refFasta,
        refFai = genome.refFai,
        refDict = genome.refDict,
        blacklistBed = genome.blacklistBed,
        image = gridssImage
    }
  }

  call mergeVcfs {
    input:
      vcfs = annotateVariants.annotatedVcf,
      outputFileName = outputFileNamePrefix + ".ann.vcf.gz",
      image = gatkImage
  }

  if (!isSomatic) {
    call linkGermlineBreakpoints {
      input:
        inputVcf = mergeVcfs.vcf,
        inputVcfIndex = mergeVcfs.vcfIndex,
        outputPrefix = outputFileNamePrefix,
        genomeModule = genome.genomeModule,
        refFasta = genome.refFasta,
        refFai = genome.refFai,
        refDict = genome.refDict,
        image = gridssImage
    }
  }

  if (isSomatic) {
    call filterSomatic {
      input:
        inputVcf = mergeVcfs.vcf,
        inputVcfIndex = mergeVcfs.vcfIndex,
        tumorSampleName = getSampleNames.tumorSampleName,
        germlineSampleName = getSampleNames.germlineSampleName,
        outputPrefix = outputFileNamePrefix,
        genomeModule = genome.genomeModule,
        refFasta = genome.refFasta,
        refFai = genome.refFai,
        refDict = genome.refDict,
        refGenomeVersion = genome.refGenomeVersion,
        ponSglFile = genome.ponSglFile,
        ponSvFile = genome.ponSvFile,
        knownHotspotFile = genome.knownHotspotFile,
        repeatMaskFile = genome.repeatMaskFile,
        image = gripssImage
    }
  }

  File callsVcf = select_first([linkGermlineBreakpoints.vcf, filterSomatic.vcf])
  File callsVcfIndex = select_first([linkGermlineBreakpoints.vcfIndex, filterSomatic.vcfIndex])

  if (symbolicVcfFormat) {
    call convertVcfFormat {
      input:
        inputVcf = callsVcf,
        inputVcfIndex = callsVcfIndex,
        outputPrefix = outputFileNamePrefix + ".symbolic",
        genomeModule = genome.genomeModule,
        refFasta = genome.refFasta,
        refFai = genome.refFai,
        refDict = genome.refDict,
        image = gridssImage
    }
  }

  output {
    File svVcf = callsVcf
    File svVcfIndex = callsVcfIndex
    File? symbolicSvVcf = convertVcfFormat.vcf
    File? symbolicSvVcfIndex = convertVcfFormat.vcfIndex
    File annotatedVcf = mergeVcfs.vcf
    File annotatedVcfIndex = mergeVcfs.vcfIndex
    File supportedAssemblyBam = addReadSupportToAssembly.bam
    File supportedAssemblyBamIndex = addReadSupportToAssembly.bamIndex
  }

  meta {
    author: "Gavin Peng"
    email: "gpeng@oicr.on.ca"
    description: "Structural variant calling for Ultima Genomics data, germline or somatic. The genome is split into shards; per shard, reads with soft clips, indels or mismatches are locally assembled into haplotypes. The haplotypes are realigned with the UA aligner and the better of the original and UA alignment is kept per haplotype, homopolymer-rich haplotypes are realigned, and every read is matched to the haplotype it supports. GRIDSS calls breakpoints from the supported assembly and annotates them shard by shard. Germline calls are then linked into events; somatic calls are filtered by GRIPSS against panels of normals. Adapted from the Ultima Genomics [structural_variant_pipeline](https://github.com/Ultimagen/healthomics-workflows/tree/main/workflows/structural_variant_pipeline) workflow (v1.37.1).\n\nSomatic calling runs when `tumorCrams` is given, with `germlineCrams` as the matched normal (optional, for tumor-only); otherwise `germlineCrams` is called as a germline sample. Read-selection thresholds and the default shard count follow the upstream germline and somatic templates. The optional pangenome (Giraffe) realignment of the upstream workflow is not included.\n\nEach task runs its tools from the upstream container images through Apptainer, so the Cromwell backend needs no container support.\n\n## Resources\n\nEverything the workflow reads besides the input alignments comes from four modules, which tasks load through their `modules` runtime attribute:\n\n* `hg38-ultima/v0` sets `HG38_ULTIMA_ROOT`: the hg38 FASTA with its `.fai` and `.dict`.\n* `hg38-ugsv-resources/1.37.1` sets `HG38_UGSV_RESOURCES_ROOT`: the calling regions and GRIDSS blacklist under `sv/`, the GRIPSS panels of normals, known fusions and repeat mask under `sv/gripss/`, laid out as in the Ultima public bucket.\n* `hg38-ua-index/v45-79372c0` sets `HG38_UA_INDEX_ROOT`: the UA aligner index and the alt-contig list.\n* `ugsv/1.37.1` sets `UGSV_ROOT` and loads apptainer: the container images, pulled by digest, under `images/`."
    dependencies: [
      {
        name: "apptainer/1.4.5",
        url: "https://apptainer.org"
      },
      {
        name: "hg38-ultima/v0",
        url: "https://github.com/Ultimagen/healthomics-workflows"
      },
      {
        name: "hg38-ugsv-resources/1.37.1",
        url: "https://github.com/Ultimagen/healthomics-workflows"
      },
      {
        name: "hg38-ua-index/v45-79372c0",
        url: "https://github.com/Ultimagen/healthomics-workflows"
      },
      {
        name: "ultimagenomics/make_examples:3.5.0",
        url: "https://hub.docker.com/r/ultimagenomics/make_examples"
      },
      {
        name: "ultimagenomics/alignment:4.1.4",
        url: "https://hub.docker.com/r/ultimagenomics/alignment"
      },
      {
        name: "ultimagenomics/rematcher:1.2.0_d0f0b0a",
        url: "https://hub.docker.com/r/ultimagenomics/rematcher"
      },
      {
        name: "ultimagenomics/gridss:1.0.2",
        url: "https://hub.docker.com/r/ultimagenomics/gridss"
      },
      {
        name: "ultimagenomics/gripss:ug_2.4.1_165b492",
        url: "https://hub.docker.com/r/ultimagenomics/gripss"
      },
      {
        name: "broadinstitute/gatk:4.6.0.0",
        url: "https://hub.docker.com/r/broadinstitute/gatk"
      }
    ]
    output_meta: {
      svVcf: {
        description: "Structural variant calls in breakend notation: linked germline calls, or GRIPSS-filtered somatic calls",
        vidarr_label: "svVcf"
      },
      svVcfIndex: {
        description: "Tabix index of svVcf",
        vidarr_label: "svVcfIndex"
      },
      symbolicSvVcf: {
        description: "svVcf with symbolic alleles (<DEL>, <INS>, <DUP>, <INV>), when symbolicVcfFormat is true",
        vidarr_label: "symbolicSvVcf"
      },
      symbolicSvVcfIndex: {
        description: "Tabix index of symbolicSvVcf",
        vidarr_label: "symbolicSvVcfIndex"
      },
      annotatedVcf: {
        description: "All GRIDSS candidate breakpoints with annotations, before linking or somatic filtering",
        vidarr_label: "annotatedVcf"
      },
      annotatedVcfIndex: {
        description: "Tabix index of annotatedVcf",
        vidarr_label: "annotatedVcfIndex"
      },
      supportedAssemblyBam: {
        description: "Realigned haplotypes annotated with the reads that support them, the assembly GRIDSS calls from",
        vidarr_label: "supportedAssemblyBam"
      },
      supportedAssemblyBamIndex: {
        description: "Index of supportedAssemblyBam",
        vidarr_label: "supportedAssemblyBamIndex"
      }
    }
  }
}

task getSampleNames {
  input {
    Array[File] tumorCrams
    Array[File] germlineCrams
    String image
    Int jobMemory = 2
    Int cores = 1
    Int timeout = 1
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    tumorCrams: "Tumor CRAM files"
    germlineCrams: "Germline or matched-normal CRAM files"
    image: "Apptainer image providing samtools"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  output {
    String tumorSampleName = read_string("tumor_sample_name.txt")
    String germlineSampleName = read_string("germline_sample_name.txt")
    Array[String] sampleNames = read_lines("sample_names.txt")
  }
}

task scatterIntervals {
  input {
    String intervalList
    Int numShards
    Int breakBandsAtMultiplesOf
    String genomeModule
    String image
    Int jobMemory = 8
    Int cores = 1
    Int timeout = 1
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    intervalList: "Picard interval list of the regions to call"
    numShards: "Number of shards to split the intervals into"
    breakBandsAtMultiplesOf: "Break intervals at multiples of this many bases"
    genomeModule: "Modules providing the reference genome and its resource files"
    image: "Apptainer image providing GATK"
    jobMemory: "Memory allocated for this job (GB); the Java heap is set 3 GB below it"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    Array[File] shardBeds = glob("beds/*.bed")
  }
}

task createAssembly {
  input {
    Array[File] tumorCrams
    Array[File] tumorCramIndexes
    Array[File] germlineCrams
    Array[File] germlineCramIndexes
    File shardBed
    String outputPrefix
    String minIndelScSize
    String minMismatchCount
    String genomeModule
    String refFasta
    String image
    Int minBase = 0
    Int minMapq = 5
    Int maxNumHaps = 40
    Int maxReadsPerPartition = 1500
    Int jobMemory = 4
    Int cores = 1
    Int timeout = 12
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    tumorCrams: "Tumor CRAM files"
    tumorCramIndexes: "Indexes of tumorCrams, in the same order"
    germlineCrams: "Germline or matched-normal CRAM files"
    germlineCramIndexes: "Indexes of germlineCrams, in the same order"
    shardBed: "BED file of the regions in this shard"
    outputPrefix: "Prefix for output file names"
    minIndelScSize: "Minimum indel or soft-clip length for a read to enter assembly; x;y gives tumor and normal values"
    minMismatchCount: "Minimum mismatch count for a read to enter assembly; x;y gives tumor and normal values"
    genomeModule: "Modules providing the reference genome and its resource files"
    refFasta: "Reference FASTA; its .fai must sit next to it"
    image: "Apptainer image providing the haplotype assembler (tool) and samtools"
    minBase: "Minimum base quality for a base to count toward assembly"
    minMapq: "Minimum mapping quality for a read to enter assembly"
    maxNumHaps: "Maximum number of haplotypes assembled per region"
    maxReadsPerPartition: "Maximum number of reads assembled per region"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File assembly = "~{outputPrefix}_hap_out_sorted.bam"
    File assemblyIndex = "~{outputPrefix}_hap_out_sorted.bam.bai"
  }
}

task mergeBams {
  input {
    Array[File] bams
    String outputPrefix
    String image
    Int jobMemory = 16
    Int cores = 8
    Int timeout = 4
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    bams: "Coordinate-sorted BAM files to merge"
    outputPrefix: "Prefix for output file names"
    image: "Apptainer image providing samtools"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job, also the samtools thread count"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  output {
    File bam = "~{outputPrefix}.bam"
    File bamIndex = "~{outputPrefix}.bam.bai"
  }
}

task alignWithUa {
  input {
    Array[File] bams
    String outputPrefix
    String genomeModule
    String uaIndex
    String refAlt
    String image
    String uaExtraArgs = "--seed-score-ratio 0.5 --vector --soft-clipping --realignment-tag re --mismatch-cost -6"
    Int jobMemory = 48
    Int cores = 16
    Int timeout = 12
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    bams: "Assembled-haplotype BAM files to realign"
    outputPrefix: "Prefix for output file names"
    genomeModule: "Modules providing the reference genome and its resource files"
    uaIndex: "UA aligner index of the reference"
    refAlt: "Alt-contig list of the reference, for ua --alt"
    image: "Apptainer image providing ua and samtools"
    uaExtraArgs: "Further ua arguments"
    jobMemory: "Memory allocated for this job (GB); ua holds the whole index in memory"
    cores: "Number of cores requested for this job, also the ua and samtools thread count"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File bam = "~{outputPrefix}.bam"
  }
}

task chooseBestAlignment {
  input {
    Array[File] bams
    String outputPrefix
    String image
    Int realignMapq = 60
    Int jobMemory = 16
    Int cores = 8
    Int timeout = 6
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    bams: "Alignments of the same haplotypes to choose between: the original, then the UA realignment"
    outputPrefix: "Prefix for output file names"
    image: "Apptainer image providing the GRIDSS helper scripts and samtools"
    realignMapq: "Minimum mapping quality for a realignment to replace the original"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job, also the samtools thread count"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  output {
    File bam = "~{outputPrefix}.bam"
    File bamIndex = "~{outputPrefix}.bam.bai"
  }
}

task alignLongHomopolymers {
  input {
    File inputBam
    File inputBamIndex
    String outputPrefix
    String genomeModule
    String refFasta
    String refFai
    String refDict
    String image
    Int homopolymerLength = 10
    Int jobMemory = 16
    Int cores = 8
    Int timeout = 6
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    inputBam: "Realigned assembled haplotypes"
    inputBamIndex: "Index of inputBam"
    outputPrefix: "Prefix for output file names"
    genomeModule: "Modules providing the reference genome and its resource files"
    refFasta: "Reference FASTA"
    refFai: "Reference FASTA index"
    refDict: "Reference sequence dictionary"
    image: "Apptainer image providing the GRIDSS helper scripts"
    homopolymerLength: "Minimum homopolymer length for a haplotype to be realigned"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File bam = "~{outputPrefix}.bam"
    File bamIndex = "~{outputPrefix}.bam.bai"
  }
}

task matchReadsToHaplotypes {
  input {
    Array[File] tumorCrams
    Array[File] tumorCramIndexes
    Array[File] germlineCrams
    Array[File] germlineCramIndexes
    File assemblyBam
    File assemblyBamIndex
    File shardBed
    String outputPrefix
    String minIndelScSize
    String minMismatchCount
    String genomeModule
    String refFasta
    String image
    Int minMapq = 5
    Int maxReadsPerWorkingArea = 50000
    Int jobMemory = 12
    Int cores = 2
    Int timeout = 12
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    tumorCrams: "Tumor CRAM files"
    tumorCramIndexes: "Indexes of tumorCrams, in the same order"
    germlineCrams: "Germline or matched-normal CRAM files"
    germlineCramIndexes: "Indexes of germlineCrams, in the same order"
    assemblyBam: "Realigned assembled haplotypes of the whole genome"
    assemblyBamIndex: "Index of assemblyBam"
    shardBed: "BED file of the regions in this shard"
    outputPrefix: "Prefix for output file names"
    minIndelScSize: "Minimum indel or soft-clip length for a read to be matched; x;y gives tumor and normal values"
    minMismatchCount: "Minimum mismatch count for a read to be matched; x;y gives tumor and normal values"
    genomeModule: "Modules providing the reference genome and its resource files"
    refFasta: "Reference FASTA; its .fai must sit next to it"
    image: "Apptainer image providing sv_rematch"
    minMapq: "Minimum mapping quality for a read to be matched"
    maxReadsPerWorkingArea: "Maximum number of reads matched per working area"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job, also the sv_rematch thread count"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File readHaplotypeAssignment = "~{outputPrefix}_rematched_hap.txt"
  }
}

task mergeMatchedReads {
  input {
    Array[File] matchedReads
    String outputPrefix
    String image
    Int jobMemory = 4
    Int cores = 1
    Int timeout = 4
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    matchedReads: "Read-to-haplotype assignments of all shards"
    outputPrefix: "Prefix for output file names"
    image: "Apptainer image providing sort and awk"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  output {
    File readHaplotypeAssignment = "~{outputPrefix}.all_matches_by_hap.txt"
  }
}

task addReadSupportToAssembly {
  input {
    File assemblyBam
    File matchedReads
    String outputPrefix
    String genomeModule
    String refFasta
    String image
    Int jobMemory = 16
    Int cores = 1
    Int timeout = 12
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    assemblyBam: "Realigned assembled haplotypes"
    matchedReads: "Read-to-haplotype assignments, ordered by haplotype"
    outputPrefix: "Prefix for output file names"
    genomeModule: "Modules providing the reference genome and its resource files"
    refFasta: "Reference FASTA; its .fai must sit next to it"
    image: "Apptainer image providing sv_rematch and samtools"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File bam = "~{outputPrefix}_assembly.support.bam"
    File bamIndex = "~{outputPrefix}_assembly.support.bam.bai"
  }
}

task identifyVariants {
  input {
    Array[String] sampleNames
    File assemblyBam
    File assemblyBamIndex
    String outputPrefix
    String genomeModule
    String refFasta
    String refFai
    String refDict
    String blacklistBed
    String image
    String gridssConfig = "variantcalling.breakpointLowQuality = 120.0 variantcalling.breakendLowQuality = 120.0 variantcalling.writeFiltered = true variantcalling.requireReadPair = false variantcalling.requireSplitRead = false variantcalling.callUnassembledBreakends = true variantcalling.callUnassembledBreakpoints = true variantcalling.breakendMaxAssemblySupportBias = 100 variantcalling.breakpointMaxAssemblySupportBias = 100 variantcalling.includeSingleAssemblyFilter = false variantcalling.requiredReadAndAssemblyBreakpointOverlap = 5"
    Int jobMemory = 16
    Int cores = 2
    Int timeout = 12
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    sampleNames: "Sample names in the order GRIDSS reports them: tumor first, then germline"
    assemblyBam: "Assembled haplotypes annotated with their supporting reads"
    assemblyBamIndex: "Index of assemblyBam"
    outputPrefix: "Prefix for output file names"
    genomeModule: "Modules providing the reference genome and its resource files"
    refFasta: "Reference FASTA"
    refFai: "Reference FASTA index"
    refDict: "Reference sequence dictionary"
    blacklistBed: "BED of regions GRIDSS ignores"
    image: "Apptainer image providing GRIDSS"
    gridssConfig: "GRIDSS configuration as space-separated key = value settings"
    jobMemory: "Memory allocated for this job (GB); the Java heap is set 2 GB below it"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File vcf = "~{outputPrefix}.vcf.gz"
    File vcfIndex = "~{outputPrefix}.vcf.gz.tbi"
  }
}

task prefilterCandidates {
  input {
    File inputVcf
    File inputVcfIndex
    String outputPrefix
    String image
    String filterExpression = "QUAL>120"
    Int jobMemory = 8
    Int cores = 1
    Int timeout = 2
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    inputVcf: "GRIDSS candidate breakpoints"
    inputVcfIndex: "Tabix index of inputVcf"
    outputPrefix: "Prefix for output file names"
    image: "Apptainer image providing bcftools"
    filterExpression: "bcftools include expression a candidate must pass to be annotated"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  output {
    File vcf = "~{outputPrefix}.vcf.gz"
    File vcfIndex = "~{outputPrefix}.vcf.gz.tbi"
  }
}

task annotateVariants {
  input {
    File inputVcf
    File inputVcfIndex
    Array[File] tumorCrams
    Array[File] tumorCramIndexes
    Array[File] germlineCrams
    Array[File] germlineCramIndexes
    File assemblyBam
    File assemblyBamIndex
    File shardBed
    String outputPrefix
    Boolean isSomatic
    String genomeModule
    String refFasta
    String refFai
    String refDict
    String blacklistBed
    String image
    String gridssConfig = "variantcalling.breakpointLowQuality = 120.0 variantcalling.breakendLowQuality = 120.0 variantcalling.writeFiltered = true variantcalling.requireReadPair = false variantcalling.requireSplitRead = false variantcalling.callUnassembledBreakends = true variantcalling.callUnassembledBreakpoints = true variantcalling.breakendMaxAssemblySupportBias = 100 variantcalling.breakpointMaxAssemblySupportBias = 100 variantcalling.includeSingleAssemblyFilter = false variantcalling.requiredReadAndAssemblyBreakpointOverlap = 5"
    Int germlineMemory = 8
    Int somaticMemory = 128
    Int cores = 1
    Int timeout = 12
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    inputVcf: "Prefiltered GRIDSS candidate breakpoints of the whole genome"
    inputVcfIndex: "Tabix index of inputVcf"
    tumorCrams: "Tumor CRAM files"
    tumorCramIndexes: "Indexes of tumorCrams, in the same order"
    germlineCrams: "Germline or matched-normal CRAM files"
    germlineCramIndexes: "Indexes of germlineCrams, in the same order"
    assemblyBam: "Assembled haplotypes annotated with their supporting reads"
    assemblyBamIndex: "Index of assemblyBam"
    shardBed: "BED file of the regions in this shard"
    outputPrefix: "Prefix for output file names"
    isSomatic: "Whether tumor reads are given, which selects somaticMemory"
    genomeModule: "Modules providing the reference genome and its resource files"
    refFasta: "Reference FASTA"
    refFai: "Reference FASTA index"
    refDict: "Reference sequence dictionary"
    blacklistBed: "BED of regions GRIDSS ignores"
    image: "Apptainer image providing GRIDSS, GATK, samtools and bcftools"
    gridssConfig: "GRIDSS configuration as space-separated key = value settings"
    germlineMemory: "Memory allocated for this job in germline calling (GB); the Java heap is set 4 GB below it"
    somaticMemory: "Memory allocated for this job in somatic calling (GB); the Java heap is set 4 GB below it"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  Int jobMemory = if isSomatic then somaticMemory else germlineMemory

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File annotatedVcf = "~{outputPrefix}.ann.vcf"
  }
}

task mergeVcfs {
  input {
    Array[File] vcfs
    String outputFileName
    String image
    Int jobMemory = 10
    Int cores = 1
    Int timeout = 2
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    vcfs: "Annotated VCFs of all shards"
    outputFileName: "Name of the merged, bgzipped VCF"
    image: "Apptainer image providing GATK"
    jobMemory: "Memory allocated for this job (GB); the Java heap is set 2 GB below it"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  output {
    File vcf = outputFileName
    File vcfIndex = outputFileName + ".tbi"
  }
}

task linkGermlineBreakpoints {
  input {
    File inputVcf
    File inputVcfIndex
    String outputPrefix
    String genomeModule
    String refFasta
    String refFai
    String refDict
    String image
    Int jobMemory = 64
    Int cores = 1
    Int timeout = 12
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    inputVcf: "Annotated GRIDSS breakpoints"
    inputVcfIndex: "Tabix index of inputVcf"
    outputPrefix: "Prefix for output file names"
    genomeModule: "Modules providing the reference genome and its resource files"
    refFasta: "Reference FASTA"
    refFai: "Reference FASTA index"
    refDict: "Reference sequence dictionary"
    image: "Apptainer image providing the GRIDSS R scripts"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File vcf = "~{outputPrefix}_linked.vcf.bgz"
    File vcfIndex = "~{outputPrefix}_linked.vcf.bgz.tbi"
  }
}

task filterSomatic {
  input {
    File inputVcf
    File inputVcfIndex
    String tumorSampleName
    String germlineSampleName
    String outputPrefix
    String genomeModule
    String refFasta
    String refFai
    String refDict
    String refGenomeVersion
    String ponSglFile
    String ponSvFile
    String knownHotspotFile
    String repeatMaskFile
    String image
    String excludeFilters = "SHORT_SR_NORMAL"
    Int? minNormalCoverage
    Int jobMemory = 64
    Int cores = 1
    Int timeout = 6
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    inputVcf: "Annotated GRIDSS breakpoints"
    inputVcfIndex: "Tabix index of inputVcf"
    tumorSampleName: "Tumor sample name in vcf"
    germlineSampleName: "Matched-normal sample name in vcf; empty for tumor-only calling"
    outputPrefix: "Prefix for output file names"
    genomeModule: "Modules providing the reference genome and its resource files"
    refFasta: "Reference FASTA"
    refFai: "Reference FASTA index"
    refDict: "Reference sequence dictionary"
    refGenomeVersion: "Reference genome version as GRIPSS names it: 37 or 38"
    ponSglFile: "GRIPSS panel of normals of single breakends"
    ponSvFile: "GRIPSS panel of normals of breakpoints"
    knownHotspotFile: "GRIPSS known fusion hotspots"
    repeatMaskFile: "GRIPSS repeat mask"
    image: "Apptainer image providing GRIPSS"
    excludeFilters: "GRIPSS filters not to apply, separated by ';'"
    minNormalCoverage: "GRIPSS minimum normal coverage; GRIPSS default when unset"
    jobMemory: "Memory allocated for this job (GB); the Java heap is set 4 GB below it"
    cores: "Number of cores requested for this job"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File vcf = "~{outputPrefix}.gripss.vcf.gz"
    File vcfIndex = "~{outputPrefix}.gripss.vcf.gz.tbi"
  }
}

task convertVcfFormat {
  input {
    File inputVcf
    File inputVcfIndex
    String outputPrefix
    String genomeModule
    String refFasta
    String refFai
    String refDict
    String image
    Int jobMemory = 64
    Int cores = 8
    Int timeout = 12
    String modules = "ugsv/1.37.1"
  }

  parameter_meta {
    inputVcf: "Structural variant calls in breakend notation"
    inputVcfIndex: "Tabix index of inputVcf"
    outputPrefix: "Prefix for output file names"
    genomeModule: "Modules providing the reference genome and its resource files"
    refFasta: "Reference FASTA"
    refFai: "Reference FASTA index"
    refDict: "Reference sequence dictionary"
    image: "Apptainer image providing the GRIDSS R scripts and bcftools"
    jobMemory: "Memory allocated for this job (GB)"
    cores: "Number of cores requested for this job, also the number of R worker processes"
    timeout: "Hours before task timeout"
    modules: "Environment modules to load; ugsv provides apptainer and the container images"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules} ~{genomeModule}"
  }

  output {
    File vcf = "~{outputPrefix}.vcf.gz"
    File vcfIndex = "~{outputPrefix}.vcf.gz.tbi"
  }
}
