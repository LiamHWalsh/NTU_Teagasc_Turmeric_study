// -------------------------------------------------------------------------
// Main Metagenomics Pipeline: Fasnicar Preprocessing + MetaPhlAn + HUMAnN3
// -------------------------------------------------------------------------

// -------------------------
// Parameters
// -------------------------
params.input_dir                = null  // Input directory with raw FASTQ files
params.output_dir               = "./results"
params.host_index_dir           = null  // Host index directory for fasnicar host decontamination
params.metaphlan_db_dir         = null  // MetaPhlAn database directory
params.metaphlan_index          = "mpa_vOct22_CHOCOPhlAnSGB_202403"
params.humann_nucleotide_db     = null  // HUMAnN3 nucleotide database (full path)
params.humann_protein_db        = null  // HUMAnN3 protein database (full path)
params.cpus                     = 8
params.email                    = null

// Container/module parameters
params.fasnicar_module          = "fasnicar/0.2.4"
params.bbmap_module             = "bbmap/39.06"
params.metaphlan_module         = "metaphlan/4.1.1"
params.humann_module            = "humann/3.8"

// HUMAnN3 pathways database
params.humann_pathways          = "metacyc"

// HUMAnN3 regrouping - set to false to skip regrouping step
params.skip_regrouping          = false

// -------------------------
// Input validation
// -------------------------
if (!params.input_dir) {
    error "❌ ERROR: --input_dir must be specified"
}
if (!params.host_index_dir) {
    error "❌ ERROR: --host_index_dir must be specified for host decontamination"
}
if (!params.metaphlan_db_dir) {
    error "❌ ERROR: --metaphlan_db_dir must be specified"
}
if (!params.humann_nucleotide_db) {
    error "❌ ERROR: --humann_nucleotide_db must be specified (full path to ChocoPhlAn database)"
}
if (!params.humann_protein_db) {
    error "❌ ERROR: --humann_protein_db must be specified (full path to UniRef database)"
}

log.info """
╔═══════════════════════════════════════════════════════════════════╗
║           METAGENOMICS PIPELINE WITH HUMANN3                      ║
╚═══════════════════════════════════════════════════════════════════╝
Input directory       : ${params.input_dir}
Output directory      : ${params.output_dir}
Host index dir        : ${params.host_index_dir}
MetaPhlAn DB dir      : ${params.metaphlan_db_dir}
MetaPhlAn index       : ${params.metaphlan_index}
HUMAnN nucleotide DB  : ${params.humann_nucleotide_db}
HUMAnN protein DB     : ${params.humann_protein_db}
CPUs per task         : ${params.cpus}
Email notifications   : ${params.email ?: 'Not set'}
═══════════════════════════════════════════════════════════════════
"""

// -------------------------
// Process 1: Find and group paired-end FASTQ files
// -------------------------
process find_fastq_pairs {
    tag "Scanning input directory"
    
    output:
    path "sample_list.txt"
    
    script:
    """
    #!/bin/bash
    # Find all R1 files and extract sample names
    find ${params.input_dir} -type f -name "*_R1_001.fastq.gz" | \
        sed 's/_R1_001.fastq.gz//' | \
        sort -u > sample_list.txt
    
    echo "Found \$(wc -l < sample_list.txt) unique samples"
    cat sample_list.txt
    """
}

// -------------------------
// Process 2: Fasnicar Preprocessing (host decontamination)
// -------------------------
process fasnicar_preprocess {
    tag "${sample_id}"
    publishDir "${params.output_dir}/01_fasnicar/${sample_id}", mode: 'copy'
    cpus params.cpus
    
    input:
    val sample_path
    
    output:
    tuple val(sample_id), path("${sample_id}/*_R1*.fastq.bz2"), path("${sample_id}/*_R2*.fastq.bz2")
    
    script:
    sample_id = sample_path.replaceAll(/.*\//, '')
    """
    #!/bin/bash
    set -e
    
    # Create output directory
    mkdir -p ${sample_id}
    
    # Find and copy the actual R1 and R2 files
    echo "Looking for files at: ${sample_path}"
    
    # Copy R1 file
    if [ -f "${sample_path}_R1_001.fastq.gz" ]; then
        cp "${sample_path}_R1_001.fastq.gz" ${sample_id}/
        echo "Copied R1 file"
    else
        echo "ERROR: R1 file not found: ${sample_path}_R1_001.fastq.gz"
        exit 1
    fi
    
    # Copy R2 file
    if [ -f "${sample_path}_R2_001.fastq.gz" ]; then
        cp "${sample_path}_R2_001.fastq.gz" ${sample_id}/
        echo "Copied R2 file"
    else
        echo "ERROR: R2 file not found: ${sample_path}_R2_001.fastq.gz"
        exit 1
    fi
    
    # Load fasnicar module and run preprocessing
    module load ${params.fasnicar_module}
    
    preprocess.new.py \
        -e .fastq.gz \
        -f _R1 \
        -r _R2 \
        -x ${params.host_index_dir} \
        --verbose \
        -n ${task.cpus} \
        -i ${sample_id}
    
    module unload ${params.fasnicar_module}
    
    echo "✅ Fasnicar preprocessing completed for ${sample_id}"
    """
}

// -------------------------
// Process 3: Interleave paired-end reads with BBMerge
// -------------------------
process interleave_reads {
    tag "${sample_id}"
    publishDir "${params.output_dir}/02_interleaved", mode: 'copy'
    cpus params.cpus
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    
    output:
    tuple val(sample_id), path("${sample_id}_interleaved.fastq")
    
    script:
    """
    #!/bin/bash
    set -e
    
    # Decompress files
    echo "Decompressing ${r1}..."
    bunzip2 -c ${r1} > r1.fastq
    
    echo "Decompressing ${r2}..."
    bunzip2 -c ${r2} > r2.fastq
    
    # Load bbmap and interleave
    module load ${params.bbmap_module}
    
    bbmerge.sh \
        in1=r1.fastq \
        in2=r2.fastq \
        out=${sample_id}_interleaved.fastq \
        threads=${task.cpus}
    
    module unload ${params.bbmap_module}
    
    echo "✅ Interleaving completed for ${sample_id}"
    """
}

// -------------------------
// Process 4: Concatenate paired-end reads for HUMAnN3
// -------------------------
process concatenate_reads {
    tag "${sample_id}"
    publishDir "${params.output_dir}/03_concatenated", mode: 'copy'
    
    input:
    tuple val(sample_id), path(r1), path(r2)
    
    output:
    tuple val(sample_id), path("${sample_id}_concat.fastq.gz")
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "Concatenating paired-end reads for ${sample_id}..."
    
    # Decompress, concatenate, and recompress
    bunzip2 -c ${r1} > r1.fastq
    bunzip2 -c ${r2} > r2.fastq
    cat r1.fastq r2.fastq | gzip > ${sample_id}_concat.fastq.gz
    
    rm r1.fastq r2.fastq
    
    echo "✅ Concatenation completed for ${sample_id}"
    """
}

// -------------------------
// Process 5: MetaPhlAn taxonomic profiling
// -------------------------
process metaphlan_profile {
    tag "${sample_id}"
    publishDir "${params.output_dir}/04_metaphlan/results", mode: 'copy', pattern: "*_metaphlan.txt"
    publishDir "${params.output_dir}/04_metaphlan/bowtie2", mode: 'copy', pattern: "*.bowtie2out.txt"
    cpus params.cpus
    
    input:
    tuple val(sample_id), path(interleaved_fastq)
    
    output:
    tuple val(sample_id), path("${sample_id}_metaphlan.txt"), emit: profile
    path "${sample_id}.bowtie2out.txt", emit: bowtie2
    
    script:
    """
    #!/bin/bash
    set -e
    
    module load ${params.metaphlan_module}
    
    metaphlan ${interleaved_fastq} \
        --input_type fastq \
        --output_file ${sample_id}_metaphlan.txt \
        --nproc ${task.cpus} \
        --add_viruses \
        --unclassified_estimation \
        --bowtie2db ${params.metaphlan_db_dir} \
        --index ${params.metaphlan_index} \
        --bowtie2out ${sample_id}.bowtie2out.txt \
        --force
    
    module unload ${params.metaphlan_module}
    
    echo "✅ MetaPhlAn profiling completed for ${sample_id}"
    """
}

// -------------------------
// Process 6: HUMAnN3 functional profiling
// -------------------------
process humann_profile {
    tag "${sample_id}"
    publishDir "${params.output_dir}/05_humann3/${sample_id}", mode: 'copy'
    cpus params.cpus
    memory '32 GB'
    
    input:
    tuple val(sample_id), path(concat_fastq), path(metaphlan_profile)
    
    output:
    tuple val(sample_id), 
          path("${sample_id}_genefamilies.tsv"),
          path("${sample_id}_pathabundance.tsv"),
          path("${sample_id}_pathcoverage.tsv"), emit: results
    path "${sample_id}_humann_temp", optional: true, emit: temp
    
    script:
    """
    #!/bin/bash
    set -e
    
    module load ${params.humann_module}
    
    echo "Running HUMAnN3 for ${sample_id}"
    echo "Input FASTQ: ${concat_fastq}"
    echo "MetaPhlAn profile: ${metaphlan_profile}"
    echo "Nucleotide DB: ${params.humann_nucleotide_db}"
    echo "Protein DB: ${params.humann_protein_db}"
    
    # Verify databases exist
    if [ ! -d "${params.humann_nucleotide_db}" ]; then
        echo "ERROR: Nucleotide database not found: ${params.humann_nucleotide_db}"
        exit 1
    fi
    echo "✓ Nucleotide database found"
    
    if [ ! -d "${params.humann_protein_db}" ]; then
        echo "ERROR: Protein database not found: ${params.humann_protein_db}"
        exit 1
    fi
    echo "✓ Protein database found"
    
    # Run HUMAnN3
    humann \
        --input ${concat_fastq} \
        --output . \
        --threads ${task.cpus} \
        --taxonomic-profile ${metaphlan_profile} \
        --pathways ${params.humann_pathways} \
        --nucleotide-database ${params.humann_nucleotide_db} \
        --protein-database ${params.humann_protein_db} \
        --remove-temp-output \
        --verbose
    
    # Rename output files to include sample ID
    # Find the actual output files (they have the input filename as prefix)
    base_name=\$(basename ${concat_fastq} .fastq.gz)
    
    if [ -f "\${base_name}_genefamilies.tsv" ]; then
        mv "\${base_name}_genefamilies.tsv" ${sample_id}_genefamilies.tsv
    fi
    
    if [ -f "\${base_name}_pathabundance.tsv" ]; then
        mv "\${base_name}_pathabundance.tsv" ${sample_id}_pathabundance.tsv
    fi
    
    if [ -f "\${base_name}_pathcoverage.tsv" ]; then
        mv "\${base_name}_pathcoverage.tsv" ${sample_id}_pathcoverage.tsv
    fi
    
    module unload ${params.humann_module}
    
    echo "✅ HUMAnN3 profiling completed for ${sample_id}"
    """
}

// -------------------------
// Process 7: Normalize HUMAnN3 gene families
// -------------------------
process humann_normalize_genefamilies {
    tag "${sample_id}"
    publishDir "${params.output_dir}/06_humann3_normalized/genefamilies", mode: 'copy'
    
    input:
    tuple val(sample_id), path(genefamilies), path(pathabundance), path(pathcoverage)
    
    output:
    tuple val(sample_id), path("${sample_id}_genefamilies_cpm.tsv")
    
    script:
    """
    #!/bin/bash
    set -e
    
    module load ${params.humann_module}
    
    # Check if input file exists
    if [[ -f "${genefamilies}" ]]; then
        echo "Processing ${sample_id}..."
        # Normalize to copies per million (CPM)
        humann_renorm_table \
            --input ${genefamilies} \
            --output ${sample_id}_genefamilies_cpm.tsv \
            --units cpm \
            --update-snames
    else
        echo "ERROR: Gene families file not found: ${genefamilies}"
        exit 1
    fi
    
    module unload ${params.humann_module}
    
    echo "✅ Gene families normalized for ${sample_id}"
    """
}

// -------------------------
// Process 8: Normalize HUMAnN3 pathway abundances
// -------------------------
process humann_normalize_pathways {
    tag "${sample_id}"
    publishDir "${params.output_dir}/06_humann3_normalized/pathways", mode: 'copy'
    
    input:
    tuple val(sample_id), path(genefamilies), path(pathabundance), path(pathcoverage)
    
    output:
    tuple val(sample_id), path("${sample_id}_pathabundance_cpm.tsv")
    
    script:
    """
    #!/bin/bash
    set -e
    
    module load ${params.humann_module}
    
    # Check if input file exists
    if [[ -f "${pathabundance}" ]]; then
        echo "Processing ${sample_id}..."
        # Normalize pathway abundance to copies per million (CPM)
        humann_renorm_table \
            --input ${pathabundance} \
            --output ${sample_id}_pathabundance_cpm.tsv \
            --units cpm \
            --update-snames
    else
        echo "ERROR: Pathway abundance file not found: ${pathabundance}"
        exit 1
    fi
    
    module unload ${params.humann_module}
    
    echo "✅ Pathway abundances normalized for ${sample_id}"
    """
}

// -------------------------
// Process 9: Regroup gene families (optional - requires utility mapping databases)
// -------------------------
process humann_regroup_genefamilies {
    tag "${sample_id}"
    publishDir "${params.output_dir}/07_humann3_regrouped", mode: 'copy'
    
    input:
    tuple val(sample_id), path(genefamilies_relab)
    
    output:
    tuple val(sample_id), 
          path("${sample_id}_*.tsv"), optional: true
    
    script:
    """
    #!/bin/bash
    set -e
    
    module load ${params.humann_module}
    
    echo "Attempting to regroup gene families for ${sample_id}..."
    
    # Try to regroup to KEGG Orthologs (may fail if database not installed)
    if humann_regroup_table --help 2>&1 | grep -q "uniref90_ko"; then
        echo "Regrouping to KEGG Orthologs..."
        humann_regroup_table \
            --input ${genefamilies_relab} \
            --output ${sample_id}_ko.tsv \
            --groups uniref90_ko
    else
        echo "⚠️  Warning: uniref90_ko database not available, skipping KO regrouping"
        touch ${sample_id}_ko.tsv.skipped
    fi
    
    # Try to regroup to EC numbers (may fail if database not installed)
    if humann_regroup_table --help 2>&1 | grep -q "uniref90_level4ec"; then
        echo "Regrouping to EC numbers..."
        humann_regroup_table \
            --input ${genefamilies_relab} \
            --output ${sample_id}_ec.tsv \
            --groups uniref90_level4ec
    else
        echo "⚠️  Warning: uniref90_level4ec database not available, skipping EC regrouping"
        touch ${sample_id}_ec.tsv.skipped
    fi
    
    # Try to regroup to GO terms (may fail if database not installed)
    if humann_regroup_table --help 2>&1 | grep -q "uniref90_go"; then
        echo "Regrouping to GO terms..."
        humann_regroup_table \
            --input ${genefamilies_relab} \
            --output ${sample_id}_go.tsv \
            --groups uniref90_go
    else
        echo "⚠️  Warning: uniref90_go database not available, skipping GO regrouping"
        touch ${sample_id}_go.tsv.skipped
    fi
    
    # Regroup to reactions (should be available by default)
    if humann_regroup_table --help 2>&1 | grep -q "uniref90_rxn"; then
        echo "Regrouping to reactions..."
        humann_regroup_table \
            --input ${genefamilies_relab} \
            --output ${sample_id}_rxn.tsv \
            --groups uniref90_rxn
    else
        echo "⚠️  Warning: uniref90_rxn database not available"
        touch ${sample_id}_rxn.tsv.skipped
    fi
    
    module unload ${params.humann_module}
    
    echo "✅ Gene family regrouping completed for ${sample_id}"
    """
}

// -------------------------
// Process 10: Merge all samples' results
// -------------------------
process merge_humann_tables {
    tag "Merging all samples"
    publishDir "${params.output_dir}/08_humann3_merged", mode: 'copy'
    
    input:
    path genefamilies_files
    path pathabundance_files
    
    output:
    path "merged_genefamilies_cpm.tsv"
    path "merged_pathabundance_cpm.tsv"
    
    script:
    """
    #!/bin/bash
    set -e
    
    module load ${params.humann_module}
    
    # Merge gene families (CPM normalized)
    humann_join_tables \
        --input . \
        --output merged_genefamilies_cpm.tsv \
        --file_name genefamilies_cpm
    
    # Merge pathway abundances (CPM normalized)
    humann_join_tables \
        --input . \
        --output merged_pathabundance_cpm.tsv \
        --file_name pathabundance_cpm
    
    module unload ${params.humann_module}
    
    echo "✅ All sample tables merged"
    """
}

// -------------------------
// Workflow
// -------------------------
workflow {
    // Step 1: Find all FASTQ pairs
    find_fastq_pairs()
    
    // Step 2: Create channel from sample list
    sample_ch = find_fastq_pairs.out
        .splitText()
        .map { it.trim() }
        .filter { it != "" }
    
    // Step 3: Fasnicar preprocessing
    // Output: [sample_id, r1_file, r2_file]
    fasnicar_out = fasnicar_preprocess(sample_ch)
    
    // Step 4a: Interleave reads for MetaPhlAn
    // Input: [sample_id, r1_file, r2_file]
    // Output: [sample_id, interleaved_fastq]
    interleave_out = interleave_reads(fasnicar_out)
    
    // Step 4b: Concatenate reads for HUMAnN3
    // Input: [sample_id, r1_file, r2_file]
    // Output: [sample_id, concat_fastq]
    concat_out = concatenate_reads(fasnicar_out)
    
    // Step 5: MetaPhlAn profiling
    // Input: [sample_id, interleaved_fastq]
    // Output (profile emit): [sample_id, metaphlan_profile]
    metaphlan_profile(interleave_out)
    
    // Step 6: Join concatenated reads with MetaPhlAn profiles for HUMAnN3
    // concat_out: [sample_id, concat_fastq]
    // metaphlan_profile.out.profile: [sample_id, metaphlan_profile]
    // After join: [sample_id, concat_fastq, metaphlan_profile]
    humann_input = concat_out.join(metaphlan_profile.out.profile)
    
    // Step 7: HUMAnN3 functional profiling
    // Input: [sample_id, concat_fastq, metaphlan_profile]
    // Output: [sample_id, genefamilies, pathabundance, pathcoverage]
    humann_profile(humann_input)
    
    // Step 8: Normalize gene families
    // Input: [sample_id, genefamilies, pathabundance, pathcoverage]
    // Output: [sample_id, genefamilies_relab]
    genefamilies_norm = humann_normalize_genefamilies(humann_profile.out.results)
    
    // Step 9: Normalize pathways
    // Input: [sample_id, genefamilies, pathabundance, pathcoverage]
    // Output: [sample_id, pathabundance_relab]
    pathways_norm = humann_normalize_pathways(humann_profile.out.results)
    
    // Step 10: Regroup gene families (optional - only if not skipped)
    if (!params.skip_regrouping) {
        humann_regroup_genefamilies(genefamilies_norm)
    }
    
    // Step 11: Merge all samples
    merge_humann_tables(
        genefamilies_norm.map { it[1] }.collect(),
        pathways_norm.map { it[1] }.collect()
    )
}

// -------------------------
// Completion message
// -------------------------
workflow.onComplete {
    log.info """
    ╔═══════════════════════════════════════════════════════════════════╗
    ║                  PIPELINE COMPLETED                               ║
    ╚═══════════════════════════════════════════════════════════════════╝
    Status      : ${workflow.success ? '✅ SUCCESS' : '❌ FAILED'}
    Duration    : ${workflow.duration}
    Output dir  : ${params.output_dir}
    
    Output directories:
    - 01_fasnicar/          : Host-decontaminated reads
    - 02_interleaved/       : Interleaved reads for MetaPhlAn
    - 03_concatenated/      : Concatenated reads for HUMAnN3
    - 04_metaphlan/         : Taxonomic profiles
    - 05_humann3/           : Raw HUMAnN3 outputs
    - 06_humann3_normalized/: Normalized abundances
    - 07_humann3_regrouped/ : Regrouped gene families (KO, EC, GO)
    - 08_humann3_merged/    : Merged tables across all samples
    ═══════════════════════════════════════════════════════════════════
    """.stripIndent()
    
    if (params.email && workflow.success) {
        sendMail(
            to: params.email,
            subject: "Pipeline completed: ${workflow.runName}",
            body: "Pipeline completed successfully!\nDuration: ${workflow.duration}\nOutput: ${params.output_dir}"
        )
    }
}