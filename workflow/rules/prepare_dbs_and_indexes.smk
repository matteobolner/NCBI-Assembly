rule bwa_mem_2_index:
    input:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        multiext(
            "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna",
            ".0123",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
        ),
    log:
        "logs/bwa-mem2_index/{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}.log",
    resources:
        mem_mb=90000,
    wrapper:
        "v3.13.8/bio/bwa-mem2/index"


rule miniprot_index:
    input:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        "{genomes_folder}/{species}/{assembly_name}/miniprot/{assembly_accession}_{assembly_name}_genomic.mpi",
    threads: 4
    shell:
        "miniprot -t{threads} -d {output} {input}"


rule get_genome_dict:
    input:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.dict",
    shell:
        "samtools dict {input} > {output}"


#########################################################################################################################################################


"""
MTGENOME MASKING - FOR NUMTS
"""


rule get_mtc_sequences_from_assembly_reports:
    input:
        reports=expand(
            "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report.tsv",
            zip,
            species=genomes["Species"],
            assembly_name=genomes["Assembly Name"],
            assembly_accession=genomes["Assembly Accession"],
        ),
        genomes="config/genomes.tsv",
    output:
        "STATS/assembly_setup/mtc_scaffolds.tsv",
        "STATS/assembly_setup/mtc_scaffolds.bed",
    script:
        "../scripts/assembled_genomes_setup/get_mtgenome_scaffolds_from_assembly_reports.py"


rule gunzip_genome_for_maskfasta:
    input:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        temp(
            "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna"
        ),
    threads: 1
    shell:
        "zcat {input} > {output}"


rule mask_mtc_sequences_in_assembly:
    input:
        fasta="{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna",
        bed="STATS/assembly_setup/mtc_scaffolds.bed",
    output:
        temp(
            "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome.fna"
        ),
    threads: 1
    shell:
        "bedtools maskfasta -fi {input.fasta} -bed {input.bed} -fo {output} "


rule bgzip_masked_genome:
    input:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome.fna",
    output:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome.fna.gz",
    threads: 4
    shell:
        "bgzip {input}"


##########################################################################################################################################################


rule makeblastdb_assembly:
    input:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        multiext(
            "{genomes_folder}/{species}/{assembly_name}/blastdb/{assembly_accession}_{assembly_name}_genomic",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto",
        ),
    params:
        species=lambda wc: wc.species,
        assembly_name=lambda wc: wc.assembly_name,
        assembly_accession=lambda wc: wc.assembly_accession,
        db_name="{genomes_folder}/{species}/{assembly_name}/blastdb/{assembly_accession}_{assembly_name}_genomic",
    shell:
        "zcat {input} | makeblastdb -dbtype nucl -in - -out {params.db_name} -title {params.assembly_name}"


rule lastdb_assembly_for_basic_alignments:
    input:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome.fna.gz",
    output:
        multiext(
            "{genomes_folder}/{species}/{assembly_name}/lastdb_basic/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome",
            ".bck",
            ".des",
            ".prj",
            ".sds",
            ".ssp",
            ".suf",
            ".tis",
        ),
    params:
        species=lambda wc: wc.species,
        assembly_name=lambda wc: wc.assembly_name,
        assembly_accession=lambda wc: wc.assembly_accession,
        db_name="{genomes_folder}/{species}/{assembly_name}/lastdb_basic/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome",
    threads: 4
    shell:
        "zcat {input} | lastdb -P {threads} {params.db_name}"


rule lastdb_assembly_for_near_orthology_with_masked_mtgenome:
    input:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome.fna.gz",
    output:
        multiext(
            "{genomes_folder}/{species}/{assembly_name}/lastdb_near/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome",
            ".bck",
            ".des",
            ".prj",
            ".sds",
            ".ssp",
            ".suf",
            ".tis",
        ),
    params:
        species=lambda wc: wc.species,
        assembly_name=lambda wc: wc.assembly_name,
        assembly_accession=lambda wc: wc.assembly_accession,
        db_name="{genomes_folder}/{species}/{assembly_name}/lastdb_near/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome",
    threads: 4
    shell:
        "zcat {input} | lastdb -P {threads} -uNEAR {params.db_name}"


rule lastdb_reference_for_distant_orthology:
    input:
        "{genomes_folder}/{species}/{ref_assembly_name}/ncbi/{ref_assembly_accession}_{ref_assembly_name}_genomic_masked_mtgenome.fna.gz",
    output:
        protected(
            multiext(
                "{genomes_folder}/{species}/{ref_assembly_name}/lastdb_distant/{ref_assembly_accession}_{ref_assembly_name}_genomic_masked_mtgenome",
                ".des",
                ".prj",
                ".sds",
                ".ssp",
                ".tis",
            )
        ),
    params:
        species=lambda wc: wc.species,
        assembly_name=lambda wc: wc.ref_assembly_name,
        assembly_accession=lambda wc: wc.ref_assembly_accession,
        db_name="{genomes_folder}/{species}/{ref_assembly_name}/lastdb_distant/{ref_assembly_accession}_{ref_assembly_name}_genomic_masked_mtgenome",
    threads: 16
    shell:
        "zcat {input} | lastdb -P {threads} -uMAM8 {params.db_name}"


rule minimap2_index_assembly:
    input:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        "{genomes_folder}/{species}/{assembly_name}/minimap2/{assembly_accession}_{assembly_name}_genomic.mni",
    shell:
        "~/work/software/minimap2-2.26_x64-linux/minimap2 -d {output} {input}"


rule unimap_index_assembly:
    input:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        "{genomes_folder}/{species}/{assembly_name}/minimap2/{assembly_accession}_{assembly_name}_genomic.umi",
    shell:
        "unimap -d {output} {input}"
