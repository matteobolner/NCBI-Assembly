import os
import scripts.assembled_genomes_setup.NCBI_Tools as NCBI_Tools


rule download_assembly:
    output:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    params:
        url=lambda wc: NCBI_Tools.get_genome_fasta_url_from_ncbi_genome_assembly_accession_and_name(
            wc.assembly_accession, wc.assembly_name
        ),
    shell:
        "rsync --copy-links --times --verbose rsync://{params.url} {output}"


rule download_assembly_report:
    output:
        temp("GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report_to_uncomment.txt")
    params:
        url=lambda wc: NCBI_Tools.get_assembly_report_from_NCBI_genome_assembly_accession_and_name(
            wc.assembly_accession, wc.assembly_name
        ),
    shell:
        "rsync --copy-links --times --verbose rsync://{params.url} {output}"

rule uncomment_assembly_report_header:
    input:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report_to_uncomment.txt"
    output:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report.txt"
    shell:
        "sed 's/# Sequence-Name/Sequence-Name/' {input} > {output}"

rule makeblastdb_assembly:
    input:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        protected(multiext(
            "GENOMES/{species}/{assembly_name}/blastdb/{assembly_accession}_{assembly_name}_genomic",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto",
        )),
    params:
        species=lambda wc: wc.species,
        assembly_name=lambda wc: wc.assembly_name,
        assembly_accession=lambda wc: wc.assembly_accession,
        db_name="GENOMES/{species}/{assembly_name}/blastdb/{assembly_accession}_{assembly_name}_genomic",
    shell:
        "zcat {input} | makeblastdb -dbtype nucl -in - -out {params.db_name} -title {params.assembly_name}"


rule lastdb_assembly:
    input:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        protected(multiext(
            "GENOMES/{species}/{assembly_name}/lastdb/{assembly_accession}_{assembly_name}_genomic",
            ".bck",
            ".des",
            ".prj",
            ".sds",
            ".ssp",
            ".suf",
            ".tis",
        )),
    params:
        species=lambda wc: wc.species,
        assembly_name=lambda wc: wc.assembly_name,
        assembly_accession=lambda wc: wc.assembly_accession,
        db_name="GENOMES/{species}/{assembly_name}/lastdb/{assembly_accession}_{assembly_name}_genomic",
    threads: 4
    shell:
        "zcat {input} | lastdb -P {threads} -uNEAR {params.db_name}"

rule lastdb_reference_for_distant_orthology:
    input:
        "GENOMES/{species}/{ref_assembly_name}/ncbi/{ref_assembly_accession}_{ref_assembly_name}_genomic.fna.gz",
    output:
        protected(multiext(
            "GENOMES/{species}/{ref_assembly_name}/lastdb_distant_orthology/{ref_assembly_accession}_{ref_assembly_name}_genomic",
            ".bck",
            ".des",
            ".prj",
            ".sds",
            ".ssp",
            ".suf",
            ".tis",
        )),
    params:
        species=lambda wc: wc.species,
        assembly_name=lambda wc: wc.ref_assembly_name,
        assembly_accession=lambda wc: wc.ref_assembly_accession,
        db_name="GENOMES/{species}/{ref_assembly_name}/lastdb_distant_orthology/{ref_assembly_accession}_{ref_assembly_name}_genomic",
    threads: 4
    shell:
        "zcat {input} | lastdb -P {threads} -uMAM8 {params.db_name}"

rule minimap2_index_assembly:
    input:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        protected("GENOMES/{species}/{assembly_name}/minimap2/{assembly_accession}_{assembly_name}_genomic.mni"),
    shell:
        "minimap2 -d {output} {input}"
