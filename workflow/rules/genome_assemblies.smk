import os
import scripts.assembled_genomes_setup.NCBI_Tools as NCBI_Tools


rule download_assembly:
    output:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    params:
        url=lambda wc: NCBI_Tools.GetGenomeFastaUrlFromNCBIGenomeAssemblyAccessionAndName(
            wc.assembly_accession, wc.assembly_name
        ),
    shell:
        "rsync --copy-links --times --verbose rsync://{params.url} {output}"


rule makeblastdb_assembly:
    input:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        multiext(
            "GENOMES/{species}/{assembly_name}/blastdb/{assembly_accession}_{assembly_name}_genomic",
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
        db_name="GENOMES/{species}/{assembly_name}/blastdb/{assembly_accession}_{assembly_name}_genomic",
    shell:
        "zcat {input} | makeblastdb -dbtype nucl -in - -out {params.db_name} -title {params.assembly_name}"


rule lastdb_assembly:
    input:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        multiext(
            "GENOMES/{species}/{assembly_name}/lastdb/{assembly_accession}_{assembly_name}_genomic",
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
        db_name="GENOMES/{species}/{assembly_name}/lastdb/{assembly_accession}_{assembly_name}_genomic",
    threads: 4
    shell:
        "zcat {input} | lastdb -P {threads} -uNEAR {params.db_name}"


rule minimap2_index_assembly:
    input:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        "GENOMES/{species}/{assembly_name}/minimap2/{assembly_accession}_{assembly_name}_genomic.mni",
    shell:
        "minimap2 -d {output} {input}"


rule last_train_reference_vs_assembly:
    input:
        ref_genome_db="GENOMES/Sus_scrofa/Sscrofa11.1/lastdb/GCA_000003025.6_Sscrofa11.1_genomic.prj",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic.fna.gz",
    output:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_train/{other_species}/{other_assembly_name}/GCA_000003025.6_{other_assembly_accession}.train",
    params:
        ref_genome_db_name="GENOMES/Sus_scrofa/Sscrofa11.1/lastdb/GCA_000003025.6_Sscrofa11.1_genomic",
        #ref_species=lambda wc: wc.ref_species,
        #ref_assembly_name=lambda wc: wc.ref_assembly_name,
        #ref_assembly_accession=lambda wc: wc.ref_assembly_accession,
        #ref_genome_db_name="GENOMES/{ref_species}/{ref_assembly_name}/lastdb/{ref_assembly_accession}_{ref_assembly_name}_genomic",
    shell:
        "zcat {input.other_genome} | last-train --revsym -E0.05 -C2 {params.ref_genome_db_name} -P 4 > {output}"

rule lastal_assembly_to_reference:
    input:
        ref_genome_db="GENOMES/Sus_scrofa/Sscrofa11.1/lastdb/GCA_000003025.6_Sscrofa11.1_genomic.prj",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic.fna.gz",
        trained_model="GENOMES/Sus_scrofa/Sscrofa11.1/last_train/{other_species}/{other_assembly_name}/GCA_000003025.6_{other_assembly_accession}.train",
    output:
        "GENOMES/Sus_scrofa/Sscrofa11.1/lastal/{other_species}/{other_assembly_name}/GCA_000003025.6_{other_assembly_accession}.maf",
    shell:
        "zcat {input.other_genome} | lastal -E0.05 -C2 --split-f=MAF+ -p {input.trained_model} {input.other_genome} > {output}"
