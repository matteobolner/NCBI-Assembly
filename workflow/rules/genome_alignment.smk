rule last_train_reference_vs_assembly:
    input:
        ref_genome_db="GENOMES/Sus_scrofa/Sscrofa11.1/lastdb/GCF_000003025.6_Sscrofa11.1_genomic.prj",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic.fna.gz",
    output:
        protected("GENOMES/Sus_scrofa/Sscrofa11.1/last_train/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.train"),
    params:
        ref_genome_db_name="GENOMES/Sus_scrofa/Sscrofa11.1/lastdb/GCF_000003025.6_Sscrofa11.1_genomic",
        #ref_species=lambda wc: wc.ref_species,
        #ref_assembly_name=lambda wc: wc.ref_assembly_name,
        #ref_assembly_accession=lambda wc: wc.ref_assembly_accession,
        #ref_genome_db_name="GENOMES/{ref_species}/{ref_assembly_name}/lastdb/{ref_assembly_accession}_{ref_assembly_name}_genomic",
    threads: 6
    shell:
        "zcat {input.other_genome} | last-train --revsym -E0.05 -C2 {params.ref_genome_db_name} -P {threads} > {output}"


rule lastal_assembly_to_reference:
    input:
        ref_genome_db="GENOMES/Sus_scrofa/Sscrofa11.1/lastdb/GCF_000003025.6_Sscrofa11.1_genomic.prj",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic.fna.gz",
        trained_model="GENOMES/Sus_scrofa/Sscrofa11.1/last_train/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.train",
    output:
        protected("GENOMES/Sus_scrofa/Sscrofa11.1/lastal/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.maf"),
    params:
        ref_genome_db="GENOMES/Sus_scrofa/Sscrofa11.1/lastdb/GCF_000003025.6_Sscrofa11.1_genomic",
    threads: 6
    shell:
        "zcat {input.other_genome} | lastal -E0.05 -C2 --split-f=MAF+ -P {threads} -p {input.trained_model} {params.ref_genome_db} > {output}"


rule last_split_assembly_and_reference:
    input:
        "GENOMES/Sus_scrofa/Sscrofa11.1/lastal/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.maf",
    output:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.maf",
    shell:
        "last-split -r -m1e-5 {input} | last-postmask > {output}"


rule convert_last_split_output_to_tab:
    input:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.maf",
    output:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.tab",
    shell:
        "maf-convert tab {input} > {output}"

rule convert_last_split_output_to_blasttab:
    input:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.maf",
    output:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.blasttab",
    shell:
        "maf-convert blasttab {input} > {output}"

rule elaborate_blasttab_output_and_get_statistics:
    input:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.blasttab",
        "GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_assembly_report.txt",
    output:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.csv",
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}_per_sequence_stats.csv",
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}_stats.csv",
    params:
        other_assembly_name=lambda wc:wc.other_assembly_name,
        other_assembly_accession=lambda wc:wc.other_assembly_accession

    script:
        "../scripts/genome_alignment/elaborate_blasttab_output_and_get_statistics.py"

# rule last_dotplot_assembly_and_reference:
#    input:
#        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.maf"
#    output:
