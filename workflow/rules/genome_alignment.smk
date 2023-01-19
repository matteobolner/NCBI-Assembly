rule last_train_all_vs_all:
    input:
        main_genome_db="GENOMES/{main_species}/{main_assembly_name}/lastdb_near/{main_assembly_accession}_{main_assembly_name}_genomic.prj",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic.fna.gz",
    output:
        protected("GENOMES/{main_species}/{main_assembly_name}/last_train/{other_species}/{other_assembly_name}/{main_assembly_accession}_{other_assembly_accession}_.train"),
    threads: 6
    shell:
        "zcat {input.other_genome} | last-train --revsym -E0.05 -C2 $(echo {input.main_genome_db} | sed 's/\.[^.]*$//') -P {threads} > {output}"

rule lastal_near_assembly_to_reference:
    input:
        ref_genome_db="GENOMES/Sus_scrofa/Sscrofa11.1/lastdb_near/GCF_000003025.6_Sscrofa11.1_genomic_masked_mtgenome.prj",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic_masked_mtgenome.fna.gz",
        trained_model="GENOMES/Sus_scrofa/Sscrofa11.1/last_train/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.train",
    output:
        protected("GENOMES/Sus_scrofa/Sscrofa11.1/lastal/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.maf"),
    params:
        ref_genome_db="GENOMES/Sus_scrofa/Sscrofa11.1/lastdb_near/GCF_000003025.6_Sscrofa11.1_genomic_masked_mtgenome",
    threads: 6
    shell:
        "zcat {input.other_genome} | lastal -E0.05 -C2 --split-f=MAF+ -P {threads} -p {input.trained_model} {params.ref_genome_db} > {output}"


rule last_split_near_assembly_and_reference:
    input:
        "GENOMES/Sus_scrofa/Sscrofa11.1/lastal/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.maf",
    output:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.maf",
    shell:
        "last-split -r -m1e-5 {input} | last-postmask > {output}"

rule convert_last_split_near_output_to_tab:
    input:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.maf",
    output:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.tab",
    shell:
        "maf-convert tab {input} > {output}"

rule convert_last_split_near_output_to_blasttab:
    input:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.maf",
    output:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.blasttab",
    shell:
        "maf-convert blasttab {input} > {output}"

rule elaborate_blasttab_near_output_and_get_statistics:
    input:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.blasttab",
        "GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_assembly_report.tsv",
    output:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/{other_assembly_accession}.csv",
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/{other_assembly_accession}_per_sequence_stats.csv",
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/{other_assembly_accession}_stats.csv",
    params:
        other_assembly_name=lambda wc:wc.other_assembly_name,
        other_assembly_accession=lambda wc:wc.other_assembly_accession

    script:
        "../scripts/genome_alignment/elaborate_blasttab_output_and_get_statistics.py"

rule merge_alignment_stats:
    input:
        alignment_stats=expand("GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/{other_assembly_accession}_stats.csv",
                zip, other_species=genomes['Species'], other_assembly_name=genomes['Assembly Name'], other_assembly_accession=genomes['Assembly Accession'])
    output:
        merged_stats="STATS/genome_alignment/tables/genomes_aligned_to_Sscrofa11.1_stats.tsv"
    run:
        stats=pd.concat([pd.read_csv(i) for i in input.alignment_stats])
        stats.to_csv(output.merged_stats)

rule last_train_distant_assembly_vs_reference:
    input:
        ref_genome_db="GENOMES/Sus_scrofa/Sscrofa11.1/lastdb_distant/GCF_000003025.6_Sscrofa11.1_genomic.prj",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic.fna.gz",
    output:
        protected("GENOMES/Sus_scrofa/Sscrofa11.1/last_train_distant/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.train"),
    shell:
        "zcat {input.other_genome} | last-train --revsym -E0.05 -C2 $(echo {input.ref_genome_db} | sed 's/\.[^.]*$//') -P {threads} > {output}"

rule lastal_distant_assembly_to_reference:
    input:
        ref_genome_db="GENOMES/Sus_scrofa/Sscrofa11.1/lastdb_distant/GCF_000003025.6_Sscrofa11.1_genomic.prj",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic.fna.gz",
        trained_model="GENOMES/Sus_scrofa/Sscrofa11.1/last_train_distant/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.train",
    output:
        protected("GENOMES/Sus_scrofa/Sscrofa11.1/lastal_distant/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.maf"),
    threads: 4
    shell:
        "zcat {input.other_genome} | lastal -E0.05 -C2 --split-f=MAF+ -m100 -P {threads} -p {input.trained_model} $(echo {input.ref_genome_db} | sed 's/\.[^.]*$//') > {output}"

rule gzip_unused_files:
    input:
        split_maf="GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.maf",
        maf="GENOMES/Sus_scrofa/Sscrofa11.1/lastal/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.maf"
    output:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.maf.gz",
        "GENOMES/Sus_scrofa/Sscrofa11.1/lastal/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.maf.gz"
    shell:
        "gzip {input.maf} && gzip {input.split_maf}"
# rule last_dotplot_assembly_and_reference:
#    input:
#        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.maf"
#    output:
