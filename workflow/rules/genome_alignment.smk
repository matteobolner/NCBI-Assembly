rule ragtag_correct:
    input:
        main_genome="GENOMES/{main_species}/{main_assembly_name}/ncbi/{main_assembly_accession}_{main_assembly_name}_genomic_masked_mtgenome.fna.gz",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic_masked_mtgenome.fna.gz",
    output:
        "ANALYSES/ragtag/{main_species}/{main_assembly_accession}_{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}/ragtag.correct.fasta"
    threads: 4
    params:
        main_species=lambda wc: wc.main_species,
        main_assembly_name=lambda wc: wc.main_assembly_name,
        main_assembly_accession=lambda wc: wc.main_assembly_accession,
        other_species=lambda wc: wc.other_species,
        other_assembly_name=lambda wc: wc.other_assembly_name,
        other_assembly_accession=lambda wc: wc.other_assembly_accession,
        outputdir="ANALYSES/ragtag/{main_species}/{main_assembly_accession}_{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}/"
    shell:
        "ragtag.py correct {input.main_genome} {input.other_genome} -t {threads} -o {params.outputdir}"

rule ragtag_scaffold:
    input:
        main_genome="GENOMES/{main_species}/{main_assembly_name}/ncbi/{main_assembly_accession}_{main_assembly_name}_genomic_masked_mtgenome.fna.gz",
        other_corrected_genome="ANALYSES/ragtag/{main_species}/{main_assembly_accession}_{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}/ragtag.correct.fasta"
    output:
        "ANALYSES/ragtag/{main_species}/{main_assembly_accession}_{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}/ragtag.scaffold.fasta"
    threads: 4
    params:
        main_species=lambda wc: wc.main_species,
        main_assembly_name=lambda wc: wc.main_assembly_name,
        main_assembly_accession=lambda wc: wc.main_assembly_accession,
        other_species=lambda wc: wc.other_species,
        other_assembly_name=lambda wc: wc.other_assembly_name,
        other_assembly_accession=lambda wc: wc.other_assembly_accession,
        outputdir="ANALYSES/ragtag/{main_species}/{main_assembly_accession}_{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}/"
    shell:
        "ragtag.py scaffold {input.main_genome} {input.other_corrected_genome} -t {threads} -r -o {params.outputdir}"

rule last_train_all_vs_all:
    input:
        main_genome_db="GENOMES/{main_species}/{main_assembly_name}/lastdb_near/{main_assembly_accession}_{main_assembly_name}_genomic_masked_mtgenome.prj",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic.fna.gz",
    output:
        "GENOMES/{main_species}/{main_assembly_name}/last_train/{other_species}/{other_assembly_name}/{main_assembly_accession}_{other_assembly_accession}.train",
    threads: 4
    shell:
        "zcat {input.other_genome} | last-train --revsym -E0.05 -C2 $(echo {input.main_genome_db} | sed 's/\.[^.]*$//') -P {threads} > {output}"

rule lastal_near_all_vs_all_masked_mtgenome:
    input:
        main_genome_db="GENOMES/{main_species}/{main_assembly_name}/lastdb_near/{main_assembly_accession}_{main_assembly_name}_genomic_masked_mtgenome.prj",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic_masked_mtgenome.fna.gz",
        trained_model="GENOMES/{main_species}/{main_assembly_name}/last_train/{other_species}/{other_assembly_name}/{main_assembly_accession}_{other_assembly_accession}.train",
    output:
        "GENOMES/{main_species}/{main_assembly_name}/lastal_near_with_masked_mtgenome/{other_species}/{other_assembly_name}/{main_assembly_accession}_{other_assembly_accession}.maf.gz",
    threads: 4
    shell:
        "zcat {input.other_genome} | lastal -E0.05 -C2 --split-f=MAF+ -P {threads} -p {input.trained_model} $(echo {input.main_genome_db} | sed 's/\.[^.]*$//') | gzip > {output}"

rule last_split_near_assembly_and_reference_masked:
    input:
        "GENOMES/{main_species}/{main_assembly_name}/lastal_near_with_masked_mtgenome/{other_species}/{other_assembly_name}/{main_assembly_accession}_{other_assembly_accession}.maf.gz",
    output:
        temp("GENOMES/{main_species}/{main_assembly_name}/last_split_near/{other_species}/{other_assembly_name}/{main_assembly_accession}_{other_assembly_accession}.split.maf")
    shell:
        "zcat {input} | last-split -r -m1e-5 - | last-postmask > {output}"

rule sort_last_split_near:
    input:
        "GENOMES/{main_species}/{main_assembly_name}/last_split_near/{other_species}/{other_assembly_name}/{main_assembly_accession}_{other_assembly_accession}.split.maf"
    output:
        "GENOMES/{main_species}/{main_assembly_name}/last_split_near/{other_species}/{other_assembly_name}/{main_assembly_accession}_{other_assembly_accession}.split.sorted.maf"
    shell:
        "maf-sort {input} > {output}"

rule maf_join_all_near_split_alignments:
    input:
        expand(
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split_near/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.sorted.maf",
        zip,
        other_species=genomes[(~genomes['Species'].isin(['Capra_hircus','Bos_taurus']))&(genomes['Assembly Name']!='Sscrofa11.1')]['Species'],
        other_assembly_name=genomes[(~genomes['Species'].isin(['Capra_hircus','Bos_taurus']))&(genomes['Assembly Name']!='Sscrofa11.1')]['Assembly Name'],
        other_assembly_accession=genomes[(~genomes['Species'].isin(['Capra_hircus','Bos_taurus']))&(genomes['Assembly Name']!='Sscrofa11.1')]['Assembly Accession'],
        )
    output:
        "ANALYSES/lastal/joined_aligmments/Sus_scrofa/Sscrofa11.1/GCF_000003025.6_vs_all_suidae.maf"
    shell:
        "maf-join {input} > {output}"


rule convert_last_split_near_output_to_tab:
    input:
        "GENOMES/{main_species}/{main_assembly_name}/last_split_near/{other_species}/{other_assembly_name}/{main_assembly_accession}_vs_{other_assembly_accession}.split.maf.gz",
    output:
        "GENOMES/{main_species}/{main_assembly_name}/last_split_near/{other_species}/{other_assembly_name}/{main_assembly_accession}_vs_{other_assembly_accession}.split.tab"
    shell:
        "zcat {input} | maf-convert tab - > {output}"

rule convert_last_split_near_output_to_blasttab:
    input:
        "GENOMES/Sus_scrofa/Sscrofa11.1/last_split/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.split.maf.gz",
    output:
        "GENOMES/{main_species}/{main_assembly_name}/last_split_near/{other_species}/{other_assembly_name}/{main_assembly_accession}_vs_{other_assembly_accession}.split.blasttab"
    shell:
        "zcat {input} | maf-convert blasttab - > {output}"

rule elaborate_blasttab_near_output_and_get_statistics:
    input:
        "GENOMES/{main_species}/{main_assembly_name}/last_split_near/{other_species}/{other_assembly_name}/{main_assembly_accession}_vs_{other_assembly_accession}.split.blasttab",
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
        ref_genome_db="GENOMES/Sus_scrofa/Sscrofa11.1/lastdb_distant/GCF_000003025.6_Sscrofa11.1_genomic_masked_mtgenome.prj",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic_masked_mtgenome.fna.gz",
    output:
        protected("GENOMES/Sus_scrofa/Sscrofa11.1/last_train_distant/{other_species}/{other_assembly_name}/GCF_000003025.6_{other_assembly_accession}.train"),
    shell:
        "zcat {input.other_genome} | last-train --revsym -E0.05 -C2 $(echo {input.ref_genome_db} | sed 's/\.[^.]*$//') -P {threads} > {output}"

rule lastal_distant_assembly_to_reference:
    input:
        ref_genome_db="GENOMES/Sus_scrofa/Sscrofa11.1/lastdb_distant/GCF_000003025.6_Sscrofa11.1_genomic_masked_mtgenome.prj",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic_masked_mtgenome.fna.gz",
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

rule minigraph_generate_graph:
    input:
        ref="GENOMES/Sus_scrofa/Sscrofa11.1/ncbi/GCF_000003025.6_Sscrofa11.1_genomic_masked_mtgenome.fna.gz",
        others=expand("GENOMES/Sus_scrofa/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic_masked_mtgenome.fna.gz",
            zip,
            other_assembly_name=genomes[(genomes['Species']!='Sus_scrofa')&(genomes['Assembly Name']!='Sscrofa11.1')]['Assembly Name'],
            other_assembly_accession=genomes[(genomes['Species']=='Sus_scrofa')&(genomes['Assembly Name']!='Sscrofa11.1')]['Assembly Accession'],
        )
    output:
        "ANALYSES/minigraph/sus_scrofa_graph.gfa"
    threads:
        4
    shell:
        "minigraph -cxggs -t{threads} {input.ref} {input.others} > {output}"

rule minigraph_call_variants:
    input:
        "ANALYSES/minigraph/sus_scrofa_graph.gfa"
    output:
        "ANALYSES/minigraph/var.bed"
    shell:
        "gfatools bubble {input} > {output}"

rule minigraph_get_paths:
    input:
        graph="ANALYSES/minigraph/sus_scrofa_graph.gfa",
        genome="GENOMES/Sus_scrofa/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic_masked_mtgenome.fna.gz"
    output:
        "ANALYSES/minigraph/genome_paths/Sus_scrofa/{other_assembly_name}/{other_assembly_accession}_{other_assembly_name}.gfa"
    shell:
        "minigraph -cxasm --call {input.graph} {input.genome} > {output}"


rule minimap_align_close_genomes:
    input:
        main_genome="GENOMES/{main_species}/{main_assembly_name}/ncbi/{main_assembly_accession}_{main_assembly_name}_genomic_masked_mtgenome.fna.gz",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic_masked_mtgenome.fna.gz",
    output:
        "ANALYSES/minimap2/intraspecies_alignment/{main_species}/{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}_vs_{main_assembly_accession}.paf"
    wildcard_constraints:
         main_assembly_accession="(GC[AF]_[0-9]{9}(\.[0-9])?)"
    threads:
        4
    shell:
        "if [ {input.main_genome} != {input.other_genome} ]; then minimap2 -cx asm5 -t{threads} --cs {input.main_genome} {input.other_genome} > {output}; else touch {output} ; fi"

rule minimap_align_less_close_genomes:
    input:
        main_genome="GENOMES/{main_species}/{main_assembly_name}/ncbi/{main_assembly_accession}_{main_assembly_name}_genomic_masked_mtgenome.fna.gz",
        other_genome="GENOMES/{other_species}/{other_assembly_name}/ncbi/{other_assembly_accession}_{other_assembly_name}_genomic_masked_mtgenome.fna.gz",
    output:
        "ANALYSES/minimap2/interspecies_alignment/{main_species}/{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}_vs_{main_assembly_accession}.paf"
    wildcard_constraints:
         main_assembly_accession="(GC[AF]_[0-9]{9}(\.[0-9])?)"
    threads:
        4
    shell:
        "if [ {input.main_genome} != {input.other_genome} ]; then minimap2 -cx asm20 -t{threads} --cs {input.main_genome} {input.other_genome} > {output}; else touch {output} ; fi"

rule sort_minimap_output_intraspecies:
    # sort by reference start coordinate
    input:
        close="ANALYSES/minimap2/intraspecies_alignment/{main_species}/{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}_vs_{main_assembly_accession}.paf",
    output:
        close="ANALYSES/minimap2/intraspecies_alignment/{main_species}/{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}_vs_{main_assembly_accession}_sorted.paf",
    shell:
        "sort -k6,6 -k8,8n {input.close} > {output.close}"

rule sort_minimap_output_interspecies:
    # sort by reference start coordinate
    input:
        less_close="ANALYSES/minimap2/interspecies_alignment/{main_species}/{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}_vs_{main_assembly_accession}.paf"
    output:
        less_close="ANALYSES/minimap2/interspecies_alignment/{main_species}/{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}_vs_{main_assembly_accession}_sorted.paf"
    shell:
        "sort -k6,6 -k8,8n {input.less_close} > {output.less_close}"

rule paftools_call_vars_intraspecies:
    input:
        ref="GENOMES/{main_species}/{main_assembly_name}/ncbi/{main_assembly_accession}_{main_assembly_name}_genomic_masked_mtgenome.fna.gz",
        paf="ANALYSES/minimap2/intraspecies_alignment/{main_species}/{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}_vs_{main_assembly_accession}_sorted.paf"
    output:
        "ANALYSES/minimap2/intraspecies_alignment/{main_species}/{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}_vs_{main_assembly_accession}.var.txt"
    shell:
        "paftools.js call {input.paf} -f {input.ref} -l 1000 -L 1000 > {output}"


rule paftools_call_vars_interspecies:
    input:
        ref="GENOMES/{main_species}/{main_assembly_name}/ncbi/{main_assembly_accession}_{main_assembly_name}_genomic_masked_mtgenome.fna.gz",
        paf="ANALYSES/minimap2/interspecies_alignment/{main_species}/{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}_vs_{main_assembly_accession}_sorted.paf"
    output:
        "ANALYSES/minimap2/interspecies_alignment/{main_species}/{main_assembly_name}/{other_species}/{other_assembly_name}/{other_assembly_accession}_vs_{main_assembly_accession}.var.txt"
    shell:
        "paftools.js call {input.paf} -f {input.ref} -l 1000 -L 1000 > {output}"


rule mumandco_test:
    input:
        ref="GENOMES/Sus_scrofa/Sscrofa11.1/ncbi/GCF_000003025.6_Sscrofa11.1_genomic_masked_mtgenome.fna",
        other="ANALYSES/ragtag/Sus_scrofa/{other_assembly_name}/{other_assembly_accession}/ragtag.scaffold.fasta",
    output:
        directory("ANALYSES/mumandco/Sus_scrofa/{other_assembly_name}/{other_assembly_accession}/"),
    threads:
        4
    shell:
        "mumandco_v3.8.sh -r {input.ref} -q {input.other} -g 2501912388 -o {output} -t {threads}"
