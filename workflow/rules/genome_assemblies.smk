import scripts.assembled_genomes_setup.NCBI_Tools as NCBI_Tools


#rule download_assembly:
#    output:
#        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
#    params:
#        url=lambda wc: NCBI_Tools.get_genome_fasta_url_from_ncbi_genome_assembly_accession_and_name(
#            wc.assembly_accession, wc.assembly_name
#        ),
#    shell:
#        "rsync --copy-links --times --verbose rsync://{params.url} {output}"

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
        temp("GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report_to_edit.txt")
    shell:
        "sed 's/# Sequence-Name/Sequence-Name/' {input} > {output}"

rule add_refseq_or_genbank_unique_column_to_report:
    input:
        report_to_edit="GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report_to_edit.txt"
    output:
        edited_report="GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report.tsv"
    params:
        species=lambda wc:wc.species
    run:
        tempreport=pd.read_table(input.report_to_edit, comment='#')
        tempreport['Sequence-ID']=tempreport['RefSeq-Accn'].where(tempreport['RefSeq-Accn']!='na', tempreport['GenBank-Accn'])
        #add annotation of sequence as mitochondrial genome
        if 'CAJOYB010001573.1' in tempreport['Sequence-ID'].tolist():
            tempindex=tempreport[tempreport['Sequence-ID']=='CAJOYB010001573.1'].index[0]
            tempreport.at[tempindex, 'Assigned-Molecule']="MT"
            tempreport.at[tempindex, 'Sequence-Name']="MT"
        tempreport['Chromosome/Scaffold Name']=tempreport['Assigned-Molecule'].where(tempreport['Assigned-Molecule']!='na', tempreport['Sequence-ID'])
        tempreport.to_csv(output.edited_report, index=False, sep='\t')

rule get_mtc_sequences_from_assembly_reports:
    input:
        reports=expand("GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report.tsv", zip, species=genomes['Species'], assembly_name=genomes['Assembly Name'], assembly_accession=genomes['Assembly Accession']),
        genomes="config/genomes.tsv"
    output:
        "STATS/assembly_setup/mtc_scaffolds.tsv",
        "STATS/assembly_setup/mtc_scaffolds.bed"
    script:
        "../scripts/assembled_genomes_setup/get_mtgenome_scaffolds_from_assembly_reports.py"

rule gunzip_genome_for_maskfasta:
    input:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz"
    output:
        temp("GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna")
    shell:
        "zcat {input} > {output}"

rule mask_mtc_sequences_in_assembly:
    input:
        fasta="GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna",
        bed="STATS/assembly_setup/mtc_scaffolds.bed"
    output:
        temp("GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome.fna")
    shell:
        "bedtools maskfasta -fi {input.fasta} -bed {input.bed} -fo {output} "

rule zip_masked_genome:
    input:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome.fna"
    output:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome.fna.gz"
    shell:
        "gzip {input}"

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


rule lastdb_assembly_for_basic_alignments:
    input:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        protected(multiext(
            "GENOMES/{species}/{assembly_name}/lastdb_basic/{assembly_accession}_{assembly_name}_genomic",
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
        db_name="GENOMES/{species}/{assembly_name}/lastdb_basic/{assembly_accession}_{assembly_name}_genomic",
    threads: 4
    shell:
        "zcat {input} | lastdb -P {threads} {params.db_name}"


rule lastdb_assembly_for_near_orthology_with_masked_mtgenome:
    input:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome.fna.gz",
    output:
        protected(multiext(
            "GENOMES/{species}/{assembly_name}/lastdb_near/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome",
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
        db_name="GENOMES/{species}/{assembly_name}/lastdb_near/{assembly_accession}_{assembly_name}_genomic_masked_mtgenome",
    threads: 4
    shell:
        "zcat {input} | lastdb -P {threads} -uNEAR {params.db_name}"


rule lastdb_reference_for_distant_orthology:
    input:
        "GENOMES/{species}/{ref_assembly_name}/ncbi/{ref_assembly_accession}_{ref_assembly_name}_genomic.fna.gz",
    output:
        protected(multiext(
            "GENOMES/{species}/{ref_assembly_name}/lastdb_distant/{ref_assembly_accession}_{ref_assembly_name}_genomic",
            ".des",
            ".prj",
            ".sds",
            ".ssp",
            ".tis",
        )),
    params:
        species=lambda wc: wc.species,
        assembly_name=lambda wc: wc.ref_assembly_name,
        assembly_accession=lambda wc: wc.ref_assembly_accession,
        db_name="GENOMES/{species}/{ref_assembly_name}/lastdb_distant/{ref_assembly_accession}_{ref_assembly_name}_genomic",
    threads: 16
    shell:
        "zcat {input} | lastdb -P {threads} -uMAM8 {params.db_name}"

rule minimap2_index_assembly:
    input:
        "GENOMES/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    output:
        protected("GENOMES/{species}/{assembly_name}/minimap2/{assembly_accession}_{assembly_name}_genomic.mni"),
    shell:
        "minimap2 -d {output} {input}"
