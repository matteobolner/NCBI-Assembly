import scripts.assembled_genomes_setup.NCBI_Tools as NCBI_Tools


rule download_assembly:
    output:
        temp("{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic_to_bgzip.fna.gz"),
    params:
        url=lambda wc: NCBI_Tools.get_genome_fasta_url_from_ncbi_genome_assembly_accession_and_name(
            wc.assembly_accession, wc.assembly_name
        ),
    shell:
        "rsync --copy-links --times --verbose rsync://{params.url} {output}"

rule bgzip_assembly:
    input:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic_to_bgzip.fna.gz",
    output:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
    shell:
        "gunzip -c {input} | bgzip > {output}"

#rule temporarily_extract_genome:
#    input:
#        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
#    output:
#        temp("{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna"),
#    shell:
#        "bgzip -d {input} -c > {output}"

rule download_assembly_report:
    output:
        temp("{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report_to_uncomment.txt")
    params:
        url=lambda wc: NCBI_Tools.get_assembly_report_from_NCBI_genome_assembly_accession_and_name(
            wc.assembly_accession, wc.assembly_name
        ),
    shell:
        "rsync --copy-links --times --verbose rsync://{params.url} {output}"

rule uncomment_assembly_report_header:
    input:
        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report_to_uncomment.txt"
    output:
        temp("{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report_to_edit.txt")
    shell:
        "sed 's/# Sequence-Name/Sequence-Name/' {input} > {output}"

rule add_refseq_or_genbank_unique_column_to_report:
    input:
        report_to_edit="{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report_to_edit.txt"
    output:
        edited_report="{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report.tsv"
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
