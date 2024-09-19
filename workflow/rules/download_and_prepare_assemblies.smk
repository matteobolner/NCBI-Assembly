rule download_assembly:
    output:
        temp("{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic_to_bgzip.fna.gz"),
    shell:
        "genomers -a {wildcards.assembly_accession} -n {wildcards.assembly_name} -g  | gunzip -c | bgzip > {output} "

#rule temporarily_extract_genome:
#    input:
#        "{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna.gz",
#    output:
#        temp("{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_genomic.fna"),
#    shell:
#        "bgzip -d {input} -c > {output}"

rule download_assembly_report_and_uncomment_header:
    output:
        temp("{genomes_folder}/{species}/{assembly_name}/ncbi/{assembly_accession}_{assembly_name}_assembly_report_to_edit.txt")
    shell:
        "genomers -a {wildcards.assembly_accession} -n {wildcards.assembly_name} -r | sed 's/# Sequence-Name/Sequence-Name/' > {output}"

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
        #if 'CAJOYB010001573.1' in tempreport['Sequence-ID'].tolist():
        #    tempindex=tempreport[tempreport['Sequence-ID']=='CAJOYB010001573.1'].index[0]
        #    tempreport.at[tempindex, 'Assigned-Molecule']="MT"
        #    tempreport.at[tempindex, 'Sequence-Name']="MT"
        tempreport['Chromosome/Scaffold Name']=tempreport['Assigned-Molecule'].where(tempreport['Assigned-Molecule']!='na', tempreport['Sequence-ID'])
        tempreport.to_csv(output.edited_report, index=False, sep='\t')
