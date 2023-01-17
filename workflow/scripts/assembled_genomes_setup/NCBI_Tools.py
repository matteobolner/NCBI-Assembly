def get_folder_url_from_ncbi_genome_assembly_accession_and_name(
    assembly_accession: str, assembly_name: str
) -> str:
    """
    to test the function:
    assembly_accession GCA_000003025.6
    assembly_name Sscrofa11.1
    """
    genbank_or_refseq_id = assembly_accession.split("_")[0]
    genome_version = assembly_accession.split(".")[1]
    genome_code = assembly_accession.split("_")[1].split(".")[0]
    genome_code_split = [genome_code[i : i + 3] for i in range(0, len(genome_code), 3)]
    url = f"ftp.ncbi.nlm.nih.gov/genomes/all/{genbank_or_refseq_id}/{genome_code_split[0]}/{genome_code_split[1]}/{genome_code_split[2]}/{assembly_accession}_{assembly_name}"
    return url


def get_genome_fasta_url_from_ncbi_genome_assembly_accession_and_name(
    assembly_accession: str, assembly_name: str
) -> str:
    """
    to test the function:
    assembly_accession GCA_000003025.6
    assembly_name Sscrofa11.1
    """
    genbank_or_refseq_id = assembly_accession.split("_")[0]
    genome_version = assembly_accession.split(".")[1]
    genome_code = assembly_accession.split("_")[1].split(".")[0]
    genome_code_split = [genome_code[i : i + 3] for i in range(0, len(genome_code), 3)]
    url = f"ftp.ncbi.nlm.nih.gov/genomes/all/{genbank_or_refseq_id}/{genome_code_split[0]}/{genome_code_split[1]}/{genome_code_split[2]}/{assembly_accession}_{assembly_name}/{assembly_accession}_{assembly_name}_genomic.fna.gz"
    return url


def get_assembly_report_from_NCBI_genome_assembly_accession_and_name(
    assembly_accession: str, assembly_name: str
) -> str:
    """
    to test the function:
    assembly_accession GCA_000003025.6
    assembly_name Sscrofa11.1
    """
    genbank_or_refseq_id = assembly_accession.split("_")[0]
    genome_version = assembly_accession.split(".")[1]
    genome_code = assembly_accession.split("_")[1].split(".")[0]
    genome_code_split = [genome_code[i : i + 3] for i in range(0, len(genome_code), 3)]
    url = f"ftp.ncbi.nlm.nih.gov/genomes/all/{genbank_or_refseq_id}/{genome_code_split[0]}/{genome_code_split[1]}/{genome_code_split[2]}/{assembly_accession}_{assembly_name}/{assembly_accession}_{assembly_name}_assembly_report.txt"
    return url
