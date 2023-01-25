import pandas as pd
import numpy as np

df = pd.read_csv(snakemake.input[0], sep="\t", header=None)

report = pd.read_table(
    snakemake.input[1],
)

df.columns = [
    "query_name",
    "reference_name",
    "identity",
    "alignment_length",
    "mismatches",
    "gap_opens",
    "query_start",
    "query_end",
    "reference_start",
    "reference_end",
]


df["query_strand"] = "+"
df["query_strand"] = df["query_strand"].where(df["query_start"] < df["query_end"], "-")
df["query_forward_start"] = df["query_start"].where(
    df["query_strand"] == "+", df["query_end"]
)
df["query_forward_end"] = df["query_end"].where(
    df["query_strand"] == "+", df["query_start"]
)
df["query_forward_start_0_based"] = df["query_forward_start"] - 1
df["query_forward_end_0_based"] = df["query_forward_end"]
df["query_aligned_length"] = (
    df["query_forward_end_0_based"] - df["query_forward_start_0_based"]
)
df["reference_start_0_based"] = df["reference_start"] - 1
df["reference_end_0_based"] = df["reference_end"]
df["reference_aligned_length"] = (
    df["reference_end_0_based"] - df["reference_start_0_based"]
)

seq_size_dict = {
    i: j for i, j in zip(report["Sequence-ID"], report["Sequence-Length"])
}

df["query_sequence_length"] = df["query_name"].apply(lambda x: seq_size_dict[x])

totalcoverage = sum(df["query_aligned_length"]) / report["Sequence-Length"].sum()
per_sequence_coverage_dict = {}
per_sequence_coverage_dict_pctg = {}

for name, group in df.groupby(by="query_name"):
    per_sequence_coverage_dict[name] = group["query_aligned_length"].sum()
    per_sequence_coverage_dict_pctg[name] = (
        group["query_aligned_length"].sum() / group["query_sequence_length"].iloc[0]
    ) * 100

for i in report["GenBank-Accn"].tolist():
    per_sequence_coverage_dict_pctg.setdefault(i, 0)
    per_sequence_coverage_dict.setdefault(i, 0)

report["Length of sequence aligned to reference"] = report["GenBank-Accn"].apply(
    lambda x: per_sequence_coverage_dict[x]
)
report["Percentage of sequence aligned to reference"] = report["GenBank-Accn"].apply(
    lambda x: per_sequence_coverage_dict_pctg[x]
)

stats_df = pd.DataFrame(
    columns=[
        "assembly_name",
        "assembly_accession",
        "total_sequence_length",
        "total_sequence_length_aligned_to_reference_genome",
        "total_sequence_length_aligned_to_reference_genome_pctg",
        "total_number_of_sequences",
        "number_of_sequences_aligned_to_reference_genome",
        "number_of_sequences_aligned_to_reference_genome_pctg",
    ]
)

assembly_name = snakemake.params.other_assembly_name
assembly_accession = snakemake.params.other_assembly_accession

total_sequence_length = report["Sequence-Length"].sum()
total_sequence_length_aligned_to_reference_genome = sum(df["query_aligned_length"])
total_sequence_length_aligned_to_reference_genome_pctg = (
    total_sequence_length_aligned_to_reference_genome / total_sequence_length
) * 100
total_number_of_sequences = len(report)
number_of_sequences_aligned_to_reference_genome = len(df["query_name"].unique())
number_of_sequences_aligned_to_reference_genome_pctg = (
    number_of_sequences_aligned_to_reference_genome / total_number_of_sequences
)

stats_df.loc[0] = [
    assembly_name,
    assembly_accession,
    total_sequence_length,
    total_sequence_length_aligned_to_reference_genome,
    total_sequence_length_aligned_to_reference_genome_pctg,
    total_number_of_sequences,
    number_of_sequences_aligned_to_reference_genome,
    number_of_sequences_aligned_to_reference_genome_pctg,
]

df.to_csv(snakemake.output[0], index=False)


number_of_reference_sequences_aligned={}
reference_sequences_aligned={}
for name,group in df.groupby(by='query_name'):
    number_of_reference_sequences_aligned[name]=(len(group['reference_name'].unique()))
    reference_sequences_aligned[name]=";".join(group['reference_name'].unique().tolist())

for i in report['Sequence-ID'].unique():
    number_of_reference_sequences_aligned.setdefault(i,0)
    reference_sequences_aligned.setdefault(i,np.nan)

report['Number of reference sequences aligned']=report['Sequence-ID'].apply(lambda x: number_of_reference_sequences_aligned[x])
report['ID of reference sequences aligned']=report['Sequence-ID'].apply(lambda x: reference_sequences_aligned[x])

report = report[
    [
        "GenBank-Accn",
        "Length of sequence aligned to reference",
        "Percentage of sequence aligned to reference",
        "Number of reference sequences aligned",
        "ID of reference sequences aligned"

    ]
]

report.to_csv(snakemake.output[1], index=False)

stats_df.to_csv(snakemake.output[2], index=False)
