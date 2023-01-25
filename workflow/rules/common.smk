import pandas as pd

configfile: "config/config.yaml"
genomes = pd.read_table(config["genomes"])

#wildcard_constraints:
#    assembly_accession="""^\b(GC[AF]_[0-9]{9}(\.[0-9])?)\b$""",
#    main_assembly_accession="""^\b(GC[AF]_[0-9]{9}(\.[0-9])?)\b$""",
#    other_assembly_accession="""^\b(GC[AF]_[0-9]{9}(\.[0-9])?)\b$""",
