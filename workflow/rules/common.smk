import pandas as pd

configfile: "config/config.yaml"
genomes = pd.read_table(config["genomes"])

wildcard_constraints:
    assembly_accession="(GC[AF]_[0-9]{9}(\.[0-9])?)",
    other_assembly_accession="(GC[AF]_[0-9]{9}(\.[0-9])?)",
    main_assembly_accession="(GC[AF]_[0-9]{9}(\.[0-9])?)",
