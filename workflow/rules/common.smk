import pandas as pd

configfile: "config/config.yaml"
genomes = pd.read_table(config["genomes"])

sus_scrofa_genomes=genomes[genomes['Species']=='Sus_scrofa']
sus_scrofa_nonref_genomes=genomes[(genomes['Species']=='Sus_scrofa')&(genomes['Assembly Name']!='Sscrofa11.1')]

wildcard_constraints:
    assembly_accession="(GC[AF]_[0-9]{9}(\.[0-9])?)",
    other_assembly_accession="(GC[AF]_[0-9]{9}(\.[0-9])?)",
    main_assembly_accession="(GC[AF]_[0-9]{9}(\.[0-9])?)",
