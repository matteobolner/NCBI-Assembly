import pandas as pd


configfile: "config/config.yaml"


genomes = pd.read_table(config["genomes"])
