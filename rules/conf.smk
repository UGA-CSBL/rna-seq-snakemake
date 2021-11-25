import pandas as pd

GENCODE_FILES = config["gencode"]["uri"]

SAMPLES = pd.read_csv(config["samples"]).set_index("Sample", drop=False)
