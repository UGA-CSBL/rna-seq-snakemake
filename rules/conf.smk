import pandas as pd

GENCODE_FILES = config["gencode"]["uri"]

SAMPLES = (
    pd.read_csv(config["samples"])
    .set_index("Sample", drop=False)
    .assign(Timepoint=lambda x: x["Sample"].str[1])
    .sort_values(["Timepoint", "SampleType"])
)
