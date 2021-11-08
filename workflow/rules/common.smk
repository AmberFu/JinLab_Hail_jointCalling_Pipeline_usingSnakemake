import pandas as pd

SAMPLES = pd.read_csv(config["samples"], sep="\t", header=None, names=['SAMPLE'], dtype={"SAMPLE": str})\
    .set_index("SAMPLE", drop=False)\
    .sort_index().SAMPLE.tolist()

INTERVALS = pd.read_csv(config['interval'], sep="\t", header=None, names=['INTERVAL'], dtype={"INTERVAL": str})\
    .set_index("INTERVAL", drop=False)\
    .sort_index().INTERVAL.tolist()
