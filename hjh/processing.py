import os
import pandas as pd

def thread_together(inside, outside):
    return (outside[0] + inside + outside[1])

def load_file(filename):
    ext = os.path.splitext(filename)[1]
    if ext == '.csv':
        return pd.read_csv(filename)
    elif ext == '.dat':
        return pd.read_table(filename)