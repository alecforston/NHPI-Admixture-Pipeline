import pandas as pd
from io import StringIO
import re

file_path = "data_2/query_clean_header.msp"

with open(file_path, 'r') as f:
    lines = f.readlines()

# Clean and split the header
header_line = lines[0].strip().replace('n snps', 'n_snps').replace('.0', '')
columns = header_line.split()

# Load the rest of the file with regex-based whitespace splitting
data_str = ''.join(lines[1:])
df = pd.read_csv(
    StringIO(data_str),
    sep=r'\s+',
    engine='python',
    names=columns,
)

# Drop columns ending in ".1"
columns_to_drop = [col for col in df.columns if re.search(r'\.1', col)]
df.drop(columns=columns_to_drop, inplace=True)

df.to_csv("data_2/query_results_cleaned_final.msp", sep='\t', index=False)
