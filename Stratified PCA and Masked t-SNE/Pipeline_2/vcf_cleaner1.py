file_path = "data_2/query_file_phased.vcf"
output_path = "data_2/cleaned_phased.vcf"

# Read file into a list of lines
with open(file_path, "r") as file:
    lines = file.readlines()

# Skip the original first line
#lines = lines[1:]

# Clean the new first line
first_line_tokens = lines[9].strip().split()
cleaned_tokens = [token.split('_')[-1] if '_' in token else token for token in first_line_tokens]
cleaned_first_line = ' '.join(cleaned_tokens)

# Replace the first line
lines[9] = cleaned_first_line + '\n'

# Write to output file
with open(output_path, "w") as out_file:
    out_file.writelines(lines)
