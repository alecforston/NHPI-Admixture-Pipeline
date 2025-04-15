input_file = "data_2/cleaned_phased.vcf"

with open(input_file, 'r') as f:
        content = f.read()
# Replace spaces with tabs
content_with_tabs = content.replace(' ', '\t')
with open(input_file, 'w') as f:
        f.write(content_with_tabs)
