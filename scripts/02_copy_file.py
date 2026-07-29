import os
import shutil

if hasattr(snakemake.input, 'source') and snakemake.input.source:
	src = snakemake.input.source
else:
	src = snakemake.input[0]
dst = snakemake.output[0]
if not os.path.exists(src):
	raise FileNotFoundError(f'Source file does not exist: {src}')
os.makedirs(os.path.dirname(dst), exist_ok=True)
shutil.copy2(src, dst)
