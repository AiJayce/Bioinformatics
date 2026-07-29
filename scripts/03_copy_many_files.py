import os
import shutil
import subprocess

outputs = list(snakemake.output)

source_dir = None
if hasattr(snakemake.input, 'source_dir'):
    source_dir = snakemake.input.source_dir
else:
    for item in list(snakemake.input):
        if os.path.isdir(item):
            source_dir = item
            break

if source_dir is None:
    raise ValueError('source_dir input is required for copy_many_files')

sources = [os.path.join(source_dir, os.path.basename(out)) for out in outputs]

for src, dst in zip(sources, outputs):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if str(dst).lower().endswith('.pdf'):
        # Normalize potentially malformed PDFs from upstream R outputs.
        gs = shutil.which('gs')
        if gs is not None:
            result = subprocess.run(
                [gs, '-o', dst, '-sDEVICE=pdfwrite', src],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            if result.returncode == 0 and os.path.exists(dst) and os.path.getsize(dst) > 0:
                continue
        shutil.copy2(src, dst)
    else:
        shutil.copy2(src, dst)
