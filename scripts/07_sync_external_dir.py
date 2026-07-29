import os
import shutil


src = snakemake.input.source_dir
done = snakemake.output.done
name = snakemake.params.name if hasattr(snakemake.params, 'name') else 'external'

if not os.path.isdir(src):
    raise FileNotFoundError(f'Source directory not found: {src}')

base_name = os.path.splitext(os.path.basename(done))[0].replace('_sync', '')
dst_root = os.path.join('results', 'external', base_name)
os.makedirs(dst_root, exist_ok=True)

n_files = 0
n_fig = 0
for root, _, files in os.walk(src):
    rel = os.path.relpath(root, src)
    out_dir = dst_root if rel == '.' else os.path.join(dst_root, rel)
    os.makedirs(out_dir, exist_ok=True)
    for fn in files:
        s = os.path.join(root, fn)
        d = os.path.join(out_dir, fn)
        shutil.copy2(s, d)
        n_files += 1
        if fn.lower().endswith(('.png', '.pdf', '.svg', '.jpg', '.jpeg', '.tif', '.tiff')):
            n_fig += 1

os.makedirs(os.path.dirname(done), exist_ok=True)
with open(done, 'w', encoding='utf-8') as f:
    f.write(f'name={name}\n')
    f.write(f'source={src}\n')
    f.write(f'destination={dst_root}\n')
    f.write(f'files={n_files}\n')
    f.write(f'figures={n_fig}\n')
