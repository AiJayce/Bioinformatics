import os
import glob

missing = []
for path in snakemake.input:
    if not os.path.exists(path):
        missing.append(path)

    if os.path.isdir(path):
        try:
            if len(os.listdir(path)) == 0:
                missing.append(path + ' (empty directory)')
        except Exception:
            missing.append(path + ' (directory unreadable)')

        if os.path.basename(os.path.normpath(path)).lower() == 'vae_model':
            model_hits = glob.glob(os.path.join(path, '**', 'vae_model.pt'), recursive=True)
            if len(model_hits) == 0:
                missing.append(path + ' (no vae_model.pt found recursively)')

if missing:
    raise FileNotFoundError('Missing required inputs: ' + ', '.join(missing))

os.makedirs(os.path.dirname(snakemake.output[0]), exist_ok=True)
with open(snakemake.output[0], 'w', encoding='utf-8') as f:
    f.write('checked\n')
