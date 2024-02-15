# Script to move files from the NCBI and the Prokka directories
import os
import sys

out_dir = os.path.abspath(snakemake.output[0])

# Make sure that the output dir exists
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
for f in snakemake.input:
    fname = os.path.split(f)[1]
    src_rel_path = os.path.relpath(f, start = out_dir)
    target = os.path.join(out_dir,fname)
    target_abs_path = os.path.abspath(target)
    os.symlink(src_rel_path, target_abs_path)
    print(f"Link {src_rel_path} to {out_dir}")
