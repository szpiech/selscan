# Selscan Verification Pipeline

Use this pipeline to **compare outputs from two different versions of selscan**. It helps identify minor bugs or differences between versions.

## Requirements
- Snakemake

## Configuration
Modify `config.yaml` to include:  
- Input filenames  
- Statistics to compute  
- Any extra parameters for `selscan`

## Workflow
The pipeline generates three sets of outputs:

- `output_sc/` – results from standard version of `selscan`   
- `output_sb/` – results from alternate version of `selscan`  
- `output/` – comparison files and figures

## Command
Inside the folder `verification_pipeline`, run:   
```
snakemake
````