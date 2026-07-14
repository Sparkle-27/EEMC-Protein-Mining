# EMCC-Protein-Mining
A scalable pipeline for large-scale protein resource construction, annotation, and enzyme discovery from extremophile microbial genomes.

## Overview

EMCC-Protein-Mining is an open-source bioinformatics workflow designed for large-scale analysis of extremophile microbial genomes (https://doi.org/10.1038/s41467-026-71145-0).

The pipeline converts genome assemblies into high-quality protein resources and enables downstream functional exploration, protein family analysis, and structure-guided enzyme discovery.

## Workflow

Genome assemblies
        |
        ↓
Prodigal gene prediction
        |
        ↓
Protein quality filtering
        |
        ↓
MMseqs2 clustering
(NR100 / NR90 / NR50)
        |
        ↓
Functional annotation
        |
        ↓
Enzyme mining
        |
        ↓
Structure prediction


## Features

- Large-scale microbial protein processing
- Protein quality control
- Non-redundant protein database construction
- Functional annotation
- Enzyme discovery framework
- HPC compatible


## Requirements

- Linux
- Python >=3.8
- Prodigal
- MMseqs2


## Installation

Clone repository:

git clone https://github.com/xxx/EMCC-Protein-Mining.git


## Usage

Example:

bash scripts/run_pipeline.sh


## Output

The pipeline generates:

- filtered protein sequences
- NR100/NR90/NR50 protein databases
- sequence mapping tables
- downstream annotation-ready datasets


## Citation

If you use this workflow, please cite:

[Publication information will be added]

## Authors

This project was developed by:

- **Huanfa Gong**  gonghuanfa@zju.edu.cn
  Pipeline development, microbial genomics analysis, computational analysis

- **Liangzhen Zheng**  astrozheng@gmail.com
   Project conception, Protein annotation, workflow optimization, computational analysis


## License

MIT License
