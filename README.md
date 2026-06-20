# FASTA Analyzer

A fast and user-friendly Python tool for analyzing FASTA files, calculating basic sequence statistics, and preparing data for downstream bioinformatics pipelines.

## Why I Built This
As a biotechnology student passionate about computational biology and the pharmaceutical domain, I needed a clean, reusable tool to quickly inspect large FASTA files (genomes, proteins, plasmids, etc.). Instead of using scattered scripts or heavy software, I built this lightweight analyzer to bridge the gap between raw sequence data and actionable insights.

## Features
- Parse single and multi-FASTA files
- Calculate GC content, sequence length, N50, total bases, etc.
- Amino acid composition and basic physicochemical properties (using Biopython + custom functions)
- Export results to CSV/JSON for further analysis
- Simple CLI + modular design for easy extension

## Technologies Used
- **Python** 3.8+
- **Biopython**
- **RDKit** (for molecular properties when sequences are translated)
- **Pandas** & **Matplotlib** for visualization
- Object-oriented design with proper error handling and logging

## Installation
```bash
git clone https://github.com/yourusername/fasta-analyzer.git
cd fasta-analyzer
pip install -r requirements.txt
