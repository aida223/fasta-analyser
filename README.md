# FASTA Analyzer Pro 🧬

A fast, bilingual (Persian/English) web tool for analyzing DNA/RNA/Protein sequences in FASTA format.

[Demo](images/demo.png)

## 🚀 Live Demo
https://fasta-analyser-8qqduidn39l6mmp3rpvkbz.streamlit.app/
(Currently may show temporary downtime due to Streamlit Cloud – try refreshing or check back soon!)
## Features
- Upload multiple FASTA files (even gzipped)
- Automatic parsing and validation
- Length and GC% calculation for each sequence
- Interactive nucleotide distribution charts (Plotly)
- Powerful motif/pattern search (case-insensitive, regex support)
- Summary table for all sequences
- Export full report as **CSV** and **PDF**
- Fully responsive — works great on mobile
- Persian/English interface

## Tech Stack
- Python
- Biopython – FASTA parsing
- Pandas – data handling
- Plotly – interactive visualizations
- Streamlit – web interface
- fpdf2 – PDF generation

## Why I Built This
I'm a 3rd-year Biotechnology student in Iran. During lab work and coursework, I often needed a quick way to check FASTA files without installing heavy software. This tool started as a personal project and evolved into a full-featured analyzer — now with motif search and report export!

My goal: make bioinformatics more accessible, especially for students in developing countries.

## Coming Soon
- Protein translation & amino acid analysis
- Basic quality control (Phred scores if FASTQ added)
- Batch processing improvements
- IUPAC ambiguity code support in motif search

## How to Run Locally
```bash
git clone https://github.com/aida233/fasta-analyser.git
cd fasta-analyser
pip install -r requirements.txt
streamlit run app.py
