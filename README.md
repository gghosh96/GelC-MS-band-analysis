\# GelC MSMS Peptide Analysis



\## Overview



This project analyzes peptide and protein quantification data derived from gel bands, integrating proteomics data with transcript information to visualize protein isoforms and peptide mappings.



---



\## Python Version



\- \*\*Python 3.10\*\* is required to run this project successfully.

\- Python versions above 3.10 (e.g., 3.11, 3.12, 3.13) may cause installation issues, especially with some dependencies like Pillow, due to compatibility problems.

\- It is recommended to create a virtual environment with Python 3.10 to ensure smooth dependency installation.



---



\## Installation



1\. Clone this repository or download the source files.



2\. Set up a Python 3.10 virtual environment (optional but recommended):



&nbsp;  ```bash

&nbsp;  python3.10 -m venv env

&nbsp;  # Windows:

&nbsp;  .\\env\\Scripts\\activate

&nbsp;  # macOS/Linux:

&nbsp;  source env/bin/activate



\## Install dependencies from the requirements.txt file:



pip install -r requirements.txt



Running the Analysis

The main script protein\_panel\_analysis.py accepts a gene name as a command-line argument and generates multiple plot panels related to that gene’s protein data.



\## USAGE:



python protein\_panel\_analysis.py <GENE\_NAME>



Replace <GENE\_NAME> with the gene symbol you want to analyze, e.g.:



python protein\_panel\_analysis.py HNRNPA2B1



\## OUTPUT



The script will generate and save the following plots as PNG files in the current directory:



<GENE\_NAME>\_panelA\_boxplot.png

Boxplot showing relative protein intensity across gel bands.



<GENE\_NAME>\_exon\_colored\_protein.png

Protein sequence colored by exon boundaries.



<GENE\_NAME>\_panelB\_heatmap.png

Heatmap of peptide intensities (tumor vs normal) across gel bands.



<GENE\_NAME>\_panelD\_isoform\_map.png

Visualization of detected peptides mapped onto isoforms.



\## Notes

The script requires internet access to collect transcript information from the UCSC Genome database.



Make sure all input files (protein quantification CSVs, peptide CSVs, and the FASTA file) are located in the correct directories as specified in the script (Gel\_Band\_Protein\_Quants/, Gel\_Band\_Peptide\_Quants/, and the FASTA file path).

