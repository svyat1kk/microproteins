This repository contains datasets and programs essential for the thesis **"Identifying Microproteins with Genetic Association Data"**
## Repository Structure

### `candidates_seeking/`
Intergration of genetic association data to determine potential involvment of microproteins in common human diseases/traits. The candidates are prioritized based on CADD score.
- **`candidate_detection.ipynb`** – main program for data exploration and selecting strong candidate
  Outputs:
  - `manhattan_plot/initial_632656.txt`: standardized SNVs from GWAS Catalog (632,656 SNVs)
  - `manhattan_plot/overlaying_2838.txt`: SNVs located in microprotein-encoding region studied in this work
  - `amino_acid_substitution/strong_candidates_14.txt`: 14 candidates that are testen for amino acid change in their sequences  

### `manhattan_plot/`
 - **`SNV_manhhattan_plot.R`** - R script to generate Manhattan plot figure.
  Outputs:
  - `Manhattan.png`

### `amino_acid_substitution/`
Test for amino acid substitution in sORFs caused by SNV's risk allele
 - **`risk_allele_consequences.py`** – Jupyter notebook for data exploration and candidate selection.  
   Outputs:
   - `mutated_sequences_14.txt`: wild-type and altered sORFs sequences