BioPython_Project:
IN-SILICO IDENTIFICATION AND FUNCTIONAL CHARACTERIZATION OF A TARGET PROTEIN 
USING SEQUENCE ANALYSIS AND HOMOLOGY-BASED ANNOTATION

PROJECT OVERVIEW

This project performs computational analysis of a protein sequence using Biopython, 
covering three major steps:

1. Sequence Quality Check (QC)
2. Homology Analysis using BLASTp
3. Functional Annotation based on BLAST results

The workflow mimics a basic bioinformatics pipeline used in research laboratories.

FOLDER STRUCTURE

- Functional_Sequence_Characterization/
  - Analysis/
    - Sequence_qc.py
    - Homology_analysis.py
  - Data/
    - Input_sequence.fasta
  - Results/
    - Qc_summary.txt
    - Blast_results.txt
    - Functional_annotation.txt

STEPS INVOLVED IN THIS PIPELINE

Step 1 – Sequence Quality Check  
The script Sequence_qc.py reads the protein from Input_sequence.fasta, calculates and reports:
- Protein ID and description (including organism and gene information)
- Full amino acid sequence in FASTA format
- Protein length (total number of amino acid residues)
- Estimated molecular weight based on standard amino acid mass values
- Amino acid composition (%) of all 20 standard residues
- Group-wise composition (Basic, Acidic, Aromatic, Polar, Non-polar)
It determines whether the sequence is suitable for further analysis and saves results in Qc_summary.txt.

Step 2 – Homology Analysis using BLASTp  
The script Homology_analysis.py submits the sequence to NCBI BLASTp (nr database), 
saves raw output in blast_result.xml, and extracts:
- Top homolog
- Percent identity
- E-value
- Alignment details
- Conserved regions
- Functional annotation
  
Results are stored in Blast_results.txt.

Step 3 – Functional Annotation
Based on the top BLAST hit, the protein is predicted to be a zinc-finger DNA-binding transcription factor, likely involved in regulating gene expression in humans.
The functional annotation is inferred using information from the UniProt database, where the top BLAST hit provides experimentally curated or computationally predicted functional data.
These results are documented in Functional_annotation.txt. 
