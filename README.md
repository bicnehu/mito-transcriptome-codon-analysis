# mito-transcriptome-codon-analysis

## 🧬 Overview

This repository contains the scripts and test data for performing comprehensive **codon usage analysis** on coding sequences (CDS) extracted from both **mitochondrial** and **transcriptomic** sequencing data.

The project is designed to automate the process of preparing sequences and calculating key metrics used in molecular evolution and gene expression studies.

---

## ✨ Key Features

The analysis pipeline focuses on:

1.  **CDS Extraction:** Programmatic extraction of Coding Sequences (CDS) from raw mitochondrial and transcriptome data.
2.  **Codon Usage Analysis:** Calculation of various codon usage indices, including:
    * **RSCU** (Relative Synonymous Codon Usage)
    * **ENC** (Effective Number of Codons)
    * **CAI** (Codon Adaptation Index)
3.  **Sequence Composition:** Calculation of overall **GC content**.

---

## 📂 Repository Contents

| File Name | Description |
| :--- | :--- |
| `mito_python_file.zip` | Compressed Python script(s) specifically for processing **mitochondrial** data and performing codon analysis. |
| `Transcriptome_python_file.zip` | Compressed Python script(s) specifically for processing **transcriptome** data and performing codon analysis. |
| `5test_mitochondria_data.zip` | Small-scale test data set (e.g., FASTA or similar format) for validating the mitochondrial analysis script. |
| `5test_transcriptome_data.zip`| Small-scale test data set (e.g., FASTA or similar format) for validating the transcriptome analysis script. |
| `.gitignore` | Standard file to specify intentionally untracked files to ignore. |

---

## 🚀 Getting Started (Inferred Workflow)

To replicate or run the analysis:

1.  **Clone the Repository:**
    ```bash
    git clone [https://github.com/bicnehu/mito-transcriptome-codon-analysis.git](https://github.com/bicnehu/mito-transcriptome-codon-analysis.git)
    cd mito-transcriptome-codon-analysis
    ```
2.  **Unzip Files:** Unzip the necessary python script(s) and data files.
    ```bash
    unzip mito_python_file.zip
    unzip Transcriptome_python_file.zip
    unzip 5test_mitochondria_data.zip
    unzip 5test_transcriptome_data.zip
    ```
3.  **Execute the Scripts:** Run the extracted Python scripts, providing the appropriate input data files (`5test_mitochondria_data` or `5test_transcriptome_data`) as input.
    * *(Note: Specific command-line arguments will depend on the implementation details within the Python files.)*

---

## 🤝 Contribution & Contact

For questions, issues, or collaboration, please open an issue on this GitHub repository.
