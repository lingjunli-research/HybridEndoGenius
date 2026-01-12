<div align="center">
<img src="https://github.com/user-attachments/assets/9c44a28d-2e48-4959-904f-ca571fc44af3">
</div>

# HybridEndoGenius

HybridEndoGenius is a **hybrid neuropeptide identification workflow** designed to improve the discovery of endogenous peptides from mass spectrometry (MS/MS) data. It integrates **database-driven** and **de novo–assisted** strategies to enable sensitive and flexible neuropeptide identification across different species and experimental settings.

The workflow is designed to be **modular, reproducible, and scalable**, and is compatible with high-throughput computing environments such as **CHTC/HTCondor**, or running via GUI.

---

## Key Features

- Hybrid peptide identification strategy combining:
  - Database search–based identification
  - De novo sequencing–assisted discovery
- Optimized for **neuropeptides**, which are typically:
  - Short
  - Low abundance
  - Poorly annotated in standard protein databases
- Flexible support for **species-specific configurations**
- Designed for **HTC execution** (HTCondor / DAGMan)

---

## Workflow Overview

At a high level, HybridEndoGenius performs the following steps:

1. **Input preparation**
   - MS/MS data (e.g. `.mzML`, `.ms2`, `.mgf`)
   - FASTA and CSV databases
2. **Database search**
   - Identification of known neuropeptides
3. **De novo sequencing**
   - Discovery of novel peptide candidates
4. **Filtering and post-processing**
   - Removal of peptides with low confidence
   - Identification of putative novel neuropeptides
5. **Result integration**
   - Database search using the FASTA generated from de novo sequencing


#### Getting started
* [User manual](https://docs.google.com/document/d/e/2PACX-1vS6e_OLpXPOV4FzfopoFB0mw024idh5BFgwWuUPmiqoiRaYqZrwkBsiOCNXqzDVY3e3VS_vP8jfgBah/pub)

#### Key references
EndoGenius:
* Fields, L.; Vu, N. Q.; Dang, T.C.; Yen, H.; Ma, M.; Wu, W.; Gray, M.; Li, L. (2024). EndoGenius: Optimized Neuropeptide Identification from Mass Spectrometry Datasets. Journal of Proteome Research. [Link](https://pubs.acs.org/doi/full/10.1021/acs.jproteome.3c00758)
* Fields, L.; Dang, T.; Gray, M.; Protya, S. S.; Li, L. (2025). EndoGenius: Enabling comprehensive identification and quantitation of endogenous peptides. bioRxiv. [Link](https://www.biorxiv.org/content/10.1101/2025.06.12.659347v1.abstract)

MotifQuest:
* Dang, T.; Fields, L.; Li, L. (2024). MotifQuest: An Automated Pipeline for Motif Database Creation to Improve Peptidomics Database Searching Programs. [Link](https://pubs.acs.org/doi/10.1021/jasms.4c00192)
---
