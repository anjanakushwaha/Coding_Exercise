# NGS Coverage Report – Congenital Myopathy & Muscular Dystrophy Panel

A Python-based bioinformatics tool for assessing next-generation sequencing (NGS) coverage across an 84-gene congenital myopathy and congenital muscular dystrophy panel.

The project analyses Sambamba coverage output, identifies genes that do not meet the minimum 30× coverage requirement, and reports failing exon coordinates.

## Click Live Results Dashboard

**[View the NGS Coverage Results Dashboard](https://anjanakushwaha.github.io/Coding_Exercise/)**

The dashboard displays the pre-generated analysis results, including:
- Gene-level coverage summary
- Genes failing the 30× coverage criterion
- Missing or unassessable genes
- Failing exon coordinates
- Overall sample summary

No file upload is required to view the published results.

## Project Objective

The aim is to assess whether every coding base across the target panel achieves at least 30× sequencing depth.

**Coverage criterion:** Every coding base must have a sequencing depth of ≥30×.

The report uses the Sambamba `percentage30` metric, which represents the percentage of bases within a region covered by at least 30 reads.

A gene is classified as failing if any reported exon has less than 100% of its bases covered at 30×.

## Workflow

1. **Input:** Sambamba per-exon coverage output.
2. **Parsing:** Read coverage data and extract gene symbols and exon coordinates.
3. **Aggregation:** Aggregate exon-level coverage metrics by gene.
4. **Quality assessment:** Identify genes and exons failing the 30× coverage criterion.
5. **Panel comparison:** Identify panel genes missing from the coverage output.
6. **Reporting:** Generate gene-level and exon-level reports.
7. **Visualisation:** Display the generated results through a static web dashboard.

## Repository Structure

```text
InterviewCHU/
│
├── coverage.py
│
├── panel_genes.txt
│
├── NGS148_34_139558_CB_CMCMD_S33_R1_001.sambamba_output.txt
│
├── output/
│   ├── gene_coverage_report.tsv
│   ├── failing_exons.tsv
│   ├── summary.txt
│   └── index.html
│
├── docs/
│   ├── index.html
│   ├── gene_coverage_report.tsv
│   ├── failing_exons.tsv
│   └── summary.txt
│
└── README.md
```

## Running the Analysis

### Requirements

- Python 3
- Python standard library

No external Python packages are required for the report-generation script.

### Command

```bash
  --python3 coverage.py \
  --input "NGS148_34_139558_CB_CMCMD_S33_R1_001.sambamba_output.txt" \
  --panel "panel_genes.txt" \
  --threshold 30 \
  --output-dir output
```

### Generated Reports

| File | Description |
|---|---|
| `gene_coverage_report.tsv` | Gene-level coverage assessment |
| `failing_exons.tsv` | Exon-level coverage failures |
| `summary.txt` | Overall analysis summary |
| `index.html` | Interactive local results viewer |

## Key Technical Considerations

- **Minimum exon coverage:** The minimum percentage of bases covered at 30× is used to determine whether a gene passes the coverage criterion.
- **Transcript handling:** Multiple transcript accessions are grouped under their corresponding gene symbol.
- **Missing genes:** Panel genes without coverage records are identified separately.
- **Data quality:** Malformed records and non-numeric coverage values are handled through validation and warnings.
- **Clinical interpretation:** This tool supports coverage assessment and reporting; it does not independently establish clinical validity or diagnostic findings.

## Technologies

- Python 3
- Next-generation sequencing (NGS)
- Sambamba coverage metrics
- TSV data processing
- HTML, CSS and JavaScript
- GitHub Pages


## Author

**Anjana Kushwaha**

Bioinformatics | Computational Biology | NGS Data Analysis | Machine Learning

---

*Developed as part of a bioinformatics coding exercise for translational genomic analysis.*
