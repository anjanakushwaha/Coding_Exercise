

import argparse
import csv
import logging
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


# CONFIGURATION 

DEFAULT_COLUMNS = {
    "chromosome": 0,
    "start": 1,
    "end": 2,
    "gene": 6,
    "percentage30": 10,
    "sample": 11,
}


# EXON CLASS DATA


@dataclass
class Exon:
    gene: str
    chromosome: str
    start: int
    end: int
    percentage: float

    @property
    def length(self):
        # Assumes 1-based inclusive genomic coordinates.
        return self.end - self.start + 1

    @property
    def coordinate(self):
        return f"{self.chromosome}:{self.start}-{self.end}"


# GENE NAME


def parse_gene_symbol(value):


    if not value:
        return None

    gene = value.split(";")[0].strip()
    #Return None for missing or placeholder symbols.
    if not gene or gene in {".", "-", "NA", "N/A", "GENE_UNKNOWN"}:
        return None

    return gene


 # COVERAGE VALUE

def parse_percentage(value):
    """Convert percentage values into numbers between 0 and 100."""

    if value is None:
        raise ValueError("Missing percentage value")

    cleaned = value.strip().replace("%", "")

    if not cleaned or cleaned.upper() in {"NA", "N/A", ".", "NULL"}:
        raise ValueError("Missing or invalid percentage")

    percentage = float(cleaned)

    if not 0 <= percentage <= 100:
        raise ValueError("Percentage must be between 0 and 100")

    return percentage


#COLUMNS VALUE


def normalize_column_name(value):
    #names easier to compare

    return re.sub(r"[^a-z0-9]", "", value.lower())


def detect_columns(header_line, threshold):
    #Detect columns from a properly tab-separated header.


    fields = header_line.rstrip("\r\n").split("\t")

    normalized = [
        normalize_column_name(field.lstrip("#"))
        for field in fields
    ]

    required = {
        "chromosome": {"chromosome", "chr"},
        "start": {"startposition", "start"},
        "end": {"endposition", "end"},
        "gene": {"genesymbolaccession", "genesymbol"},
        "percentage": {f"percentage{threshold:g}".replace(".", "")},
        "sample": {"samplename"},
    }

    detected = {}

    for name, possible_names in required.items():
        matches = [
            i for i, column in enumerate(normalized)
            if column in possible_names
        ]

        if matches:
            detected[name] = matches[0]

    if all(name in detected for name in required):
        return {
            "chromosome": detected["chromosome"],
            "start": detected["start"],
            "end": detected["end"],
            "gene": detected["gene"],
            "percentage": detected["percentage"],
            "sample": detected["sample"],
        }

    return None


#READ PANEL GENES


def read_panel_genes(panel_file):
    #Read the list of genes that can be assessed

    genes = []

    with open(panel_file, "r", encoding="utf-8-sig") as file:

        for line_number, line in enumerate(file, start=1):

            value = line.strip()

            if not value or value.startswith("#"):
                continue

            # Allow a simple one-column gene list.
            gene = parse_gene_symbol(value)

            if gene is None:
                logging.warning(
                    "Invalid panel gene at line %d: %s",
                    line_number, value
                )
                continue

            if gene not in genes:
                genes.append(gene)

    if not genes:
        raise ValueError("The panel gene list is empty.")

    return genes


# ---------------------------------------------------------
# 7. READ SAMBAMBA COVERAGE DATA
# ---------------------------------------------------------

def read_coverage_file(input_file, threshold):

    #Read Sambamba output.

    coverage_by_gene = defaultdict(dict)
    rows_seen = defaultdict(int)
    invalid_rows = defaultdict(int)

    sample_names = set()
    columns = None
    header_seen = False

    percentage_column = (
        f"percentage{threshold:g}"
    )

    with open(input_file, "r", encoding="utf-8-sig") as file:

        for line_number, line in enumerate(file, start=1):

            if not line.strip():
                continue

            # Identify the header.
            if not header_seen and "chromosome" in line.lower():
                columns = detect_columns(line, threshold)

                if columns is not None:
                    header_seen = True
                    logging.info("Detected Sambamba header.")
                    continue

                # Header may not be correctly tab-separated.
                if len(line.rstrip("\r\n").split("\t")) < 12:
                    header_seen = True
                    columns = DEFAULT_COLUMNS.copy()
                    columns["percentage"] = 10

                    if threshold != 30:
                        raise ValueError(
                            "The supplied fallback data contains "
                            "percentage30 only. A different threshold "
                            "requires the corresponding percentage "
                            "column in the input file."
                        )

                    logging.warning(
                        "Header could not be parsed. "
                        "Using the known data-row column positions."
                    )
                    continue

            if columns is None:
                columns = DEFAULT_COLUMNS.copy()
                columns["percentage"] = 10

                if threshold != 30:
                    raise ValueError(
                        f"Cannot find {percentage_column} in the "
                        "header. A corresponding percentage column "
                        "is required for this threshold."
                    )

            fields = line.rstrip("\r\n").split("\t")

            required_index = max(columns.values())

            if len(fields) <= required_index:
                logging.warning(
                    "Line %d has too few fields; skipped.",
                    line_number
                )
                continue

            gene = parse_gene_symbol(fields[columns["gene"]])

            if gene is None:
                logging.warning(
                    "Missing or invalid gene symbol at line %d.",
                    line_number
                )
                continue

            rows_seen[gene] += 1

            try:
                chromosome = fields[columns["chromosome"]].strip()
                start = int(fields[columns["start"]].strip())
                end = int(fields[columns["end"]].strip())

                if not chromosome or start < 0 or end < start:
                    raise ValueError("Invalid exon coordinates")

                percentage = parse_percentage(
                    fields[columns["percentage"]]
                )

            except (ValueError, IndexError) as error:

                invalid_rows[gene] += 1

                logging.warning(
                    "Invalid record at line %d for %s: %s",
                    line_number, gene, error
                )
                continue

            if "sample" in columns:
                sample = fields[columns["sample"]].strip()

                if sample:
                    sample_names.add(sample)

            exon = Exon(
                gene=gene,
                chromosome=chromosome,
                start=start,
                end=end,
                percentage=percentage,
            )

            # Use genomic coordinates to identify a unique exon.
            # This avoids counting the same exon repeatedly, when it appears under multiple transcript accessions.
            exon_key = (
                chromosome,
                start,
                end,
            )

            existing = coverage_by_gene[gene].get(exon_key)

            if existing is None:
                coverage_by_gene[gene][exon_key] = exon

            elif existing.percentage != percentage:
                logging.warning(
                    "Different coverage values for duplicate exon "
                    "%s in %s. Retaining the lower value.",
                    exon.coordinate, gene
                )

                if percentage < existing.percentage:
                    coverage_by_gene[gene][exon_key] = exon

    if not sample_names:
        sample_name = Path(input_file).name.split("_")[0]
    elif len(sample_names) == 1:
        sample_name = next(iter(sample_names))
    else:
        raise ValueError(
            "Multiple sample names detected in the coverage file."
        )

    return coverage_by_gene, rows_seen, invalid_rows, sample_name


#CALCULATE GENE-LEVEL COVERAGE

def summarise_gene(exons, invalid_count=0):

    
    # Calculate gene-level coverage statistics.
    # A gene fails if at least one valid exon has percentage30 below 100%.

    if not exons:
        return {
            "TotalExons": 0,
            "FailingExons": 0,
            "MinCoverage": None,
            "WorstExon": "",
            "WeightedMean": None,
            "Status": "UNASSESSABLE",
        }

    sorted_exons = sorted(
        exons,
        key=lambda exon: (
            exon.percentage,
            exon.chromosome,
            exon.start,
            exon.end,
        )
    )

    worst_exon = sorted_exons[0]

    failing_exons = [
        exon for exon in exons
        if exon.percentage < 100
    ]

    total_length = sum(exon.length for exon in exons)

    weighted_mean = (
        sum(
            exon.length * exon.percentage
            for exon in exons
        ) / total_length
    )

    if failing_exons:
        status = "FAIL"

    elif invalid_count > 0:
        status = "UNASSESSABLE"

    else:
        status = "PASS"

    return {
        "TotalExons": len(exons),
        "FailingExons": len(failing_exons),
        "MinCoverage": worst_exon.percentage,
        "WorstExon": worst_exon.coordinate,
        "WeightedMean": weighted_mean,
        "Status": status,
    }



#GENERATE GENE-LEVEL REPORT


def generate_gene_report(
    panel_genes,
    coverage_by_gene,
    invalid_rows,
    output_file,
):

    results = {}

    with open(
        output_file,
        "w",
        newline="",
        encoding="utf-8"
    ) as file:

        writer = csv.writer(file, delimiter="\t")

        writer.writerow([
            "GeneSymbol",
            "TotalExons",
            "FailingExons",
            "MinCoverage%",
            "WorstExon(chr:start-end)",
            "Status",
            "WeightedMeanCoverage%",
        ])

        for gene in panel_genes:

            exons = list(coverage_by_gene.get(gene, {}).values())

            if not exons:
                if gene not in coverage_by_gene:
                    if gene in invalid_rows:
                        status = "UNASSESSABLE"
                    else:
                        status = "MISSING"
                else:
                    status = "UNASSESSABLE"

                result = {
                    "TotalExons": 0,
                    "FailingExons": 0,
                    "MinCoverage": None,
                    "WorstExon": "",
                    "WeightedMean": None,
                    "Status": status,
                }

            else:
                result = summarise_gene(
                    exons,
                    invalid_rows.get(gene, 0)
                )

            results[gene] = result

            writer.writerow([
                gene,
                result["TotalExons"],
                result["FailingExons"],
                (
                    f"{result['MinCoverage']:.4f}"
                    if result["MinCoverage"] is not None
                    else ""
                ),
                result["WorstExon"],
                result["Status"],
                (
                    f"{result['WeightedMean']:.4f}"
                    if result["WeightedMean"] is not None
                    else ""
                ),
            ])

    return results


#GENERATE FAILING EXON REPORT

def generate_failing_exon_report(
    panel_genes,
    coverage_by_gene,
    output_file,
):

    with open(
        output_file,
        "w",
        newline="",
        encoding="utf-8"
    ) as file:

        writer = csv.writer(file, delimiter="\t")

        writer.writerow([
            "GeneSymbol",
            "Chromosome",
            "Start",
            "End",
            "Coordinate",
            "Percentage30",
            "Status",
        ])

        for gene in panel_genes:

            exons = coverage_by_gene.get(gene, {}).values()

            for exon in exons:

                if exon.percentage < 100:

                    writer.writerow([
                        gene,
                        exon.chromosome,
                        exon.start,
                        exon.end,
                        exon.coordinate,
                        f"{exon.percentage:.4f}",
                        "FAIL",
                    ])



#PRINT SAMPLE SUMMARY


def print_summary(sample_name, results, total_panel_genes):

    failing = sum(
        1 for result in results.values()
        if result["Status"] == "FAIL"
    )

    missing = sum(
        1 for result in results.values()
        if result["Status"] == "MISSING"
    )

    unassessable = sum(
        1 for result in results.values()
        if result["Status"] == "UNASSESSABLE"
    )

    passing = sum(
        1 for result in results.values()
        if result["Status"] == "PASS"
    )

    print("\n" + "=" * 55)
    print("NGS COVERAGE REPORT")
    print("=" * 55)

    print(
        f"Sample {sample_name}: "
        f"{failing}/{total_panel_genes} genes "
        "below 30x threshold"
    )

    print(f"Total panel genes: {total_panel_genes}")
    print(f"Passing genes: {passing}")
    print(f"Failing genes: {failing}")
    print(f"Missing genes: {missing}")
    print(f"Unassessable genes: {unassessable}")

    print("=" * 55)


# MAIN PROGRAM


def main():

    parser = argparse.ArgumentParser(
        description="Analyse Sambamba NGS coverage "
                    "for a clinical gene panel."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Path to Sambamba coverage output file."
    )

    parser.add_argument(
        "--panel",
        required=True,
        help="Path to panel gene list."
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=30,
        help="Read-depth threshold. Default: 30."
    )

    parser.add_argument(
        "--output-dir",
        default="output",
        help="Directory for generated reports."
    )

    parser.add_argument(
        "--log-file",
        default="coverage.log",
        help="Path to the log file."
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
        handlers=[
            logging.FileHandler(args.log_file, encoding="utf-8"),
            logging.StreamHandler(),
        ],
    )

    if args.threshold < 0:
        parser.error("--threshold must not be negative.")

    input_file = Path(args.input)
    panel_file = Path(args.panel)
    output_dir = Path(args.output_dir)

    if not input_file.is_file():
        parser.error(f"Input file not found: {input_file}")

    if not panel_file.is_file():
        parser.error(f"Panel file not found: {panel_file}")

    output_dir.mkdir(parents=True, exist_ok=True)

    try:

        logging.info("Reading panel gene list...")

        panel_genes = read_panel_genes(panel_file)

        logging.info(
            "Panel contains %d genes.",
            len(panel_genes)
        )

        logging.info("Reading Sambamba coverage file...")

        (
            coverage_by_gene,
            rows_seen,
            invalid_rows,
            sample_name,
        ) = read_coverage_file(
            input_file,
            args.threshold
        )

        logging.info(
            "Coverage records found for %d genes.",
            len(rows_seen)
        )

        gene_report = output_dir / "gene_coverage_report.tsv"

        failing_report = output_dir / "failing_exons.tsv"

        summary_report = output_dir / "summary.txt"

        logging.info("Generating gene-level report...")

        results = generate_gene_report(
            panel_genes,
            coverage_by_gene,
            invalid_rows,
            gene_report,
        )

        logging.info("Generating failing exon report...")

        generate_failing_exon_report(
            panel_genes,
            coverage_by_gene,
            failing_report,
        )

        print_summary(
            sample_name,
            results,
            len(panel_genes),
        )

        with open(
            summary_report,
            "w",
            encoding="utf-8"
        ) as file:

            for gene, result in results.items():

                if result["Status"] != "PASS":

                    file.write(
                        f"{gene}\t{result['Status']}\n"
                    )

        logging.info("Analysis completed successfully.")

        print("\nGenerated files:")
        print(gene_report)
        print(failing_report)
        print(summary_report)
        print(args.log_file)

    except (OSError, ValueError, csv.Error) as error:

        logging.exception("Analysis failed: %s", error)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
