#!/usr/bin/env python3
"""Aggregate and compare outputs produced by run_cobra_model.py."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUTPUT_CHOICES = (
    "growth-table",
    "summary-long",
    "summary-matrix",
    "summary-heatmap",
    "reaction-long",
    "reaction-matrix",
    "reaction-heatmap",
    "all",
)

ALL_OUTPUTS = tuple(
    choice for choice in OUTPUT_CHOICES
    if choice != "all"
)


def build_parser() -> argparse.ArgumentParser:
    """Build command-line argument parser."""

    parser = argparse.ArgumentParser(
        description=(
            "Combine per-model COBRA result tables into comparison "
            "tables and heat maps."
        )
    )

    parser.add_argument(
        "--results",
        type=Path,
        nargs="+",
        required=True,
        help=(
            "Result files or directories containing *.model.tsv, "
            "*.summary.tsv, and *.fluxes.tsv."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for combined tables and figures.",
    )

    parser.add_argument(
        "--outputs",
        choices=OUTPUT_CHOICES,
        nargs="+",
        default=[
            "growth-table",
            "summary-matrix",
            "summary-heatmap",
        ],
        help="Products to create. Use 'all' for every product.",
    )

    parser.add_argument(
        "--row-columns",
        nargs="+",
        default=["model_name"],
        help=(
            "Columns used to identify matrix/heat-map rows. "
            "Example: --row-columns sample_id medium"
        ),
    )

    parser.add_argument(
        "--row-separator",
        default=" | ",
        help="Separator used when combining row identifiers.",
    )

    parser.add_argument(
        "--summary-column",
        choices=(
            "metabolite",
            "canonical_metabolite",
            "modelseed_id",
            "reaction",
        ),
        default="metabolite",
        help=(
            "Column used for model.summary() matrices. "
            "canonical_metabolite compares metabolites across "
            "reconstruction systems using ModelSEED identifiers."
        ),
    )

    parser.add_argument(
        "--summary-boundary-types",
        nargs="+",
        choices=("exchange", "demand", "sink", "boundary"),
        default=["exchange"],
        help="Boundary reaction types included in summary products.",
    )

    parser.add_argument(
        "--min-abs-flux",
        type=float,
        default=1e-9,
        help="Fluxes below this absolute value are treated as zero.",
    )

    parser.add_argument(
        "--top-n",
        type=int,
        default=50,
        help=(
            "Maximum number of heat-map columns, selected by variance. "
            "Use 0 for all."
        ),
    )

    parser.add_argument(
        "--heatmap-transform",
        choices=("none", "signed-log1p"),
        default="none",
        help="Optional display-only transformation for heat maps.",
    )

    parser.add_argument(
        "--figure-format",
        choices=("png", "pdf", "svg"),
        default="png",
        help="Heat-map file format.",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Raster image resolution.",
    )

    parser.add_argument(
        "--include-nonoptimal",
        action="store_true",
        help="Include non-optimal models in flux comparisons.",
    )

    parser.add_argument(
        "--modelseed-compounds",
        type=Path,
        default=None,
        help="Path to ModelSEED Biochemistry/compounds.tsv.",
    )

    parser.add_argument(
        "--modelseed-compound-aliases",
        type=Path,
        default=None,
        help=(
            "Path to ModelSEED "
            "Biochemistry/Aliases/"
            "Unique_ModelSEED_Compound_Aliases.txt."
        ),
    )

    return parser


def read_results(
    paths: list[Path],
    suffix: str,
) -> pd.DataFrame:
    """Read and combine result tables matching a filename suffix."""

    files = []

    for path in paths:
        if not path.exists():
            raise FileNotFoundError(
                f"Result path does not exist: {path}"
            )

        if path.is_dir():
            files.extend(path.rglob(f"*{suffix}"))

        elif path.name.endswith(suffix):
            files.append(path)

    files = sorted(set(files))

    if not files:
        raise ValueError(
            f"No files ending in {suffix} were found."
        )

    return pd.concat(
        (
            pd.read_csv(path, sep="\t")
            for path in files
        ),
        ignore_index=True,
        sort=False,
    )


def read_modelseed_compounds(
    compounds_path: Path,
    aliases_path: Path,
) -> pd.DataFrame:
    """Build a lookup from ModelSEED and BiGG IDs to compound metadata."""

    compounds = pd.read_csv(
        compounds_path,
        sep="\t",
        dtype=str,
        usecols=[
            "id",
            "name",
            "abbreviation",
            "formula",
            "charge",
        ],
    ).rename(
        columns={
            "id": "modelseed_id",
            "name": "metabolite_name",
            "abbreviation": "metabolite_abbreviation",
        }
    )

    aliases = pd.read_csv(
        aliases_path,
        sep="\t",
        dtype=str,
        usecols=[
            "ModelSEED ID",
            "External ID",
            "Source",
        ],
    )

    # Keep BiGG aliases only.
    is_bigg = aliases["Source"].str.contains(
        r"(?:^|\|)BiGG1?(?:\||$)",
        regex=True,
        na=False,
    )

    aliases = (
        aliases.loc[
            is_bigg,
            ["ModelSEED ID", "External ID"],
        ]
        .rename(
            columns={
                "ModelSEED ID": "modelseed_id",
                "External ID": "lookup_id",
            }
        )
    )

    aliases["mapping_source"] = "BiGG"

    # ModelSEED IDs themselves are valid lookup IDs.
    direct = compounds[["modelseed_id"]].copy()
    direct["lookup_id"] = direct["modelseed_id"]
    direct["mapping_source"] = "ModelSEED"

    lookup = pd.concat(
        [
            direct,
            aliases[
                [
                    "modelseed_id",
                    "lookup_id",
                    "mapping_source",
                ]
            ],
        ],
        ignore_index=True,
    ).drop_duplicates()

    # Do not automatically use aliases that map to multiple
    # ModelSEED compounds.
    mapping_counts = (
        lookup.groupby("lookup_id")["modelseed_id"]
        .transform("nunique")
    )

    lookup = lookup[mapping_counts == 1]

    return lookup.merge(
        compounds,
        on="modelseed_id",
        how="left",
        validate="many_to_one",
    )


def compound_lookup_id(
    metabolite_id: str,
) -> str:
    """Remove compartment notation from a metabolite identifier."""

    identifier = re.sub(r"^M_", "", metabolite_id)

    # gapseq / ModelSEED:
    # cpd00029_c0 -> cpd00029
    if identifier.startswith("cpd"):
        return identifier.split("_", 1)[0]

    # BiGG:
    # ac_e     -> ac
    # glc__D_c -> glc__D
    return re.sub(
        r"_[a-z]\d*$",
        "",
        identifier,
    )

def metabolite_compartment(metabolite_id: str) -> str:
    """Return a normalized compartment code from a metabolite ID."""

    identifier = re.sub(r"^M_", "", metabolite_id)

    match = re.search(r"_([a-z])\d*$", identifier)

    return match.group(1) if match else ""

def annotate_metabolites(
    summary: pd.DataFrame,
    compounds: pd.DataFrame,
) -> pd.DataFrame:
    """Annotate metabolites and create a common comparison identifier."""

    annotated = summary.copy()

    annotated["lookup_id"] = (
        annotated["metabolite"]
        .astype(str)
        .map(compound_lookup_id)
    )

    annotated["compartment"] = (
        annotated["metabolite"]
        .astype(str)
        .map(metabolite_compartment)
    )

    annotated = annotated.merge(
        compounds,
        on="lookup_id",
        how="left",
        validate="many_to_one",
    )

    mapped = annotated["modelseed_id"].notna()

    annotated["canonical_metabolite"] = annotated["metabolite"]

    annotated.loc[mapped, "canonical_metabolite"] = (
        annotated.loc[mapped, "modelseed_id"]
        + annotated.loc[mapped, "compartment"].apply(
            lambda compartment: (
                f"_{compartment}"
                if compartment
                else ""
            )
        )
    )

    annotated["mapping_source"] = (
        annotated["mapping_source"]
        .fillna("unmapped")
    )

    annotation_columns = [
        "canonical_metabolite",
        "modelseed_id",
        "metabolite_name",
        "metabolite_abbreviation",
        "formula",
        "charge",
        "compartment",
        "mapping_source",
    ]

    for column in (
        "modelseed_id",
        "metabolite_name",
        "metabolite_abbreviation",
        "formula",
        "charge",
    ):
        annotated[column] = annotated[column].fillna("")

    annotated = annotated.drop(columns="lookup_id")

    columns = annotated.columns.tolist()

    for column in annotation_columns:
        columns.remove(column)

    position = columns.index("metabolite") + 1

    for column in reversed(annotation_columns):
        columns.insert(position, column)

    return annotated.loc[:, columns]

def normalize_flux(
    frame: pd.DataFrame,
    threshold: float,
) -> pd.DataFrame:
    """Convert flux to numeric and set very small values to zero."""

    result = frame.copy()

    result["flux"] = pd.to_numeric(
        result["flux"],
        errors="coerce",
    )

    result.loc[
        result["flux"].abs() < threshold,
        "flux",
    ] = 0.0

    return result


def filter_optimal(
    frame: pd.DataFrame,
    result_ids: set[str] | None,
) -> pd.DataFrame:
    """Restrict a result table to optimal model runs."""

    if result_ids is None:
        return frame

    return frame[
        frame["result_id"]
        .astype(str)
        .isin(result_ids)
    ].copy()


def make_matrix(
    frame: pd.DataFrame,
    row_columns: list[str],
    column: str,
    separator: str,
) -> pd.DataFrame:
    """Convert a long-form flux table into a model-by-feature matrix."""

    working = frame.copy()

    working["row_id"] = (
        working[row_columns]
        .fillna("")
        .astype(str)
        .agg(separator.join, axis=1)
    )

    matrix = working.pivot_table(
        index="row_id",
        columns=column,
        values="flux",
        aggfunc="sum",
        fill_value=0.0,
    )

    matrix.index.name = separator.join(row_columns)

    return matrix.sort_index()


def select_heatmap_columns(
    matrix: pd.DataFrame,
    top_n: int,
) -> pd.DataFrame:
    """Keep active columns and optionally select the most variable."""

    active = matrix.loc[
        :,
        (matrix != 0).any(axis=0),
    ]

    if (
        active.empty
        or top_n == 0
        or active.shape[1] <= top_n
    ):
        return active

    variance = (
        active.var(axis=0)
        .sort_values(ascending=False)
    )

    return active.loc[
        :,
        variance.head(top_n).index,
    ]


def write_heatmap(
    matrix: pd.DataFrame,
    output_path: Path,
    title: str,
    transform: str,
    dpi: int,
) -> None:
    """Write a flux heat map."""

    if matrix.empty:
        raise ValueError(
            f"Cannot draw {title}: no active flux columns."
        )

    if transform == "signed-log1p":
        values = (
            np.sign(matrix)
            * np.log1p(np.abs(matrix))
        )

        plotted = pd.DataFrame(
            values,
            index=matrix.index,
            columns=matrix.columns,
        )

    else:
        plotted = matrix

    figure_width = max(
        8.0,
        0.32 * plotted.shape[1] + 3.0,
    )

    figure_height = max(
        4.0,
        0.35 * plotted.shape[0] + 2.0,
    )

    figure, axis = plt.subplots(
        figsize=(figure_width, figure_height)
    )

    maximum = float(
        np.nanmax(
            np.abs(plotted.to_numpy())
        )
    )

    if not np.isfinite(maximum) or maximum == 0:
        maximum = 1.0

    image = axis.imshow(
        plotted.to_numpy(),
        aspect="auto",
        cmap="RdBu_r",
        vmin=-maximum,
        vmax=maximum,
    )

    axis.set_xticks(
        np.arange(plotted.shape[1])
    )

    axis.set_xticklabels(
        plotted.columns,
        rotation=90,
        fontsize=8,
    )

    axis.set_yticks(
        np.arange(plotted.shape[0])
    )

    axis.set_yticklabels(
        plotted.index,
        fontsize=8,
    )

    axis.set_xlabel(
        plotted.columns.name or "Feature"
    )

    axis.set_ylabel(
        plotted.index.name or "Model"
    )

    axis.set_title(title)

    colorbar = figure.colorbar(
        image,
        ax=axis,
    )

    colorbar.set_label(
        "sign(flux) × log1p(|flux|)"
        if transform == "signed-log1p"
        else "Flux"
    )

    figure.tight_layout()

    figure.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
    )

    plt.close(figure)


def main() -> int:
    """Run COBRA result aggregation."""

    parser = build_parser()
    args = parser.parse_args()

    if args.min_abs_flux < 0:
        parser.error(
            "--min-abs-flux cannot be negative."
        )

    if args.top_n < 0:
        parser.error(
            "--top-n cannot be negative."
        )

    if args.dpi <= 0:
        parser.error(
            "--dpi must be positive."
        )

    # Either supply both ModelSEED files or neither.
    if bool(args.modelseed_compounds) != bool(
        args.modelseed_compound_aliases
    ):
        parser.error(
            "--modelseed-compounds and "
            "--modelseed-compound-aliases "
            "must be supplied together."
        )

    outputs = (
        set(ALL_OUTPUTS)
        if "all" in args.outputs
        else set(args.outputs)
    )

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ------------------------------------------------------------------
    # Optional ModelSEED compound annotations
    # ------------------------------------------------------------------

    compounds = None

    if args.modelseed_compounds:
        compounds = read_modelseed_compounds(
            args.modelseed_compounds,
            args.modelseed_compound_aliases,
        )

    # ------------------------------------------------------------------
    # Model metadata
    # ------------------------------------------------------------------

    metadata = read_results(
        args.results,
        ".model.tsv",
    )

    optimal_result_ids = None

    if not args.include_nonoptimal:
        optimal_result_ids = set(
            metadata.loc[
                metadata["status"] == "optimal",
                "result_id",
            ].astype(str)
        )

    # ------------------------------------------------------------------
    # Growth table
    # ------------------------------------------------------------------

    if "growth-table" in outputs:
        growth = metadata.copy()

        sort_columns = [
            column
            for column in (
                "status",
                "maximum_biomass",
            )
            if column in growth.columns
        ]

        if sort_columns:
            growth = growth.sort_values(
                sort_columns,
                ascending=[
                    column == "status"
                    for column in sort_columns
                ],
                na_position="last",
            )

        growth.to_csv(
            args.output_dir / "growth_table.tsv",
            sep="\t",
            index=False,
        )

    # ------------------------------------------------------------------
    # Boundary / metabolite summaries
    # ------------------------------------------------------------------

    summary_outputs = {
        "summary-long",
        "summary-matrix",
        "summary-heatmap",
    }

    if outputs & summary_outputs:
        summary = read_results(
            args.results,
            ".summary.tsv",
        )

        summary = normalize_flux(
            summary,
            args.min_abs_flux,
        )

        summary = filter_optimal(
            summary,
            optimal_result_ids,
        )

        summary = summary[
            summary["boundary_type"].isin(
                args.summary_boundary_types
            )
        ].copy()

        if compounds is not None:
            summary = annotate_metabolites(
                summary,
                compounds,
            )

        if "summary-long" in outputs:
            summary.to_csv(
                args.output_dir
                / "summary_long.tsv",
                sep="\t",
                index=False,
            )

        if outputs & {
            "summary-matrix",
            "summary-heatmap",
        }:
            summary_matrix = make_matrix(
                summary,
                row_columns=args.row_columns,
                column=args.summary_column,
                separator=args.row_separator,
            )

            if "summary-matrix" in outputs:
                summary_matrix.to_csv(
                    args.output_dir
                    / "summary_flux_matrix.tsv",
                    sep="\t",
                )

            if "summary-heatmap" in outputs:
                plotted = select_heatmap_columns(
                    summary_matrix,
                    args.top_n,
                )

                write_heatmap(
                    plotted,
                    args.output_dir
                    / (
                        "summary_flux_heatmap."
                        f"{args.figure_format}"
                    ),
                    title=(
                        "Boundary-metabolite fluxes\n"
                        "positive = uptake; "
                        "negative = secretion"
                    ),
                    transform=args.heatmap_transform,
                    dpi=args.dpi,
                )

    # ------------------------------------------------------------------
    # Reaction fluxes
    # ------------------------------------------------------------------

    reaction_outputs = {
        "reaction-long",
        "reaction-matrix",
        "reaction-heatmap",
    }

    if outputs & reaction_outputs:
        fluxes = read_results(
            args.results,
            ".fluxes.tsv",
        )

        fluxes = normalize_flux(
            fluxes,
            args.min_abs_flux,
        )

        fluxes = filter_optimal(
            fluxes,
            optimal_result_ids,
        )

        if "reaction-long" in outputs:
            fluxes.to_csv(
                args.output_dir
                / "reaction_flux_long.tsv",
                sep="\t",
                index=False,
            )

        if outputs & {
            "reaction-matrix",
            "reaction-heatmap",
        }:
            reaction_matrix = make_matrix(
                fluxes,
                row_columns=args.row_columns,
                column="reaction",
                separator=args.row_separator,
            )

            if "reaction-matrix" in outputs:
                reaction_matrix.to_csv(
                    args.output_dir
                    / "reaction_flux_matrix.tsv",
                    sep="\t",
                )

            if "reaction-heatmap" in outputs:
                plotted = select_heatmap_columns(
                    reaction_matrix,
                    args.top_n,
                )

                write_heatmap(
                    plotted,
                    args.output_dir
                    / (
                        "reaction_flux_heatmap."
                        f"{args.figure_format}"
                    ),
                    title="Reaction fluxes",
                    transform=args.heatmap_transform,
                    dpi=args.dpi,
                )

    # ------------------------------------------------------------------
    # Report generated files
    # ------------------------------------------------------------------

    print("Created outputs:")

    for path in sorted(
        args.output_dir.iterdir()
    ):
        if path.is_file():
            print(path)

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())

    except KeyboardInterrupt:
        print(
            "Interrupted.",
            file=sys.stderr,
        )

        raise SystemExit(130)