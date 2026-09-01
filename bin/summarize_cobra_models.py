#!/usr/bin/env python3
"""Aggregate outputs produced by run_cobra_model.py.

The requested products are selected with --outputs, allowing a workflow to make
only the tables or figures needed for a particular run.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

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
ALL_OUTPUTS = tuple(choice for choice in OUTPUT_CHOICES if choice != "all")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Combine per-model COBRA result TSVs into selected comparison tables "
            "and heat maps."
        )
    )
    parser.add_argument(
        "--results",
        type=Path,
        nargs="+",
        required=True,
        help=(
            "Result files or directories. Directories are searched recursively for "
            "*.model.tsv, *.summary.tsv, and *.fluxes.tsv."
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
        default=["growth-table", "summary-matrix", "summary-heatmap"],
        help=(
            "Products to create. Default: growth-table summary-matrix "
            "summary-heatmap. Use 'all' for every product."
        ),
    )
    parser.add_argument(
        "--row-columns",
        nargs="+",
        default=["model_name"],
        help=(
            "Columns that identify heat-map rows. Multiple columns are joined with "
            "--row-separator. Example: --row-columns sample_id medium"
        ),
    )
    parser.add_argument(
        "--row-separator",
        default=" | ",
        help="Separator used to combine multiple row columns. Default: ' | '.",
    )
    parser.add_argument(
        "--summary-column",
        choices=("metabolite", "reaction"),
        default="metabolite",
        help="Column dimension for model.summary() matrices. Default: metabolite.",
    )
    parser.add_argument(
        "--summary-boundary-types",
        nargs="+",
        choices=("exchange", "demand", "sink", "boundary"),
        default=["exchange"],
        help="Boundary types included in summary products. Default: exchange.",
    )
    parser.add_argument(
        "--min-abs-flux",
        type=float,
        default=1e-9,
        help="Fluxes below this absolute value are treated as zero. Default: 1e-9.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=50,
        help=(
            "Maximum number of columns shown in each heat map, selected by variance. "
            "Set 0 to show all. Default: 50."
        ),
    )
    parser.add_argument(
        "--heatmap-transform",
        choices=("none", "signed-log1p"),
        default="none",
        help=(
            "Optional display-only transformation. signed-log1p preserves sign while "
            "compressing large flux ranges. Default: none."
        ),
    )
    parser.add_argument(
        "--figure-format",
        choices=("png", "pdf", "svg"),
        default="png",
        help="Heat-map file format. Default: png.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Raster image resolution. Default: 300.",
    )
    parser.add_argument(
        "--include-nonoptimal",
        action="store_true",
        help=(
            "Keep non-optimal models in matrices as all-missing rows. By default only "
            "optimal models are included in flux comparisons."
        ),
    )
    parser.add_argument(
        "--bigg-metabolites",
        type=Path,
        default=None,
        help=(
            "Optional BiGG metabolite metadata table used to annotate "
            "metabolites in the combined summary outputs."
        ),
    )
    return parser


def discover_files(paths: Iterable[Path]) -> dict[str, list[Path]]:
    discovered = {"model": [], "summary": [], "fluxes": []}

    def classify(path: Path) -> None:
        name = path.name
        if name.endswith(".model.tsv"):
            discovered["model"].append(path)
        elif name.endswith(".summary.tsv"):
            discovered["summary"].append(path)
        elif name.endswith(".fluxes.tsv"):
            discovered["fluxes"].append(path)

    for path in paths:
        if path.is_dir():
            for candidate in path.rglob("*.tsv"):
                classify(candidate)
        elif path.is_file():
            classify(path)
        else:
            raise FileNotFoundError(f"Result path does not exist: {path}")

    for key in discovered:
        discovered[key] = sorted(set(item.resolve() for item in discovered[key]))
    return discovered


def read_many(paths: list[Path], kind: str) -> pd.DataFrame:
    if not paths:
        raise ValueError(f"No {kind} result files were found.")
    frames = [pd.read_csv(path, sep="\t") for path in paths]
    return pd.concat(frames, ignore_index=True, sort=False)

def read_bigg_metabolite_names(path: Path) -> pd.DataFrame:
    """Read BiGG metabolite identifiers and names.

    The standard BiGG metabolite export contains columns including
    ``bigg_id``, ``universal_bigg_id``, and ``name``. Model-specific BiGG
    identifiers are preferred because COBRApy summary outputs normally
    include compartment suffixes such as ``_c`` or ``_e``.

    Returns
    -------
    pandas.DataFrame
        A two-column table containing ``metabolite`` and ``metabolite_name``.
    """
    if not path.is_file():
        raise FileNotFoundError(
            f"BiGG metabolite metadata file does not exist: {path}"
        )

    bigg = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    required = {"bigg_id", "name"}
    missing = required.difference(bigg.columns)
    if missing:
        raise ValueError(
            "BiGG metabolite metadata file is missing required columns: "
            + ", ".join(sorted(missing))
        )

    names = (
        bigg.loc[:, ["bigg_id", "name"]]
        .rename(
            columns={
                "bigg_id": "metabolite",
                "name": "metabolite_name",
            }
        )
        .drop_duplicates(subset="metabolite")
    )

    return names

def annotate_metabolites(
    summary: pd.DataFrame,
    metabolite_names: pd.DataFrame,
) -> pd.DataFrame:
    """Add BiGG metabolite names to a COBRApy summary table."""
    if "metabolite" not in summary.columns:
        raise KeyError(
            "Summary inputs do not contain a metabolite column."
        )

    annotated = summary.merge(
        metabolite_names,
        on="metabolite",
        how="left",
        validate="many_to_one",
    )

    annotated["metabolite_name"] = annotated["metabolite_name"].fillna("")

    columns = annotated.columns.tolist()
    columns.remove("metabolite_name")
    metabolite_position = columns.index("metabolite") + 1
    columns.insert(metabolite_position, "metabolite_name")

    return annotated.loc[:, columns]

def make_row_id(frame: pd.DataFrame, columns: list[str], separator: str) -> pd.Series:
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise KeyError(
            "Requested --row-columns are absent from input tables: "
            + ", ".join(missing)
        )
    values = frame[columns].fillna("").astype(str)
    return values.agg(separator.join, axis=1)


def normalize_flux(frame: pd.DataFrame, threshold: float) -> pd.DataFrame:
    result = frame.copy()
    result["flux"] = pd.to_numeric(result["flux"], errors="coerce")
    result.loc[result["flux"].abs() < threshold, "flux"] = 0.0
    return result


def make_matrix(
    frame: pd.DataFrame,
    row_columns: list[str],
    column: str,
    separator: str,
) -> pd.DataFrame:
    working = frame.copy()
    working["row_id"] = make_row_id(working, row_columns, separator)
    matrix = working.pivot_table(
        index="row_id",
        columns=column,
        values="flux",
        aggfunc="sum",
        fill_value=0.0,
    )
    matrix.index.name = separator.join(row_columns)
    return matrix.sort_index()


def select_heatmap_columns(matrix: pd.DataFrame, top_n: int) -> pd.DataFrame:
    active = matrix.loc[:, (matrix != 0).any(axis=0)]
    if active.empty or top_n == 0 or active.shape[1] <= top_n:
        return active
    variances = active.var(axis=0).sort_values(ascending=False)
    return active.loc[:, variances.head(top_n).index]


def transform_for_heatmap(matrix: pd.DataFrame, transform: str) -> pd.DataFrame:
    if transform == "signed-log1p":
        values = np.sign(matrix) * np.log1p(np.abs(matrix))
        return pd.DataFrame(values, index=matrix.index, columns=matrix.columns)
    return matrix


def write_heatmap(
    matrix: pd.DataFrame,
    output_path: Path,
    title: str,
    transform: str,
    dpi: int,
) -> None:
    plotted = transform_for_heatmap(matrix, transform)
    if plotted.empty:
        raise ValueError(f"Cannot draw {title}: matrix has no active flux columns.")

    figure_width = max(8.0, 0.32 * plotted.shape[1] + 3.0)
    figure_height = max(4.0, 0.35 * plotted.shape[0] + 2.0)
    figure, axis = plt.subplots(figsize=(figure_width, figure_height))

    maximum = float(np.nanmax(np.abs(plotted.to_numpy())))
    if not np.isfinite(maximum) or maximum == 0:
        maximum = 1.0

    image = axis.imshow(
        plotted.to_numpy(),
        aspect="auto",
        cmap="RdBu_r",
        vmin=-maximum,
        vmax=maximum,
    )
    axis.set_xticks(np.arange(plotted.shape[1]))
    axis.set_xticklabels(plotted.columns, rotation=90, fontsize=8)
    axis.set_yticks(np.arange(plotted.shape[0]))
    axis.set_yticklabels(plotted.index, fontsize=8)
    axis.set_xlabel(plotted.columns.name or "Feature")
    axis.set_ylabel(plotted.index.name or "Model")
    axis.set_title(title)

    colorbar = figure.colorbar(image, ax=axis)
    label = "Flux"
    if transform == "signed-log1p":
        label = "sign(flux) × log1p(|flux|)"
    colorbar.set_label(label)

    figure.tight_layout()
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def requested_outputs(values: list[str]) -> set[str]:
    selected = set(values)
    if "all" in selected:
        return set(ALL_OUTPUTS)
    return selected


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.min_abs_flux < 0:
        parser.error("--min-abs-flux cannot be negative.")
    if args.top_n < 0:
        parser.error("--top-n cannot be negative.")
    if args.dpi <= 0:
        parser.error("--dpi must be positive.")

    outputs = requested_outputs(args.outputs)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    bigg_metabolite_names: pd.DataFrame | None = None

    if args.bigg_metabolites is not None:
        try:
            bigg_metabolite_names = read_bigg_metabolite_names(
                args.bigg_metabolites
            )
        except (FileNotFoundError, ValueError, pd.errors.ParserError) as error:
            parser.error(str(error))

        try:
            files = discover_files(args.results)
        except FileNotFoundError as error:
            parser.error(str(error))

    metadata: pd.DataFrame | None = None
    summary: pd.DataFrame | None = None
    fluxes: pd.DataFrame | None = None

    needs_metadata = bool(outputs) or not args.include_nonoptimal
    if needs_metadata:
        try:
            metadata = read_many(files["model"], "model metadata")
        except ValueError as error:
            parser.error(str(error))

    optimal_result_ids: set[str] | None = None
    if metadata is not None and not args.include_nonoptimal:
        if "status" not in metadata.columns or "result_id" not in metadata.columns:
            parser.error("Model metadata must contain status and result_id columns.")
        optimal_result_ids = set(
            metadata.loc[metadata["status"] == "optimal", "result_id"].astype(str)
        )

    if "growth-table" in outputs:
        growth = metadata.copy()
        sort_columns = [column for column in ("status", "maximum_biomass") if column in growth]
        if sort_columns:
            ascending = [True if column == "status" else False for column in sort_columns]
            growth = growth.sort_values(sort_columns, ascending=ascending, na_position="last")
        growth.to_csv(args.output_dir / "growth_table.tsv", sep="\t", index=False)

    summary_outputs = {
        "summary-long",
        "summary-matrix",
        "summary-heatmap",
    }
    if outputs.intersection(summary_outputs):
        try:
            summary = read_many(files["summary"], "model summary")
        except ValueError as error:
            parser.error(str(error))
        summary = normalize_flux(summary, args.min_abs_flux)
        if optimal_result_ids is not None:
            if "result_id" not in summary.columns:
                parser.error("Summary inputs do not contain a result_id column.")
            summary = summary.loc[
                summary["result_id"].astype(str).isin(optimal_result_ids)
            ]
        if "boundary_type" not in summary.columns:
            parser.error("Summary inputs do not contain a boundary_type column.")
        summary = summary.loc[
            summary["boundary_type"].isin(args.summary_boundary_types)
        ].copy()

        if bigg_metabolite_names is not None:
            try:
                summary = annotate_metabolites(
                    summary,
                    bigg_metabolite_names,
                )
            except (KeyError, pd.errors.MergeError) as error:
                parser.error(str(error))

        if "summary-long" in outputs:
            summary.to_csv(
                args.output_dir / "summary_long.tsv",
                sep="\t",
                index=False,
            )

        if outputs.intersection({"summary-matrix", "summary-heatmap"}):
            try:
                summary_matrix = make_matrix(
                    summary,
                    row_columns=args.row_columns,
                    column=args.summary_column,
                    separator=args.row_separator,
                )
            except KeyError as error:
                parser.error(str(error))

            if "summary-matrix" in outputs:
                summary_matrix.to_csv(
                    args.output_dir / "summary_flux_matrix.tsv",
                    sep="\t",
                )

            if "summary-heatmap" in outputs:
                plotted = select_heatmap_columns(summary_matrix, args.top_n)
                try:
                    write_heatmap(
                        plotted,
                        args.output_dir / f"summary_flux_heatmap.{args.figure_format}",
                        title=(
                            "Boundary-metabolite fluxes\n"
                            "positive = uptake; negative = secretion"
                        ),
                        transform=args.heatmap_transform,
                        dpi=args.dpi,
                    )
                except ValueError as error:
                    parser.error(str(error))

    reaction_outputs = {
        "reaction-long",
        "reaction-matrix",
        "reaction-heatmap",
    }
    if outputs.intersection(reaction_outputs):
        try:
            fluxes = read_many(files["fluxes"], "reaction flux")
        except ValueError as error:
            parser.error(str(error))
        fluxes = normalize_flux(fluxes, args.min_abs_flux)
        if optimal_result_ids is not None:
            if "result_id" not in fluxes.columns:
                parser.error("Flux inputs do not contain a result_id column.")
            fluxes = fluxes.loc[
                fluxes["result_id"].astype(str).isin(optimal_result_ids)
            ]

        if "reaction-long" in outputs:
            fluxes.to_csv(
                args.output_dir / "reaction_flux_long.tsv",
                sep="\t",
                index=False,
            )

        if outputs.intersection({"reaction-matrix", "reaction-heatmap"}):
            try:
                reaction_matrix = make_matrix(
                    fluxes,
                    row_columns=args.row_columns,
                    column="reaction",
                    separator=args.row_separator,
                )
            except KeyError as error:
                parser.error(str(error))

            if "reaction-matrix" in outputs:
                reaction_matrix.to_csv(
                    args.output_dir / "reaction_flux_matrix.tsv",
                    sep="\t",
                )

            if "reaction-heatmap" in outputs:
                plotted = select_heatmap_columns(reaction_matrix, args.top_n)
                try:
                    write_heatmap(
                        plotted,
                        args.output_dir / f"reaction_flux_heatmap.{args.figure_format}",
                        title="Reaction fluxes",
                        transform=args.heatmap_transform,
                        dpi=args.dpi,
                    )
                except ValueError as error:
                    parser.error(str(error))

    print("Created outputs:")
    for path in sorted(args.output_dir.iterdir()):
        if path.is_file():
            print(path)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        raise SystemExit(130)
