#!/usr/bin/env python3
"""Optimize one COBRA model and write standardized tabular outputs.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd
from cobra.flux_analysis import pfba
from cobra.io import read_sbml_model
from cobra.util.solver import linear_reaction_coefficients


MODEL_SUFFIXES = (".xml.gz", ".sbml.gz", ".xml", ".sbml")
METADATA_COLUMNS = {
    "model_name",
    "result_id",
    "output_prefix",
    "source_model",
    "cobra_model_id",
    "status",
    "optimization_method",
    "fraction_of_optimum",
    "maximum_biomass",
    "solution_biomass",
    "objective_direction",
    "objective_reactions",
    "objective_expression",
    "solver",
    "number_reactions",
    "number_metabolites",
    "number_genes",
    "number_boundary_reactions",
    "number_exchange_reactions",
    "medium_file",
    "medium_mode",
    "missing_medium_reactions",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Load one SBML model, optimize its biomass objective, and write a "
            "model record, all reaction fluxes, and model.summary() boundary fluxes."
        )
    )
    parser.add_argument("model", type=Path, help="Input CarveMe SBML/XML model.")
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for TSV outputs.",
    )
    parser.add_argument(
        "--name",
        help=(
            "Logical model name stored in output tables. Default: the model filename "
            "without its extension (.sbml, .xml, .sbml.gz, .xml.gz, use --prefix to "
            "override)."
        ),
    )
    parser.add_argument(
        "--prefix",
        help=(
            "Filename prefix for outputs. Defaults to a filesystem-safe form of --name."
        ),
    )
    parser.add_argument(
        "--label",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help=(
            "Extra metadata copied into every output row. Repeat for sample, medium, "
            "treatment, etc. Example: --label sample_id=S1 --label medium=M9"
        ),
    )
    parser.add_argument(
        "--method",
        choices=("fba", "pfba"),
        default="pfba",
        help="Flux solution to report after maximizing biomass. Default: pfba.",
    )
    parser.add_argument(
        "--fraction-of-optimum",
        type=float,
        default=1.0,
        help="Objective fraction retained by pFBA. Default: 1.0.",
    )
    parser.add_argument(
        "--objective",
        help=(
            "Reaction ID to use as the biomass objective. If omitted, retain the "
            "objective encoded in the model."
        ),
    )
    parser.add_argument(
        "--solver",
        help="Optional COBRApy solver name, such as glpk, cplex, or gurobi.",
    )
    parser.add_argument(
        "--medium",
        type=Path,
        help=(
            "Optional CSV/TSV medium file with columns named 'exchange' and 'uptake'. "
            "Uptake values must be non-negative COBRApy medium bounds."
        ),
    )
    parser.add_argument(
        "--medium-mode",
        choices=("replace", "update"),
        default="replace",
        help=(
            "replace closes uptake for unlisted exchanges; update modifies only listed "
            "exchanges. Default: replace."
        ),
    )
    parser.add_argument(
        "--ignore-missing-medium-reactions",
        action="store_true",
        help="Ignore medium reactions absent from the model instead of failing.",
    )
    parser.add_argument(
        "--flux-threshold",
        type=float,
        default=1e-9,
        help="Fluxes below this absolute value are written as zero. Default: 1e-9.",
    )
    parser.add_argument(
        "--fail-on-nonoptimal",
        action="store_true",
        help=(
            "Exit non-zero when optimization is not optimal. Tables are still written "
            "before exiting."
        ),
    )
    return parser


def strip_model_suffix(path: Path) -> str:
    name = path.name
    for suffix in MODEL_SUFFIXES:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def safe_prefix(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
    if not cleaned:
        raise ValueError("The output prefix becomes empty after sanitization.")
    return cleaned


def parse_labels(values: Iterable[str]) -> dict[str, str]:
    labels: dict[str, str] = {}
    for item in values:
        if "=" not in item:
            raise ValueError(f"Invalid --label {item!r}; expected KEY=VALUE.")
        key, value = item.split("=", 1)
        key = key.strip()
        if not key:
            raise ValueError(f"Invalid --label {item!r}; key cannot be empty.")
        if key in METADATA_COLUMNS:
            raise ValueError(f"Label key {key!r} is reserved by the output schema.")
        if key in labels:
            raise ValueError(f"Duplicate --label key: {key!r}.")
        labels[key] = value
    return labels


def read_medium(path: Path) -> dict[str, float]:
    frame = pd.read_csv(path, sep=None, engine="python")
    required = {"exchange", "uptake"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(
            f"Medium file {path} is missing columns: {', '.join(sorted(missing))}."
        )

    if frame["exchange"].duplicated().any():
        duplicates = sorted(frame.loc[frame["exchange"].duplicated(), "exchange"].unique())
        raise ValueError(f"Duplicate exchanges in medium file: {duplicates}")

    uptake = pd.to_numeric(frame["uptake"], errors="raise")
    if (uptake < 0).any():
        bad = frame.loc[uptake < 0, "exchange"].tolist()
        raise ValueError(
            "COBRApy model.medium values must be non-negative maximum uptake rates; "
            f"negative values found for: {bad}"
        )

    return dict(zip(frame["exchange"].astype(str), uptake.astype(float), strict=True))


def apply_medium(
    model,
    medium: dict[str, float],
    mode: str,
    ignore_missing: bool,
) -> list[str]:
    reaction_ids = {reaction.id for reaction in model.reactions}
    missing = sorted(set(medium).difference(reaction_ids))
    if missing and not ignore_missing:
        raise KeyError(
            "Medium reactions absent from model: " + ", ".join(missing)
        )

    filtered = {key: value for key, value in medium.items() if key in reaction_ids}
    if mode == "replace":
        model.medium = filtered
    else:
        updated = dict(model.medium)
        updated.update(filtered)
        model.medium = updated
    return missing


def objective_details(model) -> tuple[dict[str, float], str]:
    coefficients = linear_reaction_coefficients(model)
    coefficient_map = {
        reaction.id: float(coefficient)
        for reaction, coefficient in coefficients.items()
    }
    return coefficient_map, str(model.objective.expression)


def objective_value_from_fluxes(
    fluxes: pd.Series,
    coefficients: dict[str, float],
) -> float | None:
    if not coefficients:
        return None
    return float(sum(fluxes[reaction_id] * coef for reaction_id, coef in coefficients.items()))


def zero_small_values(series: pd.Series, threshold: float) -> pd.Series:
    result = series.astype(float).copy()
    result.loc[result.abs() < threshold] = 0.0
    return result


def add_identity_columns(
    frame: pd.DataFrame,
    model_name: str,
    result_id: str,
    labels: dict[str, str],
) -> pd.DataFrame:
    result = frame.copy()
    result.insert(0, "model_name", model_name)
    result.insert(1, "result_id", result_id)
    for position, (key, value) in enumerate(labels.items(), start=2):
        result.insert(position, key, value)
    return result


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, sep="\t", index=False)
    temporary.replace(path)


def build_flux_table(model, solution, threshold: float) -> pd.DataFrame:
    records = []
    for reaction in model.reactions:
        flux = float(solution.fluxes[reaction.id])
        if abs(flux) < threshold:
            flux = 0.0
        records.append(
            {
                "reaction": reaction.id,
                "reaction_name": reaction.name or "",
                "flux": flux,
                "absolute_flux": abs(flux),
                "lower_bound": float(reaction.lower_bound),
                "upper_bound": float(reaction.upper_bound),
                "equation": reaction.reaction,
                "gene_reaction_rule": reaction.gene_reaction_rule,
            }
        )
    return pd.DataFrame.from_records(records)


def build_summary_table(model, solution, threshold: float) -> pd.DataFrame:
    summary = model.summary(solution=solution).to_frame().reset_index(drop=True)
    summary["flux"] = zero_small_values(summary["flux"], threshold)
    summary["absolute_flux"] = summary["flux"].abs()
    summary["direction"] = "inactive"
    summary.loc[summary["flux"] > 0, "direction"] = "uptake"
    summary.loc[summary["flux"] < 0, "direction"] = "secretion"

    exchange_ids = {reaction.id for reaction in model.exchanges}
    demand_ids = {reaction.id for reaction in model.demands}
    sink_ids = {reaction.id for reaction in model.sinks}

    def boundary_type(reaction_id: str) -> str:
        if reaction_id in exchange_ids:
            return "exchange"
        if reaction_id in demand_ids:
            return "demand"
        if reaction_id in sink_ids:
            return "sink"
        return "boundary"

    summary["boundary_type"] = summary["reaction"].map(boundary_type)
    return summary[
        [
            "reaction",
            "metabolite",
            "boundary_type",
            "factor",
            "flux",
            "absolute_flux",
            "direction",
        ]
    ]


def empty_flux_table() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "reaction",
            "reaction_name",
            "flux",
            "absolute_flux",
            "lower_bound",
            "upper_bound",
            "equation",
            "gene_reaction_rule",
        ]
    )


def empty_summary_table() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "reaction",
            "metabolite",
            "boundary_type",
            "factor",
            "flux",
            "absolute_flux",
            "direction",
        ]
    )


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if not 0 < args.fraction_of_optimum <= 1:
        parser.error("--fraction-of-optimum must be greater than 0 and at most 1.")
    if args.flux_threshold < 0:
        parser.error("--flux-threshold cannot be negative.")
    if not args.model.is_file():
        parser.error(f"Model file does not exist: {args.model}")
    if args.medium is not None and not args.medium.is_file():
        parser.error(f"Medium file does not exist: {args.medium}")

    try:
        labels = parse_labels(args.label)
    except ValueError as error:
        parser.error(str(error))

    model_name = args.name or strip_model_suffix(args.model)
    prefix = args.prefix or safe_prefix(model_name)
    if Path(prefix).name != prefix:
        parser.error("--prefix must be a filename prefix, not a path.")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    model_path = args.output_dir / f"{prefix}.model.tsv"
    flux_path = args.output_dir / f"{prefix}.fluxes.tsv"
    summary_path = args.output_dir / f"{prefix}.summary.tsv"

    model = read_sbml_model(str(args.model))

    if args.solver:
        model.solver = args.solver

    if args.objective:
        try:
            model.objective = model.reactions.get_by_id(args.objective)
        except KeyError as error:
            parser.error(f"Objective reaction not found: {args.objective}")
        model.objective_direction = "max"

    missing_medium_reactions: list[str] = []
    if args.medium:
        try:
            medium = read_medium(args.medium)
            missing_medium_reactions = apply_medium(
                model=model,
                medium=medium,
                mode=args.medium_mode,
                ignore_missing=args.ignore_missing_medium_reactions,
            )
        except (ValueError, KeyError) as error:
            parser.error(str(error))

    objective_reactions, objective_expression = objective_details(model)
    growth_solution = model.optimize()
    status = str(growth_solution.status)

    maximum_biomass: float | None = None
    solution_biomass: float | None = None
    flux_table = empty_flux_table()
    summary_table = empty_summary_table()

    if status == "optimal":
        maximum_biomass = float(growth_solution.objective_value)
        if args.method == "pfba":
            flux_solution = pfba(
                model,
                fraction_of_optimum=args.fraction_of_optimum,
            )
        else:
            flux_solution = growth_solution

        solution_biomass = objective_value_from_fluxes(
            flux_solution.fluxes,
            objective_reactions,
        )
        flux_table = build_flux_table(model, flux_solution, args.flux_threshold)
        summary_table = build_summary_table(model, flux_solution, args.flux_threshold)

    metadata = {
        "model_name": model_name,
        "result_id": prefix,
        **labels,
        "output_prefix": prefix,
        "source_model": str(args.model.resolve()),
        "cobra_model_id": model.id or "",
        "status": status,
        "optimization_method": args.method,
        "fraction_of_optimum": args.fraction_of_optimum,
        "maximum_biomass": maximum_biomass,
        "solution_biomass": solution_biomass,
        "objective_direction": model.objective_direction,
        "objective_reactions": json.dumps(objective_reactions, sort_keys=True),
        "objective_expression": objective_expression,
        "solver": model.solver.interface.__name__,
        "number_reactions": len(model.reactions),
        "number_metabolites": len(model.metabolites),
        "number_genes": len(model.genes),
        "number_boundary_reactions": len(model.boundary),
        "number_exchange_reactions": len(model.exchanges),
        "medium_file": str(args.medium.resolve()) if args.medium else "",
        "medium_mode": args.medium_mode if args.medium else "",
        "missing_medium_reactions": json.dumps(missing_medium_reactions),
    }

    metadata_table = pd.DataFrame([metadata])
    flux_table = add_identity_columns(flux_table, model_name, prefix, labels)
    summary_table = add_identity_columns(summary_table, model_name, prefix, labels)

    write_tsv(metadata_table, model_path)
    write_tsv(flux_table, flux_path)
    write_tsv(summary_table, summary_path)

    print(f"status={status}")
    print(f"maximum_biomass={maximum_biomass}")
    print(f"model_table={model_path}")
    print(f"flux_table={flux_path}")
    print(f"summary_table={summary_path}")

    if status != "optimal" and args.fail_on_nonoptimal:
        return 2
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted.", file=sys.stderr)
        raise SystemExit(130)
