#!/usr/bin/env python3
"""Annotate COBRA flux outputs with ModelSEED metadata and draw heatmaps."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ModelSEEDDatabase file paths relative to the root
COMPOUNDS_FILE = Path("Biochemistry/compounds.tsv")
COMPOUND_ALIASES_FILE = Path(
    "Biochemistry/Aliases/Unique_ModelSEED_Compound_Aliases.txt"
)
REACTIONS_FILE = Path("Biochemistry/reactions.tsv")
REACTION_ALIASES_FILE = Path(
    "Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt"
)
REACTION_ECS_FILE = Path(
    "Biochemistry/Aliases/Unique_ModelSEED_Reaction_ECs.txt"
)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Combine run_cobra_model.py outputs, annotate compounds and "
            "reactions with ModelSEED metadata, and draw flux heatmaps."
        )
    )
    parser.add_argument(
        "--results",
        type=Path,
        nargs="+",
        required=True,
        help="Files or directories containing *.summary.tsv and *.fluxes.tsv.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for annotated tables and heatmaps.",
    )
    parser.add_argument(
        "--modelseed-db",
        type=Path,
        required=True,
        help="Root directory of a ModelSEEDDatabase checkout.",
    )
    return parser


def required_modelseed_file(root: Path, relative_path: Path) -> Path:
    """Return a required file inside the ModelSEEDDatabase checkout."""
    path = root / relative_path
    if not path.is_file():
        raise FileNotFoundError(f"Required ModelSEED file not found: {path}")
    return path


def read_tables(paths: list[Path], suffix: str) -> pd.DataFrame:
    """Read and combine COBRA output tables matching a filename suffix."""
    files: list[Path] = []

    for path in paths:
        if not path.exists():
            raise FileNotFoundError(f"Result path does not exist: {path}")

        if path.is_dir():
            files.extend(path.rglob(f"*{suffix}"))
        elif path.name.endswith(suffix):
            files.append(path)

    files = sorted(set(files))

    if not files:
        raise ValueError(f"No files ending in {suffix} were found.")

    return pd.concat(
        (pd.read_csv(path, sep="\t") for path in files),
        ignore_index=True,
        sort=False,
    )


def join_unique(values: pd.Series) -> str:
    """Join unique non-empty strings."""
    cleaned = {
        str(value).strip()
        for value in values
        if str(value).strip()
        and str(value).strip().lower() not in {"nan", "null", "none"}
    }
    return " | ".join(sorted(cleaned))


def read_aliases(path: Path) -> pd.DataFrame:
    """Read a ModelSEED alias table."""
    aliases = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

    required = {"ModelSEED ID", "External ID", "Source"}
    missing = required.difference(aliases.columns)

    if missing:
        raise ValueError(
            f"{path} is missing required columns: {', '.join(sorted(missing))}"
        )

    return aliases


def alias_metadata(aliases: pd.DataFrame) -> pd.DataFrame:
    """Collapse external aliases and BiGG IDs by ModelSEED ID."""
    working = aliases.copy()
    working["alias"] = (
        working["Source"].astype(str)
        + ":"
        + working["External ID"].astype(str)
    )

    all_aliases = (
        working.groupby("ModelSEED ID")["alias"]
        .agg(join_unique)
        .rename("modelseed_aliases")
        .reset_index()
    )

    bigg = working[
        working["Source"].str.lower().str.startswith("bigg")
    ]

    bigg_ids = (
        bigg.groupby("ModelSEED ID")["External ID"]
        .agg(join_unique)
        .rename("bigg_ids")
        .reset_index()
    )

    return all_aliases.merge(
        bigg_ids,
        on="ModelSEED ID",
        how="left",
    )


def build_lookup(
    modelseed_ids: pd.Series,
    aliases: pd.DataFrame,
) -> tuple[pd.DataFrame, set[str]]:
    """Build lookup keys from ModelSEED IDs and BiGG aliases."""

    direct = pd.DataFrame(
        {
            "lookup_id": modelseed_ids.astype(str),
            "modelseed_id": modelseed_ids.astype(str),
            "mapping_source": "ModelSEED",
        }
    )

    bigg = aliases[
        aliases["Source"].str.lower().str.startswith("bigg")
    ][
        ["ModelSEED ID", "External ID"]
    ].rename(
        columns={
            "ModelSEED ID": "modelseed_id",
            "External ID": "lookup_id",
        }
    )

    bigg["mapping_source"] = "BiGG"

    lookup = pd.concat(
        [direct, bigg],
        ignore_index=True,
    ).drop_duplicates()

    # Find lookup IDs that refer to more than
    # one ModelSEED feature.
    counts = (
        lookup.groupby("lookup_id")["modelseed_id"]
        .nunique()
    )

    ambiguous = set(
        counts[counts > 1].index.astype(str)
    )

    # Do not automatically resolve genuinely ambiguous mappings.
    lookup = lookup[
        ~lookup["lookup_id"].isin(ambiguous)
    ].copy()

    # Multiple records that point to the SAME ModelSEED ID are
    # not ambiguous. Preserve all sources rather than choosing one.
    lookup = (
        lookup
        .groupby(
            ["lookup_id", "modelseed_id"],
            as_index=False,
        )["mapping_source"]
        .agg(join_unique)
    )

    return lookup, ambiguous


def load_compounds(
    modelseed_db: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    """Load ModelSEED compound metadata and lookup aliases."""
    compounds = pd.read_csv(
        required_modelseed_file(modelseed_db, COMPOUNDS_FILE),
        sep="\t",
        dtype=str,
        keep_default_na=False,
        usecols=["id", "name", "abbreviation", "formula", "charge"],
    ).rename(
        columns={
            "id": "modelseed_id",
            "name": "modelseed_name",
            "abbreviation": "modelseed_abbreviation",
        }
    )

    aliases = read_aliases(
        required_modelseed_file(modelseed_db, COMPOUND_ALIASES_FILE)
    )

    metadata = compounds.merge(
        alias_metadata(aliases).rename(columns={"ModelSEED ID": "modelseed_id"}),
        on="modelseed_id",
        how="left",
    )

    lookup, ambiguous = build_lookup(compounds["modelseed_id"], aliases)

    return metadata, lookup, ambiguous


def load_reactions(
    modelseed_db: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    """Load ModelSEED reaction metadata, aliases, and EC numbers."""
    reactions = pd.read_csv(
        required_modelseed_file(modelseed_db, REACTIONS_FILE),
        sep="\t",
        dtype=str,
        keep_default_na=False,
        usecols=["id", "name", "abbreviation", "definition"],
    ).rename(
        columns={
            "id": "modelseed_id",
            "name": "modelseed_name",
            "abbreviation": "modelseed_abbreviation",
        }
    )

    aliases = read_aliases(
        required_modelseed_file(modelseed_db, REACTION_ALIASES_FILE)
    )

    ecs = read_aliases(
        required_modelseed_file(modelseed_db, REACTION_ECS_FILE)
    )

    ec_numbers = (
        ecs.groupby("ModelSEED ID")["External ID"]
        .agg(join_unique)
        .rename("ec_numbers")
        .reset_index()
        .rename(columns={"ModelSEED ID": "modelseed_id"})
    )

    metadata = (
        reactions.merge(
            alias_metadata(aliases).rename(
                columns={"ModelSEED ID": "modelseed_id"}
            ),
            on="modelseed_id",
            how="left",
        )
        .merge(
            ec_numbers,
            on="modelseed_id",
            how="left",
        )
    )

    lookup, ambiguous = build_lookup(reactions["modelseed_id"], aliases)

    return metadata, lookup, ambiguous


def compound_lookup_id(identifier: str) -> str:
    """Remove model prefixes and compartment suffixes from metabolite IDs."""
    identifier = re.sub(r"^M_", "", str(identifier))

    if identifier.startswith("cpd"):
        return identifier.split("_", 1)[0]

    return re.sub(r"_[a-z]\d*$", "", identifier)


def compound_compartment(identifier: str) -> str:
    """Return a normalized compartment suffix from a metabolite ID."""
    identifier = re.sub(r"^M_", "", str(identifier))
    match = re.search(r"_([a-z])\d*$", identifier)
    return match.group(1) if match else ""


def reaction_lookup_id(identifier: str) -> str:
    """Remove model prefixes and gapseq compartment suffixes from reaction IDs."""
    identifier = re.sub(r"^R_", "", str(identifier))
    match = re.match(r"^(rxn\d{5})(?:_[a-z]\d*)?$", identifier)
    return match.group(1) if match else identifier


def reaction_compartment(identifier: str) -> str:
    """Return a normalized compartment suffix from a ModelSEED reaction ID."""
    identifier = re.sub(r"^R_", "", str(identifier))
    match = re.match(r"^rxn\d{5}_([a-z])\d*$", identifier)
    return match.group(1) if match else ""


def annotate(
    frame: pd.DataFrame,
    native_column: str,
    metadata: pd.DataFrame,
    lookup: pd.DataFrame,
    ambiguous: set[str],
    lookup_function,
    compartment_function,
) -> pd.DataFrame:
    """Join ModelSEED metadata to a COBRA output table."""
    annotated = frame.copy()
    annotated["lookup_id"] = (
        annotated[native_column].astype(str).map(lookup_function)
    )

    annotated = annotated.merge(
        lookup,
        on="lookup_id",
        how="left",
        validate="many_to_one",
    )

    annotated["mapping_source"] = annotated["mapping_source"].fillna("unmapped")
    annotated.loc[
        annotated["lookup_id"].isin(ambiguous),
        "mapping_source",
    ] = "ambiguous"

    annotated = annotated.merge(
        metadata,
        on="modelseed_id",
        how="left",
        validate="many_to_one",
    )

    compartments = (
        annotated[native_column]
        .astype(str)
        .map(compartment_function)
    )

    annotated["feature_id"] = annotated[native_column].astype(str)
    mapped = annotated["modelseed_id"].notna()

    annotated.loc[mapped, "feature_id"] = (
        annotated.loc[mapped, "modelseed_id"].astype(str)
        + compartments[mapped].map(lambda value: f"_{value}" if value else "")
    )

    return annotated.drop(columns="lookup_id")


def model_labels(frame: pd.DataFrame) -> pd.Series:
    """Return labels identifying each COBRA model/run."""
    if {"sample_id", "medium"}.issubset(frame.columns):
        return (
            frame["sample_id"].astype(str)
            + " | "
            + frame["medium"].astype(str)
        )

    if "model_name" in frame.columns:
        return frame["model_name"].astype(str)

    raise ValueError(
        "COBRA tables must contain sample_id + medium or model_name."
    )


def compound_display_labels(frame: pd.DataFrame) -> pd.Series:
    """Create readable compound labels for the heatmap."""
    labels = frame["metabolite"].astype(str).copy()
    mapped = frame["modelseed_id"].notna() & frame["modelseed_name"].notna()

    labels.loc[mapped] = (
        frame.loc[mapped, "modelseed_name"].astype(str)
        + " ["
        + frame.loc[mapped, "feature_id"].astype(str)
        + "]"
    )

    return labels


def reaction_display_labels(frame: pd.DataFrame) -> pd.Series:
    """Create readable reaction labels using ModelSEED names and EC numbers."""
    labels = frame["reaction"].astype(str).copy()
    mapped = frame["modelseed_id"].notna() & frame["modelseed_name"].notna()

    labels.loc[mapped] = frame.loc[mapped, "modelseed_name"].astype(str)

    has_ec = mapped & frame["ec_numbers"].fillna("").astype(str).ne("")
    labels.loc[has_ec] = (
        frame.loc[has_ec, "modelseed_name"].astype(str)
        + " [EC "
        + frame.loc[has_ec, "ec_numbers"].astype(str)
        + "]"
    )

    no_ec = mapped & ~has_ec
    labels.loc[no_ec] = (
        frame.loc[no_ec, "modelseed_name"].astype(str)
        + " ["
        + frame.loc[no_ec, "feature_id"].astype(str)
        + "]"
    )

    return labels


def make_matrix(
    frame: pd.DataFrame,
    display_labels: pd.Series,
) -> tuple[pd.DataFrame, list[str]]:
    """Reshape raw flux values for plotting without aggregation."""
    working = frame.copy()
    working["model_label"] = model_labels(working)
    working["display_label"] = display_labels
    working["flux"] = pd.to_numeric(working["flux"], errors="raise")

    duplicates = working.duplicated(
        subset=["model_label", "feature_id"],
        keep=False,
    )

    if duplicates.any():
        example = (
            working.loc[
                duplicates,
                ["model_label", "feature_id"],
            ]
            .drop_duplicates()
            .head(10)
        )
        raise ValueError(
            "Multiple rows map to the same ModelSEED feature within a model. "
            "Creating a heatmap would require aggregating fluxes, which this "
            "script intentionally does not do. Examples:\n"
            + example.to_string(index=False)
        )

    matrix = working.pivot(
        index="model_label",
        columns="feature_id",
        values="flux",
    )

    label_map = (
        working[["feature_id", "display_label"]]
        .drop_duplicates(subset="feature_id")
        .set_index("feature_id")["display_label"]
        .to_dict()
    )

    display_order = [
        label_map.get(feature_id, feature_id)
        for feature_id in matrix.columns
    ]

    return matrix, display_order


def write_heatmap(
    matrix: pd.DataFrame,
    labels: list[str],
    output_path: Path,
    title: str,
) -> None:
    """Draw a heatmap of raw flux values."""
    values = matrix.to_numpy(dtype=float)

    finite = np.abs(values[np.isfinite(values)])
    maximum = float(finite.max()) if finite.size else 1.0
    if maximum == 0:
        maximum = 1.0

    width = min(40.0, max(10.0, 0.25 * matrix.shape[1] + 4.0))
    height = min(30.0, max(4.0, 0.35 * matrix.shape[0] + 2.0))

    figure, axis = plt.subplots(figsize=(width, height))

    image = axis.imshow(
        values,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-maximum,
        vmax=maximum,
    )

    axis.set_xticks(np.arange(matrix.shape[1]))
    axis.set_xticklabels(labels, rotation=90, fontsize=7)

    axis.set_yticks(np.arange(matrix.shape[0]))
    axis.set_yticklabels(matrix.index, fontsize=8)

    axis.set_xlabel("")
    axis.set_ylabel("Model")
    axis.set_title(title)

    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label("Flux")

    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def main() -> int:
    """Annotate COBRA outputs and draw compound and reaction heatmaps."""
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    compound_metadata, compound_lookup, ambiguous_compounds = load_compounds(
        args.modelseed_db
    )
    reaction_metadata, reaction_lookup, ambiguous_reactions = load_reactions(
        args.modelseed_db
    )

    compound_fluxes = read_tables(args.results, ".summary.tsv")
    compound_fluxes = annotate(
        compound_fluxes,
        native_column="metabolite",
        metadata=compound_metadata,
        lookup=compound_lookup,
        ambiguous=ambiguous_compounds,
        lookup_function=compound_lookup_id,
        compartment_function=compound_compartment,
    )

    reaction_fluxes = read_tables(args.results, ".fluxes.tsv")
    reaction_fluxes = annotate(
        reaction_fluxes,
        native_column="reaction",
        metadata=reaction_metadata,
        lookup=reaction_lookup,
        ambiguous=ambiguous_reactions,
        lookup_function=reaction_lookup_id,
        compartment_function=reaction_compartment,
    )

    compound_fluxes.to_csv(
        args.output_dir / "compound_fluxes.tsv",
        sep="\t",
        index=False,
    )

    reaction_fluxes.to_csv(
        args.output_dir / "reaction_fluxes.tsv",
        sep="\t",
        index=False,
    )

    compound_matrix, compound_labels = make_matrix(
        compound_fluxes,
        compound_display_labels(compound_fluxes),
    )
    write_heatmap(
        compound_matrix,
        compound_labels,
        args.output_dir / "compound_flux_heatmap.png",
        "Compound fluxes",
    )

    reaction_matrix, reaction_labels = make_matrix(
        reaction_fluxes,
        reaction_display_labels(reaction_fluxes),
    )
    write_heatmap(
        reaction_matrix,
        reaction_labels,
        args.output_dir / "reaction_flux_heatmap.png",
        "Reaction fluxes",
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
