#!/usr/bin/env python3

"""Generate an Escher HTML visualization for a COBRA SBML model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import cobra
from escher import Builder
import libsbml


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate an Escher visualization from an SBML model."
    )

    parser.add_argument(
        "--model",
        required=True,
        type=Path,
        help="Input SBML model.",
    )

    parser.add_argument(
        "--map",
        required=True,
        type=Path,
        help="Escher map JSON used as the common layout.",
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output Escher HTML file.",
    )

    parser.add_argument(
        "--reaction-scale-preset",
        default="GaBuRd",
        help="Escher reaction color-scale preset.",
    )

    return parser.parse_args()


def get_map_reactions(map_file: Path) -> set[str]:
    """Return reaction IDs contained in an Escher map."""

    with map_file.open() as handle:
        escher_map = json.load(handle)

    reactions = escher_map[1]["reactions"]

    return {
        reaction["bigg_id"]
        for reaction in reactions.values()
        if "bigg_id" in reaction
    }


def main():
    args = parse_args()

    # Load COBRA model
    model = cobra.io.read_sbml_model(str(args.model))

    print(f"Model: {args.model}")
    print(f"Reactions: {len(model.reactions)}")
    print(f"Metabolites: {len(model.metabolites)}")
    print(f"Genes: {len(model.genes)}")

    # Check how much of the reference Escher map matches the model.
    map_reactions = get_map_reactions(args.map)
    model_reactions = {reaction.id for reaction in model.reactions}

    matching_reactions = map_reactions & model_reactions

    print(f"Escher map reactions: {len(map_reactions)}")
    print(f"Matching model/map reactions: {len(matching_reactions)}")
    print(
        f"Map coverage: "
        f"{len(matching_reactions) / len(map_reactions) * 100:.1f}%"
    )

    # Optimize model using its existing objective.
    solution = model.optimize()

    print(f"Optimization status: {solution.status}")

    if solution.status != "optimal":
        raise RuntimeError(
            f"Model optimization failed with status: {solution.status}"
        )

    print(f"Objective value: {solution.objective_value}")

    # Build the visualization.
    builder = Builder(
        map_json=str(args.map),
        model=model,
        reaction_data=solution.fluxes.to_dict(),
        metabolite_data=solution.shadow_prices.to_dict(),
        reaction_scale_preset=args.reaction_scale_preset,

        # Interactive information
        enable_tooltips=["label", "object"],
        show_gene_reaction_rules=True,
        show_metabolite_names=True,
    )



    # Write outputs (html and smbl layout)
    builder.save_html(str(args.output))

    print(f"Escher map written to: {args.output}")


if __name__ == "__main__":
    main()