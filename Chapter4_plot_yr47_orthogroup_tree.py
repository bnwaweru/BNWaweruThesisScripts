#!/usr/bin/env python3
"""Plot the orthogroup-annotated Yr47 focal NB-ARC phylogeny.

This standalone script reproduces the final A4 rectangular tree used for the
Yr47 comparative-genomics analysis. It combines three independent annotations:

* orthogroup identity controls branch and tip-marker colour;
* NLRTracker class is represented by a black symbol after each gene name;
* right-side brackets identify orthogroups that agree with phylogenetic clades.

The script also supports a biologically useful intermediate case: an
orthogroup may be non-monophyletic overall but partition cleanly into a small
number of internally pure subclades. ``OG0000019`` is handled this way in the
current dataset and is displayed as two brackets containing 5/12 and 7/12 of
its members. The bracket report records this decision explicitly.

Input formats
-------------
Tree
    Either an IQ-TREE ``.iqtree`` report containing ``Consensus tree in newick
    format:`` or a plain Newick file such as ``.treefile``.
Alignment
    Optional FASTA alignment used to verify that alignment identifiers and tree
    tips match exactly.
NLRTracker classifications
    Tab-separated table containing the columns ``cultivar``, ``seqname`` and
    ``Simple``.
Orthogroups
    Tab-separated table containing ``gene`` and ``orthogroup``. A ``genome``
    column may also be present. Empty orthogroup values are treated as
    unassigned.
Orthogroup colours
    Tab-separated table containing ``orthogroup`` and ``colour``. Colours must
    be valid Matplotlib colours, preferably hexadecimal values.

Example
-------
python plot_yr47_orthogroup_tree.py \
    --tree morex+17pGsNLR_NBARCSeqs_mafft.clipkit.iqtree \
    --alignment morex+17pGsNLR_NBARCSeqs_mafft.fa \
    --classifications all17Genomes_NBARCsWithinLOC1Intervals.tsv \
    --orthogroups yr47_orthogroup_assignments.tsv \
    --orthogroup-colours yr47_orthogroup_colours.tsv \
    --split-orthogroup OG0000019 \
    --output-prefix results/yr47_focal_orthogroup_tree

Outputs
-------
``<output-prefix>.pdf``
    Vector PDF sized to one A4 page.
``<output-prefix>.svg``
    Editable vector graphic.
``<output-prefix>.png``
    High-resolution raster image.
``<output-prefix>.orthogroup_clade_report.tsv``
    Machine-readable record of the monophyly and bracket tests.

Dependencies
------------
Python >= 3.10, Biopython 1.88 and Matplotlib 3.10.8.
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import matplotlib

# A non-interactive backend is required on most HPC compute nodes.
matplotlib.use("Agg")
import matplotlib.colors as mpl_colours
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, DrawingArea, HPacker, TextArea
from Bio import Phylo


A4_INCHES = {
    "portrait": (8.27, 11.69),
    "landscape": (11.69, 8.27),
}

# NLR type is encoded only by shape; every symbol is black in the figure.
NLR_TYPE_MARKERS = {
    "ONL": "o",  # circle
    "N": "s",  # square
    "NL": "^",  # triangle
    "ON": "D",  # diamond
    "NLO": "*",  # star
    "NO": "P",  # filled plus
}

FALLBACK_NLR_MARKER = "X"
DEFAULT_UNASSIGNED_COLOUR = "#595959"


def parse_arguments() -> argparse.Namespace:
    """Parse command-line options."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate an A4 NB-ARC phylogeny coloured by orthogroup, with "
            "black NLRTracker-class symbols after gene names."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--tree",
        required=True,
        type=Path,
        help="IQ-TREE .iqtree report or plain Newick/.treefile input.",
    )
    parser.add_argument(
        "--alignment",
        type=Path,
        help="Optional FASTA alignment for exact identifier validation.",
    )
    parser.add_argument(
        "--classifications",
        required=True,
        type=Path,
        help="NLRTracker TSV containing cultivar, seqname and Simple.",
    )
    parser.add_argument(
        "--orthogroups",
        required=True,
        type=Path,
        help="TSV containing gene and orthogroup assignments.",
    )
    parser.add_argument(
        "--orthogroup-colours",
        required=True,
        type=Path,
        help="TSV containing orthogroup and colour columns.",
    )
    parser.add_argument(
        "--output-prefix",
        required=True,
        type=Path,
        help="Path prefix for output files; do not include an extension.",
    )
    parser.add_argument(
        "--split-orthogroup",
        action="append",
        default=[],
        help=(
            "Non-monophyletic orthogroup whose maximal pure subclades should "
            "receive separate brackets. Repeat for more than one group."
        ),
    )
    parser.add_argument(
        "--rooting",
        choices=("midpoint", "as-is"),
        default="midpoint",
        help="Rooting used for rectangular display only.",
    )
    parser.add_argument(
        "--orientation",
        choices=("portrait", "landscape"),
        default="portrait",
        help="A4 page orientation.",
    )
    parser.add_argument(
        "--focal-genome",
        default="WATDE0042",
        help="Genome highlighted using bold labels and diamond tip markers.",
    )
    parser.add_argument(
        "--support-threshold",
        type=float,
        default=70.0,
        help="Minimum internal branch-support value displayed.",
    )
    parser.add_argument(
        "--font-size",
        type=float,
        default=4.95,
        help="Tip-label size in points.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="Resolution of PNG output.",
    )
    parser.add_argument(
        "--scale-bar",
        type=float,
        default=0.5,
        help="Scale-bar length in substitutions per site; set to 0 to omit.",
    )
    parser.add_argument(
        "--label-fields",
        type=int,
        choices=(1, 2, 3),
        default=2,
        help="Number of pipe-delimited fields retained in displayed tip labels.",
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("pdf", "svg", "png"),
        default=("pdf", "svg", "png"),
        help="Figure formats to generate.",
    )
    parser.add_argument(
        "--allow-tip-mismatch",
        action="store_true",
        help="Continue when alignment and tree identifiers differ.",
    )
    parser.add_argument(
        "--title",
        help="Optional figure title. Omit to maximise space on the A4 page.",
    )
    return parser.parse_args()


def check_file(path: Path | None, description: str) -> None:
    """Fail early if a supplied input path cannot be read."""
    if path is not None and not path.is_file():
        raise FileNotFoundError(f"{description} does not exist: {path}")


def require_columns(
    fieldnames: list[str] | None,
    required: set[str],
    path: Path,
) -> None:
    """Verify that a TSV header contains all required columns."""
    missing = sorted(required - set(fieldnames or []))
    if missing:
        raise ValueError(
            f"{path} is missing required column(s): {', '.join(missing)}"
        )


def extract_newick(path: Path) -> str:
    """Extract Newick from an IQ-TREE report or return a plain Newick tree."""
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    marker = "Consensus tree in newick format:"

    for index, line in enumerate(lines):
        if line.strip() == marker:
            try:
                return next(value.strip() for value in lines[index + 1 :] if value.strip())
            except StopIteration as error:
                raise ValueError(
                    f"Found '{marker}' in {path}, but no Newick tree followed it."
                ) from error

    newick = text.strip()
    if newick.startswith("(") and newick.endswith(";"):
        return newick
    raise ValueError(f"Could not recognise {path} as IQ-TREE or Newick input.")


def read_tree(path: Path, rooting: str):
    """Parse, validate, optionally midpoint-root and ladderise the tree."""
    tree = Phylo.read(io.StringIO(extract_newick(path)), "newick")
    tip_names = [tip.name for tip in tree.get_terminals()]
    if any(name is None for name in tip_names):
        raise ValueError("Every terminal node must have a sequence identifier.")
    duplicates = sorted(name for name, count in Counter(tip_names).items() if count > 1)
    if duplicates:
        raise ValueError("Duplicate tree-tip identifiers: " + ", ".join(duplicates))
    if rooting == "midpoint":
        tree.root_at_midpoint()
    # Ladderising changes display order but not topology or branch length.
    tree.ladderize(reverse=False)
    return tree


def read_fasta_ids(path: Path) -> list[str]:
    """Read the first whitespace-delimited identifier from every FASTA header."""
    identifiers: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                identifier = line[1:].strip().split(maxsplit=1)[0]
                if not identifier:
                    raise ValueError(f"Empty FASTA header in {path}")
                identifiers.append(identifier)
    if not identifiers:
        raise ValueError(f"No FASTA records found in {path}")
    return identifiers


def validate_alignment(path: Path, tree, allow_mismatch: bool) -> None:
    """Confirm that the alignment and final tree contain the same identifiers."""
    alignment_ids = read_fasta_ids(path)
    tree_ids = [tip.name for tip in tree.get_terminals()]
    duplicate_alignment_ids = sorted(
        name for name, count in Counter(alignment_ids).items() if count > 1
    )
    if duplicate_alignment_ids:
        raise ValueError(
            "Duplicate alignment identifiers: " + ", ".join(duplicate_alignment_ids)
        )

    missing_from_tree = sorted(set(alignment_ids) - set(tree_ids))
    missing_from_alignment = sorted(set(tree_ids) - set(alignment_ids))
    print(
        f"Alignment/tree validation: {len(alignment_ids)} alignment sequences; "
        f"{len(tree_ids)} tree tips."
    )
    if not missing_from_tree and not missing_from_alignment:
        print("Alignment/tree validation: all identifiers match exactly.")
        return

    parts = ["Alignment and tree identifiers differ."]
    if missing_from_tree:
        parts.append("Missing from tree: " + ", ".join(missing_from_tree))
    if missing_from_alignment:
        parts.append("Missing from alignment: " + ", ".join(missing_from_alignment))
    message = " ".join(parts)
    if allow_mismatch:
        print("WARNING: " + message, file=sys.stderr)
    else:
        raise ValueError(message)


def normalise_gene_id(identifier: str) -> str:
    """Normalise naming differences specific to the Yr47 input datasets.

    The transformations reconcile:

    * the Renan assembly prefix used in the orthogroup table;
    * shortened WATDE0042/WATDE0100 names used in tree tips;
    * terminal ``LC`` annotations and transcript suffixes such as ``.1``;
    * punctuation and case differences.

    Original files are never modified. If naming conventions change, edit this
    function and document the change in version control.
    """
    identifier = identifier.strip().upper()
    identifier = re.sub(r"^RENAN_V2(?:X1|\.1|_1)?_", "", identifier)
    identifier = identifier.replace("TRAESWAT42", "TRAESWATDE0042")
    identifier = identifier.replace("TRAESWAT100", "TRAESWATDE0100")
    identifier = re.sub(r"LC(?=\.\d+$|$)", "", identifier)
    identifier = re.sub(r"\.\d+$", "", identifier)
    return re.sub(r"[^A-Z0-9]", "", identifier)


def add_unique_mapping(
    mapping: dict[str, str],
    key: str,
    value: str,
    source: Path,
) -> None:
    """Add an annotation while rejecting contradictory duplicate rows."""
    if key in mapping and mapping[key] != value:
        raise ValueError(
            f"Conflicting annotations for normalised gene '{key}' in {source}: "
            f"'{mapping[key]}' versus '{value}'."
        )
    mapping[key] = value


def read_classifications(path: Path) -> dict[str, str]:
    """Read NLRTracker Simple classes keyed by normalised gene identifier."""
    classifications: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader.fieldnames, {"cultivar", "seqname", "Simple"}, path)
        for line_number, row in enumerate(reader, start=2):
            full_identifier = (row["seqname"] or "").strip()
            cultivar = (row["cultivar"] or "").strip()
            nlr_class = (row["Simple"] or "").strip()
            if not full_identifier or not nlr_class:
                print(
                    f"WARNING: skipped incomplete classification row {line_number} "
                    f"in {path}",
                    file=sys.stderr,
                )
                continue

            # NLRTracker prepends the cultivar name to its identifiers.
            cultivar_prefix = cultivar.upper().replace(".", "X") + "_"
            if full_identifier.upper().startswith(cultivar_prefix):
                gene_identifier = full_identifier[len(cultivar_prefix) :]
            else:
                gene_identifier = full_identifier
            add_unique_mapping(
                classifications,
                normalise_gene_id(gene_identifier),
                nlr_class,
                path,
            )
    return classifications


def read_orthogroups(path: Path) -> dict[str, str]:
    """Read non-empty orthogroup assignments keyed by normalised gene ID."""
    orthogroups: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader.fieldnames, {"gene", "orthogroup"}, path)
        for line_number, row in enumerate(reader, start=2):
            gene_identifier = (row["gene"] or "").strip()
            orthogroup = (row["orthogroup"] or "").strip()
            if not gene_identifier:
                print(
                    f"WARNING: skipped orthogroup row {line_number} without a "
                    f"gene identifier in {path}",
                    file=sys.stderr,
                )
                continue
            if not orthogroup:
                continue
            add_unique_mapping(
                orthogroups,
                normalise_gene_id(gene_identifier),
                orthogroup,
                path,
            )
    return orthogroups


def read_orthogroup_colours(path: Path) -> dict[str, str]:
    """Read and validate orthogroup colours."""
    colours: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader.fieldnames, {"orthogroup", "colour"}, path)
        for line_number, row in enumerate(reader, start=2):
            orthogroup = (row["orthogroup"] or "").strip()
            colour = (row["colour"] or "").strip()
            if not orthogroup or not colour:
                raise ValueError(f"Incomplete colour row {line_number} in {path}")
            if not mpl_colours.is_color_like(colour):
                raise ValueError(
                    f"Invalid Matplotlib colour '{colour}' on row {line_number} "
                    f"of {path}"
                )
            if orthogroup in colours and colours[orthogroup].lower() != colour.lower():
                raise ValueError(f"Conflicting colours for {orthogroup} in {path}")
            colours[orthogroup] = colour
    colours.setdefault("OG-unassigned", DEFAULT_UNASSIGNED_COLOUR)
    return colours


def tip_fields(tip_name: str) -> list[str]:
    """Split a tip expected to have Genome|Gene|Annotation fields."""
    return tip_name.split("|")


def tip_gene(tip_name: str) -> str:
    """Extract gene ID from a tree tip."""
    fields = tip_fields(tip_name)
    if len(fields) < 2:
        raise ValueError(
            f"Tree tip '{tip_name}' lacks the expected 'Genome|Gene' structure."
        )
    return fields[1]


def tip_genome(tip_name: str) -> str:
    """Extract genome ID from a tree tip."""
    return tip_fields(tip_name)[0]


def display_label(tip_name: str, number_of_fields: int) -> str:
    """Create a readable label while retaining deterministic identifiers."""
    return " | ".join(tip_fields(tip_name)[:number_of_fields])


def tip_class(tip_name: str, classifications: dict[str, str]) -> str:
    """Return the NLRTracker class or the no-symbol fallback category."""
    return classifications.get(normalise_gene_id(tip_gene(tip_name)), "Not classified")


def tip_orthogroup(tip_name: str, orthogroups: dict[str, str]) -> str | None:
    """Return a tip's orthogroup assignment when available."""
    return orthogroups.get(normalise_gene_id(tip_gene(tip_name)))


def plotted_group(tip_name: str, orthogroups: dict[str, str]) -> str:
    """Return an orthogroup or the explicit grey unassigned category."""
    return tip_orthogroup(tip_name, orthogroups) or "OG-unassigned"


def homogeneous_clade_group(clade, orthogroups: dict[str, str]) -> str | None:
    """Return an OG only when every descendant is assigned to the same OG."""
    groups = {plotted_group(tip.name, orthogroups) for tip in clade.get_terminals()}
    if len(groups) == 1 and "OG-unassigned" not in groups:
        return next(iter(groups))
    return None


def calculate_y_positions(tree) -> tuple[dict[Any, float], dict[Any, float]]:
    """Assign evenly spaced tip coordinates and centred internal coordinates."""
    terminals = tree.get_terminals()
    tip_y = {tip: float(index) for index, tip in enumerate(reversed(terminals))}
    all_y: dict[Any, float] = dict(tip_y)

    def assign(clade) -> float:
        if clade in all_y:
            return all_y[clade]
        child_positions = [assign(child) for child in clade.clades]
        all_y[clade] = (min(child_positions) + max(child_positions)) / 2.0
        return all_y[clade]

    assign(tree.root)
    return tip_y, all_y


def calculate_x_positions(tree) -> tuple[dict[Any, float], float, bool]:
    """Use cumulative branch lengths, with unit lengths only as a fallback."""
    depths = tree.depths(unit_branch_lengths=False)
    maximum = max(depths[tip] for tip in tree.get_terminals())
    used_unit_lengths = False
    if maximum <= 0:
        depths = tree.depths(unit_branch_lengths=True)
        maximum = max(depths[tip] for tip in tree.get_terminals())
        used_unit_lengths = True
        print(
            "WARNING: branch lengths were absent; plotting a unit-length cladogram.",
            file=sys.stderr,
        )
    return depths, maximum, used_unit_lengths


def analyse_global_orthogroups(tree, orthogroups: dict[str, str]) -> tuple[list[dict], list[dict]]:
    """Test each orthogroup for exact or qualified monophyly.

    Exact monophyly means the MRCA contains all and only members of the target
    orthogroup. Qualified monophyly permits additional tips only when those tips
    lack an orthogroup assignment; such brackets are dashed.
    """
    tips_by_group: dict[str, list[Any]] = defaultdict(list)
    for tip in tree.get_terminals():
        group = tip_orthogroup(tip.name, orthogroups)
        if group:
            tips_by_group[group].append(tip)

    report_rows: list[dict] = []
    brackets: list[dict] = []
    for group in sorted(tips_by_group):
        group_tips = tips_by_group[group]
        mrca = tree.common_ancestor(group_tips)
        descendants = set(mrca.get_terminals())
        group_set = set(group_tips)
        descendant_groups = {
            tip_orthogroup(tip.name, orthogroups) for tip in descendants
        }
        conflicts = sorted(
            value for value in descendant_groups if value not in (None, group)
        )
        unassigned = sum(
            tip_orthogroup(tip.name, orthogroups) is None for tip in descendants
        )

        if len(group_tips) == 1:
            status = "singleton_in_tree"
            draw = False
        elif descendants == group_set:
            status = "monophyletic_exact"
            draw = True
        elif not conflicts:
            status = "monophyletic_with_og_unassigned_tips"
            draw = True
        else:
            status = "not_monophyletic"
            draw = False

        report_rows.append(
            {
                "orthogroup": group,
                "members_in_tree": len(group_tips),
                "tips_within_global_mrca": len(descendants),
                "og_unassigned_tips_within_global_mrca": unassigned,
                "conflicting_orthogroups_within_global_mrca": ",".join(conflicts),
                "global_status": status,
                "global_bracket_drawn": "yes" if draw else "no",
                "split_subclades_drawn": "0",
                "split_subclade_sizes": "",
                "split_subclade_unassigned_counts": "",
            }
        )
        if draw:
            brackets.append(
                {
                    "group": group,
                    "clade": mrca,
                    "members": len(group_tips),
                    "total_members": len(group_tips),
                    "unassigned": unassigned,
                    "kind": "global",
                    "label": group,
                }
            )
    return report_rows, brackets


def find_split_subclades(
    tree,
    orthogroups: dict[str, str],
    target_group: str,
) -> list[dict]:
    """Find maximal pure subclades that collectively cover a split OG.

    Each accepted subclade must contain at least two target genes and no gene
    assigned to a conflicting orthogroup. Tips without an orthogroup assignment
    are allowed but recorded and shown with a dashed bracket.
    """
    parent_by_child = {
        child: parent
        for parent in tree.find_clades(order="level")
        for child in parent.clades
    }

    def is_pure_target(clade) -> bool:
        groups = [tip_orthogroup(tip.name, orthogroups) for tip in clade.get_terminals()]
        assigned_groups = {group for group in groups if group is not None}
        return groups.count(target_group) >= 2 and assigned_groups == {target_group}

    total_members = sum(
        tip_orthogroup(tip.name, orthogroups) == target_group
        for tip in tree.get_terminals()
    )
    if total_members < 2:
        raise ValueError(
            f"Cannot annotate split orthogroup {target_group}: only "
            f"{total_members} member(s) occur in the tree."
        )

    subclades: list[dict] = []
    for clade in tree.get_nonterminals(order="preorder"):
        if not is_pure_target(clade):
            continue
        parent = parent_by_child.get(clade)
        if parent is not None and is_pure_target(parent):
            continue
        members = sum(
            tip_orthogroup(tip.name, orthogroups) == target_group
            for tip in clade.get_terminals()
        )
        unassigned = sum(
            tip_orthogroup(tip.name, orthogroups) is None
            for tip in clade.get_terminals()
        )
        subclades.append(
            {
                "group": target_group,
                "clade": clade,
                "members": members,
                "total_members": total_members,
                "unassigned": unassigned,
                "kind": "split",
                "label": f"{target_group} ({members}/{total_members})",
            }
        )

    covered_members = sum(item["members"] for item in subclades)
    if len(subclades) < 2 or covered_members != total_members:
        raise ValueError(
            f"Pure subclades for {target_group} do not form a complete split: "
            f"{len(subclades)} subclade(s) cover {covered_members}/{total_members} "
            "members. Review the tree or omit --split-orthogroup for this group."
        )
    return subclades


def add_split_results_to_report(
    report_rows: list[dict],
    target_group: str,
    subclades: list[dict],
) -> None:
    """Record split-subclade decisions in the matching orthogroup report row."""
    row = next(
        (entry for entry in report_rows if entry["orthogroup"] == target_group),
        None,
    )
    if row is None:
        raise ValueError(f"No report row exists for split group {target_group}")
    row["split_subclades_drawn"] = str(len(subclades))
    row["split_subclade_sizes"] = ",".join(
        str(item["members"]) for item in subclades
    )
    row["split_subclade_unassigned_counts"] = ",".join(
        str(item["unassigned"]) for item in subclades
    )


def write_clade_report(output_prefix: Path, rows: list[dict]) -> Path:
    """Write the auditable orthogroup/clade decision table."""
    output = Path(f"{output_prefix}.orthogroup_clade_report.tsv")
    fieldnames = [
        "orthogroup",
        "members_in_tree",
        "tips_within_global_mrca",
        "og_unassigned_tips_within_global_mrca",
        "conflicting_orthogroups_within_global_mrca",
        "global_status",
        "global_bracket_drawn",
        "split_subclades_drawn",
        "split_subclade_sizes",
        "split_subclade_unassigned_counts",
    ]
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return output


def make_tip_label_box(
    label: str,
    font_size: float,
    bold: bool,
    nlr_class: str,
):
    """Pack a black gene label and its black NLRTracker symbol side by side."""
    label_area = TextArea(
        label,
        textprops={
            "color": "#222222",
            "fontsize": font_size,
            "fontweight": "bold" if bold else "normal",
        },
    )
    children = [label_area]

    if nlr_class != "Not classified":
        marker = NLR_TYPE_MARKERS.get(nlr_class, FALLBACK_NLR_MARKER)
        symbol_area = DrawingArea(7.0, 7.0, 0, 0)
        symbol = Line2D(
            [3.5],
            [3.5],
            marker=marker,
            markersize=4.5 if marker != "*" else 5.3,
            markerfacecolor="#000000",
            markeredgecolor="#000000",
            markeredgewidth=0.45,
            linestyle="none",
        )
        symbol_area.add_artist(symbol)
        children.append(symbol_area)
    return HPacker(children=children, align="center", pad=0, sep=2.0)


def plot_tree(
    tree,
    classifications: dict[str, str],
    orthogroups: dict[str, str],
    colours: dict[str, str],
    brackets: list[dict],
    args: argparse.Namespace,
):
    """Construct the rectangular publication figure."""
    tip_y, all_y = calculate_y_positions(tree)
    depths, max_depth, used_unit_lengths = calculate_x_positions(tree)
    terminals = tree.get_terminals()

    observed_groups = [plotted_group(tip.name, orthogroups) for tip in terminals]
    group_counts = Counter(observed_groups)
    nlr_counts = Counter(tip_class(tip.name, classifications) for tip in terminals)
    missing_colours = sorted(set(observed_groups) - set(colours))
    if missing_colours:
        raise ValueError(
            "No colour was supplied for observed group(s): "
            + ", ".join(missing_colours)
        )

    label_x = max_depth * 1.045
    bracket_x = max_depth * 2.08
    figure, axis = plt.subplots(figsize=A4_INCHES[args.orientation])
    # Reserve space for the orthogroup and NLR-symbol legends.
    figure.subplots_adjust(left=0.035, right=0.995, top=0.865, bottom=0.035)
    neutral_branch = "#8A8A8A"

    def draw_clade(parent, clade) -> None:
        """Draw one branch and recurse into its descendants."""
        parent_x = depths[parent]
        child_x = depths[clade]
        child_y = all_y[clade]
        if clade.is_terminal():
            group = plotted_group(clade.name, orthogroups)
        else:
            group = homogeneous_clade_group(clade, orthogroups)
        colour = colours[group] if group else neutral_branch
        width = 0.90 if clade.is_terminal() else 0.60

        axis.plot(
            [parent_x, child_x],
            [child_y, child_y],
            color=colour,
            linewidth=width,
            solid_capstyle="butt",
            zorder=2,
        )
        if not clade.is_terminal():
            child_positions = [all_y[child] for child in clade.clades]
            axis.plot(
                [child_x, child_x],
                [min(child_positions), max(child_positions)],
                color=colour,
                linewidth=0.60,
                solid_capstyle="butt",
                zorder=1,
            )
            if (
                clade.confidence is not None
                and float(clade.confidence) >= args.support_threshold
            ):
                axis.text(
                    child_x,
                    child_y + 0.30,
                    f"{float(clade.confidence):g}",
                    fontsize=max(args.font_size - 1.40, 3.0),
                    color="#4D4D4D",
                    ha="center",
                    va="bottom",
                    clip_on=False,
                )
            for child in clade.clades:
                draw_clade(clade, child)

    if tree.root.clades:
        root_y = [all_y[child] for child in tree.root.clades]
        root_group = homogeneous_clade_group(tree.root, orthogroups)
        root_colour = colours[root_group] if root_group else neutral_branch
        axis.plot(
            [0, 0],
            [min(root_y), max(root_y)],
            color=root_colour,
            linewidth=0.60,
            solid_capstyle="butt",
        )
        for child in tree.root.clades:
            draw_clade(tree.root, child)

    for tip in terminals:
        y_value = tip_y[tip]
        tip_x = depths[tip]
        group = plotted_group(tip.name, orthogroups)
        colour = colours[group]
        nlr_class = tip_class(tip.name, classifications)
        is_focal = tip_genome(tip.name).upper() == args.focal_genome.upper()

        axis.plot(
            [tip_x, label_x - max_depth * 0.006],
            [y_value, y_value],
            color=colour,
            linewidth=0.34,
            alpha=0.70,
            zorder=0,
        )
        # A fine dark edge keeps very light supplied colours visible on white.
        marker_edge = "#222222" if mpl_colours.rgb_to_hsv(mpl_colours.to_rgb(colour))[2] > 0.85 else "none"
        axis.scatter(
            [label_x],
            [y_value],
            s=4.2 if is_focal else 2.7,
            color=colour,
            marker="D" if is_focal else "o",
            edgecolors=marker_edge,
            linewidths=0.25 if marker_edge != "none" else 0,
            zorder=3,
        )

        packed_label = make_tip_label_box(
            display_label(tip.name, args.label_fields),
            font_size=args.font_size,
            bold=is_focal,
            nlr_class=nlr_class,
        )
        annotation = AnnotationBbox(
            packed_label,
            (label_x, y_value),
            xybox=(3.0, 0),
            xycoords="data",
            boxcoords="offset points",
            box_alignment=(0.0, 0.5),
            frameon=False,
            pad=0,
            zorder=4,
        )
        axis.add_artist(annotation)

    # Brackets are disjoint by construction for the present dataset. Sorting by
    # vertical position makes output deterministic across runs.
    brackets.sort(
        key=lambda item: min(tip_y[tip] for tip in item["clade"].get_terminals()),
        reverse=True,
    )
    bracket_width = max_depth * 0.030
    for item in brackets:
        y_values = [tip_y[tip] for tip in item["clade"].get_terminals()]
        y_min = min(y_values) - 0.30
        y_max = max(y_values) + 0.30
        colour = colours[item["group"]]
        line_style = "--" if item["unassigned"] else "-"
        axis.plot(
            [bracket_x, bracket_x],
            [y_min, y_max],
            color=colour,
            linewidth=1.0,
            linestyle=line_style,
            solid_capstyle="butt",
        )
        for y_end in (y_min, y_max):
            axis.plot(
                [bracket_x - bracket_width, bracket_x],
                [y_end, y_end],
                color=colour,
                linewidth=1.0,
                linestyle=line_style,
                solid_capstyle="butt",
            )
        label = item["label"] + ("*" if item["unassigned"] else "")
        axis.text(
            bracket_x + max_depth * 0.025,
            (y_min + y_max) / 2,
            label,
            fontsize=args.font_size + 0.10,
            fontweight="bold",
            color="#222222",
            ha="left",
            va="center",
        )
    axis.text(
        bracket_x,
        len(terminals) - 0.05,
        "Orthogroup-supported clades",
        fontsize=args.font_size + 0.20,
        fontweight="bold",
        color="#252525",
        ha="center",
        va="bottom",
    )

    # Only orthogroups represented in this tree are included in the legend.
    colour_order = [group for group in colours if group in group_counts]
    group_handles = [
        Line2D(
            [0],
            [0],
            color=colours[group],
            marker="s",
            markeredgecolor="#222222",
            markeredgewidth=0.25,
            markersize=3.8,
            linewidth=1.2,
            label=f"{group} (n={group_counts[group]})",
        )
        for group in colour_order
    ]
    figure.legend(
        handles=group_handles,
        title="Orthogroup colour",
        loc="upper left",
        bbox_to_anchor=(0.035, 0.995),
        ncol=5,
        frameon=False,
        fontsize=4.6,
        title_fontsize=5.1,
        handlelength=1.35,
        handletextpad=0.35,
        columnspacing=0.9,
        borderaxespad=0,
    )

    classified_types = sorted(
        category for category in nlr_counts if category != "Not classified"
    )
    preferred_type_order = ["ONL", "N", "NL", "ON", "NLO", "NO"]
    type_order = [value for value in preferred_type_order if value in classified_types]
    type_order.extend(value for value in classified_types if value not in type_order)
    nlr_handles = [
        Line2D(
            [0],
            [0],
            color="#000000",
            marker=NLR_TYPE_MARKERS.get(nlr_class, FALLBACK_NLR_MARKER),
            markerfacecolor="#000000",
            markeredgecolor="#000000",
            markersize=(
                4.8
                if NLR_TYPE_MARKERS.get(nlr_class, FALLBACK_NLR_MARKER) == "*"
                else 4.0
            ),
            linestyle="none",
            label=f"{nlr_class} (n={nlr_counts[nlr_class]})",
        )
        for nlr_class in type_order
    ]
    nlr_handles.append(
        Line2D(
            [0],
            [0],
            color="none",
            marker=None,
            linestyle="none",
            label=(
                "No symbol = not classified "
                f"(n={nlr_counts.get('Not classified', 0)})"
            ),
        )
    )
    figure.legend(
        handles=nlr_handles,
        title="NLRTracker class: black symbol after gene name",
        loc="upper left",
        bbox_to_anchor=(0.035, 0.905),
        ncol=6,
        frameon=False,
        fontsize=4.7,
        title_fontsize=5.1,
        handlelength=1.0,
        handletextpad=0.30,
        columnspacing=0.9,
        borderaxespad=0,
    )

    scale_y = -1.65
    if args.scale_bar > 0:
        scale = args.scale_bar
        axis.plot([0, scale], [scale_y, scale_y], color="#333333", linewidth=0.8)
        axis.plot(
            [0, 0],
            [scale_y - 0.15, scale_y + 0.15],
            color="#333333",
            linewidth=0.8,
        )
        axis.plot(
            [scale, scale],
            [scale_y - 0.15, scale_y + 0.15],
            color="#333333",
            linewidth=0.8,
        )
        scale_text = (
            f"{scale:g} branch-length units"
            if used_unit_lengths
            else f"{scale:g} substitutions per site"
        )
        axis.text(
            scale / 2,
            scale_y - 0.42,
            scale_text,
            fontsize=max(args.font_size - 0.15, 4.0),
            color="#333333",
            ha="center",
            va="top",
        )

    axis.text(
        label_x,
        scale_y,
        f"Bold diamonds: {args.focal_genome}",
        fontsize=max(args.font_size - 0.25, 4.0),
        color="#333333",
        ha="left",
        va="center",
    )

    split_notes = []
    for target_group in args.split_orthogroup:
        target_brackets = [
            item for item in brackets if item["group"] == target_group and item["kind"] == "split"
        ]
        sizes = ", ".join(
            f'{item["members"]}/{item["total_members"]}' for item in target_brackets
        )
        split_notes.append(f"{target_group}: {len(target_brackets)} pure clades ({sizes})")
    footer_note = "; ".join(split_notes)
    if any(item["unassigned"] for item in brackets):
        footer_note += ("; " if footer_note else "") + "dashed*: includes OG-unassigned genes"
    if footer_note:
        axis.text(
            bracket_x - max_depth * 0.16,
            scale_y,
            footer_note,
            fontsize=max(args.font_size - 0.90, 3.8),
            color="#333333",
            ha="left",
            va="center",
        )

    if args.title:
        figure.suptitle(args.title, fontsize=8.0, fontweight="normal", y=0.999)

    axis.set_xlim(-max_depth * 0.01, max_depth * 2.75)
    axis.set_ylim(-2.5, len(terminals) - 0.15)
    axis.axis("off")
    return figure, group_counts, nlr_counts


def save_figure(
    figure,
    output_prefix: Path,
    formats: list[str],
    dpi: int,
) -> list[Path]:
    """Save requested outputs using metadata accepted by each backend."""
    title = "Orthogroup-coloured phylogeny of Yr47-interval NB-ARC sequences"
    outputs: list[Path] = []
    for file_format in formats:
        output = Path(f"{output_prefix}.{file_format}")
        if file_format == "pdf":
            figure.savefig(
                output,
                metadata={
                    "Title": title,
                    "Subject": "Midpoint-rooted IQ-TREE consensus phylogeny",
                    "Keywords": "Yr47, NB-ARC, NLR, orthogroup, phylogeny",
                },
            )
        elif file_format == "svg":
            figure.savefig(
                output,
                metadata={
                    "Title": title,
                    "Description": "Rectangular IQ-TREE consensus phylogeny",
                    "Keywords": "Yr47, NB-ARC, NLR, orthogroup, phylogeny",
                },
            )
        elif file_format == "png":
            figure.savefig(
                output,
                dpi=dpi,
                facecolor="white",
                metadata={"Title": title},
            )
        outputs.append(output)
    return outputs


def main() -> None:
    """Validate inputs, annotate the tree, render outputs and report decisions."""
    args = parse_arguments()
    for path, description in (
        (args.tree, "Tree input"),
        (args.alignment, "Alignment input"),
        (args.classifications, "Classification input"),
        (args.orthogroups, "Orthogroup input"),
        (args.orthogroup_colours, "Orthogroup-colour input"),
    ):
        check_file(path, description)
    if args.font_size <= 0:
        raise ValueError("--font-size must be greater than zero")
    if args.dpi <= 0:
        raise ValueError("--dpi must be greater than zero")

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    tree = read_tree(args.tree, args.rooting)
    if args.alignment:
        validate_alignment(args.alignment, tree, args.allow_tip_mismatch)

    classifications = read_classifications(args.classifications)
    orthogroups = read_orthogroups(args.orthogroups)
    colours = read_orthogroup_colours(args.orthogroup_colours)

    report_rows, brackets = analyse_global_orthogroups(tree, orthogroups)
    globally_bracketed = {item["group"] for item in brackets}
    for target_group in args.split_orthogroup:
        if target_group in globally_bracketed:
            raise ValueError(
                f"{target_group} already forms a globally bracketed clade and "
                "should not also be requested with --split-orthogroup."
            )
        if target_group not in colours:
            raise ValueError(f"No colour is defined for {target_group}")
        subclades = find_split_subclades(tree, orthogroups, target_group)
        brackets.extend(subclades)
        add_split_results_to_report(report_rows, target_group, subclades)

    report_path = write_clade_report(args.output_prefix, report_rows)
    figure, group_counts, nlr_counts = plot_tree(
        tree,
        classifications,
        orthogroups,
        colours,
        brackets,
        args,
    )
    output_paths = save_figure(figure, args.output_prefix, args.formats, args.dpi)
    plt.close(figure)

    print(f"Tree tips plotted: {len(tree.get_terminals())}")
    print(
        "Orthogroup counts: "
        + ", ".join(f"{key}={value}" for key, value in sorted(group_counts.items()))
    )
    print(
        "NLRTracker classes: "
        + ", ".join(f"{key}={value}" for key, value in sorted(nlr_counts.items()))
    )
    for output in output_paths:
        print(f"Wrote: {output}")
    print(f"Wrote: {report_path}")


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2) from error
