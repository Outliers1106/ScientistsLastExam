"""Alignment-blind negative controls for the withdrawn instance family."""


def _fold(items):
    tree = items[0]
    for item in items[1:]:
        tree = f"({tree},{item})"
    return tree


def _blocks(taxa):
    return _fold([_fold(taxa[i:i+5]) for i in range(0, len(taxa), 5)]) + ";"


def build_tree(problem):
    return _blocks(sorted(problem["taxa"], key=lambda name: int(name[1:])))


def row_blocks(problem):
    return _blocks(problem["taxa"])
