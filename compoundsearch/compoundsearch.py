import os
from collections import defaultdict

from pymatgen import Composition

__author__ = 'Anubhav Jain <ajain@lbl.gov>'

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_LIT_FILE = os.path.join(module_dir, "compounds_list", "thermoelectrics.txt")

# TODO: make a proper command line interface


def get_group_formula(composition):
    """
    Returns the formula string, but with elements replaced by their group num
    in the Periodic Table.

    e.g., Y2SiO5 becomes (3)2-(14)1-(16)5

    Args:
        composition: (Composition)

    Returns: (str)

    """
    new_comp_dict = defaultdict(float)

    for e in composition.elements:
        new_comp_dict[e.group] += composition[e]

    new_comp = Composition(new_comp_dict).reduced_composition

    form = ""
    sorted_elements = sorted(new_comp.elements, key=lambda x: x.Z)

    for x in sorted_elements:
        amt = new_comp[x] if new_comp[x] != int(new_comp[x]) else int(new_comp[x])
        l = x.Z
        form += "({}){}-".format(l, amt)

    return form[0:-1]


def get_fuzzy_formula(composition):
    """
    Returns the formula string, but with elements anonymized by atom type

    e.g., Y2SiO5 becomes Tm1Y3

    Args:
        composition: (Composition)

    Returns: (str)

    """
    x = defaultdict(float)

    for e in composition.elements:
        if e.is_lanthanoid or e.is_actinoid:
            x["Re"] += composition[e]
        elif e.is_alkali or e.is_alkaline:
            x["A"] += composition[e]
        elif e.is_transition_metal:
            x["Tm"] += composition[e]
        elif e.Z in [13, 31, 49, 50, 81, 82, 83]:
            x["B"] += composition[e]
        elif e.is_metalloid or e.Z in [6, 7, 8, 15, 16, 34] or e.is_halogen:
            x["Y"] += composition[e]
        elif e.is_noble_gas:
            x["Nb"] += composition[e]
        else:
            x["U"] += composition[e]

    c = Composition(x).get_reduced_composition_and_factor()[0]

    my_form=""
    for label in ["A", "B", "Tm", "Re", "Y", "Nb", "U"]:
        if c[label]:
            amt = c[label] if c[label] != int(c[label]) else int(c[label])
            my_form += "{}{}".format(label, amt)

    return my_form


def get_fuzzy_formula_strict(composition):
    """
    Returns the formula string, but with elements anonymized by atom type using
    a greater number of anonymizations than the normal fuzzy_formula method.

    e.g., Y2SiO5 becomes Tm2X1Y5

    Args:
        composition: (Composition)

    Returns: (str)

    """

    x = defaultdict(float)

    for e in composition.elements:
        if e.is_lanthanoid:
            x["Ln"] += composition[e]
        elif e.is_actinoid:
            x["Ac"] += composition[e]
        elif e.is_alkali:
            x["A"] += composition[e]
        elif e.is_alkaline:
            x["B"] += composition[e]
        elif e.is_transition_metal:
            x["Tm"] += composition[e]
        elif e.is_metalloid:
            x["X"] += composition[e]
        elif e.Z in [13, 31, 49, 50, 81, 82, 83]:
            x["C"] += composition[e]
        elif e.Z in [6, 7, 8, 15, 16, 34]:
            x["Y"] += composition[e]
        elif e.is_halogen:
            x["Z"] += composition[e]
        elif e.is_noble_gas:
            x["Nb"] += composition[e]
        else:
            x["U"] += composition[e]

    c = Composition(x).get_reduced_composition_and_factor()[0]

    my_form=""
    for label in ["A", "B", "C", "Tm", "Ln", "Ac", "X", "Y", "Z", "Nb", "U"]:
        if c[label]:
            amt = c[label] if c[label] != int(c[label]) else int(c[label])
            my_form += "{}{}".format(label, amt)

    return my_form


# TODO: be able to customize the scoring system
# TODO: tune the scoring system a bit better
def match_compositions(formula1, formula2):
    c1 = Composition(formula1).get_reduced_composition_and_factor()[0]
    c2 = Composition(formula2).get_reduced_composition_and_factor()[0]

    score = 0
    matches = []

    if c1 == c2:
        return 100, ["exact formula"]

    if c1.anonymized_formula == c2.anonymized_formula:
        score += 15
        matches.append("anonymized formula")

    if set(c1.elements) == set(c2.elements):
        score += 25
        matches.append("chemical system")

    if get_group_formula(c1) == get_group_formula(c2):
        score += 30
        matches.append("group formula")

    if get_fuzzy_formula_strict(c1) == get_fuzzy_formula_strict(c2):
        score += 20
        matches.append("fuzzy formula (strict)")

    elif get_fuzzy_formula(c1) == get_fuzzy_formula(c2):
        score += 15
        matches.append("fuzzy formula")

    return score, matches


class CompoundSearch:
    # TODO: be able to search more than one lit file
    def __init__(self, lit_file=DEFAULT_LIT_FILE):

        self.line_data = {}  #  (line no.) -> {"formula_pretty", "ref1", "ref2", "ref3"}

        if not os.path.exists(lit_file):
            raise ValueError("Cannot find lit file: {}".format(lit_file))

        with open(lit_file) as f:
            line_no = 1
            ref1 = ""  # most recent level-1 tag
            ref2 = ""  # most recent level-2 tag
            ref3 = ""  # most recent level-3 tag

            for line in f:
                line = line.strip()
                if line.startswith("###"):
                    ref3 = line[3:]
                elif line.startswith("##"):
                    ref2 = line[2:]
                    ref3 = ""
                elif line.startswith("#"):
                    ref1 = line[1:]
                    ref2 = ""
                    ref3 = ""
                elif line:
                    c = Composition(line)
                    self.line_data[line_no] = \
                        {"formula_pretty": c.get_reduced_formula_and_factor()[0],
                         "ref1": ref1, "ref2": ref2, "ref3": ref3}
                line_no += 1

    def search(self, target_formula):
        target_data = {}
        # add "score" and "matches" keys to target_data
        for line_no in self.line_data:
            score, matches = match_compositions(self.line_data[line_no]["formula_pretty"], target_formula)
            target_data[line_no] = {}
            target_data[line_no]["score"] = score
            target_data[line_no]["matches"] = matches

        # TODO: make this a Pandas dataframe so that you can sort however you want
        # print header
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format("formula_pretty", "score", "matches", "line_no", "ref1", "ref2", "ref3"))

        for line_no in self.line_data:
            ld = self.line_data[line_no]
            td = target_data[line_no]
            if td["score"] > 0:
                print "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
                    ld["formula_pretty"], td["score"], td["matches"], line_no,
                    ld["ref1"], ld["ref2"], ld["ref3"])

if __name__ == "__main__":
    cs = CompoundSearch()
    cs.search("InCl")