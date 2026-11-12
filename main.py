"""
Stoichiometry Assistant - CS1 Final Project
Author: Teo Valianatos

A Drafter web app that:
  * balances chemical equations automatically,
  * finds the limiting reactant, leftover reactant masses, and theoretical
    yield from the masses you start with,
  * decides whether the reaction is exothermic or endothermic and computes
    the heat flow per mole of each product,
  * calculates percent yield from the mass of product you actually harvested.

Based on OpenStax Chemistry 2e, Chapter 4 (Chemical Reactions and
Stoichiometry) and Chapter 5 (Thermochemistry).

All chemistry is written in plain Python (no ChemPy / PubChem) so that the
site still works after it is deployed to the browser.
"""
from dataclasses import dataclass
from drafter import *
from bakery import assert_equal


# ===========================================================================
# GLOBAL CONSTANTS (constants are allowed; there are no global variables)
# ===========================================================================

# (symbol, name, molar mass in g/mol) in order of atomic number.
# Radioactive elements use the mass number of their most stable isotope.
ELEMENTS = [
    ("H", "Hydrogen", 1.008), ("He", "Helium", 4.0026),
    ("Li", "Lithium", 6.94), ("Be", "Beryllium", 9.0122),
    ("B", "Boron", 10.81), ("C", "Carbon", 12.011),
    ("N", "Nitrogen", 14.007), ("O", "Oxygen", 15.999),
    ("F", "Fluorine", 18.998), ("Ne", "Neon", 20.180),
    ("Na", "Sodium", 22.990), ("Mg", "Magnesium", 24.305),
    ("Al", "Aluminum", 26.982), ("Si", "Silicon", 28.085),
    ("P", "Phosphorus", 30.974), ("S", "Sulfur", 32.06),
    ("Cl", "Chlorine", 35.45), ("Ar", "Argon", 39.948),
    ("K", "Potassium", 39.098), ("Ca", "Calcium", 40.078),
    ("Sc", "Scandium", 44.956), ("Ti", "Titanium", 47.867),
    ("V", "Vanadium", 50.942), ("Cr", "Chromium", 51.996),
    ("Mn", "Manganese", 54.938), ("Fe", "Iron", 55.845),
    ("Co", "Cobalt", 58.933), ("Ni", "Nickel", 58.693),
    ("Cu", "Copper", 63.546), ("Zn", "Zinc", 65.38),
    ("Ga", "Gallium", 69.723), ("Ge", "Germanium", 72.630),
    ("As", "Arsenic", 74.922), ("Se", "Selenium", 78.971),
    ("Br", "Bromine", 79.904), ("Kr", "Krypton", 83.798),
    ("Rb", "Rubidium", 85.468), ("Sr", "Strontium", 87.62),
    ("Y", "Yttrium", 88.906), ("Zr", "Zirconium", 91.224),
    ("Nb", "Niobium", 92.906), ("Mo", "Molybdenum", 95.95),
    ("Tc", "Technetium", 98.0), ("Ru", "Ruthenium", 101.07),
    ("Rh", "Rhodium", 102.91), ("Pd", "Palladium", 106.42),
    ("Ag", "Silver", 107.87), ("Cd", "Cadmium", 112.41),
    ("In", "Indium", 114.82), ("Sn", "Tin", 118.71),
    ("Sb", "Antimony", 121.76), ("Te", "Tellurium", 127.60),
    ("I", "Iodine", 126.90), ("Xe", "Xenon", 131.29),
    ("Cs", "Cesium", 132.91), ("Ba", "Barium", 137.33),
    ("La", "Lanthanum", 138.91), ("Ce", "Cerium", 140.12),
    ("Pr", "Praseodymium", 140.91), ("Nd", "Neodymium", 144.24),
    ("Pm", "Promethium", 145.0), ("Sm", "Samarium", 150.36),
    ("Eu", "Europium", 151.96), ("Gd", "Gadolinium", 157.25),
    ("Tb", "Terbium", 158.93), ("Dy", "Dysprosium", 162.50),
    ("Ho", "Holmium", 164.93), ("Er", "Erbium", 167.26),
    ("Tm", "Thulium", 168.93), ("Yb", "Ytterbium", 173.05),
    ("Lu", "Lutetium", 174.97), ("Hf", "Hafnium", 178.49),
    ("Ta", "Tantalum", 180.95), ("W", "Tungsten", 183.84),
    ("Re", "Rhenium", 186.21), ("Os", "Osmium", 190.23),
    ("Ir", "Iridium", 192.22), ("Pt", "Platinum", 195.08),
    ("Au", "Gold", 196.97), ("Hg", "Mercury", 200.59),
    ("Tl", "Thallium", 204.38), ("Pb", "Lead", 207.2),
    ("Bi", "Bismuth", 208.98), ("Po", "Polonium", 209.0),
    ("At", "Astatine", 210.0), ("Rn", "Radon", 222.0),
    ("Fr", "Francium", 223.0), ("Ra", "Radium", 226.0),
    ("Ac", "Actinium", 227.0), ("Th", "Thorium", 232.04),
    ("Pa", "Protactinium", 231.04), ("U", "Uranium", 238.03),
    ("Np", "Neptunium", 237.0), ("Pu", "Plutonium", 244.0),
    ("Am", "Americium", 243.0), ("Cm", "Curium", 247.0),
    ("Bk", "Berkelium", 247.0), ("Cf", "Californium", 251.0),
    ("Es", "Einsteinium", 252.0), ("Fm", "Fermium", 257.0),
    ("Md", "Mendelevium", 258.0), ("No", "Nobelium", 259.0),
    ("Lr", "Lawrencium", 266.0), ("Rf", "Rutherfordium", 267.0),
    ("Db", "Dubnium", 268.0), ("Sg", "Seaborgium", 269.0),
    ("Bh", "Bohrium", 270.0), ("Hs", "Hassium", 269.0),
    ("Mt", "Meitnerium", 278.0), ("Ds", "Darmstadtium", 281.0),
    ("Rg", "Roentgenium", 282.0), ("Cn", "Copernicium", 285.0),
    ("Nh", "Nihonium", 286.0), ("Fl", "Flerovium", 289.0),
    ("Mc", "Moscovium", 290.0), ("Lv", "Livermorium", 293.0),
    ("Ts", "Tennessine", 294.0), ("Og", "Oganesson", 294.0),
]

# (formula, common name) shown on the Common Molar Masses page.
COMMON_COMPOUNDS = [
    ("H2O", "Water"), ("CO2", "Carbon dioxide"), ("O2", "Oxygen gas"),
    ("H2", "Hydrogen gas"), ("N2", "Nitrogen gas"), ("NH3", "Ammonia"),
    ("CH4", "Methane"), ("C3H8", "Propane"), ("C8H18", "Octane"),
    ("C6H12O6", "Glucose"), ("C2H5OH", "Ethanol"), ("NaCl", "Table salt"),
    ("HCl", "Hydrochloric acid"), ("H2SO4", "Sulfuric acid"),
    ("NaOH", "Sodium hydroxide"), ("Ca(OH)2", "Calcium hydroxide"),
    ("CaCO3", "Calcium carbonate"), ("NaHCO3", "Baking soda"),
    ("Fe2O3", "Iron(III) oxide (rust)"), ("Al2O3", "Aluminum oxide"),
    ("KMnO4", "Potassium permanganate"), ("CuSO4", "Copper(II) sulfate"),
]

# Standard enthalpies of formation at 298 K in kJ/mol, from OpenStax
# Chemistry 2e Appendix G. Each formula is assumed to be in its usual
# state at room temperature (for example, H2O is liquid water).
# Elements in their standard state (O2, H2, Fe, ...) are 0 and are handled
# in get_enthalpy_of_formation, so they are not listed here.
ENTHALPIES_OF_FORMATION = [
    ("H2O", -285.8), ("CO2", -393.5), ("CO", -110.5), ("CH4", -74.6),
    ("C2H6", -84.0), ("C3H8", -103.8), ("C4H10", -125.6),
    ("C8H18", -250.1), ("C2H5OH", -277.6), ("CH3OH", -239.2),
    ("C6H12O6", -1273.3), ("NH3", -45.9), ("NO", 91.3), ("NO2", 33.2),
    ("N2O4", 11.1), ("SO2", -296.8), ("SO3", -395.7), ("H2S", -20.6),
    ("HCl", -92.3), ("NaCl", -411.2), ("CaCO3", -1206.9),
    ("CaO", -634.9), ("MgO", -601.6), ("Fe2O3", -824.2),
    ("Al2O3", -1675.7), ("H2O2", -187.8), ("O3", 142.7),
]

# ===========================================================================
# DATACLASSES
# ===========================================================================

@dataclass
class Species:
    """
    One chemical (element or compound) in a reaction.

    formula: chemical formula, such as "H2O" or "Ca(OH)2"
    coefficient: number in front of the formula in the balanced equation
    molar_mass: grams per mole (0.0 means the formula was not valid)
    mass: for a reactant, the grams the user started with;
          for a product, the grams the user actually harvested
          (0.0 means not entered yet)
    """
    formula: str
    coefficient: int
    molar_mass: float
    mass: float


@dataclass
class State:
    """
    Everything the site remembers between pages.

    equation: the current balanced equation, like "2 H2 + O2 → 2 H2O"
    reactants: the reactants of the equation on the calculation page
    products: the products of the equation on the calculation page
    total_reactant_mass: sum of all reactant masses entered (grams)
    enthalpy: heat of reaction in kJ per mole of reaction as written
    enthalpy_source: where the enthalpy came from ("" = unknown)
    """
    equation: str
    reactants: list[Species]
    products: list[Species]
    total_reactant_mass: float
    enthalpy: float
    enthalpy_source: str


def make_new_state() -> State:
    """Return a brand-new, empty State for starting the site."""
    return State("", [], [], 0.0, 0.0, "")

# ===========================================================================
# SMALL HELPERS
# ===========================================================================

def format_number(value: float) -> str:
    """Round a number to 3 decimal places and turn it into text."""
    return str(round(value, 3))


def is_number(text: str) -> bool:
    """
    Check whether text is a (possibly negative) decimal number like
    "12", "-571.6" or ".5", without using try/except.
    """
    cleaned = text.strip()
    if cleaned.startswith("-"):
        cleaned = cleaned[1:]
    digit_count = 0
    dot_count = 0
    for character in cleaned:
        if character.isdigit():
            digit_count += 1
        elif character == ".":
            dot_count += 1
        else:
            return False
    return digit_count > 0 and dot_count <= 1


def greatest_common_divisor(first: int, second: int) -> int:
    """Return the greatest common divisor using Euclid's algorithm."""
    first = abs(first)
    second = abs(second)
    while second != 0:
        first, second = second, first % second
    return first


def least_common_multiple(first: int, second: int) -> int:
    """Return the least common multiple of two positive integers."""
    return first * second // greatest_common_divisor(first, second)


# ===========================================================================
# ELEMENTS AND FORMULAS
# ===========================================================================

def get_atomic_mass(symbol: str) -> float:
    """Return an element's molar mass, or 0.0 if it is not an element."""
    for element in ELEMENTS:
        if element[0] == symbol:
            return element[2]
    return 0.0


def read_number(text: str, index: int) -> list[int]:
    """
    Read the whole number that starts at text[index].
    Returns [number, index right after the number]. If there is no number
    there, the number is 1 (like the hidden 1 in "H2O" after the O).
    """
    digits = ""
    while index < len(text) and text[index].isdigit():
        digits += text[index]
        index += 1
    if digits == "":
        return [1, index]
    return [int(digits), index]


def add_atoms(counts: dict[str, int], symbol: str, amount: int) -> None:
    """Add amount atoms of symbol to the counts dictionary."""
    if symbol in counts:
        counts[symbol] += amount
    else:
        counts[symbol] = amount


def count_atoms(formula: str) -> dict[str, int]:
    """
    Count the atoms of each element in a formula.
    count_atoms("Ca(OH)2") -> {"Ca": 1, "O": 2, "H": 2}
    Returns an empty dictionary if the formula is not valid.

    ADVANCED TECHNIQUE - a "stack" (list used as a pile of dictionaries):
    Every "(" puts a new empty dictionary on top of the pile. Atoms are
    always added to the dictionary on top. At ")" the top dictionary is
    taken off, multiplied by the number after the ")", and added into the
    dictionary underneath. This is how nested groups like Ca(OH)2 or
    Fe2(SO4)3 get counted correctly.
    """
    stack = [{}]
    index = 0
    while index < len(formula):
        character = formula[index]
        if character == "(":
            stack.append({})
            index += 1
        elif character == ")":
            if len(stack) == 1:
                return {}
            group = stack.pop()
            number_and_index = read_number(formula, index + 1)
            index = number_and_index[1]
            for symbol in group:
                add_atoms(stack[-1], symbol,
                          group[symbol] * number_and_index[0])
        elif character.isupper():
            symbol = character
            index += 1
            if index < len(formula) and formula[index].islower():
                symbol += formula[index]
                index += 1
            if get_atomic_mass(symbol) == 0.0:
                return {}
            number_and_index = read_number(formula, index)
            index = number_and_index[1]
            add_atoms(stack[-1], symbol, number_and_index[0])
        else:
            return {}
    if len(stack) != 1:
        return {}
    return stack[0]


def calculate_molar_mass(formula: str) -> float:
    """Return the molar mass of a formula in g/mol, or 0.0 if invalid."""
    counts = count_atoms(formula)
    total = 0.0
    for symbol in counts:
        total += get_atomic_mass(symbol) * counts[symbol]
    return total


# ===========================================================================
# READING EQUATIONS TYPED BY THE USER
# ===========================================================================

def split_side(side_text: str) -> list[str]:
    """Split "2 H2 + O2" into the terms ["2 H2", "O2"]."""
    terms = []
    for piece in side_text.split("+"):
        cleaned = piece.strip()
        if cleaned != "":
            terms.append(cleaned)
    return terms


def make_species(term: str) -> Species:
    """Turn one term like "2 H2O" or "2H2O" into a Species."""
    number_and_index = read_number(term, 0)
    formula = term[number_and_index[1]:].strip()
    return Species(formula, number_and_index[0],
                   calculate_molar_mass(formula), 0.0)


def parse_side(side_text: str) -> list[Species]:
    """Turn one side of an equation into a list of Species."""
    species_list = []
    for term in split_side(side_text):
        species_list.append(make_species(term))
    return species_list


def check_sides(reactants: list[Species], products: list[Species]) -> str:
    """
    Make sure both sides exist and contain real formulas.
    Returns an error message, or "" if everything is valid.
    """
    if len(reactants) == 0 or len(products) == 0:
        return "Both the reactant side and the product side need a formula."
    for species in reactants + products:
        if species.molar_mass == 0.0 or species.coefficient == 0:
            return ("'" + species.formula + "' is not a valid formula. "
                    "Check capital letters, e.g. CO (carbon monoxide) "
                    "vs Co (cobalt).")
    return ""


def side_to_text(species_list: list[Species]) -> str:
    """Turn a list of Species into text like "2 H2 + O2"."""
    terms = []
    for species in species_list:
        if species.coefficient == 1:
            terms.append(species.formula)
        else:
            terms.append(str(species.coefficient) + " " + species.formula)
    return " + ".join(terms)


def equation_to_text(reactants: list[Species],
                     products: list[Species]) -> str:
    """Build the full equation text, like "2 H2 + O2 → 2 H2O"."""
    return side_to_text(reactants) + " → " + side_to_text(products)


def parse_equation(equation: str) -> list[list[Species]]:
    """
    Split a stored equation at the arrow.
    Returns [reactants, products]; both are empty if there is no arrow.
    """
    sides = equation.split("→")
    if len(sides) != 2:
        return [[], []]
    return [parse_side(sides[0]), parse_side(sides[1])]


# ===========================================================================
# CHECKING AND BALANCING
# ===========================================================================

def total_atoms(species_list: list[Species]) -> dict[str, int]:
    """Count every atom on one side, including the coefficients."""
    totals = {}
    for species in species_list:
        counts = count_atoms(species.formula)
        for symbol in counts:
            add_atoms(totals, symbol, counts[symbol] * species.coefficient)
    return totals


def is_balanced(reactants: list[Species], products: list[Species]) -> bool:
    """Return True if both sides have the same number of each atom."""
    return total_atoms(reactants) == total_atoms(products)


def list_elements(species_list: list[Species]) -> list[str]:
    """List each element that appears, in order, without repeats."""
    elements = []
    for species in species_list:
        for symbol in count_atoms(species.formula):
            if symbol not in elements:
                elements.append(symbol)
    return elements


def build_atom_matrix(reactants: list[Species],
                      products: list[Species]) -> list[list[int]]:
    """
    Make one row per element and one column per species. Each number is
    how many atoms of that element the species has; products are negative
    because they are on the other side of the arrow.
    """
    matrix = []
    for symbol in list_elements(reactants + products):
        row = []
        for species in reactants:
            row.append(count_atoms(species.formula).get(symbol, 0))
        for species in products:
            row.append(-count_atoms(species.formula).get(symbol, 0))
        matrix.append(row)
    return matrix


def eliminate(target: list[int], pivot: list[int], column: int) -> list[int]:
    """
    Make target[column] zero by combining it with the pivot row, using only
    whole numbers (no fractions), then shrink the row by its common factor.
    """
    new_row = []
    for position in range(len(target)):
        new_row.append(target[position] * pivot[column]
                       - pivot[position] * target[column])
    common = 0
    for value in new_row:
        common = greatest_common_divisor(common, value)
    if common > 1:
        for position in range(len(new_row)):
            new_row[position] = new_row[position] // common
    return new_row


def solve_coefficients(matrix: list[list[int]], column_count: int) -> list[int]:
    """
    Find the smallest positive whole numbers x so that matrix * x = 0.
    Returns [] if there is no single answer.

    ADVANCED TECHNIQUE - Gauss-Jordan elimination (linear algebra):
    Balancing is the same as solving one equation per element, e.g. for
    a H2 + b O2 -> c H2O:  H: 2a - 2c = 0  and  O: 2b - c = 0.
    We clear out each column one at a time (like solving a system of
    equations by elimination in algebra class). Exactly one column should
    be left "free". We give that species a value, read every other
    coefficient off its row, then divide by the common factor so the
    coefficients are as small as possible.
    """
    rows = []
    for row in matrix:
        rows.append(row[:])
    pivot_columns = []
    for column in range(column_count):
        pivot_row = len(pivot_columns)
        found = -1
        for row_index in range(pivot_row, len(rows)):
            if found == -1 and rows[row_index][column] != 0:
                found = row_index
        if found != -1:
            rows[pivot_row], rows[found] = rows[found], rows[pivot_row]
            for row_index in range(len(rows)):
                if row_index != pivot_row and rows[row_index][column] != 0:
                    rows[row_index] = eliminate(rows[row_index],
                                                rows[pivot_row], column)
            pivot_columns.append(column)
    free_columns = []
    for column in range(column_count):
        if column not in pivot_columns:
            free_columns.append(column)
    if len(free_columns) != 1:
        return []
    return read_solution(rows, pivot_columns, free_columns[0], column_count)


def read_solution(rows: list[list[int]], pivot_columns: list[int],
                  free_column: int, column_count: int) -> list[int]:
    """
    Helper for solve_coefficients: after elimination, row i says
    pivot * x[pivot column] + value * x[free column] = 0.
    Choose x[free column] so every coefficient is a whole number.
    """
    free_value = 1
    for row_index in range(len(pivot_columns)):
        pivot = abs(rows[row_index][pivot_columns[row_index]])
        free_value = least_common_multiple(free_value, pivot)
    solution = [0] * column_count
    solution[free_column] = free_value
    for row_index in range(len(pivot_columns)):
        column = pivot_columns[row_index]
        solution[column] = (-rows[row_index][free_column] * free_value
                            // rows[row_index][column])
    common = 0
    for value in solution:
        common = greatest_common_divisor(common, value)
    for position in range(column_count):
        solution[position] = solution[position] // common
        if solution[position] <= 0:
            return []
    return solution


def balance(reactants: list[Species], products: list[Species]) -> bool:
    """
    Set the coefficient of every species so the equation is balanced.
    Returns False (and changes nothing) if it cannot be balanced.
    """
    matrix = build_atom_matrix(reactants, products)
    coefficients = solve_coefficients(matrix, len(reactants) + len(products))
    if len(coefficients) == 0:
        return False
    position = 0
    for species in reactants + products:
        species.coefficient = coefficients[position]
        position += 1
    return True


# ===========================================================================
# STOICHIOMETRY CALCULATIONS
# ===========================================================================

def reaction_extent(state: State) -> float:
    """
    Moles of reaction that can happen: the smallest (moles / coefficient)
    among reactants whose mass was entered. 0.0 if no mass was entered.
    """
    smallest = -1.0
    for reactant in state.reactants:
        if reactant.mass > 0:
            ratio = reactant.mass / reactant.molar_mass / reactant.coefficient
            if smallest < 0 or ratio < smallest:
                smallest = ratio
    if smallest < 0:
        return 0.0
    return smallest


def find_limiting_reactant(state: State) -> str:
    """Return the formula of the limiting reactant ("" if none)."""
    extent = reaction_extent(state)
    for reactant in state.reactants:
        if reactant.mass > 0:
            ratio = reactant.mass / reactant.molar_mass / reactant.coefficient
            if ratio == extent:
                return reactant.formula
    return ""


def theoretical_yield(state: State, product: Species) -> float:
    """Grams of product made if the limiting reactant is fully used."""
    return reaction_extent(state) * product.coefficient * product.molar_mass


def leftover_mass(state: State, reactant: Species) -> float:
    """Grams of a reactant left over once the reaction stops."""
    used = reaction_extent(state) * reactant.coefficient * reactant.molar_mass
    return max(0.0, reactant.mass - used)


def get_enthalpy_of_formation(formula: str) -> float:
    """Return a formula's standard enthalpy of formation in kJ/mol."""
    for entry in ENTHALPIES_OF_FORMATION:
        if entry[0] == formula:
            return entry[1]
    return 0.0


def has_enthalpy_data(formula: str) -> bool:
    """
    True if the formula is in the enthalpy table or is a pure element
    (an element in its standard state has an enthalpy of formation of 0).
    """
    for entry in ENTHALPIES_OF_FORMATION:
        if entry[0] == formula:
            return True
    return len(count_atoms(formula)) == 1


def can_use_enthalpy_table(state: State) -> bool:
    """True if every species in the reaction has enthalpy data."""
    for species in state.reactants + state.products:
        if not has_enthalpy_data(species.formula):
            return False
    return len(state.reactants) > 0


def enthalpy_from_table(state: State) -> float:
    """
    Hess's law: ΔH(reaction) = sum of n·ΔHf(products)
                               - sum of n·ΔHf(reactants)
    """
    total = 0.0
    for product in state.products:
        total += product.coefficient * get_enthalpy_of_formation(
            product.formula)
    for reactant in state.reactants:
        total -= reactant.coefficient * get_enthalpy_of_formation(
            reactant.formula)
    return total


def describe_heat(enthalpy: float) -> str:
    """Say whether a reaction with this ΔH is exothermic or endothermic."""
    if enthalpy < 0:
        return "Exothermic: ΔH is negative, so the reaction releases heat."
    elif enthalpy > 0:
        return "Endothermic: ΔH is positive, so the reaction absorbs heat."
    return "Thermoneutral: ΔH is zero, so no net heat flows."
def percent_yield(state: State, product: Species) -> float:
    """
    Percent yield = harvested mass / theoretical yield * 100.
    Returns 0.0 if there is no theoretical yield to compare against.
    """
    theoretical = theoretical_yield(state, product)
    if theoretical == 0.0:
        return 0.0
    return product.mass / theoretical * 100


# ===========================================================================
# PAGE BUILDERS (helper functions that make Pages; not routes)
# ===========================================================================

def navigation() -> list:
    """The button at the bottom of every page that goes back home."""
    return [HorizontalRule(), Button("Back to Home", index)]


def error_lines(error: str) -> list:
    """Show an error message only if there is one."""
    if error == "":
        return []
    return ["⚠ " + error]


def render_balance_page(state: State, error: str) -> Page:
    """The page where the user types an unbalanced equation."""
    return Page(state, [Header("Balance an Equation")]
                + error_lines(error)
                + ["Type the formulas on each side, separated by +. "
                   "Example: H2 + O2 → H2O",
                   "Reactants:", TextBox("reactant_side", ""),
                   "→ Products:", TextBox("product_side", ""),
                   Button("Balance", balance_reaction)]
                + navigation())


def atom_count_rows(reactants: list[Species],
                    products: list[Species]) -> list[list[str]]:
    """Rows comparing atom counts on each side (for 'showing the work')."""
    left = total_atoms(reactants)
    right = total_atoms(products)
    rows = [["Element", "Reactant side", "Product side"]]
    for symbol in list_elements(reactants + products):
        rows.append([symbol, str(left.get(symbol, 0)),
                     str(right.get(symbol, 0))])
    return rows


def render_balance_result(state: State) -> Page:
    """
    Show the balanced equation, then let the user choose between going
    home or calculating the stats of this reaction.
    """
    sides = parse_equation(state.equation)
    return Page(state, [Header("Balanced Equation"),
                        state.equation,
                        "Atom check (every element must match):",
                        Table(atom_count_rows(sides[0], sides[1])),
                        "What would you like to do next?",
                        Button("Calculate Reaction Stats", move_equation),
                        Button("Back to Home", index)])


def species_rows(state: State) -> list[list[str]]:
    """Rows describing every reactant and product on the calculation page."""
    rows = [["Formula", "Role", "Coefficient", "Molar mass (g/mol)",
             "Mass entered (g)"]]
    for reactant in state.reactants:
        rows.append([reactant.formula, "Reactant (starting mass)",
                     str(reactant.coefficient),
                     format_number(reactant.molar_mass),
                     format_number(reactant.mass)])
    for product in state.products:
        rows.append([product.formula, "Product (harvested mass)",
                     str(product.coefficient),
                     format_number(product.molar_mass),
                     format_number(product.mass)])
    return rows


def formulas_of(species_list: list[Species]) -> list[str]:
    """List just the formulas, for SelectBox options."""
    formulas = []
    for species in species_list:
        formulas.append(species.formula)
    return formulas


def mass_entry_section(state: State) -> list:
    """Steps 2-4 of the calculation page (shown once an equation exists)."""
    if can_use_enthalpy_table(state):
        heat_note = ("ΔH can be found from the built-in enthalpy table, "
                     "but a value you type here will be used instead.")
    else:
        heat_note = ("Some formulas are not in the built-in enthalpy "
                     "table, so type ΔH here if you want heat results.")
    all_species = state.reactants + state.products
    return [Header("Current Equation", 2), state.equation,
            Table(species_rows(state)),
            Header("Step 2: Enter the masses", 2),
            "Pick a reactant and enter the grams you started with. "
            "Then pick a product and enter the grams you harvested. "
            "Repeat for each one.",
            "Chemical:",
            SelectBox("species_choice", formulas_of(all_species),
                      all_species[0].formula),
            "Mass in grams:", TextBox("mass_text", ""),
            Button("Set Mass", input_mass),
            "Total reactant mass: "
            + format_number(state.total_reactant_mass) + " g",
            Header("Step 3: Heat of reaction (optional)", 2),
            heat_note,
            "ΔH in kJ per mole of reaction:", TextBox("enthalpy_text", ""),
            Header("Step 4: Calculate", 2),
            Button("Calculate", calculate)]


def render_stats_page(state: State, error: str) -> Page:
    """The reaction statistics (calculation) page."""
    content = [Header("Calculate Reaction Stats")] + error_lines(error)
    content += [Header("Step 1: Enter a balanced equation", 2),
                "Reactants:",
                TextBox("reactant_side", side_to_text(state.reactants)),
                "→ Products:",
                TextBox("product_side", side_to_text(state.products)),
                Button("Load Equation", load_equation)]
    if len(state.reactants) > 0:
        content += mass_entry_section(state)
    return Page(state, content + navigation())


def mole_rows(state: State) -> list[list[str]]:
    """Grams → moles → moles ÷ coefficient for each reactant."""
    rows = [["Reactant", "Mass (g)", "Molar mass (g/mol)", "Moles",
             "Moles ÷ coefficient"]]
    for reactant in state.reactants:
        if reactant.mass > 0:
            moles = reactant.mass / reactant.molar_mass
            rows.append([reactant.formula, format_number(reactant.mass),
                         format_number(reactant.molar_mass),
                         format_number(moles),
                         format_number(moles / reactant.coefficient)])
    return rows


def leftover_rows(state: State) -> list[list[str]]:
    """How much of each reactant is left over."""
    rows = [["Reactant", "Started with (g)", "Left over (g)"]]
    for reactant in state.reactants:
        if reactant.mass > 0:
            rows.append([reactant.formula, format_number(reactant.mass),
                         format_number(leftover_mass(state, reactant))])
        else:
            rows.append([reactant.formula, "not entered",
                         "assumed in excess"])
    return rows


def yield_rows(state: State) -> list[list[str]]:
    """Theoretical yield, harvested mass and percent yield of each product."""
    rows = [["Product", "Theoretical yield (g)", "Harvested (g)",
             "Percent yield"]]
    for product in state.products:
        theoretical = format_number(theoretical_yield(state, product))
        if product.mass > 0:
            rows.append([product.formula, theoretical,
                         format_number(product.mass),
                         format_number(percent_yield(state, product))
                         + " %"])
        else:
            rows.append([product.formula, theoretical, "not entered", "-"])
    return rows


def yield_warnings(state: State) -> list:
    """Warn about any percent yield over 100%."""
    for product in state.products:
        if percent_yield(state, product) > 100:
            return ["A percent yield over 100% usually means the product "
                    "was still wet or impure, or a mass was measured "
                    "wrong."]
    return []


def heat_lines(state: State) -> list:
    """Exothermic/endothermic and heat flow per mole of product."""
    if state.enthalpy_source == "":
        return ["No ΔH available for this reaction. Go back and type one "
                "in to see heat results."]
    lines = ["ΔH = " + format_number(state.enthalpy)
             + " kJ per mole of reaction (" + state.enthalpy_source + ")",
             describe_heat(state.enthalpy)]
    rows = [["Product", "Heat flow per mole of product (kJ/mol)"]]
    for product in state.products:
        rows.append([product.formula,
                     format_number(state.enthalpy / product.coefficient)])
    total_heat = state.enthalpy * reaction_extent(state)
    return lines + [Table(rows),
                    "Total heat flow for your amounts: "
                    + format_number(total_heat) + " kJ"]


def render_results_page(state: State) -> Page:
    """Show the full breakdown of the reaction, step by step."""
    limiting = find_limiting_reactant(state)
    return Page(state, [Header("Reaction Results"), state.equation,
                        Header("Limiting Reactant", 2),
                        Table(mole_rows(state)),
                        "Limiting reactant: " + limiting
                        + " (it has the smallest moles ÷ coefficient)",
                        Table(leftover_rows(state)),
                        Header("Theoretical and Percent Yield", 2),
                        Table(yield_rows(state))]
                + yield_warnings(state)
                + [Header("Heat of Reaction", 2)]
                + heat_lines(state)
                + [Button("Change Masses", stats_page)]
                + navigation())


# ===========================================================================
# ROUTES
# ===========================================================================

@route
def index(state: State) -> Page:
    """Home page with exactly five options."""
    return Page(state, [Header("Stoichiometry Assistant"),
                        "Choose an option:",
                        Button("Balance an Equation", balance_page),
                        Button("Calculate Reaction Stats", stats_page),
                        Button("Periodic Table", periodic_table),
                        Button("Common Enthalpies", common_enthalpies),
                        Button("Molar Masses of Common Compounds",
                               common_molar_masses)])


@route
def balance_page(state: State) -> Page:
    """Show the balancing page."""
    return render_balance_page(state, "")


@route
def balance_reaction(state: State, reactant_side: str,
                     product_side: str) -> Page:
    """Balance the equation the user typed and save it in the State."""
    reactants = parse_side(reactant_side)
    products = parse_side(product_side)
    error = check_sides(reactants, products)
    if error != "":
        return render_balance_page(state, error)
    if not balance(reactants, products):
        return render_balance_page(state, "That reaction can't be balanced "
                                   "in a single way. Check the formulas.")
    state.equation = equation_to_text(reactants, products)
    return render_balance_result(state)


@route
def move_equation(state: State) -> Page:
    """Move the balanced equation over to the calculation page."""
    sides = parse_equation(state.equation)
    if len(sides[0]) == 0:
        return render_balance_page(state, "Balance an equation first.")
    state.reactants = sides[0]
    state.products = sides[1]
    state.total_reactant_mass = 0.0
    state.enthalpy_source = ""
    return render_stats_page(state, "")


@route
def stats_page(state: State) -> Page:
    """Show the calculation page."""
    return render_stats_page(state, "")


@route
def load_equation(state: State, reactant_side: str,
                  product_side: str) -> Page:
    """Load a balanced equation typed straight into the calculation page."""
    reactants = parse_side(reactant_side)
    products = parse_side(product_side)
    error = check_sides(reactants, products)
    if error != "":
        return render_stats_page(state, error)
    if not is_balanced(reactants, products):
        return render_stats_page(state, "That equation is not balanced. "
                                 "Use Balance an Equation on the home page "
                                 "to fix it first.")
    state.reactants = reactants
    state.products = products
    state.equation = equation_to_text(reactants, products)
    state.total_reactant_mass = 0.0
    state.enthalpy_source = ""
    return render_stats_page(state, "")


@route
def input_mass(state: State, species_choice: str, mass_text: str) -> Page:
    """
    Save a mass: the starting mass for a reactant, or the harvested mass
    for a product.
    """
    if not is_number(mass_text) or float(mass_text) <= 0:
        return render_stats_page(state, "The mass must be a positive "
                                 "number of grams.")
    for product in state.products:
        if product.formula == species_choice:
            product.mass = float(mass_text)
    total = 0.0
    for reactant in state.reactants:
        if reactant.formula == species_choice:
            reactant.mass = float(mass_text)
        total += reactant.mass
    state.total_reactant_mass = total
    return render_stats_page(state, "")


@route
def calculate(state: State, enthalpy_text: str) -> Page:
    """Run every calculation and show the results page."""
    if reaction_extent(state) == 0.0:
        return render_stats_page(state, "Enter the mass of at least one "
                                 "reactant first.")
    if enthalpy_text.strip() != "":
        if not is_number(enthalpy_text):
            return render_stats_page(state, "ΔH must be a number, "
                                     "like -571.6")
        state.enthalpy = float(enthalpy_text)
        state.enthalpy_source = "entered by you"
    elif can_use_enthalpy_table(state):
        state.enthalpy = enthalpy_from_table(state)
        state.enthalpy_source = ("from standard enthalpies of formation, "
                                 "OpenStax Appendix G")
    else:
        state.enthalpy_source = ""
    return render_results_page(state)


@route
def periodic_table(state: State) -> Page:
    """List every element with its molar mass."""
    rows = [["#", "Symbol", "Name", "Molar mass (g/mol)"]]
    for number in range(len(ELEMENTS)):
        element = ELEMENTS[number]
        rows.append([str(number + 1), element[0], element[1],
                     str(element[2])])
    return Page(state, [Header("Periodic Table of Molar Masses"),
                        Table(rows)] + navigation())


@route
def common_enthalpies(state: State) -> Page:
    """List the standard enthalpies of formation the site knows."""
    rows = [["Formula", "ΔHf° (kJ/mol)"]]
    for entry in ENTHALPIES_OF_FORMATION:
        rows.append([entry[0], str(entry[1])])
    return Page(state, [Header("Common Enthalpies of Formation"),
                        "Standard enthalpies of formation at 298 K "
                        "(OpenStax Chemistry 2e, Appendix G). Each compound "
                        "is in its usual state at room temperature, so H2O "
                        "is liquid water. Pure elements such as O2, H2 and "
                        "Fe are 0.",
                        Table(rows)] + navigation())


@route
def common_molar_masses(state: State) -> Page:
    """List molar masses of common compounds (calculated, not typed in)."""
    rows = [["Formula", "Name", "Molar mass (g/mol)"]]
    for compound in COMMON_COMPOUNDS:
        rows.append([compound[0], compound[1],
                     format_number(calculate_molar_mass(compound[0]))])
    return Page(state, [Header("Molar Masses of Common Compounds"),
                        Table(rows)] + navigation())



# ===========================================================================
# UNIT TESTS
# ===========================================================================

def water_state() -> State:
    """Test helper: 2 H2 + O2 → 2 H2O with 10 g H2 and 32 g O2 entered."""
    return State("2 H2 + O2 → 2 H2O",
                 [Species("H2", 2, 2.016, 10.0),
                  Species("O2", 1, 31.998, 32.0)],
                 [Species("H2O", 2, 18.015, 0.0)],
                 42.0, 0.0, "")


def harvested_state() -> State:
    """Test helper: water_state, plus 30 g of H2O harvested."""
    state = water_state()
    state.products[0].mass = 30.0
    return state


# --- helpers ---------------------------------------------------------------
assert_equal(is_number("-571.6"), True)
assert_equal(is_number("12"), True)
assert_equal(is_number("abc"), False)
assert_equal(is_number("1.2.3"), False)
assert_equal(is_number(""), False)
assert_equal(greatest_common_divisor(12, 18), 6)
assert_equal(least_common_multiple(4, 6), 12)
assert_equal(format_number(18.01500001), "18.015")

# --- formulas --------------------------------------------------------------
assert_equal(get_atomic_mass("O"), 15.999)
assert_equal(get_atomic_mass("Xx"), 0.0)
assert_equal(count_atoms("H2O"), {"H": 2, "O": 1})
assert_equal(count_atoms("Ca(OH)2"), {"Ca": 1, "O": 2, "H": 2})
assert_equal(count_atoms("Fe2(SO4)3"), {"Fe": 2, "S": 3, "O": 12})
assert_equal(count_atoms("Xy2"), {})
assert_equal(count_atoms("H2O)"), {})
assert_equal(count_atoms("co2"), {})
assert_equal(calculate_molar_mass("H2O"), 18.015)
assert_equal(calculate_molar_mass("C6H12O6"), 180.156)
assert_equal(calculate_molar_mass("Qq"), 0.0)
assert_equal(make_species("2 H2O"), Species("H2O", 2, 18.015, 0.0))
assert_equal(make_species("3O2"), Species("O2", 3, 31.998, 0.0))
assert_equal(split_side(" 2 H2 + O2 "), ["2 H2", "O2"])
assert_equal(len(parse_side("")), 0)
assert_equal(check_sides(parse_side("H2"), parse_side("")),
             "Both the reactant side and the product side need a formula.")
assert_equal(check_sides(parse_side("H2"), parse_side("H2O")), "")
assert_equal(side_to_text([Species("H2", 2, 2.016, 0.0),
                           Species("O2", 1, 31.998, 0.0)]), "2 H2 + O2")
assert_equal(formulas_of(water_state().reactants), ["H2", "O2"])
assert_equal(parse_equation("no arrow"), [[], []])
assert_equal(len(parse_equation("2 H2 + O2 → 2 H2O")[0]), 2)

# --- balancing -------------------------------------------------------------
assert_equal(is_balanced(parse_side("2 H2 + O2"), parse_side("2 H2O")), True)
assert_equal(is_balanced(parse_side("H2 + O2"), parse_side("H2O")), False)
assert_equal(list_elements(parse_side("CH4 + O2")), ["C", "H", "O"])
assert_equal(build_atom_matrix(parse_side("H2 + O2"), parse_side("H2O")),
             [[2, 0, -2], [0, 2, -1]])
assert_equal(solve_coefficients([[2, 0, -2], [0, 2, -1]], 3), [2, 1, 2])
assert_equal(solve_coefficients([[1, -1]], 3), [])


def balanced_text(left: str, right: str) -> str:
    """Test helper: balance an equation and return its text."""
    reactants = parse_side(left)
    products = parse_side(right)
    if not balance(reactants, products):
        return "cannot balance"
    return equation_to_text(reactants, products)


assert_equal(balanced_text("H2 + O2", "H2O"), "2 H2 + O2 → 2 H2O")
assert_equal(balanced_text("CH4 + O2", "CO2 + H2O"),
             "CH4 + 2 O2 → CO2 + 2 H2O")
assert_equal(balanced_text("C8H18 + O2", "CO2 + H2O"),
             "2 C8H18 + 25 O2 → 16 CO2 + 18 H2O")
assert_equal(balanced_text("Fe + O2", "Fe2O3"), "4 Fe + 3 O2 → 2 Fe2O3")
assert_equal(balanced_text("Al + CuSO4", "Al2(SO4)3 + Cu"),
             "2 Al + 3 CuSO4 → Al2(SO4)3 + 3 Cu")
assert_equal(balanced_text("H2", "O2"), "cannot balance")

# --- stoichiometry ---------------------------------------------------------
assert_equal(reaction_extent(water_state()), 32.0 / 31.998)
assert_equal(reaction_extent(make_new_state()), 0.0)
assert_equal(find_limiting_reactant(water_state()), "O2")
assert_equal(theoretical_yield(water_state(), water_state().products[0]),
             2 * 18.015 * 32.0 / 31.998)
assert_equal(leftover_mass(water_state(), water_state().reactants[1]), 0.0)
assert_equal(leftover_mass(water_state(), water_state().reactants[0]),
             10.0 - 2 * 2.016 * 32.0 / 31.998)
assert_equal(get_enthalpy_of_formation("H2O"), -285.8)
assert_equal(has_enthalpy_data("O2"), True)
assert_equal(has_enthalpy_data("KMnO4"), False)
assert_equal(can_use_enthalpy_table(water_state()), True)
assert_equal(can_use_enthalpy_table(make_new_state()), False)
assert_equal(enthalpy_from_table(water_state()), -571.6)
assert_equal(describe_heat(-571.6),
             "Exothermic: ΔH is negative, so the reaction releases heat.")
assert_equal(describe_heat(180.0),
             "Endothermic: ΔH is positive, so the reaction absorbs heat.")
assert_equal(describe_heat(0.0),
             "Thermoneutral: ΔH is zero, so no net heat flows.")
assert_equal(mole_rows(water_state())[2],
             ["O2", "32.0", "31.998", "1.0", "1.0"])
assert_equal(leftover_rows(water_state())[1], ["H2", "10.0", "5.968"])
assert_equal(percent_yield(water_state(), water_state().products[0]), 0.0)
assert_equal(percent_yield(harvested_state(), harvested_state().products[0]),
             30.0 / (2 * 18.015 * 32.0 / 31.998) * 100)
assert_equal(percent_yield(make_new_state(), Species("H2O", 2, 18.015, 5.0)),
             0.0)
assert_equal(mole_rows(water_state())[2],
             ["O2", "32.0", "31.998", "1.0", "1.0"])
assert_equal(leftover_rows(water_state())[1], ["H2", "10.0", "5.968"])
assert_equal(yield_rows(water_state())[1],
             ["H2O", "36.032", "not entered", "-"])
assert_equal(yield_rows(harvested_state())[1],
             ["H2O", "36.032", "30.0", "83.259 %"])
assert_equal(yield_warnings(harvested_state()), [])
assert_equal(len(yield_warnings(State("2 H2 + O2 → 2 H2O",
                                      [Species("H2", 2, 2.016, 10.0)],
                                      [Species("H2O", 2, 18.015, 99.0)],
                                      10.0, 0.0, ""))), 1)
assert_equal(error_lines(""), [])
assert_equal(error_lines("Oops"), ["⚠ Oops"])

# --- routes (every route has at least one test) ----------------------------
assert_equal(index(make_new_state()),
             Page(make_new_state(),
                  [Header("Stoichiometry Assistant"),
                   "Choose an option:",
                   Button("Balance an Equation", balance_page),
                   Button("Calculate Reaction Stats", stats_page),
                   Button("Periodic Table", periodic_table),
                   Button("Common Enthalpies", common_enthalpies),
                   Button("Molar Masses of Common Compounds",
                          common_molar_masses)]))

assert_equal(balance_page(make_new_state()),
             Page(make_new_state(),
                  [Header("Balance an Equation"),
                   "Type the formulas on each side, separated by +. "
                   "Example: H2 + O2 → H2O",
                   "Reactants:", TextBox("reactant_side", ""),
                   "→ Products:", TextBox("product_side", ""),
                   Button("Balance", balance_reaction),
                   HorizontalRule(), Button("Back to Home", index)]))

assert_equal(balance_reaction(make_new_state(), "H2 + O2", "H2O").state
             .equation, "2 H2 + O2 → 2 H2O")
assert_equal(balance_reaction(make_new_state(), "H2 + O2", "H2O").content,
             [Header("Balanced Equation"), "2 H2 + O2 → 2 H2O",
              "Atom check (every element must match):",
              Table([["Element", "Reactant side", "Product side"],
                     ["H", "4", "4"], ["O", "2", "2"]]),
              "What would you like to do next?",
              Button("Calculate Reaction Stats", move_equation),
              Button("Back to Home", index)])
assert_equal(balance_reaction(make_new_state(), "H2 + Oo", "H2O").content[1],
             "⚠ 'Oo' is not a valid formula. Check capital letters, e.g. "
             "CO (carbon monoxide) vs Co (cobalt).")
assert_equal(balance_reaction(make_new_state(), "H2", "O2").content[1],
             "⚠ That reaction can't be balanced in a single way. "
             "Check the formulas.")

assert_equal(move_equation(State("2 H2 + O2 → 2 H2O", [], [], 0.0, 0.0,
                                 "")).state.reactants,
             [Species("H2", 2, 2.016, 0.0), Species("O2", 1, 31.998, 0.0)])
assert_equal(move_equation(make_new_state()).content[1],
             "⚠ Balance an equation first.")

assert_equal(stats_page(make_new_state()).content[0],
             Header("Calculate Reaction Stats"))
assert_equal(stats_page(water_state()).content[8], "2 H2 + O2 → 2 H2O")
assert_equal(stats_page(water_state()).content[13],
             SelectBox("species_choice", ["H2", "O2", "H2O"], "H2"))

assert_equal(load_equation(make_new_state(), "2 H2 + O2", "2 H2O").state
             .products, [Species("H2O", 2, 18.015, 0.0)])
assert_equal(load_equation(make_new_state(), "H2 + O2", "H2O").content[1],
             "⚠ That equation is not balanced. Use Balance an Equation on "
             "the home page to fix it first.")

assert_equal(input_mass(load_equation(make_new_state(), "2 H2 + O2",
                                      "2 H2O").state, "O2", "16")
             .state.reactants[1].mass, 16.0)
assert_equal(input_mass(water_state(), "H2", "4").state.total_reactant_mass,
             36.0)
assert_equal(input_mass(water_state(), "H2O", "30").state.products[0].mass,
             30.0)
assert_equal(input_mass(water_state(), "H2O", "30").state
             .total_reactant_mass, 42.0)
assert_equal(input_mass(water_state(), "H2", "-4").content[1],
             "⚠ The mass must be a positive number of grams.")

assert_equal(calculate(water_state(), "").state.enthalpy, -571.6)
assert_equal(calculate(water_state(), "-500").state.enthalpy_source,
             "entered by you")
assert_equal(calculate(water_state(), "hot").content[1],
             "⚠ ΔH must be a number, like -571.6")
assert_equal(calculate(move_equation(State("2 H2 + O2 → 2 H2O", [], [], 0.0,
                                           0.0, "")).state, "")
             .content[1], "⚠ Enter the mass of at least one reactant first.")
assert_equal(calculate(water_state(), "").content[4],
             "Limiting reactant: O2 (it has the smallest moles ÷ "
             "coefficient)")
assert_equal(calculate(harvested_state(), "").content[7],
             Table([["Product", "Theoretical yield (g)", "Harvested (g)",
                     "Percent yield"],
                    ["H2O", "36.032", "30.0", "83.259 %"]]))
assert_equal(calculate(harvested_state(), "").content[9],
             "ΔH = -571.6 kJ per mole of reaction (from standard enthalpies "
             "of formation, OpenStax Appendix G)")
assert_equal(calculate(harvested_state(), "").content[10],
             "Exothermic: ΔH is negative, so the reaction releases heat.")
assert_equal(calculate(harvested_state(), "").content[11],
             Table([["Product", "Heat flow per mole of product (kJ/mol)"],
                    ["H2O", "-285.8"]]))

assert_equal(periodic_table(make_new_state()).content[1].rows[1],
             ["1", "H", "Hydrogen", "1.008"])
assert_equal(len(periodic_table(make_new_state()).content[1].rows), 119)

assert_equal(common_enthalpies(make_new_state()).content[2].rows[1],
             ["H2O", "-285.8"])
assert_equal(len(common_enthalpies(make_new_state()).content[2].rows),
             len(ENTHALPIES_OF_FORMATION) + 1)

assert_equal(common_molar_masses(make_new_state()).content[1].rows[1],
             ["H2O", "Water", "18.015"])



# ===========================================================================
# SITE SETUP
# ===========================================================================

set_site_information(
    author="teov@udel.edu",
    description=("A stoichiometry assistant that balances equations and "
                 "finds the limiting reactant, leftovers, theoretical and "
                 "percent yield, and heat of reaction."),
    sources=["OpenStax Chemistry 2e, Chapter 4: Chemical Reactions and "
             "Stoichiometry",
             "OpenStax Chemistry 2e, Appendix G: Standard Thermodynamic "
             "Properties",
             "IUPAC standard atomic weights"],
    planning=["CS1_Final_Project_0__Pick_a_Problem.pdf",
              "CS1_Final_Project_1__Design.pdf"],
    links=["https://github.com/YOUR-USERNAME/YOUR-REPO"]
)
set_website_title("Stoichiometry Assistant")
start_server(make_new_state())

