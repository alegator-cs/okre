import re
import math
import itertools

#############################################
# Preprocessing: parse the factored equation #
#############################################

def preprocess_equation(eq_str):
    """
    Given an input string such as "5*x0*x1+2*x1*x2=9",
    this function does two things:
      (a) It splits the LHS into terms.
      (b) For each term, it extracts the coefficient and the list of x–variables (in order).
          It then assigns a new variable name y{i} for that term.
    Returns:
      - new_eq_str: a string of the form "5*y0+2*y1=9"
      - term_data: a list of dictionaries, one per term. Each dictionary has:
            'y'    : the new variable name (e.g. "y0")
            'coeff': the numeric coefficient (as an int)
            'xvars': a list of the x–variable names (e.g. ["x0", "x1"])
      - rhs: the right‐hand side as an integer.
    Assumptions: The input equation is of the form
         <term> (+ <term>)* = <number>
    and each term is of the form (optional coefficient)* product-of-x–variables,
    where an x–variable looks like "x<number>".
    """
    # Remove spaces
    eq_str = eq_str.replace(" ", "")
    # Split LHS and RHS
    parts = eq_str.split("=")
    if len(parts) != 2:
        raise ValueError("Equation must contain exactly one '=' sign.")
    lhs_str, rhs_str = parts
    try:
        rhs = int(rhs_str)
    except ValueError:
        raise ValueError("Right-hand side must be an integer.")
    
    # Split the LHS by '+' signs (we assume no minus signs for simplicity)
    term_strs = re.split(r'\+', lhs_str)
    term_data = []
    new_terms = []
    for idx, term in enumerate(term_strs):
        # Each term is expected to be like "5*x0*x1" or maybe "x0*x1" (implicit coefficient 1)
        # First, split by '*' signs.
        factors = term.split('*')
        # Look at the first factor: if it is purely numeric, that's the coefficient.
        if re.fullmatch(r'\d+', factors[0]):
            coeff = int(factors[0])
            xvars = factors[1:]
        else:
            coeff = 1
            xvars = factors
        # We require at least one xvar in the product.
        if not xvars:
            raise ValueError("Each term must contain at least one x-variable.")
        y_name = f"y{idx}"
        term_data.append({'y': y_name, 'coeff': coeff, 'xvars': xvars})
        new_terms.append(f"{coeff}*{y_name}")
    new_eq_str = "+".join(new_terms) + "=" + str(rhs)
    return new_eq_str, term_data, rhs

#############################################
# Solving the linear (preprocessed) equation #
#############################################

def solve_linear_equation(new_eq_str, term_data, rhs):
    """
    The new equation is of the form sum_i (coeff_i * y_i) = rhs.
    For each term i, we use a trivial bound for y_i:
         0 <= y_i <= floor(rhs / coeff_i)
    (Assuming coeff_i > 0.)
    Since the numbers are small, we enumerate over all possibilities.
    Returns a list of solutions in y-space.
    Each solution is a tuple (y0, y1, ...), in the order of term_data.
    """
    bounds = []
    for term in term_data:
        c = term['coeff']
        # trivial bound: 0 <= y <= floor(rhs/c)
        ub = math.floor(rhs / c)
        bounds.append((0, ub))
    # Enumerate over the Cartesian product of bounds.
    ranges = [range(lb, ub+1) for (lb, ub) in bounds]
    solutions = []
    for candidate in itertools.product(*ranges):
        total = sum(term_data[i]['coeff'] * candidate[i] for i in range(len(term_data)))
        if total == rhs:
            solutions.append(candidate)
    return solutions

#############################################
# Factorization: factor n into exactly r factors #
#############################################

def factorize_n(n, r):
    """
    Return all r-tuples of positive integers (a1, ..., ar) such that
          a1 * a2 * ... * ar == n.
    For simplicity, we assume n > 0 and r >= 1.
    (If n == 1, then the only factorization is (1,1,...,1).)
    """
    if r == 1:
        return [(n,)]
    solutions = []
    # Loop over possible first factor from 1 to n.
    for a in range(1, n+1):
        if n % a == 0:
            for tail in factorize_n(n // a, r-1):
                solutions.append((a,) + tail)
    return solutions

#############################################
# Enforcing equivalence constraints among x-vars #
#############################################

def filter_by_equivalences(solution_factorizations, term_data):
    """
    Given:
      - solution_factorizations: a list with one element per term.
        For term i, solution_factorizations[i] is a list of factorizations,
        where each factorization is a tuple of integers (the factors corresponding
        to the x variables in that term).
      - term_data: a list (one per term) with a key 'xvars' that is the ordered list
        of original x-variable names for that term.
    Some x–variables appear in more than one term; we must enforce that their
    factorizations agree.
    
    This function builds a dictionary mapping each x–variable to the list of positions
    where it appears: each position is (term_index, factor_index).
    Then it filters the Cartesian product over the factorizations for each term, keeping only
    those combinations for which every x–variable that appears multiple times is assigned the same value.
    
    Returns a list of combined factorizations. Each combined factorization is a tuple
    with one element per term, where each element is a factorization tuple.
    """
    # Build mapping: x_var -> list of (term index, factor index)
    x_positions = {}
    for i, term in enumerate(term_data):
        for j, xvar in enumerate(term['xvars']):
            x_positions.setdefault(xvar, []).append((i, j))
    
    # For each term, we have a list of possible factorizations.
    all_combos = itertools.product(*solution_factorizations)
    valid = []
    for combo in all_combos:
        consistent = True
        # For each x variable that appears in more than one term,
        # check that the factors in those positions are equal.
        for xvar, pos_list in x_positions.items():
            vals = []
            for (ti, fi) in pos_list:
                # For term ti, the chosen factorization is combo[ti]
                vals.append(combo[ti][fi])
            if len(set(vals)) > 1:
                consistent = False
                break
        if consistent:
            valid.append(combo)
    return valid

#############################################
# Postprocessing: combine with factorization lookup #
#############################################

def postprocess_solution(y_solution, term_data):
    """
    Given a solution in y-space (a tuple, one entry per term) and term_data,
    for each term i, factorize y_solution[i] into exactly r factors, where
    r = len(term_data[i]['xvars']).
    Then, take the Cartesian product over terms and filter by equivalence
    constraints (i.e. if the same x appears in multiple terms, the corresponding factors must agree).
    Returns the list of combined factorizations.
    Each combined factorization is a tuple with one element per term,
    and each element is a tuple of factors (corresponding to the x–variables in that term).
    """
    per_term_factorizations = []
    for i, y_val in enumerate(y_solution):
        r = len(term_data[i]['xvars'])
        # For y=0, factorization is not well defined (or one might allow a unique factorization)
        # Here we assume y > 0.
        if y_val <= 0:
            # For simplicity, if y_val == 0, we require all factors to be 0.
            # (This is ad hoc; adjust as needed.)
            fac = (0,) * r
            per_term_factorizations.append([fac])
        else:
            facs = factorize_n(y_val, r)
            per_term_factorizations.append(facs)
    # Now filter by equivalences.
    valid = filter_by_equivalences(per_term_factorizations, term_data)
    return valid

#############################################
# Main driver
#############################################

def main():
    # Get input equation string.
    eq_str = input("Enter a factored linear Diophantine equation (e.g. '5*x0*x1+2*x1*x2=9'): ")
    # Preprocess: convert product terms to new y variables.
    new_eq_str, term_data, rhs = preprocess_equation(eq_str)
    print("Preprocessed equation:", new_eq_str)
    print("Term data:")
    for term in term_data:
        print(f"  {term['y']}: coefficient {term['coeff']}, factors {term['xvars']}")
    
    # Solve the linear equation in the y variables.
    # Use trivial bounds: for each term, 0 <= y <= floor(rhs / coeff)
    y_solutions = solve_linear_equation(new_eq_str, term_data, rhs)
    if not y_solutions:
        print("No solutions in y-space found.")
        return
    print("Solutions in y-space:")
    for sol in y_solutions:
        print(sol)
    
    # For each y-solution, perform postprocessing: factorize each y into the required number of factors,
    # then filter by the equivalence constraints.
    overall_factorizations = []
    for y_sol in y_solutions:
        factorizations = postprocess_solution(y_sol, term_data)
        overall_factorizations.append(factorizations)
    
    # Print results.
    # Each overall_factorizations element corresponds to one solution in y-space.
    # For example, for y_sol = (1,2) one might obtain something like:
    #   [ (((1,1)), ((1,2))) ]
    print("\nFinal factorizations (each item is a tuple with one factorization per term):")
    for i, sol_fac in enumerate(overall_factorizations):
        print(f"Solution in y-space {y_solutions[i]} yields factorizations:")
        for combo in sol_fac:
            # combo is a tuple, one entry per term (each entry is a tuple of factors)
            print(combo)

if __name__ == '__main__':
    main()
