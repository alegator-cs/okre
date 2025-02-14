#!/usr/bin/env python3
import sys
import re
import math
import itertools

#############################################
# Preprocessing: parse the factored equation #
#############################################

def preprocess_equation(eq_str):
    """
    Given an input string such as "5*x0*x1+2*x1*x2+3=9",
    this function splits the left-hand side into terms, extracts the coefficient
    and the ordered list of x-variable names for each term, and assigns a new
    y-variable name for each nonconstant term.
    
    Constant-only terms (with no x-variable) are not converted;
    instead, their numeric value is subtracted from the rhs.
    
    Returns:
      - new_eq_str: a string like "5*y0+2*y1=<adjusted_rhs>"
      - term_data: a list of dicts, one per nonconstant term, each with keys:
            'y'    : the new variable name (e.g. "y0")
            'coeff': the numeric coefficient (as an int)
            'xvars': a list of the x-variable names (e.g. ["x0", "x1"])
      - rhs: the adjusted right-hand side as an integer.
    """
    # Remove spaces.
    eq_str = eq_str.replace(" ", "")
    parts = eq_str.split("=")
    if len(parts) != 2:
        raise ValueError("Equation must contain exactly one '=' sign.")
    lhs_str, rhs_str = parts
    try:
        rhs = int(rhs_str)
    except ValueError:
        raise ValueError("Right-hand side must be an integer.")
    
    term_strs = re.split(r'\+', lhs_str)
    term_data = []
    new_terms = []
    constant_sum = 0
    for term in term_strs:
        # Split the term by '*' characters.
        factors = term.split('*')
        # If there's only one factor and it is entirely numeric, it's a constant.
        if len(factors) == 1 and re.fullmatch(r'\d+', factors[0]):
            constant_sum += int(factors[0])
        else:
            # Otherwise, if the first factor is numeric, that's the coefficient.
            if re.fullmatch(r'\d+', factors[0]):
                coeff = int(factors[0])
                xvars = factors[1:]
            else:
                coeff = 1
                xvars = factors
            # If no x-variable is present, treat it as a constant.
            if not xvars:
                constant_sum += coeff
            else:
                # Generate a new y-variable name using the current count.
                y_name = f"y{len(term_data)}"
                term_data.append({'y': y_name, 'coeff': coeff, 'xvars': xvars})
                new_terms.append(f"{coeff}*{y_name}")
    new_rhs = rhs - constant_sum
    new_eq_str = "+".join(new_terms) + "=" + str(new_rhs)
    return new_eq_str, term_data, new_rhs

#############################################
# Solve the linear (preprocessed) equation  #
#############################################

def solve_linear_equation(new_eq_str, term_data, rhs):
    """
    The new equation is of the form sum_i (coeff_i * y_i) = rhs.
    For each nonconstant term, we impose the trivial bound:
         0 <= y <= floor(rhs/coeff)
    Then we enumerate over the Cartesian product of these ranges.
    Returns a list of solutions in y-space, each as a tuple.
    """
    bounds = []
    for term in term_data:
        c = term['coeff']
        ub = math.floor(rhs / c)
        bounds.append((0, ub))
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
    Returns all r-tuples of positive integers (a1, ..., ar) such that
          a1 * a2 * ... * ar == n.
    Assumes n > 0 and r >= 1.
    """
    if r == 1:
        return [(n,)]
    solutions = []
    for a in range(1, n+1):
        if n % a == 0:
            for tail in factorize_n(n // a, r-1):
                solutions.append((a,) + tail)
    return solutions

#############################################
# Enforce equivalence constraints among x-vars #
#############################################

def filter_by_equivalences(solution_factorizations, term_data):
    """
    Given a list of factorizations per term (solution_factorizations) and the
    term_data (which contains the ordered x-variable names for each term),
    ensure that if a particular x-variable appears in more than one term, the
    corresponding factor is identical.
    
    Returns a list of combined factorizations (one tuple per overall solution).
    """
    # Build mapping: xvar -> list of (term index, factor index)
    x_positions = {}
    for i, term in enumerate(term_data):
        for j, xvar in enumerate(term['xvars']):
            x_positions.setdefault(xvar, []).append((i, j))
    all_combos = itertools.product(*solution_factorizations)
    valid = []
    for combo in all_combos:
        consistent = True
        for xvar, pos_list in x_positions.items():
            vals = []
            for (ti, fi) in pos_list:
                vals.append(combo[ti][fi])
            if len(set(vals)) > 1:
                consistent = False
                break
        if consistent:
            valid.append(combo)
    return valid

#############################################
# Postprocessing: factorize each y term         #
#############################################

def postprocess_solution(y_solution, term_data):
    """
    For each y_solution (a tuple, one per nonconstant term) and corresponding term_data,
    factorize y into exactly r factors (where r = number of x-vars in that term),
    then filter the overall Cartesian product by equivalence constraints.
    Returns a list of combined factorizations.
    """
    per_term_factorizations = []
    for i, y_val in enumerate(y_solution):
        r = len(term_data[i]['xvars'])
        if y_val <= 0:
            fac = (0,) * r
            per_term_factorizations.append([fac])
        else:
            facs = factorize_n(y_val, r)
            per_term_factorizations.append(facs)
    valid = filter_by_equivalences(per_term_factorizations, term_data)
    return valid

#############################################
# Main driver
#############################################

def main():
    # Use command-line argument if provided; otherwise, prompt the user.
    if len(sys.argv) > 1:
        eq_str = sys.argv[1]
    else:
        eq_str = input("Enter a factored linear Diophantine equation (e.g. '5*x0*x1+2*x1*x2+3=9'): ")
    
    try:
        new_eq_str, term_data, new_rhs = preprocess_equation(eq_str)
    except ValueError as ve:
        print("Error in preprocessing:", ve)
        return
    print("Preprocessed equation:", new_eq_str)
    print("Term data (nonconstant terms only):")
    for term in term_data:
        print(f"  {term['y']}: coefficient {term['coeff']}, factors {term['xvars']}")
    print("Adjusted RHS:", new_rhs)
    
    y_solutions = solve_linear_equation(new_eq_str, term_data, new_rhs)
    if not y_solutions:
        print("No solutions in y-space found.")
        return
    print("Solutions in y-space:")
    for sol in y_solutions:
        print(sol)
    
    overall_factorizations = []
    for y_sol in y_solutions:
        factorizations = postprocess_solution(y_sol, term_data)
        overall_factorizations.append(factorizations)
    
    print("\nFinal factorizations (each item is a tuple with one factorization per nonconstant term):")
    for i, sol_fac in enumerate(overall_factorizations):
        print(f"Solution in y-space {y_solutions[i]} yields factorizations:")
        for combo in sol_fac:
            print(combo)

if __name__ == '__main__':
    main()

