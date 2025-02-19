#!/usr/bin/env python3
import sys
import re
import math
import itertools
import numpy as np
from sympy import symbols, diophantine, Eq, sympify, lambdify
from scipy.optimize import linprog

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
        factors = term.split('*')
        # If the term is a constant (only one factor and numeric), add it.
        if len(factors) == 1 and re.fullmatch(r'\d+', factors[0]):
            constant_sum += int(factors[0])
        else:
            if re.fullmatch(r'\d+', factors[0]):
                coeff = int(factors[0])
                xvars = factors[1:]
            else:
                coeff = 1
                xvars = factors
            if not xvars:
                constant_sum += coeff
            else:
                y_name = f"y{len(term_data)}"
                term_data.append({'y': y_name, 'coeff': coeff, 'xvars': xvars})
                new_terms.append(f"{coeff}*{y_name}")
    new_rhs = rhs - constant_sum
    new_eq_str = "+".join(new_terms) + "=" + str(new_rhs)
    return new_eq_str, term_data, new_rhs

#############################################
# Helper to build a "human-friendly" constraint string
#############################################

def constraint_to_string(a_list, c_val, free_params, bound_val, is_lower=True):
    sign = ">=" if is_lower else "<="
    parts = []
    for coeff, p in zip(a_list, free_params):
        if abs(coeff) < 1e-9:
            continue
        if coeff == 1:
            parts.append(f"{p}")
        elif coeff == -1:
            parts.append(f"-{p}")
        else:
            parts.append(f"{coeff:g}*{p}")
    left_side = " + ".join(parts) if parts else "0"
    if abs(c_val) > 1e-9:
        if left_side == "0":
            left_side = f"{c_val:g}"
        else:
            if c_val > 0:
                left_side += f" + {c_val:g}"
            else:
                left_side += f" - {abs(c_val):g}"
    else:
        if not parts:
            left_side = "0"
    return f"{left_side} {sign} {bound_val}"

#############################################
# Solve the linear equation in terms of free parameters
#############################################

def solve_linear_equation(new_eq_str, term_data, rhs):
    y_vars = [symbols(term['y'], integer=True) for term in term_data]
    eq_sym = Eq(sum(t['coeff']*v for t,v in zip(term_data,y_vars)), rhs)
    
    sol_set = diophantine(eq_sym)
    if not sol_set:
        return []
    sol = next(iter(sol_set))
    if not isinstance(sol, (tuple,list)):
        sol = (sol,)
    
    print("General solution for y variables:")
    for i, expr in enumerate(sol):
        print(f"  {y_vars[i]} = {expr}")
    
    # Identify free params
    free_params = set()
    for expr in sol:
        free_params |= expr.free_symbols
    free_params = list(free_params)
    print("Free parameters:", [str(p) for p in free_params])
    
    # Build lambdas
    f_funcs = [lambdify(free_params, expr, "math") for expr in sol]
    
    # trivial bounds
    trivial_bounds = []
    for term in term_data:
        trivial_bounds.append((0, math.floor(rhs / term['coeff'])))
    
    # Build constraints
    coeff_constraints = []
    for expr in sol:
        d = expr.as_coefficients_dict()
        c_val = d.get(1, 0)
        a_list = [float(d.get(p, 0)) for p in free_params]
        coeff_constraints.append((a_list, float(c_val)))
    
    print("\nConstraints for each y_i:")
    for i, (a_list, c_val) in enumerate(coeff_constraints):
        L_i, U_i = trivial_bounds[i]
        cstr_lower = constraint_to_string(a_list, c_val, free_params, L_i, True)
        cstr_upper = constraint_to_string(a_list, c_val, free_params, U_i, False)
        print(f"  y{i} >= {L_i} => {cstr_lower}")
        print(f"  y{i} <= {U_i} => {cstr_upper}")
    
    A_ub = []
    b_ub = []
    for i, (a, c_val) in enumerate(coeff_constraints):
        L_i, U_i = trivial_bounds[i]
        # lower: y_i >= L_i => a·t + c >= L => -a·t <= c - L
        A_ub.append([-x for x in a])
        b_ub.append(c_val - L_i)
        # upper: y_i <= U_i => a·t + c <= U => a·t <= U - c
        A_ub.append(a)
        b_ub.append(U_i - c_val)
    
    A_ub = np.array(A_ub)
    b_ub = np.array(b_ub)
    
    # Solve LP for each free param with no fallback bounding
    free_bounds = [(None, None)] * len(free_params)
    lower_free = []
    upper_free = []
    
    for j in range(len(free_params)):
        c_obj = [0]*len(free_params)
        c_obj[j] = 1
        res_min = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, bounds=free_bounds, method="highs")
        if not res_min.success or res_min.status in (3,4):
            print(f"No feasible/min bound for {free_params[j]} => region unbounded or infeasible.")
            return []
        min_val = math.ceil(res_min.fun)
        
        c_obj[j] = -1
        res_max = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, bounds=free_bounds, method="highs")
        if not res_max.success or res_max.status in (3,4):
            print(f"No feasible/max bound for {free_params[j]} => region unbounded or infeasible.")
            return []
        max_val = math.floor(-res_max.fun)
        
        if min_val > max_val:
            print(f"Contradictory bounds for {free_params[j]}: [{min_val}, {max_val}]. No solutions.")
            return []
        
        lower_free.append(min_val)
        upper_free.append(max_val)
    
    print("\nBounds for free parameters (inferred from constraints):")
    for p, lo, hi in zip(free_params, lower_free, upper_free):
        print(f"  {p} in [{lo}, {hi}]")
    
    # Enumerate
    free_ranges = [range(lo, hi+1) for lo, hi in zip(lower_free,upper_free)]
    candidate_y_solutions = []
    for vals in itertools.product(*free_ranges):
        candidate = tuple(int(round(f(*vals))) for f in f_funcs)
        valid = True
        for i, y_val in enumerate(candidate):
            L_i, U_i = trivial_bounds[i]
            if y_val < L_i or y_val > U_i:
                valid = False
                break
        if valid:
            candidate_y_solutions.append(candidate)
    
    return candidate_y_solutions

#############################################
# Factorization, filtering, etc.
#############################################

def factorize_n(n, r):
    if r == 1:
        return [(n,)]
    solutions = []
    for a in range(1, n+1):
        if n % a == 0:
            for tail in factorize_n(n // a, r-1):
                solutions.append((a,) + tail)
    return solutions

def filter_by_equivalences(solution_factorizations, term_data):
    x_positions = {}
    for i, term in enumerate(term_data):
        for j, xvar in enumerate(term['xvars']):
            x_positions.setdefault(xvar, []).append((i, j))
    all_combos = itertools.product(*solution_factorizations)
    valid = []
    for combo in all_combos:
        consistent = True
        for xvar, pos_list in x_positions.items():
            vals = [combo[ti][fi] for (ti, fi) in pos_list]
            if len(set(vals)) > 1:
                consistent = False
                break
        if consistent:
            valid.append(combo)
    return valid

def postprocess_solution(y_solution, term_data):
    per_term_factorizations = []
    for i, y_val in enumerate(y_solution):
        r = len(term_data[i]['xvars'])
        if y_val <= 0:
            fac = (0,)*r
            per_term_factorizations.append([fac])
        else:
            facs = factorize_n(y_val, r)
            per_term_factorizations.append(facs)
    return filter_by_equivalences(per_term_factorizations, term_data)

def solutions_as_dicts(overall_factorizations, term_data):
    sol_dicts = []
    for combo in overall_factorizations:
        d = {}
        valid = True
        for i, factorization in enumerate(combo):
            for j, xvar in enumerate(term_data[i]['xvars']):
                if xvar in d:
                    if d[xvar] != factorization[j]:
                        valid = False
                        break
                else:
                    d[xvar] = factorization[j]
            if not valid:
                break
        if valid:
            sol_dicts.append(d)
    return sol_dicts

#############################################
# Parse the original eq for verification
#############################################

def parse_original_equation(eq_str):
    """
    Parse the original eq_str (like '4*x0+1*x1*x2=20') into a Sympy Eq object
    using the x-variables found. We'll use the same approach as a minimal parse.
    """
    from sympy import symbols, Eq, sympify
    eq_str_nospace = eq_str.replace(" ", "")
    parts = eq_str_nospace.split("=")
    if len(parts) != 2:
        raise ValueError("Original eq must have '='.")
    lhs_str, rhs_str = parts
    # find x-variables
    var_names = sorted(set(re.findall(r'x\d+', lhs_str)), key=lambda s: int(s[1:]))
    var_dict = {name: symbols(name, integer=True) for name in var_names}
    lhs_expr = sympify(lhs_str, locals=var_dict)
    rhs_expr = sympify(rhs_str, locals=var_dict)
    eq_orig = Eq(lhs_expr, rhs_expr)
    return eq_orig, var_dict

#############################################
# Main driver
#############################################

def main():
    if len(sys.argv) > 1:
        eq_str = sys.argv[1]
    else:
        eq_str = input("Enter a factored linear Diophantine equation: ")
    
    try:
        new_eq_str, term_data, new_rhs = preprocess_equation(eq_str)
    except ValueError as e:
        print("Error in preprocessing:", e)
        return
    
    print("Preprocessed equation:", new_eq_str)
    print("Term data (nonconstant terms only):")
    for t in term_data:
        print(f"  {t['y']}: coeff {t['coeff']}, xvars {t['xvars']}")
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
        combos = postprocess_solution(y_sol, term_data)
        overall_factorizations.extend(combos)
    
    sol_dicts = solutions_as_dicts(overall_factorizations, term_data)
    print("\nFinal candidate solutions (dict mapping x-variable to value):")
    
    # Now parse the *original* eq_str for verification:
    try:
        eq_orig, var_dict = parse_original_equation(eq_str)
    except ValueError:
        eq_orig = None
        var_dict = {}
        print("(Couldn't parse original eq for verification.)")
    
    for d in sol_dicts:
        # Print the dictionary
        print(d, end='')
        # If we have eq_orig, let's verify
        if eq_orig is not None and var_dict:
            # Build a subs map from x0->..., x1->... in 'd' 
            # but eq_orig might have a subset of those variables
            subs_map = {}
            for k,v in d.items():
                if k in var_dict:
                    subs_map[var_dict[k]] = v
            lhs_val = eq_orig.lhs.subs(subs_map)
            rhs_val = eq_orig.rhs.subs(subs_map)
            if lhs_val == rhs_val:
                print(" (success)")
            else:
                print(" (error)")
        else:
            print("")

if __name__ == "__main__":
    main()

