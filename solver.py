#!/usr/bin/env python3
import sys
import re
import math
import itertools
import numpy as np
from sympy import symbols, diophantine, Eq, lambdify
from scipy.optimize import linprog

#############################################
# Helper: parse an x-variable specification #
#############################################
def parse_xvar_spec(spec):
    """
    Parse an x-variable specification string.

    Valid formats:
      "x0"         -> (0, sys.maxsize)
      "x0:3"       -> (3, 3)
      "x0:3,"      -> (3, sys.maxsize)
      "x0:,5"      -> (0, 5)
      "x0:3,5"     -> (3, 5)
    """
    global_max = sys.maxsize
    if ':' not in spec:
        # No colon => no bounds => (0, max)
        return (0, global_max)
    # Split on colon
    parts = spec.split(':', 1)
    bounds_spec = parts[1]
    if ',' not in bounds_spec:
        # Single integer => exact
        try:
            val = int(bounds_spec)
            return (val, val)
        except:
            return (0, global_max)
    else:
        lower_str, upper_str = bounds_spec.split(',', 1)
        lower = int(lower_str) if lower_str else 0
        upper = int(upper_str) if upper_str else global_max
        return (lower, upper)

#############################################
# Preprocessing: parse the factored equation #
#############################################
def preprocess_equation(eq_str):
    """
    E.g. "5*x0*x1+2*x1:x2+3=9" -> parse LHS into terms, subtract constant from RHS,
    record nonconstant terms in term_data with 'y', 'coeff', 'xvars'.
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
        if len(factors) == 1 and re.fullmatch(r'\d+', factors[0]):
            # pure constant
            constant_sum += int(factors[0])
        else:
            # maybe leading numeric coefficient
            if re.fullmatch(r'\d+', factors[0]):
                coeff = int(factors[0])
                xvars = factors[1:]
            else:
                coeff = 1
                xvars = factors
            if not xvars:
                # no xvars => pure constant
                constant_sum += coeff
            else:
                y_name = f"y{len(term_data)}"
                term_data.append({'y': y_name, 'coeff': coeff, 'xvars': xvars})
                new_terms.append(f"{coeff}*{y_name}")
    new_rhs = rhs - constant_sum
    new_eq_str = "+".join(new_terms) + "=" + str(new_rhs)
    return new_eq_str, term_data, new_rhs

#############################################
# Helper: build a "human-friendly" constraint
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
    """
    We build an equation sum(coeff_i * y_i) = rhs, call diophantine,
    parse the parametric solution, and then set trivial bounds for each y_i
    based on xvar specs. We solve param bounds with linprog, then enumerate.
    """
    # Create sympy symbols
    from sympy import Eq, diophantine, symbols, lambdify
    y_vars = [symbols(term['y'], integer=True) for term in term_data]
    eq_sym = Eq(sum(t['coeff'] * v for t,v in zip(term_data, y_vars)), rhs)
    
    sol_set = diophantine(eq_sym)
    if not sol_set:
        return []
    sol = next(iter(sol_set))
    if not isinstance(sol, (tuple, list)):
        sol = (sol,)
    
    print("General solution for y variables:")
    for i, expr in enumerate(sol):
        print(f"  {y_vars[i]} = {expr}")
    
    # Identify free parameters
    free_params = set()
    for expr in sol:
        free_params |= expr.free_symbols
    free_params = list(free_params)
    print("Free parameters:", [str(p) for p in free_params])
    
    # Build lambdas
    f_funcs = [lambdify(free_params, expr, "math") for expr in sol]
    
    # For each y_i, compute trivial bounds
    # from xvar specs: product of (lb..ub)
    # then intersect with [0..floor(rhs/coeff)].
    trivial_bounds = []
    for term in term_data:
        coeff = term['coeff']
        xvars = term['xvars']
        prod_lower = 1
        prod_upper = 1
        debraced_vars = []
        for xspec in xvars:
            if xspec.startswith("{") and xspec.endswith("}"):
                xspec = xspec[1:-1]  # remove outer braces
            debraced_vars.append(xspec)
        for xspec in debraced_vars:
            lb, ub = parse_xvar_spec(xspec)
            prod_lower *= lb
            prod_upper *= ub
        # y_i's lower bound = max(0, coeff * prod_lower)
        # y_i's upper bound = min(floor(rhs/coeff), coeff * prod_upper)
        low_candidate = coeff * prod_lower
        high_candidate = coeff * prod_upper
        global_upper = math.floor(rhs / coeff)
        L_i = max(0, low_candidate)
        U_i = min(global_upper, high_candidate)
        trivial_bounds.append((L_i, U_i))
    
    print("\nConstraints for each y_i (using xvar specs):")
    # Build the coefficient list
    from sympy import Add
    sol_coeffs = []
    for expr in sol:
        d = expr.as_coefficients_dict()
        c_val = float(d.get(1, 0))
        a_list = []
        for p in free_params:
            a_list.append(float(d.get(p, 0)))
        sol_coeffs.append((a_list, c_val))
    
    for i, (a_list, c_val) in enumerate(sol_coeffs):
        L_i, U_i = trivial_bounds[i]
        low_str = constraint_to_string(a_list, c_val, free_params, L_i, True)
        up_str = constraint_to_string(a_list, c_val, free_params, U_i, False)
        print(f"  y{i} => {term_data[i]['xvars']} : {term_data[i]['coeff']}")
        print(f"     Lower: {low_str}")
        print(f"     Upper: {up_str}")
    
    # Build A_ub, b_ub
    A_ub = []
    b_ub = []
    for i, (a_list, c_val) in enumerate(sol_coeffs):
        L_i, U_i = trivial_bounds[i]
        # y_i >= L_i => -a·t <= c_val - L_i
        A_ub.append([-a for a in a_list])
        b_ub.append(c_val - L_i)
        # y_i <= U_i => a·t <= U_i - c_val
        A_ub.append(a_list)
        b_ub.append(U_i - c_val)
    
    A_ub = np.array(A_ub, dtype=float)
    b_ub = np.array(b_ub, dtype=float)
    
    # Solve LP for each free param with no fallback bounding
    free_bounds_lp = [(None, None)] * len(free_params)
    lower_free = []
    upper_free = []
    
    for j in range(len(free_params)):
        c_obj = [0]*len(free_params)
        c_obj[j] = 1
        res_min = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, bounds=free_bounds_lp, method="highs")
        if not res_min.success:
            print(f"No feasible minimum bound for {free_params[j]}.")
            return []
        min_val = math.ceil(res_min.fun)
        
        c_obj[j] = -1
        res_max = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, bounds=free_bounds_lp, method="highs")
        if not res_max.success:
            print(f"No feasible maximum bound for {free_params[j]}.")
            return []
        max_val = math.floor(-res_max.fun)
        
        if min_val > max_val:
            print(f"Contradictory bounds for {free_params[j]}: [{min_val}, {max_val}]. No solutions.")
            return []
        
        lower_free.append(min_val)
        upper_free.append(max_val)
    
    print("\nFinal bounds for free parameters:")
    for p, lo, hi in zip(free_params, lower_free, upper_free):
        print(f"  {p} in [{lo}, {hi}]")
    
    # Enumerate integer values
    free_ranges = [range(lo, hi+1) for lo, hi in zip(lower_free, upper_free)]
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
# factorize_n, filter_by_equivalences, etc.
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
            fac = (0,) * r
            per_term_factorizations.append([fac])
        else:
            facs = factorize_n(y_val, r)
            per_term_factorizations.append(facs)
    valid = filter_by_equivalences(per_term_factorizations, term_data)
    return valid

def solutions_as_dicts(overall_factorizations, term_data):
    """
    Convert each combined factorization (one factorization tuple per nonconstant term)
    into a dictionary mapping just the base x-variable name (e.g. 'x0') to its integer value.

    For example, if xvar is 'x0:1,' or 'x0:2,5', we strip off everything after the colon
    and store only 'x0' as the dictionary key.
    """
    sol_dicts = []
    for combo in overall_factorizations:
        d = {}
        valid = True
        for i, factorization in enumerate(combo):
            for j, xvar in enumerate(term_data[i]['xvars']):
                # Extract the base variable name (e.g. 'x0') by ignoring any suffix after ':'
                colon_pos = xvar.find(':')
                if xvar.startswith("{") and xvar.endswith("}"):
                    xvar = xvar[1:-1]  # remove the outer braces&
                base_var = xvar[:colon_pos - 1] if colon_pos != -1 else xvar

                if base_var in d:
                    # If we've already assigned a value to this base_var, check consistency
                    if d[base_var] != factorization[j]:
                        valid = False
                        break
                else:
                    d[base_var] = factorization[j]
            if not valid:
                break
        if valid:
            sol_dicts.append(d)
    return sol_dicts

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
    for d in sol_dicts:
        print(d)

if __name__ == "__main__":
    main()

