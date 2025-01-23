from sympy import symbols, sympify

# Helper function to evaluate expressions
def parse_expression(base_vars ,expr):
    if isinstance(expr, list):
        return [parse_expression(e) for e in expr]
    elif isinstance(expr, str):
        return float(sympify(expr, base_vars))
            
