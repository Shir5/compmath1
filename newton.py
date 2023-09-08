from util import calculate_derivative


def newton_method(f, df, x0, a, b, tol=1e-6, max_iter=100):
    """
    Newton's method for finding the root of a nonlinear equation.

    :param f: The function for which we want to find the root.
    :param df: The derivative of the function f.
    :param x0: Initial guess for the root.
    :param tol: Tolerance (stop when |f(x)| < tol).
    :param max_iter: Maximum number of iterations.
    :return: List of tuples (approximate root, f(x)) and a flag indicating convergence.
    """
    if (x0 == None):
        if f(a) * calculate_derivative(calculate_derivative((f(a), a)), a) > 0:
            x0 = a
        elif f(b) * calculate_derivative(calculate_derivative((f(b), b)), b) < 0:
            x0 = b
    x = x0
    root_and_values = [(x, f(x))]
    for i in range(max_iter):
        # Check if the derivative is zero (division by zero is not allowed)
        if abs(df(x)) < 1e-12:
            return root_and_values, False

        # Compute the next approximation using Newton's formula
        x_next = x - f(x) / df(x)

        # Check if the solution has converged
        if abs(f(x_next)) < tol:
            root_and_values.append((x_next, f(x_next)))
            return root_and_values, True

        root_and_values.append((x_next, f(x_next)))
        x = x_next

    return root_and_values, False


# Equation: f(x) = x^3 - 2x^2 + 2 = 0
# Derivative: f'(x) = 3x^2 - 4x

def my_function(x):
    return x ** 3 - 2 * x ** 2 + 2


def my_derivative(x):
    return 3 * x ** 2 - 4 * x


# Initial guess

# Find the roots using Newton's method
# root_and_values, converged = newton_method(my_function, my_derivative, initial_guess)

# if converged:
#     print("Approximate roots and corresponding f(x) values (iterations):")
# for i, (root, fx) in enumerate(root_and_values):
#     print(f"Iteration {i + 1}: Root = {root:.3f}, f(x) = {fx:.3f}")
# else:
#     print("Newton's method did not converge after {} iterations.".format(len(root_and_values) - 1))

