import sympy as sp

from util import count_accuracy


def secant_method(f, x0, x1, eps):
    """
        Аргументы функции:

            f: функция, корень которой необходимо найти.
            x0: начальная точка.
            x1: начальная точка.
            eps: заданная точность.

        Функция возвращает найденное значение корня уравнения.
        """
    count = 0
    while abs(f(x1)) > eps:
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        x0 = x1
        x1 = x2
        count += 1
    print(f"Количество итераций: {count}")
    print(f"f(x) = {f(x1)}")
    return x1