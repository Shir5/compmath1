def bisection_method(f, a, b, eps):
    """
        Данная функция реализует метод бисекции для нахождения корня уравнения func(x) = 0 на интервале [a, b] с заданной точностью acc.

        Аргументы функции:
            func: функция, корень которой необходимо найти.
            a: левая граница интервала.
            b: правая граница интервала.
            acc: заданная точность.

        Функция возвращает найденное значение корня уравнения.
    """
    count = 0
    while (b - a) / 2 > eps:
        c = (a + b) / 2
        if f(c) == 0:
            return c
        elif f(c) * f(a) < 0:
            b = c
        else:
            a = c
        count += 1
    print(f"Количество итераций: {count}")
    print(f"f(x) = {f((a + b) / 2)}")
    return (a + b) / 2