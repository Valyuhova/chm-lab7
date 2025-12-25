import math

def f(x, y):
    return math.cos(y) / (x + 2.0) + 0.3 * y * y

def rk4_step(x, y, h):
    k1 = f(x, y)
    k2 = f(x + h/2.0, y + h*k1/2.0)
    k3 = f(x + h/2.0, y + h*k2/2.0)
    k4 = f(x + h,     y + h*k3)
    return y + (h/6.0) * (k1 + 2*k2 + 2*k3 + k4)

def solve_cauchy_rk4(x0=0.0, y0=0.0, x_end=1.0, h0=0.1, eps=1e-3, h_min=1e-12):
    x = x0
    y = y0
    h = h0

    rows = []

    n = 0
    rows.append((n, x, y, h, 0.0))

    while x < x_end:
        if x + h > x_end:
            h = x_end - x

        while True:
            if h < h_min:
                raise RuntimeError("Крок став надто малим. Перевірте параметри або поведінку рівняння.")

            y_h  = rk4_step(x, y, h)
            y_h2 = rk4_step(x, y, h/2.0)
            y_h2 = rk4_step(x + h/2.0, y_h2, h/2.0)

            err = abs(y_h2 - y_h) / 15.0

            if err <= eps:
                x = x + h
                y = y_h2
                n += 1
                rows.append((n, x, y, h, err))
                break
            else:
                h = h / 2.0

        if err < eps / 16.0:
            h = h * 2.0

    return rows

if __name__ == "__main__":
    x0 = 0.0
    y0 = 0.0
    x_end = 1.0
    h0 = 0.1
    eps = 0.001

    table = solve_cauchy_rk4(x0, y0, x_end, h0, eps)

    print(f"Задача Коші: y' = cos(y)/(x+2) + 0.3*y^2,  y({x0}) = {y0}")
    print(f"Точність: eps = {eps}")
    print(f"Інтервал: [{x0}, {x_end}]")
    print("\n  n        x            y              h           оцінка_похибки")
    print("---------------------------------------------------------------------")
    for n, x, y, h, err in table:
        print(f"{n:3d}  {x:10.6f}  {y:14.10f}  {h:10.6f}  {err:14.10f}")

    n_last, x_last, y_last, h_last, err_last = table[-1]
    print("\nПідсумок:")
    print(f"y({x_last}) ≈ {y_last} (останній крок h={h_last}, оцінка похибки={err_last})")