import sys

from methods.gauss_seidel import GaussSeidel
from methods.jacobi import Jacobi

from methods.method import Window
from methods.method import print_grid

def main():
    window = Window((0.0, 0.0), (1.0, 1.0), 32,
                    lambda x, y: (6 * (x ** 4) * y + 12 * (x ** 2) * (y ** 3)),
                    lambda x, y: (x ** 4) * (y ** 3),
                    lambda x, y: (x ** 4) * (y ** 3),
                    lambda x, y: (x ** 4) * (y ** 3),
                    lambda x, y: (x ** 4) * (y ** 3))

    eps = 0.01

    # method = Jacobi()
    # result = method.solve(window, eps)
    #
    # grid = result[0]
    # iters = result[1]
    #
    # print_grid(grid)
    # print(iters)

    # method = GaussSeidel()
    # result = method.solve(window, eps)
    #
    # grid = result[0]
    # iters = result[1]
    #
    # print_grid(grid)
    # print(iters)


if __name__ == "__main__":
    sys.exit(main())
