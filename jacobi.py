from methods.method import Method


class Jacobi(Method):
    def proceed_rhs(self, u_grid, win):
        updated = [[0.0 for i in range(win.num + 1)] for j in range(win.num + 1)]
        for row in range(win.num + 1):
            for col in range(win.num + 1):
                x = win.dx * col + win.lowLeft[0]
                y = win.dy * row + win.lowLeft[1]
                if row == 0 or col == 0 or row == win.num or col == win.num:
                    updated[row][col] = u_grid[row][col]
                else:
                    val = (u_grid[row - 1][col] + u_grid[row + 1][col] + u_grid[row][col - 1] + u_grid[row][col + 1] +
                           win.dx * win.dy * win.f(x, y)) / 4.0
                    updated[row][col] = val
        return updated
