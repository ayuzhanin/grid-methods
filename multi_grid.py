from methods.gauss_seidel import GaussSeidel
import copy


class MultiGrid(GaussSeidel):
    def rougher(self, some_grid):
        return [[some_grid[j][i] for i in range(0, len(some_grid), 2)]
                for j in range(0, len(some_grid), 2)]

    def proceed_grid(self, rhs_grid, win):
        update = [[0.0 for i in range(win.num + 1)] for j in range(win.num + 1)]
        for row in range(win.num + 1):
            for col in range(win.num + 1):
                x = win.dx * col + win.lowLeft[0]
                y = win.dy * row + win.lowLeft[1]
                if row == 0:
                    update[row][col] = win.low(x, y)
                elif row == win.num:
                    update[row][col] = win.up(x, y)
                elif col == 0:
                    update[row][col] = win.left(x, y)
                elif col == win.num:
                    update[row][col] = win.right(x, y)
                else:
                    val = (update[row - 1][col] + update[row][col - 1] - win.dx * win.dy * rhs_grid[row][col]) / 4.0
                    update[row][col] = val
        return update

    def solve(self, win, eps):
        win_h = win
        u = self.init_grid(win_h)
        u = self.proceed_rhs(u, win_h)
        r_h = [[0.0 for i in range(win_h.num + 1)] for j in range(win_h.num + 1)]

        for row in range(win_h.num + 1):
            for col in range(win_h.num + 1):
                x = win_h.dx * col + win_h.lowLeft[0]
                y = win_h.dy * row + win_h.lowLeft[1]
                r_h[u][row] = - win_h.f(x, y) - u[row][col]

        win_2h = copy.deepcopy(win_h)
        win_2h.dx *= 2
        win_2h.dy *= 2
        win_2h.num /= 2
        r_2h = self.rougher(r_h)
        nu_2h = self.proceed_grid(r_2h, win_2h)

        rr_2h = [[0.0 for i in range(win_2h.num + 1)] for j in range(win_2h.num + 1)]
        for row in range(win_2h.num + 1):
            for col in range(win_2h.num + 1):
                rr_2h[row][col] = r_2h[row][col] - nu_2h[row][col]

        win_4h = copy.deepcopy(win_2h)
        win_4h.dx *= 2
        win_4h.dy *= 2
        win_4h.num /= 2
        r_4h = self.rougher(r_2h)
        nu_4h = self.proceed_grid(r_4h, win_4h)

        rr_4h = [[0.0 for i in range(win_4h.num + 1)] for j in range(win_4h.num + 1)]
        for row in range(win_4h.num + 1):
            for col in range(win_4h.num + 1):
                rr_4h[row][col] = r_4h[row][col] - nu_4h[row][col]
