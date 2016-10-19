from methods.method import Method


class Jacobi(Method):
    def proceed_rhs(self, prev, win):
        u_upd = [[0.0 for i in range(win.num + 1)] for j in range(win.num + 1)]
        for row in range(win.num + 1):
            for col in range(win.num + 1):
                x = win.dx * col + win.lowLeft[0]
                y = win.dy * row + win.lowLeft[1]
                if row == 0 or col == 0 or row == win.num or col == win.num:
                    u_upd[row][col] = prev[row][col]
                else:
                    val = (prev[row - 1][col] + prev[row + 1][col] + prev[row][col - 1] + prev[row][col + 1] +
                           win.dx * win.dy * win.f(x, y)) / 4.0
                    u_upd[row][col] = val
        return u_upd
