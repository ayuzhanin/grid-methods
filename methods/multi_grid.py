from methods.method import Method


class MultiGrid(Method):
    def rougher(self, grid):
        return [[grid[j][i] for i in range(0, len(grid), 2)]
                for j in range(0, len(grid), 2)]

    def interpolate(self, grid):
        def f(i, j, val):
            if i % 2 == 0 and j % 2 == 0:
                return val
            else:
                return 0.0

        g = [[f(i, j, grid[int(j / 2)][int(i / 2)]) for i in range(len(grid[0]) * 2)] for j in range(len(grid) * 2)]
        for i in range(1, len(g) - 2, 2):
            for j in range(1, len(g[0]) - 2, 2):
                g[i][j] = (g[i - 1][j - 1] + g[i - 1][j + 1] + g[i + 1][j - 1] + g[i + 1][j + 1]) / 4.0

        for i in range(1, len(g) - 2, 1):
            for j in range(1, len(g[0]) - 2, 1):
                if (j % 2 == 0 and i % 2 == 1) or (i % 2 == 0 and j % 2 == 1):
                    g[i][j] = (g[i - 1][j] + g[i + 1][j] + g[i][j - 1] + g[i][j + 1]) / 4.0

        for i in range(1, len(g) - 2):
            if not i % 2 == 0:
                g[i][0] = (g[i - 1][0] + g[i][1] + g[i + 1][0]) / 3.0
                g[i][len(g) - 2] = (g[i - 1][len(g) - 2] + g[i][len(g) - 3] + g[i + 1][len(g) - 2]) / 3.0

        for i in range(1, len(g[0]) - 2):
            if not i % 2 == 0:
                g[0][i] = (g[0][i - 1] + g[1][i] + g[0][i + 1]) / 3.0
                g[len(g[0]) - 2][i] = (g[len(g[0]) - 2][i - 1] + g[len(g[0]) - 3][i] + g[len(g[0]) - 2][i + 1]) / 3.0

        for i in range(0, len(g)):
            g[i][len(g) - 1] = (g[i][len(g) - 2] + g[i][len(g) - 3]) / 2.0

        for i in range(0, len(g[0])):
            g[len(g[0]) - 1][i] = (g[len(g[0]) - 2][i] + g[len(g[0]) - 3][i]) / 2.0

        g[0][0] = (g[1][0] + g[0][1] + g[1][1]) / 3.0

        return g

    def proceed_grid(self, rhs_grid, win):
        """Computing discrepancy (updated) for L_h nu_h = r_h equation

        Initializing grid (updated) with zeros
        And then obtaining discrepancy as updated = - L_h nu_h - r_h
        Использовать модуль для r_2h или нет?
        """
        updated = [[0.0 for i in range(win.num + 1)] for j in range(win.num + 1)]
        for row in range(win.num + 1):
            for col in range(win.num + 1):
                x = win.dx * col + win.lowLeft[0]
                y = win.dy * row + win.lowLeft[1]
                if row == 0:
                    updated[row][col] = win.low(x, y)
                elif row == win.num:
                    updated[row][col] = win.up(x, y)
                elif col == 0:
                    updated[row][col] = win.left(x, y)
                elif col == win.num:
                    updated[row][col] = win.right(x, y)
                else:
                    val = (updated[row - 1][col] + updated[row][col - 1] + win.dx * win.dy * rhs_grid[row][col]) / 4.0
                    updated[row][col] = val
        return updated

    def rr(self, nu, r, win):
        """Computing discrepancy (updated) for L_h nu_h = r_h equation

        Initializing grid (updated) with zeros
        And then obtaining discrepancy as updated = - L_h nu_h - r_h
        Использовать модуль для r_2h или нет?
        """
        updated = [[0.0 for i in range(win.num + 1)] for j in range(win.num + 1)]
        for row in range(win.num + 1):
            for col in range(win.num + 1):
                updated[row][col] = - ((nu[row - 1][col] - 2 * nu[row][col] + nu[row + 1][col] +
                                        nu[row][col - 1] - 2 * nu[row][col] + nu[row][col + 1])
                                       / (win.dx * win.dy)) - r[row][col]
        return updated

    def sum(self, a, b):
        """Computing sum of two grids (updated)

        Initializing grid (updated) with zeros
        And then obtaining sum as updated = a + b
        """
        updated = [[0.0 for i in range(len(a))] for j in range(len(a[0]))]
        for row in range(len(a)):
            for col in range(len(a[0])):
                updated[row][col] = a[row][col] + b[row][col]
        return updated

    def solve(self, win, eps):
        """"Solving equation L_h u_h = f_h

         Initializing grid: fill grid with zeros and values at borders using conditions and self.init_grid method
         And then perform one V loop iteration using self.iterate
        """

        win_h = win

        prev = self.init_grid(win)
        cur = self.iterarate(prev, win)
        dis_prev = self.discrepancy(prev, win)
        dis_cur = self.discrepancy(cur, win)
        diff = abs(dis_cur - dis_prev)
        counter = 1
        while diff > eps:
            prev = cur
            cur = self.iterarate(prev, win)
            dis_prev = self.discrepancy(prev, win)
            dis_cur = self.discrepancy(cur, win)
            diff = abs(dis_prev - dis_cur)
            counter += 1
        return cur, counter

    def iterarate(self, grid, win):
        win_h = win
        u = self.proceed_rhs(grid, win_h)

        # Computing discrepancy (r_h) for L_h u_h = f_h equation
        # Initializing grid (r_h) with zeros
        # And then obtaining discrepancy as r_h = - L_h u_h - f_h
        # Использовать модуль для r_h или нет?
        r_h = [[0.0 for i in range(win_h.num + 1)] for j in range(win_h.num + 1)]
        for row in range(win_h.num + 1):
            for col in range(win_h.num + 1):
                x = win_h.dx * col + win_h.lowLeft[0]
                y = win_h.dy * row + win_h.lowLeft[1]
                r_h[u][row] = - ((u[row - 1][col] - 2 * u[row][col] + u[row + 1][col] +
                                  u[row][col - 1] - 2 * u[row][col] + u[row][col + 1])
                                 / (win.dx * win.dy)) - win.f(x, y)

        # Create a copy of window at previous step
        # Halving number of intervals, and updating intervals' lengths
        # Then roughing grid obtained at previous step by taking one of every two nearly located points
        # Using self.rougher for that
        # Proceeding Gauss-Seidel iteration for the case when in instead of function
        # grid represented at the right hand side of equation instead. Using self.proceed_grid for that
        import copy
        win_2h = copy.deepcopy(win_h)
        win_2h.num /= 2
        win_2h.update()
        r_2h = self.rougher(r_h)
        nu_2h = self.proceed_grid(r_2h, win_2h)

        # Computing discrepancy (rr_2h) for L_2h nu_2h = r_2h equation
        rr_2h = self.rr(nu_2h, r_2h, win_2h)

        # Create a copy of window at previous step
        # Halving number of intervals, and updating intervals' lengths
        # Then roughing grid obtained at previous step by taking one of every two nearly located points
        # Using self.rougher for that
        # Proceeding Gauss-Seidel iteration for the case when in instead of function
        # grid represented at the right hand side of equation instead. Using self.proceed_grid for that
        win_4h = copy.deepcopy(win_2h)
        win_4h.num /= 2
        win_4h.update()
        r_4h = self.rougher(rr_2h)
        nu_4h = self.proceed_grid(r_4h, win_4h)

        # Computing discrepancy (rr_4h) for L_4h nu_4h = r_4h equation
        rr_4h = self.rr(nu_4h, r_4h, win_4h)

        # Повторили то же самое для еще более грубой сетки
        win_8h = copy.deepcopy(win_4h)
        win_8h.num /= 2
        win_8h.update()
        r_8h = self.rougher(rr_4h)
        nu_8h = self.proceed_grid(r_8h, win_8h)

        # Теперь нужно двигаться в сторону уменьшения сетки (т. е. в обратную сторону)

        nu_4h = self.sum(nu_4h, self.interpolate(nu_8h))
        method = GaussSeidel()
        nu_4h = method.proceed_rhs(nu_4h, win_4h)

        nu_2h = self.sum(nu_2h, self.interpolate(nu_4h))
        nu_2h = method.proceed_rhs(nu_2h, win_2h)

        u = self.sum(u, self.interpolate(nu_2h))

        return u


