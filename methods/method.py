class Window:
    def __init__(self, low_left, up_right,
                 num, f,
                 up, low, left, right):
        self.lowLeft = low_left
        self.upRight = up_right
        self.num = num
        self.f = f
        self.up = up
        self.low = low
        self.left = left
        self.right = right
        self.dx = (self.upRight[1] - self.lowLeft[1]) / self.num
        self.dy = (self.upRight[0] - self.lowLeft[0]) / self.num

    def update(self):
        self.dx = (self.upRight[1] - self.lowLeft[1]) / self.num
        self.dy = (self.upRight[0] - self.lowLeft[0]) / self.num


def print_grid(grid):
    for row in range(len(grid)):
        for col in range(len(grid[0])):
            print(grid[row][col], end=" ")
        print(" ")


class Method:
    def init_grid(self, win):
        grid = [[0.0 for i in range(win.num + 1)] for j in range(win.num + 1)]
        for row in range(win.num + 1):
            for col in range(win.num + 1):
                x = win.dx * col + win.lowLeft[0]
                y = win.dy * row + win.lowLeft[1]
                if row == 0:
                    grid[row][col] = win.low(x, y)
                elif row == win.num:
                    grid[row][col] = win.up(x, y)
                elif col == 0:
                    grid[row][col] = win.left(x, y)
                elif col == win.num:
                    grid[row][col] = win.right(x, y)
                else:
                    grid[row][col] = 0
        return grid

    def discrepancy(self, grid, win):
        import sys
        maxima = sys.float_info.min
        for row in (range(len(grid))[1:-1]):
            for col in (range(len(grid))[1:-1]):
                x = win.dx * col + win.lowLeft[0]
                y = win.dy * row + win.lowLeft[1]
                lh_uh = (grid[row - 1][col] - 2 * grid[row][col] + grid[row + 1][col]) / (win.dx * win.dy) + \
                        (grid[row][col - 1] - 2 * grid[row][col] + grid[row][col + 1]) / (win.dx * win.dy)
                dis = abs(- lh_uh - win.f(x, y))
                maxima = max(dis, maxima)
        return maxima

    def proceed_rhs(self, u_grid, win):
        return [[0.0 for i in range(win.num + 1)] for j in range(win.num + 1)]

    def solve(self, win, eps):
        prev = self.init_grid(win)
        cur = self.proceed_rhs(prev, win)
        dis_prev = self.discrepancy(prev, win)
        dis_cur = self.discrepancy(cur, win)
        diff = abs(dis_cur - dis_prev)
        counter = 1
        while diff > eps:
            prev = cur
            cur = self.proceed_rhs(prev, win)
            dis_prev = self.discrepancy(prev, win)
            dis_cur = self.discrepancy(cur, win)
            diff = abs(dis_prev - dis_cur)
            counter += 1
        return cur, counter
