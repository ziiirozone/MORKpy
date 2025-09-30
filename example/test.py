from MORKpy import *
from matplotlib import pyplot as plt
from math import sqrt,sin,cos

def ENDMORK1_weight_function(_N):
    return [[0.], [1.]]

def ENDMORK1_nodes():
    return [0., 1.]

def ENDMORK1_weight_graph():
    return [[False, False], [True, False]]

ENDMORK1 = NDMORKPy(
    ENDMORK1_weight_function,
    ENDMORK1_nodes(),
    ENDMORK1_weight_graph(),
)

def INDMORK1_weight_function(_N):
    return [[1.], [1.]]

def INDMORK1_nodes():
    return [1., 1.]

def INDMORK1_weight_graph():
    return [[True, False], [True, False]]

INDMORK1 = NDMORKPy(
    INDMORK1_weight_function,
    INDMORK1_nodes(),
    INDMORK1_weight_graph(),
)

def ENDMORK2_weight_function(N):
    return [
        [0., 0.],
        [2**(-N), 0.],
        [1. - 2. / (1. + N), 2. / (1. + N)],
    ]

def ENDMORK2_nodes():
    return [0., 0.5, 1.]

def ENDMORK2_weight_graph():
    return [
        [False, False, False],
        [True, False, False],
        [True, True, False],
    ]

ENDMORK2 = NDMORKPy(
    ENDMORK2_weight_function,
    ENDMORK2_nodes(),
    ENDMORK2_weight_graph(),
)

def INDMORK2_weight_function(N):
    return [[2**(-N)], [1.]]

def INDMORK2_nodes():
    return [0.5, 1.]

def INDMORK2_weight_graph():
    return [[True, False], [True, False]]

INDMORK2 = NDMORKPy(
    INDMORK2_weight_function,
    INDMORK2_nodes(),
    INDMORK2_weight_graph(),
)

def ENDMORK3_weight_function(N):
    return [
        [0., 0., 0.],
        [3**(-N), 0., 0.],
        [
            (2. / 3)**N * (N - 1.) / (1. + N),
            (2. / 3)**N * 2. / (1. + N),
            0.,
        ],
        [
            1. - 9. * N / (2. * (1. + N) * (2. + N)),
            6. * (N - 1.) / ((1. + N) * (2. + N)),
            3. * (4. - N) / (2. * (1. + N) * (2. + N)),
        ],
    ]

def ENDMORK3_nodes():
    return [0., 1. / 3., 2. / 3., 1.]

def ENDMORK3_weight_graph():
    return [
        [False, False, False, False],
        [True, False, False, False],
        [True, True, False, False],
        [True, True, True, False],
    ]

ENDMORK3 = NDMORKPy(
    ENDMORK3_weight_function,
    ENDMORK3_nodes(),
    ENDMORK3_weight_graph(),
)

def INDMORK3_weight_function(_N):
    return [[0., 0.], [0., 1.], [0.5, 0.5]]

def INDMORK3_nodes():
    return [0., 1., 1.]

def INDMORK3_weight_graph():
    return [
        [False, False, False],
        [False, True, False],
        [True, True, False],
    ]

INDMORK3 = NDMORKPy(
    INDMORK3_weight_function,
    INDMORK3_nodes(),
    INDMORK3_weight_graph(),
)

def INDMORK3_1_weight_function(_N):
    return [[0., 0.], [0.5, 0.5], [0.5, 0.5]]

def INDMORK3_1_nodes():
    return [0., 1., 1.]

def INDMORK3_1_weight_graph():
    return [
        [False, False, False],
        [True, True, False],
        [True, True, False],
    ]

INDMORK3_1 = NDMORKPy(
    INDMORK3_1_weight_function,
    INDMORK3_1_nodes(),
    INDMORK3_1_weight_graph(),
)

def ENDMORK4_1_weight_function(N):
    return [
        [0., 0., 0., 0.],
        [2**(-N), 0., 0., 0.],
        [
            2**(-N) * (N - 1.) / (1. + N),
            2**(1 - N) / (1. + N),
            0.,
            0.,
        ],
        [(N - 1.) / (1. + N), (1. - N) / (1. + N), 1., 0.],
        [
            N**(2) / ((1. + N) * (2. + N)),
            2. * N / ((1. + N) * (2. + N)),
            2. * N / ((1. + N) * (2. + N)),
            (2. - N) / ((1. + N) * (2. + N)),
        ],
    ]

def ENDMORK4_1_nodes():
    return [0., 0.5, 0.5, 1., 1.]

def ENDMORK4_1_weight_graph():
    return [
        [False, False, False, False, False],
        [True, False, False, False, False],
        [True, True, False, False, False],
        [True, True, True, False, False],
        [True, True, True, True, False],
    ]

ENDMORK4_1 = NDMORKPy(
    ENDMORK4_1_weight_function,
    ENDMORK4_1_nodes(),
    ENDMORK4_1_weight_graph(),
)

def ENDMORK4_2_weight_function(N):
    return [
        [0., 0., 0., 0.],
        [2**(-N), 0., 0., 0.],
        [
            2**(-N) * N / (1. + N),
            2**(-N) / (1. + N),
            0.,
            0.,
        ],
        [
            (N - 1.) / (1. + N),
            2. * (N - 2.) / (1. + N),
            2. * (3. - N) / (1. + N),
            0.,
        ],
        [
            N**2 / ((1. + N) * (2. + N)),
            0.,
            4. * N / ((1. + N) * (2. + N)),
            (2. - N) / ((1. + N) * (2. + N)),
        ],
    ]

def ENDMORK4_2_nodes():
    return [0., 0.5, 0.5, 1., 1.]

def ENDMORK4_2_weight_graph():
    return [
        [False, False, False, False, False],
        [True, False, False, False, False],
        [True, True, False, False, False],
        [True, True, True, False, False],
        [True, True, True, True, False],
    ]

ENDMORK4_2 = NDMORKPy(
    ENDMORK4_2_weight_function,
    ENDMORK4_2_nodes(),
    ENDMORK4_2_weight_graph(),
)

def INDMORK4_weight_function(N):
    sqrt3 = sqrt(3)
    no1 = 0.5 - sqrt3 / 6.
    no2 = 0.5 + sqrt3 / 6.
    return [
        [
            no1**N / (1. + N) * (1. + N / 2. * (1. + sqrt3)),
            -sqrt3 * N / (1. + N) * no1**(N + 1),
        ],
        [
            sqrt3 * N / (1. + N) * no2**(N + 1),
            no2**N / (1. + N) * (1. + N / 2. * (1. - sqrt3)),
        ],
        [
            0.5 + sqrt3 * (N - 1.) / (2. * (1. + N)),
            0.5 - sqrt3 * (N - 1.) / (2. * (1. + N)),
        ],
    ]

def INDMORK4_nodes():
    return [0.5 - sqrt(3) / 6., 0.5 + sqrt(3) / 6., 1.]

def INDMORK4_weight_graph():
    return [
        [True, True, False],
        [True, True, False],
        [True, True, False],
    ]

INDMORK4 = NDMORKPy(
    INDMORK4_weight_function,
    INDMORK4_nodes(),
    INDMORK4_weight_graph(),
)

def ERK1_weights():
    return [[0.], [1.]]

def ERK1_nodes():
    return [0., 1.]

ERK1 = RKPy(ERK1_weights(), ERK1_nodes())


def IRK1_weights():
    return [[1.], [1.]]

def IRK1_nodes():
    return [1., 1.]

IRK1 = RKPy(IRK1_weights(), IRK1_nodes())

def ERK2_weights():
    return [[0., 0.], [0.5, 0.], [0., 1.]]

def ERK2_nodes():
    return [0., 0.5, 1.]

ERK2 = RKPy(ERK2_weights(), ERK2_nodes())

def IRK2_weights():
    return [[0.5], [1.]]

def IRK2_nodes():
    return [0.5, 1.]

IRK2 = RKPy(IRK2_weights(), IRK2_nodes())

def IRK3_weights():
    return [[0., 0.], [0., 1.], [0.5, 0.5]]

def IRK3_nodes():
    return [0., 1., 1.]

IRK3 = RKPy(IRK3_weights(), IRK3_nodes())

def IRK3_1_weights():
    return [[0., 0.], [0.5, 0.5], [0.5, 0.5]]

def IRK3_1_nodes():
    return [0., 1., 1.]

IRK3_1 = RKPy(IRK3_1_weights(), IRK3_1_nodes())

def ERK3_weights():
    return [
        [0., 0., 0.],
        [1. / 3., 0., 0.],
        [0., 2. / 3., 0.],
        [1. / 4., 0., 3. / 4.],
    ]

def ERK3_nodes():
    return [0., 1. / 3., 2. / 3., 1.]

ERK3 = RKPy(ERK3_weights(), ERK3_nodes())

def ERK4_1_weights():
    return [
        [0., 0., 0., 0.],
        [0.5, 0., 0., 0.],
        [0., 0.5, 0., 0.],
        [0., 0., 1., 0.],
        [1. / 6., 1. / 3., 1. / 3., 1. / 6.],
    ]

def ERK4_1_nodes():
    return [0., 0.5, 0.5, 1., 1.]

ERK4_1 = RKPy(ERK4_1_weights(), ERK4_1_nodes())

def ERK4_2_weights():
    return [
        [0., 0., 0., 0.],
        [0.5, 0., 0., 0.],
        [0.25, 0.25, 0., 0.],
        [0., -1., 2., 0.],
        [1. / 6., 0., 2. / 3., 1. / 6.],
    ]

def ERK4_2_nodes():
    return [0., 0.5, 0.5, 1., 1.]

ERK4_2 = RKPy(ERK4_2_weights(), ERK4_2_nodes())

def IRK4_weights():
    sqrt3 = sqrt(3)
    return [
        [1. / 4., 1. / 4. - sqrt3 / 6.],
        [1. / 4. + sqrt3 / 6., 1. / 4.],
        [0.5, 0.5],
    ]

def IRK4_nodes():
    sqrt3 = sqrt(3)
    return [0.5 - sqrt3 / 6., 0.5 + sqrt3 / 6., 1.]

IRK4 = RKPy(IRK4_weights(), IRK4_nodes())

t0 = 0
y0 = [[-1,-1,1,1]]
f = lambda t,y : [-y[0][1]]
solution = lambda t: [sin(t) - cos(t),- sin(t) - cos(t),- sin(t) + cos(t),sin(t) + cos(t)]

h_left = 0.0001
h_right = 1
N = 20
h = [h_left * (h_right/h_left)**(i/(N-1)) for i in range(N)]
error = []

method = IRK4

for i in range(N):
    y = method.approximate(t0,h[i],f,y0)
    sol = solution(t0+h[i])
    error.append([abs(y[0][0] - sol[0]),abs(y[0][1] - sol[1]),abs(y[0][2] - sol[2]),abs(y[0][3] - sol[3])])
plt.xscale('log')
plt.yscale('log')
plt.plot(h,error)
plt.grid()
plt.show()

'''
h = 1
N = 20
t = [t0 + k * h for k in range(N+1)]
y = [y0 for _ in range(N+1)]
for k in range(N):
    y[k+1] = method.approximate(t[k],h,f,y[k])
extracted_y = [i[0] for i in y]
plt.plot(t,extracted_y)

plt.grid()
plt.show()

method = INDMORK4

t = [t0 + k * h for k in range(N+1)]
y = [y0 for _ in range(N+1)]
for k in range(N):
    y[k+1] = method.approximate(t[k],h,f,y[k])
extracted_y1 = [i[0] for i in y]
plt.plot(t,extracted_y1)

plt.grid()
plt.show()
'''
