import numpy


def get_knutepunkt(filename):
    array = les_input(filename)
    return array[0]


def get_elementer(filename):
    array = les_input(filename)
    return array[1]


def get_laster(filename):
    array = les_input(filename)
    return array[2]


def get_Iprofil(int):
    return 1


def get_areal(param):
    return 0.1


def les_input(filename):
    info = []
    with open(filename) as f:
        global n_punkt
        n_punkt = int(f.readline())
        info = f.readlines()
        global n_element
        knutepunkt = []
        elementer = []
        laster = []
        total_array = []
        n_element = 1
        for i in range(0, len(info)):
            arr = list(map(int, info[i].split()))
            if i < n_punkt:
                knutepunkt.append(arr)
            if i > n_punkt and len(arr) > 1:
                elementer.append(arr)
                n_element = n_element + 1
            if i > (n_punkt + n_element):
                laster.append(arr)
        total_array.append(knutepunkt)
        total_array.append(elementer)
        total_array.append(laster)
        return total_array


def sysMatrise_setup():
    global SysMatrise
    SysMatrise = numpy.zeros((3*n_punkt, 3*n_punkt))


def sysMatrise_adder(ki, number, ende1, ende2):
    if number == 0:
        for i in range(3):
            for j in range(3):
                SysMatrise[ende1 + i][ende1 + j] = ki[i][j]
    if number == 1:
        for i in range(3):
            for j in range(3):
                SysMatrise[ende2 + i][ende2 + j] = ki[i][j]
    if number == 2:
        for i in range(3):
            for j in range(3):
                SysMatrise[ende2 + i][ende1 + j] = ki[i][j]
    if number == 3:
        for i in range(3):
            for j in range(3):
                SysMatrise[ende2 + i][ende2 + j] = ki[i][j]


def get_Tmatrise(theta):
    s = numpy.sin(theta)
    c = numpy.cos(theta)
    T_matrise = numpy.array([[c, -s, 0, 0, 0, 0], [s, c, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],
                             [0, 0, 0, c, -s, 0], [0, 0, 0, s, c, 0], [0, 0, 0, 0, 0, 1]])
    return T_matrise


def kMatrise(filename):
    elementer = get_elementer('rammedata.txt')
    knutepunkt = get_knutepunkt('rammedata.txt')
    laster = get_laster('rammedata.txt')

    for i in range(0, len(elementer)):
        ende1 = elementer[i][0]
        ende2 = elementer[i][1]

        x_ende1 = knutepunkt[ende1][0]
        y_ende1 = knutepunkt[ende1][1]
        x_ende2 = knutepunkt[ende2][0]
        y_ende2 = knutepunkt[ende2][1]
        L = numpy.sqrt((x_ende1 - x_ende2) ** 2 + (y_ende1 - y_ende2) ** 2)
        E = elementer[i][2]
        I = get_Iprofil(elementer[i][3])
        A = get_areal(elementer[i][3])
        theta = numpy.arcsin((y_ende2 - y_ende1) / L)
        k = numpy.array([[E * A / L, 0., 0., -E * A / L, 0., 0.],
                         [0., 12 * E * I / (L ** 3), -6 * E * I / (L ** 2), 0., -12 * E * I / (L ** 3),
                          -6 * E * I / (L ** 2)],
                         [0., -6 * E * I / (L ** 2), 4 * E * I / L, 0., 6 * E * I / (L ** 2), 2 * E * I / L],
                         [E * A / L, 0., 0., E * A / L, 0., 0.],
                         [0., -12 * E * I / (L ** 3), 6 * E * I / (L ** 2), 0., 12 * E * I / (L ** 3), 6 * E * I / L],
                         [0., -6 * E * I / (L ** 2), 2 * E * I / L, 0., 6 * E * I / L, 4 * E * I / L]])

        T = get_Tmatrise(theta)
        kg_temp = numpy.dot(T, k)
        T_inv = numpy.linalg.inv(T)
        kg = numpy.dot(kg_temp, T_inv)

        k1 = [[kg[0][0], kg[0][1], kg[0][2]], [kg[1][0], kg[1][1], kg[1][2]], [kg[2][0], kg[2][1], kg[2][2]]]
        k2 = [[kg[0][3], kg[0][4], kg[0][5]], [kg[1][3], kg[1][4], kg[1][5]], [kg[2][3], kg[2][4], kg[2][5]]]
        k3 = [[kg[3][0], kg[3][1], kg[3][2]], [kg[4][0], kg[4][1], kg[4][2]], [kg[5][0], kg[5][1], kg[5][2]]]
        k4 = [[kg[3][3], kg[3][4], kg[3][5]], [kg[4][3], kg[4][4], kg[4][5]], [kg[5][3], kg[5][4], kg[5][5]]]
        k_mini = numpy.array([k1, k2, k3, k4])
        for l in range(4):
            sysMatrise_adder(k_mini[l], l, ende1, ende2)

if __name__ == '__main__':
    les_input('rammedata.txt')
    sysMatrise_setup()
    kMatrise('rammedata.txt')
    print(SysMatrise)
