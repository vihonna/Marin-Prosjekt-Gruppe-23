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
    return 1000000

def get_areal(param):
    return 10000


def les_input(filename):
    info = []
    with open(filename) as f:
        n_punkt = int(f.readline())
        info = f.readlines()

        knutepunkt = []
        elementer  = []
        laster     = []
        total_array= []
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
        return(total_array)


def local_k(filename):

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
        L = numpy.sqrt((x_ende1 - x_ende2)**2 + (y_ende1 - y_ende2)**2)
        E = elementer[i][2]
        I = get_Iprofil(elementer[i][3])
        A = get_areal(elementer[i][3])

        k = numpy.array([[E * A / L, 0., 0., E * A / L, 0., 0.],
                         [0., 12 * E * I / (L ** 3), -6 * E * I / (L ** 2), 0., -12 * E * I / (L ** 3),
                          -6 * E * I / (L ** 2)],
                         [0., -6 * E * I / (L ** 2), 4 * E * I / L, 0., 6 * E * I / (L ** 2), 2 * E * I / L],
                         [E * A / L, 0., 0., E * A / L, 0., 0.],
                         [0., -12 * E * I / (L ** 3), 6 * E * I / (L ** 2), 0., 12 * E * I / (L ** 3), 6 * E * I / L],
                         [0., -6 * E * I / (L ** 2), 2 * E * I / L, 0., 6 * E * I / L, 4 * E * I / L]])





    #Use a breakpoint in the code line below to debug your script.

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    get_knutepunkt('rammedata.txt')
    get_elementer('rammedata.txt')
    get_laster('rammedata.txt')
    local_k('rammedata.txt')
