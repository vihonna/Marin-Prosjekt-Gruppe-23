# coding=utf-8
# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

# les_input er en funksjon som tar inn data fra en fil og skriver ut formatet

def get_knutepunkt(filename):
    array = les_input(filename)
    return(array[0])

def get_elementer(filename):
    array = les_input(filename)
    return(array[1])

def get_laster(filename):
    array = les_input(filename)
    return(array[2])

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
        for i in range(1, len(info)):
            arr = list(map(int, info[i].split()))
            if i < n_punkt:
                knutepunkt.append(arr)
            if i > n_punkt and len(arr) > 1:
                elementer.append(arr)
                n_element = n_element + 1
            if i > (n_punkt + n_element):
                laster.append(arr)
        print(laster)
        total_array.append(knutepunkt)
        total_array.append(elementer)
        total_array.append(laster)
        return(total_array)

    # Use a breakpoint in the code line below to debug your script.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    get_knutepunkt('rammedata.txt')
    get_elementer('rammedata.txt')
    get_laster('rammedata.txt')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
