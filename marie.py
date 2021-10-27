from main import get_knutepunkt
from main import get_laster
from main import get_elementer
from main import les_input

import numpy as np


def get_theta():
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

        theta_1 = np.arctan(x_ende1/y_ende1)
        theta_2 = np.arctan(x_ende2/y_ende2)

    

get_theta()




# def get_transformasjonsmatrise(): #en funksjon som lager transformasjonsmatrisen
#    #t = 0
#    # T = [np.cos(t), np.cos(t), np.sin(t), np.cos(t) ]"
#     matrix = []
#     row=[]
#     for i in range(6):
#         a=[]
#         for j in range(6):
#             a.append('i')
#         matrix.append(a)

#     for j in range(6):
#         for k in range(6):
#             print(matrix[j][k], end = ' ')
#         print()

# get_transformasjonsmatrise()

t = 2

T = np.array( [ [ np.cos(t), -np.sin(t), 0. ,0., 0., 0.] , 
[ np.sin(t), np.cos(t), 0., 0., 0., 0.] ,
[0., 0., 1, 0., 0., 0.] , 
[0., 0., 0., np.cos(t), -np.sin(t), 0.], 
[0., 0., 0., np.sin(t), np.cos(t), 0.], 
[0., 0., 0., 0., 0., 1] ])