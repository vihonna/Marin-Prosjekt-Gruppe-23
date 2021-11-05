import numpy


def les_input(filename):
    info = []
    with open(filename) as f:

        # Disse variablene vil brukes i veldig mange andre deler av oppgaven, så for å forenkle lagres dem som globale
        global n_punkt
        global n_element
        global n_laster

        # Dette er arrays som hentes ut fra inputfilen, og vil brukes i flere funksjoner
        global knutepunkt
        global elementer
        global laster

        n_punkt = int(f.readline())
        info = f.readlines()

        # Dersom disse listene enda er uferdige legges en 'temp' til for at de ikke skrives til de globale variablene med det første
        knutepunkt_temp = []
        elementer_temp = []
        laster_temp = []

        a = 0
        for i in range(0, len(info)):
            arr = list(map(int, info[i].split()))
            if (len(arr) < 2) and (a == 0):
                n_element = int(arr[0])
                a = 1
            elif (len(arr) < 2) and (a == 1):
                n_laster = int(arr[0])
                a = 2
            elif a == 0:
                arr[0] = arr[0]*0.001
                arr[1] = arr[1]*0.001
                knutepunkt_temp.append(arr)
            elif a == 1:
                arr[2] = arr[2]*10**6
                elementer_temp.append(arr)
            else:
                arr[2] = arr[2]*1000
                arr[3] = arr[3]*1000
                laster_temp.append(arr)

        knutepunkt = numpy.array(knutepunkt_temp)
        elementer = numpy.array(elementer_temp)
        laster = numpy.array(laster_temp)


def matriser_setup(filename):
    les_input(filename)

    # Her lager programmet to tomme svære matriser, hvor n_punkt er antall knutepunkt
    global SysMatrise
    SysMatrise = numpy.zeros((3 * n_punkt, 3 * n_punkt))
    global R_matrise
    R_matrise = numpy.zeros((3 * n_punkt, 1))

    # De to funksjonene som står for alt av utregning av matrisene, de returnerer ingenting da de bare skriver til globale arrays/matriser
    fyll_SystemMatrise()
    fyll_Rmatrise()


def fyll_SystemMatrise():
    # elementer og knutepunkt er 2D arrays hvor formatet er slik:
    # elementer[i] = [lokal ende 1, lokal ende 2, Emodul, I-profil nummer]
    # knutepunkt[i] = [x-koordinat, y-koordinat, fast innspent "ja(1)/nei(0)"]

    # Ettersom at lengde og vinkel er nyttig i lastberegingene lagres alle lengder og vinkler globalt
    global L_array
    global theta_array
    global k_array
    L_array = numpy.zeros(n_element)
    theta_array = numpy.zeros(n_element)
    k_array = []

    for i in range(0, len(elementer)):

        # Generell informasjon om hver bjelke
        ende1 = elementer[i][0]
        ende2 = elementer[i][1]
        x_ende1 = knutepunkt[ende1][0]
        y_ende1 = knutepunkt[ende1][1]
        x_ende2 = knutepunkt[ende2][0]
        y_ende2 = knutepunkt[ende2][1]

        # Henter nyttig info om tverrsnitt
        Tverrsnittsdata = get_Tverrsnittsdata(elementer[i][3])

        # Finner verdier som trengs for lokal k matrise for hver enkelt bjelke
        L = numpy.sqrt((x_ende1 - x_ende2) ** 2 + (y_ende1 - y_ende2) ** 2)
        #print('L for bjelke ', i, ' er ', L)
        E = elementer[i][2]
        I = Tverrsnittsdata[0]
        A = Tverrsnittsdata[1]
        theta = numpy.arcsin((y_ende2 - y_ende1) / L)

        # Hvis en ende er fastinnspent skal stivheten økes kraftig, det gjøres ved å legge til en variabel a1/a2
        arr = fastInnspentSjekk(ende1, ende2)
        a1 = arr[0]
        a2 = arr[1]

        # Lagrer lengde og vinklene til bjelkene i globale arrays deklarert i begynnelsen av funksjonen
        L_array[i] = L
        theta_array[i] = theta

        # Setter opp den generelle stivhetsmatrisen for et element
        k = numpy.array([[a1*E*A/L, 0., 0., -a2*E*A/L, 0., 0.],
                         [0., a1*12*E*I/(L**3), a1*(-6)*E*I/(L**2), 0., a2*(-12)*E*I/(L**3), a2*(-6)*E*I/(L ** 2)],
                         [0., a1*(-6)*E*I/(L**2), a1*4*E*I/L, 0., a2*6*E*I/(L**2), a2*2*E*I/L],
                         [-a1*E*A/L, 0., 0., a2*E*A/L, 0., 0.],
                         [0., a1*(-12)*E*I/(L**3), a1*6*E*I/(L**2), 0., a2*12*E*I/(L ** 3), a2*6*E*I/L],
                         [0., a1*(-6)*E*I/(L**2), a1*2*E*I/L, 0., a2*6*E*I/(L**2), a2*4*E*I/L]])

        #Lagrer lokale k matriser i en liste som senere brukes til å regne ut hvert enkelt endemoment osv
        k_array.append(k)

        # Omgjør elementets lokale stivhetsmatrise til global stivhetsmatrise
        T = get_Tmatrise(theta)
        kg_temp = numpy.dot(T, k)
        T_inv = numpy.matrix.transpose(T)
        kg = numpy.dot(kg_temp, T_inv)

        #print('kg for bjelke ', i, ' er: ', kg.astype(int))

        # Deler opp 6x6 matrisen i fire 3x3 matriser, som seperat skal legges inn i system matrisen
        k1 = [[kg[0][0], kg[0][1], kg[0][2]], [kg[1][0], kg[1][1], kg[1][2]], [kg[2][0], kg[2][1], kg[2][2]]]
        k2 = [[kg[0][3], kg[0][4], kg[0][5]], [kg[1][3], kg[1][4], kg[1][5]], [kg[2][3], kg[2][4], kg[2][5]]]
        k3 = [[kg[3][0], kg[3][1], kg[3][2]], [kg[4][0], kg[4][1], kg[4][2]], [kg[5][0], kg[5][1], kg[5][2]]]
        k4 = [[kg[3][3], kg[3][4], kg[3][5]], [kg[4][3], kg[4][4], kg[4][5]], [kg[5][3], kg[5][4], kg[5][5]]]

        #print('k4 for bjelke ', i, ' er: ', k4)
        #print('for bjelke ', i, ' er den lok stivhetsmatrisen: ')
        #print(kg.astype(int))
        # Adderer inn hver 3x3 del av k matrisen, hvor adder funksjon også trenger å vite node (endei)
        k_mini = numpy.array([k1, k2, k3, k4])
        for l in range(4):
            sysMatrise_adder(k_mini[l], l, ende1, ende2)
        #print('SysMatrise nummer ', i, ': ', SysMatrise.astype(int))

def fyll_Rmatrise():
    # Hver laster[i] er et array med lengde 5, som består av noder, intensitet og bjelke den ligger på
    # Laster[i] = [node1, node2, intensitet node 1, intensitet node2, elementnummer]

    # Lastene vil også brukes til å finne lokale endekrefter/moment og lagres derfor i et eget array
    global R_array
    R_array = [0]*n_element

    #Lager en variabel som brukes til å lage momentdiagram ved å vite hvilken type last det er
    a = 0

    for i in range(n_laster):


        # Ettersom at L_array og theta_array er globale variabler kan funksjonen kalle på dem direkte
        L = L_array[laster[i][4]]
        theta = theta_array[laster[i][4]]

        # Lastene vil også brukes til å finne lokale endekrefter/moment og lagres derfor i et eget array

        # For å vite hvor store m og q er trengs det å vite hvilken type last det er
        if laster[i][0] == laster[i][1]:  # punktlast, nr = 1
            P = laster[i][2]
            m1 = 0
            m2 = 0
            q1 = P
            q2 = 0
            a = 1
        elif laster[i][2] == 0:  # Trekantlast hvor lokal ende 1 er 0
            P = laster[i][3]
            m1 = -P * L * L / 30
            m2 = P * L * L / 20
            q1 = 3 * P * L / 20
            q2 = 7 * P * L / 20
        elif laster[i][3] == 0:  # Trekantlast hvor lokal ende 2 er 0
            P = laster[i][2]
            m1 = -P * L * L / 20
            m2 = P * L * L / 30
            q1 = 7 * P * L / 20
            q2 = 3 * P * L / 20
        elif laster[i][2] == laster[i][3]:
            P = laster[i][2]
            m1 = -P * L * L / 12
            m2 = P * L * L / 12
            q1 = P * L / 2
            q2 = P * L / 2
        elif (laster[i][2] > laster[i][3]) or (laster[i][3] > laster[i][2]):  # Firkant + trekantlast
            Pmin = min(laster[i][2], laster[i][3])
            Pmaks = max(laster[i][2], laster[i][3])
            Pdiff = Pmaks - Pmin

            # Her benyttes superposisjon, hvor de summeres til slutt
            m1 = -Pmin * L * L / 12
            m2 = Pmin * L * L / 12
            q1 = Pmin * L / 6
            q2 = Pmin * L / 3
            if laster[i][2] > laster[i][3]:
                m1 = m1 + Pdiff * L * L / 20
                m2 = m2 - Pdiff * L * L / 30
                q1 = q1 + 7 * Pdiff * L / 20
                q2 = q2 + 3* Pdiff * L / 20
            if laster[i][3] > laster[i][2]:
                m1 = m1 - Pdiff * L * L / 30
                m2 = m2 + Pdiff * L * L / 20
                q1 = q1 + 3 * Pdiff * L / 20
                q2 = q2 + 7 * Pdiff * L / 20

        # Programmet legger til mi og vi til begge de lokale endene
        Rmatrise_adder(laster[i][0], m1, q1, 0, theta)
        Rmatrise_adder(laster[i][1], m2, q2, 0, theta)
        R_array[i] = [0, m1, q1, 0, m2, q2, a]


def sysMatrise_adder(ki, number, ende1, ende2):
    # Variabelen number er den som forteller hvilken av de 4 3x3 matrisene som skal legges til
    # Denne metoden å appendere til Systemmatrisen er hentet fra den første forelesningen i oktober
    # For feks element mellom node 1 og 2 vil:
    # SysMatrise[3*ende1 + i][3*ende1 + j] = ki[i][j] ---> Sysmatrise[3 + 0][3 + 0] = [EA/L]
    #print('number: ', number, 'ki: ', ki)
    if number == 0:
        for i in range(3):
            for j in range(3):
                SysMatrise[3 * ende1 + i][3 * ende1 + j] = SysMatrise[3 * ende1 + i][3 * ende1 + j] + ki[i][j]
    if number == 1:
        for i in range(3):
            for j in range(3):
                SysMatrise[3 * ende1 + i][3 * ende2 + j] = SysMatrise[3 * ende2 + i][3 * ende2 + j] + ki[i][j]
    if number == 2:
        for i in range(3):
            for j in range(3):
                SysMatrise[3 * ende2 + i][3 * ende1 + j] = SysMatrise[3 * ende2 + i][3 * ende1 + j] + ki[i][j]
    if number == 3:
        for i in range(3):
            for j in range(3):
                SysMatrise[3 * ende2 + i][3 * ende2 + j] = SysMatrise[3 * ende2 + i][3 * ende2 + j] + ki[i][j]


def Rmatrise_adder(node, mi, qi, ni, theta):
    arr = numpy.array([qi*numpy.sin(theta), qi*numpy.cos(theta), mi])
    for i in range(3):
        R_matrise[3*node + i] = R_matrise[3*node + i] + arr[i]


def get_Tmatrise(theta):
    # Lager variablene som skal inn, theta er vinkelen bjelke[i] har på det globale koordinatsystemet
    s = numpy.sin(theta)
    c = numpy.cos(theta)

    # Setter opp den generelle Tmatrisen, som kun består av sinus og cosinus som variabler
    T_matrise = numpy.array([[c, -s, 0, 0, 0, 0], [s, c, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],
                             [0, 0, 0, c, -s, 0], [0, 0, 0, s, c, 0], [0, 0, 0, 0, 0, 1]])
    return T_matrise


def fastInnspentSjekk(ende1, ende2):
    a1 = 1
    a2 = 1

    if knutepunkt[ende1][2] == 1:
        a1 = 10 ** 6
    if knutepunkt[ende2][2] == 1:
        a2 = 10 ** 6
    return [a1, a2]


def get_Tverrsnittsdata(profil):
    return [1.29*10**(2), 2]


def finn_deformasjoner():
    SysMatrise_invertert = numpy.linalg.inv(SysMatrise)
    deformasjoner = numpy.linalg.solve(SysMatrise, R_matrise)
    #print('Deformasjoner: ')
    #print(deformasjoner)
    return deformasjoner


def kalkuler_krefter():
    for i in range(0, n_element):
        r = finn_deformasjoner()
        ende1 = elementer[i][0]
        ende2 = elementer[i][1]
        r_i = numpy.zeros((6, 1))
        for j in range(0, 6):
            if j < 3:
                r_i[j] = (r[ende1*3 + j])
            else:
                r_i[j] = (r[ende2*3 + (j-3)])

        if R_array[i] != 0:
            R_i = numpy.zeros((6, 1))
            for j in range(0, 6):
                R_i[j] = R_array[i][j]

        Tmatrise = get_Tmatrise(theta_array[i])
        Tinv = numpy.linalg.inv(Tmatrise)
        v = numpy.dot(Tinv, r_i)
        k_i = numpy.array(k_array[i])
        krefter = numpy.dot(k_i, v) - R_i
        print('for bjelke ', i, 'blir kreftene: ')
        print(krefter)

        #Lager momentdiagram:




if __name__ == '__main__':
    matriser_setup('rammedata.txt')
    print(SysMatrise.astype(int))
    #print(SysMatrise)
    print(finn_deformasjoner())
    kalkuler_krefter()
    print(R_matrise)
    # print(int_SysMatrise)
    #print(R_matrise.astype(int))
