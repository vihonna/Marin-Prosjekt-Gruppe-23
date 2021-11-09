import math
import numpy


def les_input(filename):
    info = []
    with open(filename) as f:

        # Disse variablene vil brukes i veldig mange andre deler av oppgaven, så for å forenkle lagres dem som globale
        global n_punkt
        global n_element
        global n_laster
        global n_tverrsnitt

        # Dette er arrays som hentes ut fra inputfilen, og vil brukes i flere funksjoner
        global knutepunkt
        global elementer
        global laster
        global tverrsnitt

        # Det er satt krav til inputfilen om at første tall skal være npunktet, så første linje programmet leser vil
        # da inneholde inputfil
        n_punkt = int(f.readline())

        info = f.readlines()

        # Dersom disse listene enda er uferdige legges en 'temp' til for at de ikke skrives til de globale variablene med det første
        knutepunkt_temp = []
        elementer_temp = []
        laster_temp = []
        tverrsnitt_temp = []

        a = 0

        # Her skilles det mellom om hva en gitt linje av tekstfilen inneholder.
        for i in range(0, len(info)):  # Går gjennom filen linje for linje
            arr = list(map(int, info[i].split()))  # Lager en tekstlinje om til et array som inneholder heltall
            if (len(arr) < 2) and (a == 0):  # Linjen inneholder tallet som tilsvarer n_element
                n_element = int(arr[0])
                a = 1
            elif (len(arr) < 2) and (a == 1):  # Linjen inneholder tallet som tilsvarer antall elementer
                n_laster = int(arr[0])
                a = 2
            elif (len(arr) < 2) and (a == 2):
                n_tverrsnitt = int(arr[0])
                a = 3
            elif a == 0:  # Linjen består av et array som appenderer til knutepunkt array
                arr[0] = arr[0]*0.001
                arr[1] = arr[1]*0.001
                knutepunkt_temp.append(arr)
            elif a == 1: # Linjen består av et array som appenderer til elementer array
                arr[2] = arr[2]*10**6
                elementer_temp.append(arr)
            elif a == 2:  # Linjen bestpr av et array som appenderer til laster array
                arr[2] = arr[2]*1000
                arr[3] = arr[3]*1000
                laster_temp.append(arr)
            elif a == 3:  # Linjen bestpr av et array som appenderer til tverrsnitt array
                arr = numpy.dot(0.001, arr)
                arr[0] = 1000
                tverrsnitt_temp.append(arr)

        knutepunkt = numpy.array(knutepunkt_temp)
        elementer = numpy.array(elementer_temp)
        laster = numpy.array(laster_temp)
        tverrsnitt = numpy.array(tverrsnitt_temp)


def matriser_setup(filename):
    les_input(filename)

    # Her lager programmet to tomme svære matriser, hvor n_punkt er antall knutepunkt
    global SysMatrise
    SysMatrise = numpy.zeros((3 * n_punkt, 3 * n_punkt))
    global R_matrise
    R_matrise = numpy.zeros((3 * n_punkt, 1))


def get_Tverrsnittsdata(profil):
    pi = math.pi
    fy = 0  # Initialverdi

    # parameteren som tas inn i denne funksjonen er et tall som beskriver hvilken utregning av I og A som tilhører
    # hvilken profil. Funksjonen her har kun som oppgave å ta inn geometri mål og regne ut, ikke å også bestemme disse
    # parametrene
    # Tverrsnittet skrives ulikt avhengig av om det regner ut et rør eller Iprofil.
    # Hvis rør    : Tverrsnitt[i] = [1, R, r, 0, 0]
    # Hvis Iprofil: Tverrsnitt[i] = [0, B, H, b, h]
    # B = bredde flens, H = tykkelse flens, b = tykkelse steg, h = høyde steg

    if profil == 0:
        R = tverrsnitt[0][1] / 2
        t = tverrsnitt[0][2]
        r = R - t
        I = pi * (R**4 - r**4) / 4
        A = pi * (R**2 - r**2)
        z = R

    elif profil == 1:
        B = tverrsnitt[2][1]
        H = tverrsnitt[2][2]
        b = tverrsnitt[2][3]
        h = tverrsnitt[2][4]
        I = 2 * ((B * H * H * H) / 12 + B * B * H * H * (H + h * 0.5)) + (b * h * h * h) / 12
        A = 2 * B * H + b * h
        z = (2*H + h)/2

    elif profil == 2:
        B = tverrsnitt[2][1]
        H = tverrsnitt[2][2]
        b = tverrsnitt[2][3]
        h = tverrsnitt[2][4]
        I = 2 * ((B * H * H * H) / 12 + B * B * H * H * (H + h * 0.5)) + (b * h * h * h) / 12
        A = 2 * B * H + b * h
        z = (2 * H + h) / 2

    elif profil == 4:  # Testprofil (alt blir bare 1)
        I = pi * (0.001**4)/4
        A = pi * (0.001**2)
    else:
        print('Valid cross section not found')

    # Finner flytspenning, første element i array tilsvarer om det er stål eller alu
    if tverrsnitt[profil][0] == 0:  # Stål
        fy = 300000
    if tverrsnitt[profil][0] == 1:  # Aluminium
        fy = 100000
    return [I, A, z, fy]


def get_Tmatrise(theta):
    # Lager variablene som skal inn, theta er vinkelen bjelke[i] har på det globale koordinatsystemet
    s = numpy.sin(theta)
    c = numpy.cos(theta)

    # Setter opp den generelle Tmatrisen, som kun består av sinus og cosinus som variabler
    T_matrise = numpy.array([[c, -s, 0, 0, 0, 0], [s, c, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],
                             [0, 0, 0, c, -s, 0], [0, 0, 0, s, c, 0], [0, 0, 0, 0, 0, 1]])
    return T_matrise


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

        # Henter ut de to lokale bjelkene
        ende1 = elementer[i][0]
        ende2 = elementer[i][1]

        # Koordinatene til de lokale endene til bjelke nummer i
        x_ende1 = knutepunkt[ende1][0]
        y_ende1 = knutepunkt[ende1][1]
        x_ende2 = knutepunkt[ende2][0]
        y_ende2 = knutepunkt[ende2][1]

        # Henter nyttig info om tverrsnitt
        Tverrsnittsdata = get_Tverrsnittsdata(elementer[i][3])

        # Finner verdier som trengs for lokal k matrise for hver enkelt bjelke
        L = numpy.sqrt((x_ende1 - x_ende2) ** 2 + (y_ende1 - y_ende2) ** 2)
        E = elementer[i][2]
        I = Tverrsnittsdata[0]
        A = Tverrsnittsdata[1]

        # Finner vinkelen ved hjelp av en sinusfunksjon, arcsin fant vi ut at kun henter ut verdier mellom 90 og -90
        # grader, så i oppbygningen av inputfil må man være veldig konsekvent med valg av lokal ende 1 og 2 for at
        # vinkel mellom alltid vil mindre enn 90 grader og større enn -90.
        theta = numpy.arcsin((y_ende2 - y_ende1)/L)

        # Lagrer lengde og vinklene til bjelkene i globale arrays deklarert i begynnelsen av funksjonen
        L_array[i] = L
        theta_array[i] = theta

        # Setter opp den generelle stivhetsmatrisen for et element
        k = numpy.array([[E*A/L, 0., 0., -E*A/L, 0., 0.],
                         [0., 12*E*I/(L**3), (-6)*E*I/(L**2), 0., (-12)*E*I/(L**3), (-6)*E*I/(L ** 2)],
                         [0., (-6)*E*I/(L**2), 4*E*I/L, 0., 6*E*I/(L**2), 2*E*I/L],
                         [-E*A/L, 0., 0., E*A/L, 0., 0.],
                         [0., (-12)*E*I/(L**3), 6*E*I/(L**2), 0., 12*E*I/(L ** 3), 6*E*I/(L**2)],
                         [0., (-6)*E*I/(L**2), 2*E*I/L, 0., 6*E*I/(L**2), 4*E*I/L]])

        # Lagrer lokale k matriser i en liste som senere brukes til å regne ut kreftene på hver bjelke
        k_array.append(k)

        # Henter transformasjonsmatrisen T for bjelke i, som transponeres og bruker disse til å globalisere k matrisen
        T = get_Tmatrise(theta)
        T_inv = numpy.matrix.transpose(T)
        kg = T.dot(k).dot(T_inv)

        p1 = ende1*3  # Global posisjon til lokal ende 1 i stivhetsmatrisen
        p2 = ende2*3  # Global posisjon til lokal ende 2 i stivhetsmatrisen

        # Deler opp 6x6 matrisen i fire 3x3 matriser, som seperat skal legges inn i system matrisen
        SysMatrise[p1:p1+3, p1:p1+3] += kg[0:3, 0:3]  # kii
        SysMatrise[p2:p2+3, p2:p2+3] += kg[3:6, 3:6]  # kjj
        SysMatrise[p1:p1+3, p2:p2+3] += kg[0:3, 3:6]  # kij
        SysMatrise[p2:p2+3, p1:p1+3] += kg[3:9, 0:3]  # kji


def fyll_Rmatrise():
    # Hver laster[i] er et array med lengde 5, som består av noder, intensitet og bjelke den ligger på
    # Laster[i] = [node1, node2, intensitet node 1, intensitet node2, elementnummer, fortegn]

    # Lastene vil også brukes til å finne lokale endekrefter/moment og lagres derfor i et eget array
    global R_matrise_c
    R_matrise_c = numpy.zeros((6 * n_element, 1))

    for i in range(n_laster):

        # Definerer endenodene for hver last
        node1 = laster[i][0]
        node2 = laster[i][1]

        # Ettersom at L_array og theta_array er globale variabler kan funksjonen kalle på dem direkte
        L = L_array[laster[i][4]]
        theta = theta_array[laster[i][4]]

        # Fortegn beskriver om lasten peker på bjelken eller fra bjelken, peker den fra vil denne verdien være -1
        fortegn = laster[i][5]

        # Tall som beskriver om en gitt last skal tas med i beregning av S
        s_incl = 1

        # For å vite hvor store m og q er trengs det å vite hvilken type last det er
        if laster[i][0] == laster[i][1]:  # Punktlast
            P = laster[i][2]
            m1 = 0
            m2 = 0
            q1 = P * fortegn
            q2 = 0
            s_incl = 0
        else:  # Jevnt fordelt last
            P1 = laster[i][2]
            P2 = laster[i][3]
            m1 = ( -P1*L*L/30 - P2*L*L/20 ) * fortegn
            m2 = ( P1*L*L/20 + P2*L*L/30  ) * fortegn
            q1 = ( (7*P1*L + 3*P2*L) / 20 ) * fortegn
            q2 = ( (3*P1*L + 7*P2*L) / 20 ) * fortegn

        # Programmet legger inn i lastvektor, samt at gitte last omdannes til å passe i globalt koordinatysystem
        # Legger inn for endenode 1:
        R_matrise[3 * node1 + 0][0] += -q1 * numpy.sin(theta)
        R_matrise[3 * node1 + 1][0] += q1 * numpy.cos(theta)
        R_matrise[3 * node1 + 2][0] += m1

        # Legger inn for endenode 2:
        R_matrise[3 * node2 + 0][0] += -q2 * numpy.sin(theta)
        R_matrise[3 * node2 + 1][0] += q2 * numpy.cos(theta)
        R_matrise[3 * node2 + 2][0] += m2

        element_i = laster[i][4]
        # Trengs også å beholde verdiene i lokalt system
        R_matrise_c[6 * element_i  + 0][0] += 0
        R_matrise_c[6 * element_i + 1][0] += q1
        R_matrise_c[6 * element_i + 2][0] += m1
        R_matrise_c[6 * element_i + 3][0] += 0
        R_matrise_c[6 * element_i + 4][0] += q2
        R_matrise_c[6 * element_i + 5][0] += m2


        print('Globalisert')
        print('node nummer: ', laster[i][0], ' qi: ', q1*numpy.cos(theta), ' mi: ', m1, ' n1: ', -q1*numpy.sin(theta))
        print('node nummer: ', laster[i][1], ' qi: ', q2*numpy.cos(theta), ' mi: ', m2, ' n1: ', -q2*numpy.sin(theta))
        print(' ')
        print(' ')


def fastInnspentSjekk():
    # Funksjonen iterer gjennom alle knutepunkt, og dem som er fast innspent er markert med 1. Adderer dermed en
    # fjær med svært høy stivhet til knutepunktets diagonal
    for node in range(0, n_punkt):
        if knutepunkt[node][2] == 1:
            SysMatrise[3*node][3*node] += 10**20
            SysMatrise[3*node + 1][3*node + 1] += 10**20
            SysMatrise[3*node + 2][3*node + 2] += 10**20


def finn_deformasjoner():

    # Legger fast innspent fjærer til nodene som er fast innspent.
    fastInnspentSjekk()

    # Regner ut deformasjoner i det globale koordinatsystemet
    deformasjoner = numpy.linalg.solve(SysMatrise, -R_matrise)
    return deformasjoner


def kalkuler_krefter():

    #Iterer gjennom hver bjelke
    for i in range(0, n_element):

        # Henter de lokale endene
        ende1 = elementer[i][0]
        ende2 = elementer[i][1]

        # r er vektoren for deformasjoner i globalt koordinatsystem
        r = finn_deformasjoner()

        # r_i er en 6x1 vektor som inneholde globale deformasjoner i hver ende av bjelken
        r_i = numpy.zeros((6, 1))
        for j in range(0, 6):
            if j < 3:
                r_i[j] = r[3*ende1 + j]
            else:
                r_i[j] = r[3*ende2 + (j - 3)]

        # Omdannes r_i til sitt opprinnelige globale koordinatsystem
        T = get_Tmatrise(theta_array[i])
        T_trans = numpy.matrix.transpose(T)
        v = numpy.matmul(T_trans, r_i)

        # R_matrise_c er matrisen som inneholdt innfestningskreftene i hver ende av bjelken
        # Den er altså på størrelsesorden n_element x 6
        R_i = numpy.zeros((6, 1))
        for j in range(0, 6):
            R_i[j] = R_matrise_c[i*6 + j]

        # Henter den lokale stivhetsmatrisen for det gitte elementet
        k_i = numpy.array(k_array[i])
        print('k_i: ')
        print(k_i)

        # Løser ligningene for kreftene i S
        S = numpy.matmul(k_i, v)

        # Printer kreftene ut oversiktilig
        print(' ')
        print('Last: ')
        print(R_i)
        print('for bjelke ', i, 'blir kreftene: ')
        print(S)
        midtmoment = midtmomenter(S, i)
        sjekk_spenning(i, midtmoment)


def midtmomenter(krefter, element):

    # Legger initialverdiene
    L = L_array[element]
    p1 = 0
    p2 = 0

    # Hvis bjelken er påført last appenderes verdien i hver node til p1 og p2
    for i in range(0, n_laster):
        if laster[i][4] == element:
            p1 = laster[i][2]
            p2 = laster[i][3]

    # Bidrag fra endemomenter:
    endemoment_bidrag = ( krefter[2][0] + krefter[5][0]) / 2

    # Bidrag fra last
    last_bidrag = (p1 * L)/16 + (p2 * L)/16

    midtmoment = endemoment_bidrag + last_bidrag

    print('Midtmoment på bjelke ', element, ' er:', midtmoment)

def sjekk_spenning(Mmax, i):
    tverrsnitt = get_Tverrsnittsdata(elementer[i][3])
    z = tverrsnitt[2]
    I = tverrsnitt[1]
    fy = tverrsnitt[0]
    spenning = (Mmax*z)/I
    if spenning > fy:
        print('Underdimensjonert for element: ', i)


if __name__ == '__main__':
    # Initialiserer, og danner SysMatrise og lastvektor
    matriser_setup('rammedata.txt')

    # Fyller inn systemmatrise og lastvektor
    fyll_SystemMatrise()
    fyll_Rmatrise()

    # Regner ut deformasjoner
    r1 = finn_deformasjoner()

    #Printer ut deformasjoner i mm
    print('Deformasjoner i mm')
    print(numpy.around(numpy.dot(1000, r1)))
    kalkuler_krefter()
