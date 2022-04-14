import copy
import math
import pathlib
import random
import json
import matplotlib.pyplot as plt

path = str(pathlib.Path(__file__).parent.resolve())

with open(path + "\input.json") as fjsonin:
    fin = json.load(fjsonin)
fout = open(path + "\output.txt", "w")

# Date care vor fi citite din input
# Nr de indivizi din populatie
n = 0
# Domeniul de definitie al functiei
a, b = 0, 0
# Constantele functiei
A, B, C = 0, 0, 0
# Precizia (cate zecimale pe discretizare)
precizie = 0
# Probabilitatea de recombinare
probRecomb = 0
# Probablitatea de mutatie
probMutatie = 0
# Nr de etape dupa care oprim algoritmul
nrEtape = 0

# Date determinate ulterior
# Populatia curenta
populatie = []

# Fitness populatie
fitnessPopulatie = []

# Indivizii trecuti mai departe
selectPopulatie = []

# Populatia pe care se vor aplica mutatii
nouaPopulatie = []

# Individul elitist
elitist = []

# Lungimea unui cromozom
l = 0

# MinFitness
minFitness = float('inf')

# Vector de max pe functie
maxFunction = []


# Variabile pt plot
maxHistory = []
avgHistory = []


# Citire input
def readInput():
    global n, a, b, A, B, C, precizie, probRecomb, probMutatie, nrEtape
    n = fin['n']
    a, b = float(fin['a']), float(fin['b'])
    A, B, C = float(fin['A']), float(fin['B']), float(fin['C'])
    precizie = int(fin['precizie'])
    probRecomb = float(fin['probRecomb'])
    probMutatie = float(fin['probMutatie'])
    nrEtape = fin['nrEtape']


# Afisare in fisierul de output a ce am citit
def parseInput():
    fout.write(f"Marimea populatiei: {n}\n")
    fout.write(f"f(x) = {A}x^2 + {B}x + {C}\n")
    fout.write(f"Definita pe [{a}, {b})\n")
    fout.write(f"Precizia discretizarii: {precizie} zecimale\n")
    fout.write(f"Probabilitatea de recombinare: {probRecomb}\n")
    fout.write(f"Probabilitatea de mutatie: {probMutatie}\n")
    fout.write(f"Numarul de etape: {nrEtape}\n")


# Codificarea din curs
def codificare():
    global l
    aux = (b - a) * (10 ** precizie)
    aux = math.log(aux, 2)
    aux = math.ceil(aux)
    return aux


# Transformare din array in string
def parseIndividToString(individ):
    result = ''.join(str(i) for i in individ)
    return result


# Transformare din array in int
def parseIndividToInt(individ):
    result = parseIndividToString(individ)
    result = int(result, 2)
    result = (b - a) / (pow(2, l) - 1) * result + a
    return round(result, precizie)


# Generarea de populatie random de inceput
def setPopulatie():
    global populatie
    populatie = [[random.getrandbits(1) for i in range(l)] for j in range(n)]


# Generarea fitnessului pt fiecare individ si normalizare sa fie pozitive
def generateFitness():
    global populatie, fitnessPopulatie, minFitness
    fitnessPopulatie.clear()
    minFitness = float('inf')
    for individ in populatie:
        individInt = parseIndividToInt(individ)
        fitnessIndivid = functie(individInt)
        fitnessPopulatie.append(fitnessIndivid)
        minFitness = min(minFitness, fitnessIndivid)
    if minFitness < 0:
        for i in range(n):
            fitnessPopulatie[i] -= minFitness


# Scrierea populatiei in fisier
def parsePopulatie(populatie):
    for i in range(len(populatie)):
        individString = parseIndividToString(populatie[i])
        individInt = parseIndividToInt(populatie[i])
        individFitness = fitnessFunction(individInt)
        fout.write(f"Individ {i+1}: {individString} cu int = {individInt} cu f = {individFitness}\n")
    fout.write('\n')


# Generare max 
def generateFunction():
    global populatie, maxFunction
    maxFunction.clear()
    for individ in populatie:
        individInt = parseIndividToInt(individ)
        functieIndivid = functie(individInt)
        maxFunction.append(functieIndivid)


# Functia
def functie(x):
    return A * x ** 2 + B * x + C


# Fitness
def fitnessFunction(x):
    global minFitness
    if minFitness < 0:
        return functie(x) - minFitness
    return functie(x)


# Cauta intr-o lista list locul lui x prin cautare binara.
def cautare(list, x, s, d):
    if x<=list[s]:
        return s
    elif x>=list[d]:
        return d + 1
    elif s < d:
        mid = int((s + d) / 2)
        if list[mid + 1] > x and list[mid] <= x:
            return mid+1
        elif list[mid]<=x and list[mid+1]<=x:
            return cautare(list, x, mid+1, d)
        else:
            return cautare(list, x, s, mid-1)


# Selectia naturala
def selectie(afis = False):
    global n, a, b, precizie, populatie, fitnessPopulatie, selectPopulatie, elitist

    # Probabilitatile ca indivizii sa supravietuiasca
    probPopulatie = []

    # Intervalele pentru metoda ruletei
    intervalePopulatie = []

    if afis:
        fout.write("\nSituatia etapei 1 de selectie:\n\n")

    generateFitness()

    # Calculam fitnesul total
    fitnessTotal = sum(fitnessPopulatie)

    # Se gaseste elitistul
    maxFitness = max(fitnessPopulatie)
    indexElitist = fitnessPopulatie.index(maxFitness)
    
    # Se elimina elitistul din calcule, el e trecut deja
    elitist = copy.copy(populatie[indexElitist])
    populatie.pop(indexElitist)
    fitnessPopulatie.pop(indexElitist)
    
    if afis:
        fout.write(f"Elitistul generatiei este: {elitist}\n")

    m = n - 1

    # Se fac probabilitatile pt populatie
    for i in range(m):
        probPopulatie.append(fitnessPopulatie[i] / fitnessTotal)
    
    if afis:
        fout.write(f"Fitness total: {fitnessTotal}\n\nProbabilitati:\n")
        for i in range(m):
            fout.write(f"Individ {i+1}: cu probabilitate = {probPopulatie[i]}\n")
    

    # Metoda ruletei - generare de intervale
    curr = 0
    intervalePopulatie.append(0)
    for i in range(m):
        curr += probPopulatie[i]
        intervalePopulatie.append(curr)
    intervalePopulatie[m] = float(1)
    
    if afis:
        fout.write("\nIntervale:\n")
        for i in range(m):
            fout.write(f"Individ {i+1}: cu intervalul = [{intervalePopulatie[i]}, {intervalePopulatie[i+1]})\n")
        fout.write("\nGenerare de numere:\n")
    
    # Invartirea ruletei
    nrIndivizi = 0
    while nrIndivizi < m:
        nr = random.uniform(0, 1)
        individ = cautare(intervalePopulatie, nr, 0, m)
        if afis:
            fout.write(f"S-a generat {nr} si deci am ales individul {individ}\n")
        selectPopulatie.append(copy.deepcopy(populatie[individ - 1]))
        nrIndivizi += 1

    if afis:
        fout.write("\nPopulatia dupa selectie:\n")
        parsePopulatie(selectPopulatie)
        fout.write('\n')


# Combinare de 2 indivizi
def comb2(individ1, individ2, punctRupere):
    individ1[punctRupere:], individ2[punctRupere:] = individ2[punctRupere:], individ1[punctRupere:]
    cpIndivid1, cpIndivid2 = individ1, individ2
    return cpIndivid1, cpIndivid2


# Combinare de 3 indivizi
def comb3(individ1, individ2, individ3, punctRupere):
    individ1[punctRupere:], individ2[punctRupere:], individ3[punctRupere:] = individ2[punctRupere:], individ3[punctRupere:], individ1[punctRupere:]
    cpIndivid1, cpIndivid2, cpIndivid3 = individ1, individ2, individ3
    return cpIndivid1, cpIndivid2, cpIndivid3


# Etapa de recombinare a genelor
def recombinare(afis = False):
    global n, a, b, precizie, probRecomb, selectPopulatie, nouaPopulatie
    participantiPopulatie = []
    m = n - 1

    # Facem un shuffle sa se poate combina oricine cu oricine mereu
    random.shuffle(selectPopulatie)

    if afis:
        fout.write("Etapa de recombinare:\n\nPopulatia dupa shuffle:\n")
        parsePopulatie(selectPopulatie)
        fout.write("Se fac alegerile:\n")

    # Se decide luand in calcul probabilitatea cine se recombina
    for i in range(m):
        nr = random.uniform(0, 1)
        if nr < probRecomb:
            participantiPopulatie.append(copy.deepcopy(selectPopulatie[i]))
            if afis:
                fout.write(f"A fost ales individul {i+1}.\n")
        else:
            nouaPopulatie.append(copy.deepcopy(selectPopulatie[i]))
    
    if afis:
        fout.write("\nDeci populatia selectata pentru recombinare este:\n")
        parsePopulatie(participantiPopulatie)

    par = len(participantiPopulatie)

    # Ii combinam
    if par == 0:
        pass

    elif par == 1:
        nouaPopulatie.append(copy.deepcopy(participantiPopulatie[0]))
        if afis:
            fout.write(f"A trecut mai departe individul {participantiPopulatie[0]}\n")

    elif par % 2 == 0:
        for i in range(1, par, 2):
            individ1 = copy.deepcopy(participantiPopulatie[i-1])
            individ2 = copy.deepcopy(participantiPopulatie[i])
            punctRupere = random.randint(0, l)

            if afis:
                fout.write(f"Se combina individul {i}: {parseIndividToString(individ1)} cu individul {i+1}: {parseIndividToString(individ2)}\n")
                fout.write(f"Punct de rupere: {punctRupere}\n")

            noulIndivid1, noulIndivid2 = copy.deepcopy(comb2(individ1, individ2, punctRupere))
            
            if afis:
                fout.write(f"Ies indivizii: {parseIndividToString(noulIndivid1)} si {parseIndividToString(noulIndivid2)}\n\n")
            
            nouaPopulatie.append(copy.deepcopy(noulIndivid1))
            nouaPopulatie.append(copy.deepcopy(noulIndivid2))
    
    else:
        individ1 = copy.deepcopy(participantiPopulatie[0])
        individ2 = copy.deepcopy(participantiPopulatie[1])
        individ3 = copy.deepcopy(participantiPopulatie[2])
        punctRupere = random.randint(0, l)

        if afis:
            fout.write(f"Se combina individul 1: {parseIndividToString(individ1)} cu 2: {parseIndividToString(individ2)} si 3: {parseIndividToString(individ3)}\n")
            fout.write(f"Punct de rupere: {punctRupere}\n")

        noulIndivid1, noulIndivid2, noulIndivid3 = copy.deepcopy(comb3(individ1, individ2, individ3, punctRupere))

        if afis:
            fout.write(f"Ies indivizii: {parseIndividToString(noulIndivid1)} si {parseIndividToString(noulIndivid2)} si {parseIndividToString(noulIndivid3)}\n\n")

        nouaPopulatie.append(copy.deepcopy(noulIndivid1))
        nouaPopulatie.append(copy.deepcopy(noulIndivid2))
        nouaPopulatie.append(copy.deepcopy(noulIndivid3))

        for i in range(4, par, 2):
            individ1 = copy.deepcopy(participantiPopulatie[i-1])
            individ2 = copy.deepcopy(participantiPopulatie[i])
            punctRupere = random.randint(0, l)

            if afis:
                fout.write(f"Se combina individul {i}: {parseIndividToString(individ1)} cu individul {i+1}: {parseIndividToString(individ2)}\n")
                fout.write(f"Punct de rupere: {punctRupere}\n")

            noulIndivid1, noulIndivid2 = copy.deepcopy(comb2(individ1, individ2, punctRupere))
            
            if afis:
                fout.write(f"Ies indivizii: {parseIndividToString(noulIndivid1)} si {parseIndividToString(noulIndivid2)}\n\n")
            
            nouaPopulatie.append(copy.deepcopy(noulIndivid1))
            nouaPopulatie.append(copy.deepcopy(noulIndivid2))

    if afis:
        fout.write(f"Noua populatie dupa recombinare (ultimii {par} au fost recombinati):\n")
        parsePopulatie(nouaPopulatie)


# Mutatie populatie
def mutatie(afis = False):
    global n, a, b, precizie, probMutatie, selectPopulatie, nouaPopulatie

    if afis:
        fout.write("\nEtapa de mutatie:\n")
    
    m = n - 1
    mutatie = False

    # Facem mutatia
    for i in range(m):
        nr = random.uniform(0, 1)
        if nr < probMutatie:
            alea = random.randint(0, l - 1)
            
            if afis:
                mutatie = True
                fout.write(f"Sa mutat gena {alea+1} a individului {i}:\n")
                fout.write(f"S-a transformat din {parseIndividToString(nouaPopulatie[i])} ")

            nouaPopulatie[i][alea] = int(not nouaPopulatie[i][alea])

            if afis:
                fout.write(f"in {parseIndividToString(nouaPopulatie[i])}\n\n")
    
    # Adaugam inapoi elitistul
    nouaPopulatie.append(copy.deepcopy(elitist))

    if afis:
        if not mutatie:
            fout.write("Nu a avut loc nicio mutatie\n\n")
        fout.write("Populatia dupa mutatii, care va fi pasata la pasul urmator:\n")
        parsePopulatie(nouaPopulatie)


# Wrapper la tot
def generatie(etapa, afis = False):
    global n, a, b, probMutatie, maxHistory, avgHistory, maxFunction, probRecomb, precizie, populatie, fitnessPopulatie, nouaPopulatie, selectPopulatie, l, elitist, minFitness
    generateFitness()
    generateFunction()
    if afis:
        parsePopulatie(populatie)
        fout.write(f"Inainte de etapa {etapa} average fitness: {sum(fitnessPopulatie) / n} si average max {sum(maxFunction) / n} cu max {max(maxFunction)}\n")
        maxHistory.append(max(maxFunction))
        avgHistory.append(sum(maxFunction) / n)
    
    selectie(afis)
    recombinare(afis)
    mutatie(afis)

    populatie = copy.deepcopy(nouaPopulatie)
    generateFitness()
    generateFunction()
    fout.write(f"Dupa etapa {etapa} average fitness: {sum(fitnessPopulatie) / n} si average max {sum(maxFunction) / n} cu max {max(maxFunction)}\n")

    maxHistory.append(max(maxFunction))
    avgHistory.append(sum(maxFunction) / n)

    nouaPopulatie.clear()
    fitnessPopulatie.clear()
    selectPopulatie.clear()
    elitist.clear()



if __name__ == "__main__":
    readInput()
    parseInput()

    l = codificare()

    setPopulatie()
    fout.write("\nPopulatia de inceput:\n")

    etapa = 1
    while etapa <= nrEtape:
        if etapa == 1:
            generatie(etapa, True)
        else:
            generatie(etapa)
        etapa += 1

    etape = [i for i in range(0, nrEtape + 1)]
    plt.plot(etape, maxHistory)
    plt.plot(etape, avgHistory)
    plt.xlabel('Nr etapa')
    plt.ylabel('Albastru maxim\nPortocaliu medie')
    plt.show()

    fout.close()