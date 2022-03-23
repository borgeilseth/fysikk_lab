# TFY41xx Fysikk vaaren 2021.
#
# Programmet tar utgangspunkt i hoeyden til de 8 festepunktene.
# Deretter beregnes baneformen y(x) ved hjelp av 7 tredjegradspolynomer,
# et for hvert intervall mellom to festepunkter, slik at baade banen y,
# dens stigningstall y' = dy/dx og dens andrederiverte
# y'' = d2y/dx2 er kontinuerlige i de 6 indre festepunktene.
# I tillegg velges null krumning (andrederivert)
# i banens to ytterste festepunkter (med bc_type='natural' nedenfor).
# Dette gir i alt 28 ligninger som fastlegger de 28 koeffisientene
# i de i alt 7 tredjegradspolynomene.

# De ulike banene er satt opp med tanke paa at kula skal
# (1) fullfoere hele banen selv om den taper noe mekanisk energi underveis;
# (2) rulle rent, uten aa gli ("slure").

# Vi importerer noedvendige biblioteker:
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

# Horisontal avstand mellom festepunktene er 0.200 m
h = 0.200
xfast = np.asarray([0, h, 2 * h, 3 * h, 4 * h, 5 * h, 6 * h, 7 * h])


# Vi begrenser starthøyden (og samtidig den maksimale høyden) til
# å ligge mellom 250 og 300 mm
ymax = 300
# yfast: tabell med 8 heltall mellom 50 og 300 (mm); representerer
# høyden i de 8 festepunktene
yfast = np.array([0.297, 0.237, 0.198, 0.165, 0.125, 0.158, 0.171, 0.103])


# inttan: tabell med 7 verdier for (yfast[n+1]-yfast[n])/h (n=0..7); dvs
# banens stigningstall beregnet med utgangspunkt i de 8 festepunktene.
inttan = np.diff(yfast) / h
attempts = 1

file1 = open("tracker data fil", "r")
Lines = file1.readlines()

tExp = []
xExp = []
yExp = []
vExp = []

count = 0
for line in Lines:
    count += 1
    if count <= 2:
        continue

    numbers = line.split(" ")
    tExp.append(float(numbers[0]))
    xExp.append(float(numbers[1]))
    yExp.append(float(numbers[2]) - 0.10)
    vExp.append(float(numbers[3]))


vExp = np.asarray(vExp)
tExp = np.asarray(tExp)
xExp = np.asarray(xExp)
yExp = np.asarray(yExp)

# Når programmet her har avsluttet while-løkka, betyr det at
# tallverdiene i tabellen yfast vil resultere i en tilfredsstillende bane.

# Programmet beregner deretter de 7 tredjegradspolynomene, et
# for hvert intervall mellom to nabofestepunkter.


# Med scipy.interpolate-funksjonen CubicSpline:
cs = CubicSpline(xfast, yfast, bc_type="natural")

xmin = 0.000
xmax = 1.401
dx = 0.001

x = np.arange(xmin, xmax, dx)


# Differentiating to calculate experimental results
xDiff = np.array(xExp)
yDiff = np.array(yExp)

dyDiff = np.zeros(len(xDiff) - 2)

for i in range(len(dyDiff)):
    dyDiff[i] = (yDiff[i + 2] - yDiff[i]) / (xDiff[i + 2] - xDiff[i])

d2yDiff = np.zeros(len(xDiff) - 4)

for i in range(len(d2yDiff)):
    d2yDiff[i] = (dyDiff[i + 2] - dyDiff[i]) / (xDiff[i + 3] - xDiff[i + 1])

# Merk: De deriverte har nå to færre elementer! De mangler den første og siste verdien
# For å bruke disse i beregninger må vi derfor slice arraysene.

# Eksempel:
# A = dy/x[1:len(x)-1:]


# funksjonen arange returnerer verdier paa det "halvaapne" intervallet
# [xmin,xmax), dvs slik at xmin er med mens xmax ikke er med. Her blir
# dermed x[0]=xmin=0.000, x[1]=xmin+1*dx=0.001, ..., x[1400]=xmax-dx=1.400,
# dvs x blir en tabell med 1401 elementer

Nx = len(x)
y = cs(x)  # y=tabell med 1401 verdier for y(x)
dy = cs(x, 1)  # dy=tabell med 1401 verdier for y'(x)
d2y = cs(x, 2)  # d2y=tabell med 1401 verdier for y''(x)

# Eksempel: Plotter banens form y(x)
baneform = plt.figure("y(x)", figsize=(12, 6))
plt.plot(x, y)
plt.plot(xExp, yExp)
plt.title("Banens form")
plt.xlabel("$x$ (m)", fontsize=20)
plt.ylabel("$y(x)$ (m)", fontsize=20)
plt.ylim(0.0, 0.40)
plt.grid()
plt.show()

v = np.sqrt((2 * 9.81 * (y[0] - y)) / (1 + (2 / 5)))

kulefart = plt.figure("y(x)", figsize=(12, 6))
plt.plot(x, v)
plt.plot(xExp, vExp)
plt.title("Kulas fart")
plt.xlabel("$x$ (m)", fontsize=20)
plt.ylabel("$v(x)$ (m/s)", fontsize=20)
plt.grid()
plt.show()

# k = np.(d2y / (1 + dy**2)**(3/2))

b = np.arctan(dy)
bExp = np.arctan(dyDiff)

helling = plt.figure("y(x)", figsize=(12, 6))
plt.plot(x, b)
plt.plot(xExp[1:-1], bExp)
plt.title("Helningsvinkel")
plt.xlabel("$x$ (m)", fontsize=20)
plt.ylabel("$vinekl$ (rad)", fontsize=20)
plt.grid()
plt.show()


vx = np.multiply(v, np.cos(b))
vxExp = np.multiply(vExp[1:-1], np.cos(bExp))

t = np.zeros(Nx)
for i in range(1, Nx):
    mean = (vx[i] + vx[i - 1]) / 2
    t[i] = t[i - 1] + 0.001 / mean


tidsutvikling = plt.figure("y(x)", figsize=(12, 6))
plt.plot(t, x)
plt.plot(tExp, xExp)
plt.title("Tidsutvikling")
plt.xlabel("$t$ (s)", fontsize=20)
plt.ylabel("$x$ (m)", fontsize=20)
plt.grid()
plt.show()


k = np.divide(d2y, (1 + dy ** 2) ** (3 / 2))
kEpx = np.divide(d2yDiff, (1 + dyDiff[1:-1] ** 2) ** (3 / 2))

krumning = plt.figure("y(x)", figsize=(12, 6))
plt.plot(x, k)
plt.plot(xExp[2:-2], kEpx)
plt.title("Krumning")
plt.xlabel("$x$ (m)", fontsize=20)
plt.ylabel("$r$ (m)", fontsize=20)
plt.ylim(-10, 10)
plt.grid()
plt.show()

aExp = np.multiply(vExp[2:-2], vExp[2:-2], kEpx)
a = np.multiply(v, v, k)

nExp = np.array(0.031 * (9.81 * np.cos(bExp[1:-1]) + aExp))
n = np.array(0.031 * (9.81 * np.cos(b) + a))

normalkraft = plt.figure("y(x)", figsize=(12, 6))
plt.plot(x, n)
plt.plot(xExp[2:-2], nExp)
plt.title("Normalkraft")
plt.xlabel("$x$ (m)", fontsize=20)
plt.ylabel("$N/MG$", fontsize=20)
plt.grid()
plt.show()

# Figurer kan lagres i det formatet du foretrekker:
baneform.savefig("baneform.png", bbox_inches="tight")
kulefart.savefig("kulefart.png", bbox_inches="tight")
helling.savefig("helling.png", bbox_inches="tight")
tidsutvikling.savefig("tidsutvikling.png", bbox_inches="tight")
krumning.savefig("krumning.png", bbox_inches="tight")
normalkraft.savefig("normalkraft.png", bbox_inches="tight")

# baneform.savefig("baneform.png", bbox_inches='tight')
# baneform.savefig("baneform.eps", bbox_inches='tight')


print("Festepunkthøyder (m)", yfast)
print("Banens høyeste punkt (m)", np.max(y))
print("Kulas sluttfart er (m/s)", v[-1])

