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
from cmath import cos
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

# Horisontal avstand mellom festepunktene er 0.200 m
h = 0.200
xfast=np.asarray([0,h,2*h,3*h,4*h,5*h,6*h,7*h])


# Vi begrenser starthøyden (og samtidig den maksimale høyden) til
# å ligge mellom 250 og 300 mm
ymax = 300
# yfast: tabell med 8 heltall mellom 50 og 300 (mm); representerer
# høyden i de 8 festepunktene
yfast = np.array([0.297, 0.237, 0.198, 0.165, 0.125, 0.158, 0.171, 0.103])


# inttan: tabell med 7 verdier for (yfast[n+1]-yfast[n])/h (n=0..7); dvs
# banens stigningstall beregnet med utgangspunkt i de 8 festepunktene.
inttan = np.diff(yfast)/h
attempts=1




# Når programmet her har avsluttet while-løkka, betyr det at
# tallverdiene i tabellen yfast vil resultere i en tilfredsstillende bane. 

#Programmet beregner deretter de 7 tredjegradspolynomene, et
#for hvert intervall mellom to nabofestepunkter.


#Med scipy.interpolate-funksjonen CubicSpline:
cs = CubicSpline(xfast, yfast, bc_type='natural')

xmin = 0.000
xmax = 1.401
dx = 0.001

x = np.arange(xmin, xmax, dx)   

#funksjonen arange returnerer verdier paa det "halvaapne" intervallet
#[xmin,xmax), dvs slik at xmin er med mens xmax ikke er med. Her blir
#dermed x[0]=xmin=0.000, x[1]=xmin+1*dx=0.001, ..., x[1400]=xmax-dx=1.400, 
#dvs x blir en tabell med 1401 elementer

Nx = len(x)
y = cs(x)       #y=tabell med 1401 verdier for y(x)
dy = cs(x,1)    #dy=tabell med 1401 verdier for y'(x)
d2y = cs(x,2)   #d2y=tabell med 1401 verdier for y''(x)

#Eksempel: Plotter banens form y(x)
baneform = plt.figure('y(x)',figsize=(12,6))
plt.plot(x,y, xfast, yfast, "*")
plt.title('Banens form')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$y(x)$ (m)',fontsize=20)
plt.ylim(0.0,0.40)
plt.grid()
plt.show()

v = np.sqrt((2 * 9.81 * (y[0] - y))/(1 + (2/5)))

kulefart = plt.figure('y(x)',figsize=(12,6))
plt.plot(x,v)
plt.title('Kulas fart')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$v(x)$ (m/s)',fontsize=20)
plt.grid()
plt.show()

# k = np.(d2y / (1 + dy**2)**(3/2))

b = np.arctan(dy)

helling = plt.figure('y(x)',figsize=(12,6))
plt.plot(x,b)
plt.title('Helningsvinkel')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$vinekl$ (rad)',fontsize=20)
plt.grid()
plt.show()


vx = np.multiply(v, np.cos(b))


t = np.zeros(Nx)
for i in range(1, Nx):
   mean = (vx[i] + vx[i- 1])/2
   t[i] = t[i - 1] + 0.001 / mean


tidsutvikling = plt.figure('y(x)',figsize=(12,6))
plt.plot(t, x)
plt.title('Tidsutvikling')
plt.xlabel('$t$ (s)',fontsize=20)
plt.ylabel('$x$ (m)',fontsize=20)
plt.grid()
plt.show()


k = np.divide(d2y, (1 + dy**2)**(3/2))

krumning = plt.figure('y(x)',figsize=(12,6))
plt.plot(x, k)
plt.title('Krumning')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$r$ (m)',fontsize=20)
plt.grid()
plt.show()


a = np.multiply(v, v, k)
print(len(a))

n = np.array( 0.031*(9.81 * np.cos(b) + a))

normalkraft = plt.figure('y(x)',figsize=(12,6))
plt.plot(x, n)
plt.title('Normalkraft')
plt.xlabel('$x$ (m)',fontsize=20)
plt.ylabel('$N/MG$',fontsize=20)
plt.grid()
plt.show()

#Figurer kan lagres i det formatet du foretrekker:
baneform.savefig("baneform.png", bbox_inches='tight')
kulefart.savefig("kulefart.png", bbox_inches='tight')
helling.savefig("helling.png", bbox_inches='tight')
tidsutvikling.savefig("tidsutvikling.png", bbox_inches='tight')
krumning.savefig("krumning.png", bbox_inches='tight')
normalkraft.savefig("normalkraft.png", bbox_inches='tight')

#baneform.savefig("baneform.png", bbox_inches='tight')
#baneform.savefig("baneform.eps", bbox_inches='tight')


print('Antall forsøk',attempts)
print('Festepunkthøyder (m)',yfast)
print('Banens høyeste punkt (m)',np.max(y))

print('NB: SKRIV NED festepunkthøydene når du/dere er fornøyd med banen.')
print('Eller kjør programmet på nytt inntil en attraktiv baneform vises.')