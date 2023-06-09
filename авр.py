from scipy.fft import fft, ifft
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import numpy as np
import random as r


def mapping(values_x, a, b): 
    y = []
    for i in range(len(values_x)):
        y.append(a*values_x[i] + b)
    return y

def correlogram(N3):
    c = (ifft(np.abs(X)**2).real/N)[0:N3]
    return c

def smoothed_periodogram(N3, c, k):    
    Cm = []
    Wm = []

    for j in range(N3):
        Wm = np.append(Wm, (1 -2*a) + 2*a*np.cos(np.pi*j/N3))
        Cm = np.append(Cm, c[j]*Wm[j])

    for i in range(N2-len(c)):
        Cm = np.append(Cm, 0)

    Cm_fft = fft(Cm)
    Dj = (2*Cm_fft.real/N3 - Cm_fft[0]/N3)[:N1+1]

    axes[k][1].plot(nu, Dj, color = "crimson", linewidth = 1.3,)
    axes[k][1].set_title(f"Сглаженная периодограмма (N* = N/{int(N/N3)}, a = 0.25)")


# дано:
N = 230; q = 0.01
dt = 1              # шаг временной дискретизации
dnu = 2/N/dt        # ширина главного лепестка спектрального окна
nu_c = 1/2/dt       # частота Найквиста

t = [ k*dt for k in range(N) ]

X1 = 9.0; A1 = 1; nu1 = 0.1
phi1 = 0; gamma = 0.50
alpha = 0.1; beta = 0.05
sigma_n = np.sqrt((A1**2)/(2*gamma))
eps = [r.normalvariate(0,1) for i in range(N)]



# создадим модельный ряд уже с 2 частотами, nu1 = 0.1,
# nu2 --- отличается от nu1 на четное количество частоты Найквиста:
x = [alpha + beta*t[k] + A1*np.cos(2*np.pi*(4*nu_c + nu1)*t[k]-phi1) + \
     A1*np.cos(2*np.pi*nu1*t[k]-phi1) + sigma_n*eps[k] for k in range(N)]



# шаг временной дискретизации (для того чтобы продемонстрировать
#  aliasing  сделаем его поменьше)
dt1 = 1/4
t1 = [ k*dt1 for k in range(N) ]


# создадим модельный ряд c измененным dt, c то есть частота найквиста для 
# такого ряда больше, следовательно:
x1 = [alpha + beta*t1[k] + A1*np.cos(2*np.pi*(2*nu_c + nu1)*t1[k]-phi1) +\
       A1*np.cos(2*np.pi*nu1*t1[k]-phi1) + sigma_n*eps[k] for k in range(N)]



# + A1*np.cos(2*np.pi*(4*nu_c + nu1)*t[k]-phi1)
# x = [alpha + beta*t[k] + A1*np.cos(2*np.pi*nu1*t[k]-phi1) + sigma_n*eps[k] for k in range(N)]


print('delta nu = ', dnu)
print(' частота Найквиста = ', nu_c)

# plt.plot(t, x)
# plt.plot(t, x4)



# исключение тренда и центрирование ряда
popt, _ = curve_fit(mapping, t, x)
a, b = popt

y = mapping(t, a, b)
y1 = mapping(t1, a, b)
x_new = [x[i] - y[i] for i in range(N)]
x1_new = [x1[i] - y1[i] for i in range(N)]
print(len(x1_new))
# x1_new = x1

#5 оценивание дисперсии
sigma2 = 0
for i in range(len(x_new)):
    sigma2 = sigma2 + x_new[i]**2
sigma2 = sigma2/(N-1)

print(sigma2)



# 4 вычисление периодограммы
N1 = 512
N2 = 2*N1

for i in range(N2-N):
    x_new = np.append(x_new, 0)
    x1_new= np.append(x1_new, 0)
    # x2= np.append(x2, 0)
print(len(x1_new))

X = fft(x_new)
X1 = fft(x1_new)
# X2 = fft(x2)
print(len(X))

D = []
# D1 = []
# D2 = []
nu = []
nu1 = []
dnu = 1/(N2*dt)
dnu1 = 1/(N2*dt1)
for j in range(N1 + 1):
    D = np.append(D, (np.abs(X[j])/N)**2)
    # D1 = np.append(D1, (np.abs(X1[j])/N)**2)
    # D2 = np.append(D2, (np.abs(X2[j])/N)**2)
    nu = np.append(nu, j*dnu)
    # nu1 = np.append(nu1, j*dnu1)

q = sigma2*X1/N


# 7 коррелограмма
N3 = int(N/10)
N4 = int(N/2)
a = 0.25
c1 = correlogram(N3)
c2 = correlogram(N4)

# строим графики:

fig, axes = plt.subplots(3, 2, figsize = (15, 15))
axes[0][0].plot(t, x, color = "crimson", linewidth = 1.3)
axes[0][0].set_title("График модельного ряда")
axes[1][0].plot(t, x_new[0:230], color = "crimson", linewidth = 1.3,)
axes[1][0].set_title("График центрированного ряда")
axes[2][0].plot(t[:N4], c2, color = "crimson", linewidth = 1.3,)
axes[2][0].set_title("График смещенной функции коррелограммы")
axes[2][0].set_xlabel("Time, s")

axes[0][1].set_title("Периодограмма Шустера модельного ряда и 99%-й порог обнаружения сигнала в шумах")
axes[0][1].plot(nu, D, color = "crimson", linewidth = 1.3,)
# axes[0][1].plot(nu, D1[:513], linewidth = 1.3,)

# axes[0][1].axhline(y = q)
axes[2][1].set_xlabel("Frequency, 1/s")

# 9-11 взвешенная коррелограмма и сглаженная периодограмма
smoothed_periodogram(N3, c1, 1)
smoothed_periodogram(N4, c2, 2)



# axes[0].plot(nu, D, color = "crimson", linewidth = 1.3,)
# plt.axhline(y = q)
# axes[1].plot(nu1, D1[:513], linewidth = 1.3, )
# axes[1].set_xlabel("Frequency, 1/s")




plt.savefig("image1.png")
plt.show()