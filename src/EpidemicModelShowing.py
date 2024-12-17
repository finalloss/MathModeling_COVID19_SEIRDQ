import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

#总人数
N = 10000000
#传染率
E_Infectious_Rate = 0.6
I_Infectious_Rate = 0.7
#致死率
Lethallty_Rate = 0.1

#潜伏期
Incubation = 14
#无症状感染者时空接触平均人群数
E_Contact = 10
#感染者（已发病）接触平均人群数
I_Contact = 2
#治疗率
Recover_Rate = 0.2
#抗体持续时间
Antibody_duration = 90
#疫苗覆盖率
Vaccine_Coverage = 0.4
day = 0
#隔离人群比
Quaranitine_Rate = 0.8   
#隔离自愈概率
Self_recovery = 0.1

def SEIR(t, SEIR):
    #易感人群，无症状感染者，感染者（已发病），康复者， 病死患者，隔离人员
    Susceptible, Exposed, Infectious, Recovered = SEIR
    dSusceptible_today = - E_Infectious_Rate * Exposed * Susceptible/N + 1/Antibody_duration * Recovered
    dExposed_today = E_Infectious_Rate * Exposed * Susceptible/N - \
                    1/Incubation * Exposed
    dInfectious_today = Exposed * 1/Incubation - Infectious * Recover_Rate 
    dRecovered_today = Infectious * Recover_Rate - 1/Antibody_duration * Recovered
    return np.array([dSusceptible_today,dExposed_today,dInfectious_today,dRecovered_today])


def SEIRDQ(t, SEIRDQ):
    #易感人群，无症状感染者，感染者（已发病），康复者， 病死患者，隔离人员
    Susceptible, Exposed, Infectious, Recovered, Death, Quaranitine = SEIRDQ
    global I_Infectious_Rate,E_Infectious_Rate,Lethallty_Rate
    day = int(t)
    dSusceptible_today = - E_Contact * E_Infectious_Rate * Exposed * Susceptible/(N-Death) - \
                        I_Contact * I_Infectious_Rate * Exposed * Susceptible/(N-Death) + \
                        1/Antibody_duration * Recovered     
    dExposed_today = E_Contact * E_Infectious_Rate * Exposed * Susceptible/(N-Death) + \
                    I_Contact * I_Infectious_Rate*Exposed*Susceptible/(N-Death) - \
                    1/Incubation * Exposed - Exposed * Quaranitine_Rate
    dInfectious_today = Exposed * 1/Incubation - Infectious * Recover_Rate - Infectious * Lethallty_Rate - Infectious * Quaranitine_Rate
    dRecovered_today = Infectious * Recover_Rate - 1/Antibody_duration * Recovered + Self_recovery * Quaranitine 
    dDeath_today = Infectious * Lethallty_Rate
    dQuaranitine_today = (Exposed + Infectious) * Quaranitine_Rate - Self_recovery * Quaranitine 

    return np.array([dSusceptible_today,dExposed_today,dInfectious_today,dRecovered_today,dDeath_today,dQuaranitine_today])

""" tmax, n = 40,40
Susceptible0, Exposed0, Infectious0, Recovered0, Death0, Quaranitine0 = (1-Vaccine_Coverage)*N-1,1,0,Vaccine_Coverage*N,0,0
t = np.linspace(0, tmax, n)
ht = solve_ivp(SEIRDQ, [0, tmax], [Susceptible0, Exposed0, Infectious0, Recovered0, Death0, Quaranitine0], t_eval = t)

Susceptible = ht.y[0]
Exposed = ht.y[1]
Infectious = ht.y[2]
Recovered = ht.y[3]
Death = ht.y[4]
Quaranitine = ht.y[5]

plt.figure()
plt.plot(t, Susceptible, 'blue', label = 'Susceptible')
plt.plot(t, Exposed, 'yellow', label = 'Exposed')
plt.plot(t, Recovered, 'green', label = 'Recovered')
plt.plot(t, Infectious, 'r-', label = 'Infectious')
plt.plot(t, Death, 'black', label = 'Death')
plt.plot(t, Quaranitine, 'purple', label = 'Quaranitine')
plt.title("Time evolution SEIRDQ")

plt.legend(loc=1)
# 添加模型参数
plt.savefig(f"./img/SEIRDQ Model with params Vaccine_Coverage={Vaccine_Coverage},day={tmax},Antibody_duration={Antibody_duration}.png")
plt.show() """


tmax, n = 90, 90
Susceptible0, Exposed0, Infectious0, Recovered0 = (1-Vaccine_Coverage)*N-1,1,0,Vaccine_Coverage*N
t = np.linspace(0, tmax, n)
ht = solve_ivp(SEIR, [0, tmax], [Susceptible0, Exposed0, Infectious0, Recovered0], t_eval = t)

Susceptible = ht.y[0]
Exposed = ht.y[1]
Infectious = ht.y[2]
Recovered = ht.y[3]

plt.figure()
plt.plot(t, Susceptible, 'blue', label = 'Susceptible')
plt.plot(t, Exposed, 'yellow', label = 'Exposed')
plt.plot(t, Recovered, 'green', label = 'Recovered')
plt.plot(t, Infectious, 'r-', label = 'Infectious')
plt.title("Time evolution SEIR")
plt.legend(loc=1)
plt.savefig(f"./img/SEIR Model with params Vaccine_Coverage={Vaccine_Coverage},day={tmax},Antibody_duration={Antibody_duration}.png")
plt.show()
plt.close()
