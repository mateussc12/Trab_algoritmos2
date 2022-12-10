import funcoes as fc
from sympy import *
from matplotlib import pyplot as plt
import numpy as np

# -------------------------------------------------------------------- TAREFA 1
fc.Grafico_FMP()
fc.Grafico_FDA()

# -------------------------------------------------------------------- TAREFA 2

P = [1 / 5, 1 / 2, 4 / 5]
n = 1000
fc.meu_geornd(P, n)

# -------------------------------------------------------------------- TAREFA 3 e TAREFA 4

x, k, lamba, fa, fb, fm = symbols('x k gama f(a) f(b) f(m)')
init_printing(use_unicode=True)

plt.style.use('ggplot')

p_weibull = (k / lamba) * pow((x / lamba), k - 1) * exp(-1 * pow(x / lamba, k))

f_weibull = 1 - exp(-1 * pow(x / lamba, k))

lambas = [1, 1, 1, 1, 2, 2]
ks = [0.5, 1, 1.5, 5, 1, 5]
xs0 = np.linspace(0, 2.5, 100)

fs_weibull = []
ps_weibull = []
for i in range(len(lambas)):
    fs_weibull.append(f_weibull.subs([(lamba, lambas[i]), (k, ks[i])]))
    ps_weibull.append(p_weibull.subs([(lamba, lambas[i]), (k, ks[i])]))

f1_weibull = []
f2_weibull = []
f3_weibull = []
f4_weibull = []
f5_weibull = []
f6_weibull = []

p1_weibull = []
p2_weibull = []
p3_weibull = []
p4_weibull = []
p5_weibull = []
p6_weibull = []

for i in xs0:
    f1_weibull.append(fs_weibull[0].subs(x, i))
    f2_weibull.append(fs_weibull[1].subs(x, i))
    f3_weibull.append(fs_weibull[2].subs(x, i))
    f4_weibull.append(fs_weibull[3].subs(x, i))
    f5_weibull.append(fs_weibull[4].subs(x, i))
    f6_weibull.append(fs_weibull[5].subs(x, i))

    p1_weibull.append(ps_weibull[0].subs(x, i))
    p2_weibull.append(ps_weibull[1].subs(x, i))
    p3_weibull.append(ps_weibull[2].subs(x, i))
    p4_weibull.append(ps_weibull[3].subs(x, i))
    p5_weibull.append(ps_weibull[4].subs(x, i))
    p6_weibull.append(ps_weibull[5].subs(x, i))

func, f = plt.subplots()
f.plot(xs0, f1_weibull, label="lambda = 1, k = 0.5")
f.plot(xs0, f2_weibull, label="lambda = 1, k = 1")
f.plot(xs0, f3_weibull, label="lambda = 1, k = 1.5")
f.plot(xs0, f4_weibull, label="lambda = 1, k = 5")
f.plot(xs0, f5_weibull, label="lambda = 2, k = 1")
f.plot(xs0, f6_weibull, label="lambda = 2, k = 5")

f.legend()
f.set_title('FDA da Distribuição Weibull')

prob, p = plt.subplots()

p.plot(xs0, p1_weibull, label="lambda = 1, k = 0.5")
p.plot(xs0, p2_weibull, label="lambda = 1, k = 1")
p.plot(xs0, p3_weibull, label="lambda = 1, k = 1.5")
p.plot(xs0, p4_weibull, label="lambda = 1, k = 5")
p.plot(xs0, p5_weibull, label="lambda = 2, k = 1")
p.plot(xs0, p6_weibull, label="lambda = 2, k = 5")

p.legend()
p.set_title('Função Densidade de Probabilidade da Distribuição Wiebull')


# comparação integração
a, b = 0, 1
# Valor analítico
v_anali = abs(float(fs_weibull[3].subs(x, a) - fs_weibull[3].subs(x, b)))

num_traps = [1, 4, 6, 10]
num_intervalos = [2, 4, 6, 10]
valores_trap = []
valores_sims = []

for i in range(len(num_traps)):
    valores_trap.append(float(fc.regra_trap_composto(a, b, num_traps[i] + 1, ps_weibull[3])[1]))
    valores_sims.append(float(fc.regra_de_simpson_um_terco(a, b, num_intervalos[i] + 1, ps_weibull[3])[1]))

print('')
print('')
print('|================================= TABELA DE RESULTADOS ==============================|')
print(f'| Integral de {ps_weibull[3]} com os intervalos indo de a = {a} e b = {b}               |')
print('|-------------------------------------------------------------------------------------|')
print(f'| Valor exato em x ={b} : {v_anali}                                            |')
print('|--------------|----------------------------------|-----------------------------------|')
print('|   Método     |            Valor aproximado      |              Erro                 | ')
print('|--------------|----------------------------------|-----------------------------------|')
print(f'| I_TR_1       |       {valores_trap[0]:.16f}         |        {abs(v_anali - valores_trap[0]):.16f}         |')
print(f'| I_TR_4       |       {valores_trap[1]:.16f}         |        {abs(v_anali - valores_trap[1]):.16f}         |')
print(f'| I_TR_6       |       {valores_trap[2]:.16f}         |        {abs(v_anali - valores_trap[2]):.16f}         |')
print(f'| I_TR_10      |       {valores_trap[3]:.16f}         |        {abs(v_anali - valores_trap[3]):.16f}         |')
print(f'| I_SR_2       |       {valores_sims[0]:.16f}         |        {abs(v_anali - valores_sims[0]):.16f}         |')
print(f'| I_TR_4       |       {valores_sims[1]:.16f}         |        {abs(v_anali - valores_sims[1]):.16f}         |')
print(f'| I_TR_6       |       {valores_sims[2]:.16f}         |        {abs(v_anali - valores_sims[2]):.16f}         |')
print(f'| I_TR_10      |       {valores_sims[3]:.16f}         |        {abs(v_anali - valores_sims[3]):.16f}         |')
print('|--------------|----------------------------------|-----------------------------------|')

#gerar gráfico da integral da distribuição de weibull em relação a x

trap, t = plt.subplots()

n = 1
A = []
i = 0
start = 0
A_ant = 0

for j in range(0, n):
    for x2 in xs0[(j)*round(len(xs0)/n):(j+1)*round(len(xs0)/n)]:
       # print(p4_weibull[i])
        A.append(A_ant + (p4_weibull[(j)*round(len(xs0)/n)] + p4_weibull[i])*x2/2)
       # print(A[i])
        i+=1
      #  print(f'{i}---------')
    A_ant = A[len(A)-1]

t.plot(xs0, A, label='ITR_1 lambda=1 k=5')

t.legend()
t.set_title('Integração por trapezios da FDP de Weibull')

# -------------------------------------------------------------------- TAREFA 5
k_s = [5., 1., 0.5]
lamba_s = [1, 2, 5]

x_zes = np.arange(1, 100.)/50.

boll = True
i = 0

fc.weib(x_zes, lamba_s[i], k_s[i])

while boll:
    cont = len(k_s)

    if i < cont:
        fc.graf_comp(x_zes, lamba_s[i], k_s[i])
        
        i += 1
    else:
        boll = False

r = fc.reestima()
#print(r)

plt.tight_layout()
plt.show()
