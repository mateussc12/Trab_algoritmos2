from matplotlib import pyplot as plt
from sympy import *
import numpy as np
import random

x = symbols('x')
init_printing(use_unicode=True)

plt.style.use('ggplot')


def Grafico_FMP():
    Y = []
    repeticao = 10
    enis = []

    def calc_prob_insuss(p, num_repet):
        prob_y_ser_n = []
        prob = [1 / 5, 1 / 2, 4 / 5]

        for i in range(num_repet + 1):
            op_prob = ((1 - prob[p]) ** i) * prob[p]
            prob_y_ser_n.append(op_prob)

        return prob_y_ser_n

    for ques in range(3):
        Y.append(calc_prob_insuss(ques, repeticao))

    for i in range(repeticao + 1):
        enis.append(i)

    fmp1, fmp = plt.subplots()
    fmp.plot(enis, Y[0], label='P=0.2', linewidth=0.1, alpha=1, marker='s', color='blue')
    fmp.plot(enis, Y[1], label='P=0.5', linewidth=0.1, alpha=1, marker='>', color='red')
    fmp.plot(enis, Y[2], label='P=0.8', linewidth=0.1, alpha=1, marker='*', color='orange')
    fmp.legend(fontsize='medium')
    fmp.set_title("Distribuição Geometrica: FMP ")
    fmp.set_ylabel('P(Y=n)')
    fmp.set_xlabel('n')

    grafi, grafico = plt.subplots()
    grafico.set_yscale('log')

    grafico.plot(enis, Y[0], label='P=0.2', linewidth=0.1, alpha=1, marker='s', color='blue')
    grafico.plot(enis, Y[1], label='P=0.5', linewidth=0.1, alpha=1, marker='>', color='red')
    grafico.plot(enis, Y[2], label='P=0.8', linewidth=0.1, alpha=1, marker='*', color='orange')
    grafico.legend(fontsize='medium')
    grafico.set_title("Distribuição Geometrica: FMP com Ln P(Y=n)")
    grafico.set_ylabel('ln P(Y=n)')
    grafico.set_xlabel('n')


def Grafico_FDA():
    Y = []
    repeticao = 10
    enis = []

    def calc(p, num_repet):
        prob = []
        Pb = [1 / 5, 1 / 2, 4 / 5]
        acum = 0

        for i in range(num_repet + 1):
            op_prob = acum + Pb[p] * pow(1 - Pb[p], i)
            prob.append(op_prob)
            acum = op_prob

        return prob

    for ques in range(3):
        Y.append(calc(ques, repeticao))

    for i in range(repeticao + 1):
        enis.append(i)

    fmp1, fmp = plt.subplots()

    fmp.plot(enis, Y[0], label='P=0.2', linewidth=0.1, alpha=1, marker='s', color='blue')
    fmp.plot(enis, Y[1], label='P=0.5', linewidth=0.1, alpha=1, marker='>', color='red')
    fmp.plot(enis, Y[2], label='P=0.8', linewidth=0.1, alpha=1, marker='*', color='orange')
    fmp.legend(fontsize='medium')
    fmp.set_title("Distribuição Geometrica: FDA ")
    fmp.set_ylabel('P(Y=n)')
    fmp.set_xlabel('n')

    grafi, grafico = plt.subplots()
    grafico.set_yscale('log')

    grafico.plot(enis, Y[0], label='P=0.2', linewidth=0.1, alpha=1, marker='s', color='blue')
    grafico.plot(enis, Y[1], label='P=0.5', linewidth=0.1, alpha=1, marker='>', color='red')
    grafico.plot(enis, Y[2], label='P=0.8', linewidth=0.1, alpha=1, marker='*', color='orange')
    grafico.legend(fontsize='medium')
    grafico.set_title("Distribuição Geometrica: FDA com Ln P(Y=n)")
    grafico.set_ylabel('ln P(Y=n)')
    grafico.set_xlabel('n')


def regra_trap_composto(a, b, num_pontos, func):
    num_trap = num_pontos - 1
    h = abs(b - a) / num_trap

    pontos = [a]
    for i in range(1, num_pontos):
        pontos.append(pontos[i - 1] + h)

    somatorio = 0
    for i in range(num_trap):
        somatorio += (func.subs(x, pontos[i]) + func.subs(x, pontos[i + 1]))

    regra_trap_composta = (h / 2) * somatorio

    valores = [pontos, regra_trap_composta]

    return valores


def regra_de_simpson_um_terco(a, b, num_pontos, func):
    num_trap = num_pontos - 1
    h = abs(b - a) / num_trap

    pontos = [a]
    for i in range(1, num_pontos):
        pontos.append(pontos[i - 1] + h)

    somatorio_1 = 0
    for i in range(1, num_pontos, 2):
        somatorio_1 += func.subs(x, pontos[i])

    somatorio_2 = 0
    for i in range(2, num_pontos - 1, 2):
        somatorio_2 += func.subs(x, pontos[i])

    re_simpson = ((b - a) / (3 * num_trap)) * (func.subs(x, a) + 4 * somatorio_1 + 2 * somatorio_2 + func.subs(x, b))

    valores = [pontos, re_simpson]

    return valores


def meu_geornd(P, n):
    P_aprox = []
    erro = []
    for i in range(3):
        if 0 < P[i] < 1:
            k = np.floor(- np.random.exponential(1, n) / np.log(1 - P[i]))
        elif P[i] == 0:
            k = np.inf
        elif P[i] == 1:
            k = 0
        elif P[i] < 0 or P[i] > 1:
            k = np.nan
        P_aprox.append(n / (n + sum(k)))
        erro.append(abs(P[i] - P_aprox[i]))
    print('')
    print('|---------------------------------- MEU GEORND -----------------------------------|')
    for i in range(3):
        print(f'|   Valor Real = {P[i]}       Valor Aproximado = {P_aprox[i]:.5f}          Erro = {erro[i]:.5f}     |')
    print('|---------------------------------------------------------------------------------|')


def weib(x, lamba, k):
    return (k / lamba) * (x / lamba) ** (k - 1) * np.exp(-(x / lamba) ** k)

def reestima():
    k = symbols('k')

    k0 = 1
    lamba = 5
    x1 = []
    
    random.seed(1234567890)
    
    '''for i  in range(1000):
        x1.append(weib(random.random(), lamba, k0))'''
    
    x = np.random.weibull(k0, 1000)
    
    '''lista = []
    
    for i in x:
        lista.append(float(i))
    
    x = np.array(lista)'''
    
    _1_phi = (np.mean((x**k)*np.log(x)))/(np.mean(x**k)) - np.mean(np.log(x))
    
    phi = 1/_1_phi
    
    k_j1 = phi.subs(k, k0)
    
    #print(k_j1)
    
    lamba_ = (np.mean(x**k_j1))**(1/k_j1)
    
    #print(lamba_)
    
    return [k_j1, lamba_]

def graf_comp(x_zes, lamba, k):
    
    fig6, compa = plt.subplots()
    fig6 = plt.figure(figsize=(100,100))

    count, bins, ignored = compa.hist(np.random.weibull(k, 1000), alpha=0.3, color='#F51720',
                                      label='Números gerados por método de Wielbull')

    scale = count.max() / weib(x_zes, lamba, k).max()

    compa.plot(x_zes, weib(x_zes, lamba, k) * scale, label=f'Função densidade de probabilidade de Wielbull com k = {k} e '
                                                         f'lambda = {lamba}', color='#F51720')

    compa.set_title("1000 números gerados por método de Wielbull comparados com a respectiva FDP")

    compa.legend(fontsize='xx-small', loc=0)
    