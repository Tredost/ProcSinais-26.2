"""Testes numéricos e figuras dos nove filtros do A2."""
from pathlib import Path
import csv

import matplotlib.pyplot as plt
import numpy as np

import A2_filtros as f

ROOT = Path(__file__).resolve().parent
FIR = "#176c78"
IIR = "#a7422c"
ENTRADA = "#697485"


def plano_z(ax, polos, zeros):
    # Mesmo círculo e convenções de plot_z_plane do material da aula.
    angulo = np.linspace(0, 2 * np.pi, 401)
    ax.plot(np.cos(angulo), np.sin(angulo), "--", color=ENTRADA, linewidth=1)
    ax.axhline(0, color="#b3bac1", linewidth=.8)
    ax.axvline(0, color="#b3bac1", linewidth=.8)
    for pontos, marcador, cor, rotulo in [(zeros, "o", FIR, "zero"), (polos, "x", IIR, "polo")]:
        unicos, contagens = np.unique(np.round(pontos, 8), return_counts=True)
        ax.scatter(unicos.real, unicos.imag, marker=marcador,
                   facecolors="none" if marcador == "o" else cor,
                   edgecolors=cor if marcador == "o" else None, label=rotulo)
        for ponto, quantidade in zip(unicos, contagens):
            if quantidade > 1:
                ax.annotate(str(quantidade), (ponto.real, ponto.imag),
                            xytext=(6, 6), textcoords="offset points", color=cor)
    limite = min(1.85, max(1.22, 1.1 * max(np.abs(np.r_[polos, zeros]), default=0)))
    fora = np.count_nonzero((abs(zeros.real) > limite) | (abs(zeros.imag) > limite))
    if fora:
        ax.text(.02, .04, f"{fora} zero(s) fora da escala", transform=ax.transAxes,
                fontsize=8, color=ENTRADA)
    ax.set(xlabel="Re(z)", ylabel="Im(z)", xlim=(-limite, limite), ylim=(-limite, limite))
    ax.set_aspect("equal")
    ax.legend(loc="upper right", frameon=False)


def espectro(x, fs):
    janela = np.hanning(len(x))
    frequencias = np.fft.rfftfreq(len(x), 1 / fs)
    amplitude = 2 * np.abs(np.fft.rfft(x * janela)) / janela.sum()
    amplitude[0] /= 2
    if len(x) % 2 == 0:
        amplitude[-1] /= 2
    return frequencias, amplitude


def amplitude_tom(x, fs, hz):
    n = np.arange(len(x))
    return 2 * abs(np.mean(x * np.exp(-2j * np.pi * hz * n / fs)))


def casos():
    fs = 22050
    tons = (600, 2500, 5500)
    yield "01_pb_fir", "Passa-baixa FIR", fs, tons, (0,), f.passa_baixa_fir(fs, 2000, 48)
    yield "02_pb_iir", "Passa-baixa IIR", fs, tons, (0,), f.passa_baixa_iir(fs, 1400, 4)
    yield "03_pa_fir", "Passa-alta FIR", fs, tons, (2,), f.passa_alta_fir(fs, 3900, 48)
    yield "04_pa_iir", "Passa-alta IIR", fs, tons, (2,), f.passa_alta_iir(fs, 3900, 4)
    yield "05_bp_fir", "Passa-banda FIR", fs, tons, (1,), f.passa_banda_fir(fs, 1500, 3500, 48)
    yield "06_bp_iir", "Passa-banda IIR", fs, tons, (1,), f.passa_banda_iir(fs, 2500, 3)
    yield "07_rb_fir", "Rejeita-banda FIR", fs, tons, (0, 2), f.rejeita_banda_fir(fs, 1500, 3500, 48)
    yield "08_rb_iir", "Rejeita-banda IIR", fs, tons, (0, 2), f.rejeita_banda_iir(fs, 2500, 0.9)
    yield "09_notch", "Notch IIR de 60 Hz", 1000, (25, 60, 170), (0, 2), f.notch_60_hz(1000)


def relatorio(destino, linhas):
    cabecalho = """# A2 · Biblioteca de filtros digitais

A biblioteca calcula os coeficientes com NumPy. Matplotlib gera as figuras. A
função `filtrar` aplica a equação de diferenças, amostra a amostra, com repouso inicial.

Em todos os casos, a função de transferência é
`H(z) = (b₀ + b₁z⁻¹ + … + bₘz⁻ᵐ) / (1 + a₁z⁻¹ + … + aₙz⁻ⁿ)`.
A ordem é o maior grau em z⁻¹. Os vetores abaixo estão na ordem dos atrasos. O
CSV traz os coeficientes em precisão de trabalho, raízes, ordens e os ganhos
observados. O número ao lado de uma raiz no plano Z indica multiplicidade. O plano
Z amplia a vizinhança do círculo unitário. O CSV traz as raízes que ficam
fora do recorte e as aproximações numéricas de raízes múltiplas.

## Como os coeficientes são calculados

- **FIR**. A sinc representa um passa-baixa ideal. A janela de Hamming limita
  a resposta a uma quantidade finita de amostras. Para o passa-alta, subtraio
  a sinc do impulso central. Para o passa-banda, subtraio duas sinc. Para o
  rejeita-banda, subtraio esse passa-banda do impulso central.
- **IIR Butterworth**. Os polos do protótipo analógico ficam espaçados no
  semicírculo esquerdo. A frequência é pré-distorcida com
  Ωc = 2fs·tan(πfc/fs). A transformação bilinear leva os polos ao plano Z.
  Os zeros ficam em −1 no passa-baixa e em +1 no passa-alta.
- **IIR de banda**. Uso uma seção biquadrática com frequência central e Q.
  O passa-banda tem zeros em DC e Nyquist. O rejeita-banda tem zeros no
  círculo unitário no ângulo da frequência central. O notch é esse mesmo
  rejeita-banda com Q = 20 para concentrar a rejeição em torno de 60 Hz.

## Escolha da ordem

Os FIR do banco usam ordem 48 (49 coeficientes) para separar os tons de
teste com margem, com 24 amostras de atraso. Mantendo janela e bordas,
encurtar o filtro alarga a transição e reduz custo e atraso. A janela de
Hamming controla lóbulos laterais em troca de uma transição mais larga
que a da janela retangular de mesmo comprimento. Os Butterworth de
ordem 4 têm passagem monotônica. Aumentar a ordem acentua a transição,
mas custa mais operações e pode aumentar a variação de atraso.
Nos biquads, uma seção de ordem 2 já permite obter um par de zeros e
ajustar a largura pelo Q. Os testes verificam as frequências abaixo.

Cada sinal combina três senoides de amplitude 0,3 e ruído gaussiano de
desvio padrão 0,012 com semente 9. As medidas usam os últimos quatro
segundos, após três segundos de acomodação. A FFT usa janela de Hann e
a mesma escala para entrada e saída. O ganho é medido por projeção no
tom. O piso visual em dB evita log de zero. A fase fica indefinida nos
zeros da resposta.

| Filtro | Ordem | Ganho medido nos três tons, em dB |
|---|---:|---|
"""
    for linha in linhas:
        cabecalho += f'| {linha["nome"]} | {linha["ordem"]} | {linha["ganhos"]} |\n'
    cabecalho += """
## Comparação FIR e IIR sob a mesma especificação

A tabela de ordem, ripple, rejeição e custo obtida por execução está em
`../A3_resultados/A3_comparativo.csv`. O A3 aplica a mesma biblioteca
e compara os dois passa-baixa com uma única máscara de magnitude.

## Resultados por filtro
"""
    for linha in linhas:
        cabecalho += f"""
### {linha["nome"]}

- b = `{linha["b"]}`
- a = `{linha["a"]}`

![Plano Z, magnitude, fase, tempo e espectro]({linha["id"]}.png)
"""
    (destino / "A2_relatorio.md").write_text(cabecalho, encoding="utf-8")


def executar(destino=ROOT / "A2_resultados"):
    destino.mkdir(parents=True, exist_ok=True)
    linhas = []
    for codigo, nome, fs, tons, passam, (b, a) in casos():
        t = np.arange(7 * fs) / fs
        x = sum(0.3 * np.sin(2 * np.pi * hz * t) for hz in tons)
        x += np.random.default_rng(9).normal(0, 0.012, t.size)
        y = f.filtrar(x, b, a)
        intervalo = slice(3 * fs, None)
        ganhos = [float(f.db(amplitude_tom(y[intervalo], fs, hz) /
                             amplitude_tom(x[intervalo], fs, hz))) for hz in tons]
        polos, zeros = f.polos_zeros(b, a)
        assert np.all(np.abs(polos) < 1) and np.isfinite(y).all(), nome
        for i, ganho in enumerate(ganhos):
            assert abs(ganho) < 2 if i in passam else ganho < -10, (nome, tons[i], ganho)
        if codigo == "09_notch":
            assert ganhos[1] < -50, ganhos
        frequencias = np.linspace(0, fs / 2, 4097)
        h = f.resposta(b, a, fs, frequencias)
        fase = np.unwrap(np.angle(h))
        fase[np.abs(h) < 1e-6] = np.nan

        fig, axs = plt.subplots(2, 3, figsize=(12.5, 7), layout="constrained")
        fig.suptitle(nome, fontsize=18)
        plano_z(axs[0, 0], polos, zeros)
        axs[0, 0].set_title("Plano Z")
        axs[0, 1].plot(frequencias, np.maximum(f.db(h), -90), color=FIR)
        axs[0, 1].set(title="Magnitude", xlabel="Hz", ylabel="dB", ylim=(-93, 4))
        axs[0, 2].plot(frequencias, fase, color=IIR)
        axs[0, 2].set(title="Fase desenrolada", xlabel="Hz", ylabel="rad")
        duracao = .18 if fs == 1000 else .013
        janela = slice(3 * fs, 3 * fs + int(duracao * fs))
        axs[1, 0].plot(t[janela] - 3, x[janela], label="entrada", color=ENTRADA)
        axs[1, 0].plot(t[janela] - 3, y[janela], label="saída", color=IIR)
        axs[1, 0].set(title="Antes e depois no tempo", xlabel="s", ylabel="amplitude")
        axs[1, 0].legend(frameon=False)
        for sinal, legenda, cor in [(x, "entrada", ENTRADA), (y, "saída", IIR)]:
            freq, amp = espectro(sinal[intervalo], fs)
            axs[1, 1].plot(freq, f.db(amp), color=cor, alpha=.85, label=legenda)
        axs[1, 1].set(title="Espectro antes e depois", xlabel="Hz", ylabel="dB ref. 1",
                      ylim=(-110, 0), xlim=(0, min(fs / 2, max(tons) * 1.3)))
        axs[1, 1].legend(frameon=False)
        axs[1, 2].axis("off")
        ordem = max(len(b), len(a)) - 1
        axs[1, 2].text(.05, .75, f"Ordem {ordem}\n\n" + "\n".join(
            f"{hz:g} Hz: {ganho:+.2f} dB" for hz, ganho in zip(tons, ganhos)),
            va="top", fontsize=14, transform=axs[1, 2].transAxes)
        fig.savefig(destino / f"{codigo}.png", dpi=150)
        plt.close(fig)
        linhas.append(dict(id=codigo, nome=nome, fs=fs, ordem=ordem, b=b.tolist(), a=a.tolist(),
                           polos=polos.tolist(), zeros=zeros.tolist(),
                           ganhos=" / ".join(f"{hz:g} Hz {g:+.2f}" for hz, g in zip(tons, ganhos))))
        print(nome, linhas[-1]["ganhos"], flush=True)
    with (destino / "A2_medidas.csv").open("w", encoding="utf-8-sig", newline="") as arquivo:
        writer = csv.DictWriter(arquivo, fieldnames=list(linhas[0]), delimiter=";")
        writer.writeheader()
        writer.writerows(linhas)
    relatorio(destino, linhas)


if __name__ == "__main__":
    executar()
