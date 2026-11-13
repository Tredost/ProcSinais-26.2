"""Tema 9: compara dois passa-baixa sobre uma mesma gravação de fala."""
import argparse
import csv
import json
from pathlib import Path
import wave

import matplotlib.pyplot as plt
import numpy as np

import A2_filtros as f
from A2_banco_testes import ENTRADA, FIR, IIR, amplitude_tom, espectro, plano_z

ROOT = Path(__file__).resolve().parent
FP, FR = 1500, 4000
AP, AR = 1, 30
TOM = 5000


def medidas_filtro(b, a, fs):
    passagem = np.linspace(0, FP, 20001)
    rejeicao = np.linspace(FR, fs / 2, 20001)
    h = f.resposta(b, a, fs, passagem)
    ganho = f.db(h)
    fase = np.unwrap(np.angle(h))
    w = 2 * np.pi * passagem / fs
    atraso = -np.gradient(fase, w, edge_order=2)
    coef_reta = np.polyfit(w, fase, 1)
    residuo = fase - np.polyval(coef_reta, w)
    polos, _ = f.polos_zeros(b, a)
    return dict(
        ordem=max(len(b), len(a)) - 1,
        multiplicacoes=len(b) + len(a) - 1,
        somas=len(b) + len(a) - 2,
        desvio_max_db=float(max(abs(ganho))),
        ripple_pp_db=float(np.ptp(ganho)),
        rejeicao_min_db=float(-max(f.db(f.resposta(b, a, fs, rejeicao)))),
        atraso_min=float(min(atraso)),
        atraso_max=float(max(atraso)),
        residuo_fase_max_rad=float(max(abs(residuo))),
        raio_max=float(max(abs(polos), default=0)),
    )


def projetar(fs):
    historico = []
    fc_fir = (FP + FR) / 2
    for ordem in range(2, 160, 2):
        b_fir, a_fir = f.passa_baixa_fir(fs, fc_fir, ordem)
        m_fir = medidas_filtro(b_fir, a_fir, fs)
        historico.append(dict(familia="FIR", fc_hz=fc_fir, **m_fir))
        if atende(m_fir):
            break
    else:
        raise ValueError("Nenhuma ordem FIR testada atingiu a máscara.")

    # Butterworth: menor corte que permite perda de até AP em FP.
    # Qualquer corte maior reduziria a atenuação em FR para a mesma ordem.
    for ordem in range(1, 13):
        t = np.tan(np.pi * FP / fs) / (10 ** (AP / 10) - 1) ** (1 / (2 * ordem))
        fc_iir = fs / np.pi * np.arctan(t)
        b_iir, a_iir = f.passa_baixa_iir(fs, fc_iir, ordem)
        m_iir = medidas_filtro(b_iir, a_iir, fs)
        historico.append(dict(familia="IIR", fc_hz=fc_iir, **m_iir))
        if atende(m_iir):
            break
    else:
        raise ValueError("Nenhuma ordem Butterworth testada atingiu a máscara.")
    return (b_fir, a_fir), (b_iir, a_iir), historico


def atende(medidas):
    return medidas["desvio_max_db"] <= AP + 1e-8 and medidas["rejeicao_min_db"] >= AR


def ler_wav(caminho):
    with wave.open(str(caminho), "rb") as arquivo:
        if arquivo.getsampwidth() != 2:
            raise ValueError("Use WAV PCM de 16 bits.")
        fs = arquivo.getframerate()
        canais = arquivo.getnchannels()
        x = np.frombuffer(arquivo.readframes(8 * fs), dtype="<i2")
    x = x.reshape(-1, canais).mean(axis=1)
    if fs != 22050 or len(x) < 2 * fs or not np.any(x):
        raise ValueError("Use pelo menos dois segundos de voz em 22.050 Hz.")
    return fs, x / np.max(np.abs(x)) * .65


def gravar_wav(caminho, fs, x):
    assert np.max(np.abs(x)) < 1, "A comparação deve preservar a escala sem saturação."
    with wave.open(str(caminho), "wb") as arquivo:
        arquivo.setparams((1, 2, fs, 0, "NONE", "not compressed"))
        arquivo.writeframes(np.round(x * 32767).astype("<i2").tobytes())


def escrever_csv(caminho, linhas):
    with caminho.open("w", encoding="utf-8-sig", newline="") as arquivo:
        writer = csv.DictWriter(arquivo, fieldnames=list(linhas[0]), delimiter=";")
        writer.writeheader()
        writer.writerows(linhas)


def graficos(destino, fs, t, voz, sinais, filtros):
    fig, ax = plt.subplots(figsize=(10, 3), layout="constrained")
    ax.plot(t, voz, color=FIR, linewidth=.8)
    ax.set(title="Gravação original", xlabel="Tempo (s)", ylabel="Amplitude", xlim=(0, t[-1]))
    fig.savefig(destino / "01_voz_original.png", dpi=180)
    plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(8.5, 4.2), layout="constrained")
    for ax, nome, (b, a) in zip(axs, ("FIR", "IIR"), filtros):
        polos, zeros = f.polos_zeros(b, a)
        if nome == "IIR":
            zeros = np.full(len(b) - 1, -1.0)  # B(z) é proporcional a (z+1)^4.
        plano_z(ax, polos, zeros)
        ax.set_title(nome)
    fig.savefig(destino / "02_plano_z.png", dpi=180)
    plt.close(fig)

    freq = np.linspace(0, fs / 2, 1103)
    fp = np.linspace(0, FP, 151)
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4), layout="constrained")
    curvas = []
    for nome, cor, (b, a) in zip(("FIR", "IIR"), (FIR, IIR), filtros):
        h = f.resposta(b, a, fs, freq)
        fase = np.unwrap(np.angle(f.resposta(b, a, fs, fp)))
        atraso = -np.gradient(fase, 2 * np.pi * fp / fs, edge_order=2)
        mag = np.maximum(f.db(h), -75)
        axs[0].plot(freq / 1000, mag, color=cor, label=nome)
        axs[1].plot(fp / 1000, fase, color=cor)
        axs[2].plot(fp / 1000, atraso, color=cor)
        curvas.append(dict(nome=nome, mag=mag.tolist(), fase=fase.tolist(), atraso=atraso.tolist()))
    axs[0].axvspan(FP / 1000, FR / 1000, color="#edf0f2")
    axs[0].axhline(-AR, color=ENTRADA, linestyle=":", linewidth=1)
    axs[0].set(title="Magnitude", ylabel="dB", xlabel="kHz", xlim=(0, 8), ylim=(-75, 3))
    axs[0].legend(frameon=False)
    axs[1].set(title="Fase na passagem", ylabel="rad", xlabel="kHz")
    axs[2].set(title="Atraso de grupo na passagem", ylabel="amostras", xlabel="kHz")
    fig.savefig(destino / "03_respostas.png", dpi=180)
    plt.close(fig)

    # Seleciona automaticamente uma janela com fala forte no primeiro segundo.
    meio = int(np.argmax(np.abs(voz[fs // 4:fs]))) + fs // 4
    trecho = slice(meio - int(.003 * fs), meio + int(.005 * fs))
    tm = 1000 * (t[trecho] - t[trecho][0])
    freq_fft, _ = espectro(sinais[0], fs)
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.2), layout="constrained")
    demo = []
    for x, nome, cor in zip(sinais, ("Entrada com tom", "Saída FIR", "Saída IIR"), (ENTRADA, FIR, IIR)):
        axs[0].plot(tm, x[trecho], label=nome, color=cor, linewidth=1)
        _, amp = espectro(x, fs)
        mag = f.db(amp)
        axs[1].plot(freq_fft / 1000, mag, label=nome, color=cor, alpha=.8, linewidth=.8)
        # Envelope de pico em bins de 20 Hz mantém visível o pico estreito.
        passo = max(1, round(20 / (fs / len(x))))
        faixa = freq_fft < 8000
        indices = np.arange(0, np.count_nonzero(faixa), passo)
        pico = np.maximum.reduceat(mag[faixa], indices)
        centro = freq_fft[indices] / 1000
        demo.append(dict(nome=nome, tempo=x[trecho].tolist(), espectro=np.maximum(pico, -100).tolist()))
    axs[0].set(title="Trecho de fala com interferência", xlabel="Tempo no trecho (ms)", ylabel="Amplitude")
    axs[1].axvline(TOM / 1000, color=ENTRADA, linestyle=":", linewidth=1)
    axs[1].set(title="Espectro de amplitude", xlabel="Frequência (kHz)", ylabel="dB ref. 1",
               xlim=(0, 8), ylim=(-100, -8))
    axs[0].legend(loc="upper right", frameon=False, fontsize=8)
    fig.savefig(destino / "04_antes_depois.png", dpi=180)
    plt.close(fig)

    passo_voz = max(1, len(voz) // 420)
    inicios = np.arange(0, len(voz), passo_voz)
    vmax = np.maximum.reduceat(voz, inicios)
    vmin = np.minimum.reduceat(voz, inicios)
    return dict(voz_t=t[inicios].tolist(), voz_max=vmax.tolist(), voz_min=vmin.tolist(),
                freq_khz=(freq / 1000).tolist(), passagem_khz=(fp / 1000).tolist(),
                curvas=curvas, demo=demo, tempo_ms=tm.tolist(), fft_khz=centro.tolist())


def executar(arquivo, destino):
    destino.mkdir(parents=True, exist_ok=True)
    fs, voz = ler_wav(arquivo)
    t = np.arange(len(voz)) / fs
    interferencia = .24 * np.sin(2 * np.pi * TOM * t)
    entrada = voz + interferencia
    fir, iir, historico = projetar(fs)
    sinais = [entrada, f.filtrar(entrada, *fir), f.filtrar(entrada, *iir)]
    resultados = []
    coeficientes = {}
    for nome, (b, a), y in zip(("FIR", "IIR"), (fir, iir), sinais[1:]):
        # Isola a contribuição do tom por linearidade e descarta o transiente.
        tom_saida = f.filtrar(interferencia, b, a)
        ganho_tom = float(f.db(amplitude_tom(tom_saida[fs:], fs, TOM) /
                               amplitude_tom(interferencia[fs:], fs, TOM)))
        ganho_teorico = float(f.db(f.resposta(b, a, fs, [TOM])[0]))
        assert abs(ganho_tom - ganho_teorico) < .01
        m = medidas_filtro(b, a, fs)
        assert atende(m) and m["raio_max"] < 1 and np.isfinite(y).all()
        resultados.append(dict(familia=nome, **m, ganho_tom_db=ganho_tom))
        polos, zeros = f.polos_zeros(b, a)
        if nome == "IIR":
            zeros = np.full(len(b) - 1, -1.0)
        coeficientes[nome] = dict(b=b.tolist(), a=a.tolist(),
                                  polos=[[float(p.real), float(p.imag)] for p in polos],
                                  zeros=[[float(z.real), float(z.imag)] for z in zeros])
    ordens_escolhidas = {m["familia"]: m["ordem"] for m in resultados}
    for tentativa in historico:
        assert atende(tentativa) == (tentativa["ordem"] == ordens_escolhidas[tentativa["familia"]])
    escrever_csv(destino / "A3_comparativo.csv", resultados)
    escrever_csv(destino / "A3_busca_ordens.csv", historico)
    for nome, x in zip(("entrada_com_tom", "saida_FIR", "saida_IIR"), sinais):
        gravar_wav(destino / (nome + ".wav"), fs, x)
    dados = graficos(destino, fs, t, voz, sinais, (fir, iir))
    dados.update(medidas=resultados, coeficientes=coeficientes, historico=historico)
    (destino / "A3_dados.json").write_text(json.dumps(dados), encoding="utf-8")
    for m in resultados:
        print(m)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--entrada", type=Path, default=ROOT / "audio" / "voz_original.wav")
    parser.add_argument("--destino", type=Path, default=ROOT / "A3_resultados")
    args = parser.parse_args()
    executar(args.entrada, args.destino)
