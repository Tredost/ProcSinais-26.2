"""Coeficientes e aplicação de filtros causais com condições iniciais nulas."""
import numpy as np


def _validar_frequencias(fs, *frequencias):
    if not (0 < frequencias[0] and frequencias[-1] < fs / 2
            and all(a < b for a, b in zip(frequencias, frequencias[1:]))):
        raise ValueError("Use frequências crescentes entre zero e Nyquist.")


def _impulso_ideal(fc, fs, m):
    # np.sinc(v) = sen(pi*v)/(pi*v), incluindo o limite em v = 0.
    return (2 * fc / fs) * np.sinc(2 * fc / fs * m)


def _fir(fs, ordem, tipo, *bordas):
    _validar_frequencias(fs, *bordas)
    if ordem < 2 or ordem % 2:
        raise ValueError("O FIR usa ordem par e comprimento ímpar.")
    n = np.arange(ordem + 1)
    m = n - ordem // 2
    delta = (m == 0).astype(float)
    baixa = _impulso_ideal(bordas[0], fs, m)
    if tipo == "passa-baixa":
        ideal = baixa
        referencia = 0
    elif tipo == "passa-alta":
        ideal = delta - baixa
        referencia = fs / 2
    else:
        banda = _impulso_ideal(bordas[1], fs, m) - baixa
        ideal = banda if tipo == "passa-banda" else delta - banda
        referencia = np.mean(bordas) if tipo == "passa-banda" else 0
    janela = 0.54 - 0.46 * np.cos(2 * np.pi * n / ordem)
    b = ideal * janela
    b /= abs(resposta(b, [1], fs, [referencia])[0])
    return b, np.array([1.0])


def passa_baixa_fir(fs, fc, ordem):
    return _fir(fs, ordem, "passa-baixa", fc)


def passa_alta_fir(fs, fc, ordem):
    return _fir(fs, ordem, "passa-alta", fc)


def passa_banda_fir(fs, f1, f2, ordem):
    return _fir(fs, ordem, "passa-banda", f1, f2)


def rejeita_banda_fir(fs, f1, f2, ordem):
    return _fir(fs, ordem, "rejeita-banda", f1, f2)


def _pre_distorcer(f, fs):
    return 2 * fs * np.tan(np.pi * f / fs)


def _butterworth(fs, fc, ordem, passa_alta=False):
    _validar_frequencias(fs, fc)
    if ordem < 1:
        raise ValueError("A ordem do Butterworth deve ser positiva.")
    wc = _pre_distorcer(fc, fs)
    # Polos do protótipo analógico no semiplano esquerdo.
    k = np.arange(ordem)
    prototipo = np.exp(1j * np.pi * (2 * k + ordem + 1) / (2 * ordem))
    polos_s = wc * prototipo
    polos_z = (2 * fs + polos_s) / (2 * fs - polos_s)
    a = np.poly(polos_z).real
    zero = 1 if passa_alta else -1
    b = np.poly([zero] * ordem).astype(float)
    q = -1 if passa_alta else 1
    # No ponto de referência, H = B(q)/A(q) = 1, com q = z⁻¹.
    b *= np.polynomial.polynomial.polyval(q, a) / np.polynomial.polynomial.polyval(q, b)
    return b, a


def passa_baixa_iir(fs, fc, ordem):
    return _butterworth(fs, fc, ordem)


def passa_alta_iir(fs, fc, ordem):
    # A transformação s -> wc/s inverte o conjunto dos polos protótipo.
    # Como |p|=1 e o conjunto é conjugado, resultam os mesmos polos.
    return _butterworth(fs, fc, ordem, passa_alta=True)


def _biquad_banda(fs, f0, q, rejeita):
    _validar_frequencias(fs, f0)
    if q <= 0:
        raise ValueError("Q deve ser positivo.")
    c = np.cos(2 * np.pi * f0 / fs)
    alfa = np.sin(2 * np.pi * f0 / fs) / (2 * q)
    # Bilinear de um ressonador de ordem 2. A forma rejeita tem
    # zeros exp(±j*w0). Na forma passa, os zeros ficam em ±1.
    a = np.array([1 + alfa, -2 * c, 1 - alfa])
    b = np.array([1, -2 * c, 1]) if rejeita else np.array([alfa, 0, -alfa])
    return b / a[0], a / a[0]


def passa_banda_iir(fs, f0, q):
    return _biquad_banda(fs, f0, q, rejeita=False)


def rejeita_banda_iir(fs, f0, q):
    return _biquad_banda(fs, f0, q, rejeita=True)


def notch_60_hz(fs, q=20):
    return rejeita_banda_iir(fs, 60, q)


def resposta(b, a, fs, frequencias):
    """Avalia H no círculo unitário sem recorrer a uma função de projeto."""
    q = np.exp(-2j * np.pi * np.asarray(frequencias) / fs)
    return np.polynomial.polynomial.polyval(q, b) / np.polynomial.polynomial.polyval(q, a)


def polos_zeros(b, a):
    # Multiplicar H por z^ordem no numerador e no denominador
    # inclui os polos ou zeros na origem devidos aos atrasos.
    comprimento = max(len(b), len(a))
    zeros = np.roots(np.pad(b, (0, comprimento - len(b))))
    polos = np.roots(np.pad(a, (0, comprimento - len(a))))
    return polos, zeros


def filtrar(x, b, a):
    """Equação de diferenças. Os projetos entregam a[0] = 1."""
    entrada = np.asarray(x, dtype=float).tolist()
    b = list(map(float, b))
    realimentacao = list(map(float, a[1:]))
    saida = [0.0] * len(entrada)
    for n, valor in enumerate(entrada):
        soma = b[0] * valor
        for k in range(1, min(len(b), n + 1)):
            soma += b[k] * entrada[n - k]
        for k in range(1, min(len(realimentacao) + 1, n + 1)):
            soma -= realimentacao[k - 1] * saida[n - k]
        saida[n] = soma
    return np.asarray(saida)


def db(magnitude):
    return 20 * np.log10(np.maximum(np.abs(magnitude), 1e-12))
