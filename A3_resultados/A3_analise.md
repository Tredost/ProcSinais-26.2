# A3 · Fala com interferência de 5 kHz

O experimento usa os primeiros 8 segundos da gravação de Garth Comira
para *The Ashiel Mystery, capítulo 2*, exemplo `libri1` do corpus LibriSpeech.
A gravação foi convertida de OGG para WAV PCM, em mono, a 22.050 Hz.
No ensaio, a voz tem pico ajustado a 0,65 e recebe uma senoide de 5 kHz
com amplitude 0,24. A interferência adicionada permite medir a rejeição
separadamente por superposição. Voz e interferência passam juntas pelo filtro
na demonstração.

[Identificação da gravação](https://librosa.org/doc/main/recordings.html)
· [Arquivo de origem](https://librosa.org/data/audio/5703-47212-0000.ogg)
· [Corpus LibriSpeech, CC BY 4.0](https://www.openslr.org/12/)
· [Licença](https://creativecommons.org/licenses/by/4.0/)

## Uma máscara para os dois projetos

- De 0 a 1.500 Hz, desvio máximo de 1 dB em relação a 0 dB.
- De 4.000 Hz a Nyquist (11.025 Hz), atenuação mínima de 30 dB.
- O intervalo entre as duas bordas é a faixa de transição. Os cortes
  internos de projeto podem ser diferentes.

O FIR usa uma sinc com corte de projeto em 2.750 Hz, janela de Hamming e
coeficientes simétricos. Testo ordens pares crescentes com esse mesmo corte.
A ordem 22 rejeita 29,04 dB no pior ponto da banda e falha. A ordem 24
rejeita 33,73 dB e é a primeira aprovada nessa busca.

O IIR é Butterworth com transformação bilinear. Para cada ordem N, calculo
o menor corte que mantém a perda na borda passante em até 1 dB. Na fórmula
analógica de magnitude, a perda é
`10 log10(1 + (Ω/Ωc)^(2N))`, com `Ω = 2fs·tan(πf/fs)`.
Igualando a perda a 1 dB em 1.500 Hz, obtenho
`Ωc = Ωp/(10^(1/10) − 1)^(1/(2N))`. Subir o corte reduziria a rejeição.
A ordem 3 atinge só 22,37 dB. A ordem 4 usa corte de 1.765,35 Hz e atinge
31,75 dB. A busca usa 20.001 pontos em cada banda, incluindo as bordas.

## Resultado

| Medida | FIR Hamming | IIR Butterworth |
|---|---:|---:|
| Ordem | 24 | 4 |
| Multiplicações por amostra | 25 | 9 |
| Somas ou subtrações por amostra | 24 | 8 |
| Variação de magnitude na passagem, pico a pico | 0,19 dB | 1,00 dB |
| Pior atenuação na rejeição | 33,73 dB | 31,75 dB |
| Atraso de grupo na passagem | 12,00 amostras | 5,08–7,91 amostras |
| Maior desvio de fase em relação à reta ajustada na passagem | ≈ 0 rad | 0,129 rad |
| Ganho medido sobre o tom de 5 kHz | −60,26 dB | −42,12 dB |
| Raio do polo mais externo | 0 | 0,830 |

O custo conta as operações da forma direta em regime. Usar a simetria
do FIR para somar pares antes da multiplicação reduziria as 25
multiplicações para 13. As duas respostas foram medidas na mesma escala,
sem normalizar as saídas separadamente.

O FIR tem 24 polos na origem na representação em potências de z. Sua ROC
causal é `|z| > 0` e a resposta ao impulso termina após 25 amostras.
A simetria dos coeficientes torna a fase linear na banda passante, com atraso
fixo de 12 amostras (0,54 ms). Os zeros próximos aos arcos de alta frequência
reduzem a resposta nessas regiões. A rejeição maior em 5 kHz é pontual.

O IIR tem quatro zeros em `z = −1`. Esse ponto é Nyquist e a resposta ali
é zero. Em 5 kHz, que fica em outro ângulo no círculo, a atenuação é
parcial. Os dois pares de polos internos moldam a transição e fazem a
resposta ao impulso decair. A ROC causal `|z| > 0,830` (raio arredondado)
inclui o círculo unitário. Ambos são estáveis. O IIR apresenta atrasos
diferentes para componentes de frequências diferentes na banda passante.

**Escolha.** Usaria o FIR para preservar a forma relativa das componentes
da banda passante. Usaria o IIR se o orçamento de operações fosse
menor e a fase não linear fosse aceitável. A comparação vale para estes
projetos e para esta máscara. A ordem 24 é o mínimo na busca FIR com
Hamming, corte fixo e ordem par. O ensaio usa uma gravação e um único
tom injetado, sem teste de inteligibilidade ou medição em microcontrolador.

## Resultados para conferência

- [Forma de onda original](01_voz_original.png)
- [Polos e zeros](02_plano_z.png)
- [Magnitude, fase e atraso de grupo](03_respostas.png)
- [Antes e depois no tempo e no espectro](04_antes_depois.png)
- [Entrada com interferência](entrada_com_tom.wav) · [Saída FIR](saida_FIR.wav) · [Saída IIR](saida_IIR.wav)
- [Medidas dos dois filtros](A3_comparativo.csv) · [Ordens testadas](A3_busca_ordens.csv)
