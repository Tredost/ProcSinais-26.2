# AP2 · Processamento de Sinais

**Ian Esteves Amaral Francisco · Tema 9 · Comparativo FIR e IIR**

| Entrega | Arquivo ou pasta |
|---|---|
| A1 · Transformada Z | [Caderno Jupyter](A1_Caderno_Transformada_Z.ipynb) com cálculos e plano Z já executado |
| A2 · Filtros digitais | [Biblioteca](A2_filtros.py), [banco de testes](A2_banco_testes.py) e [relatório com nove gráficos](A2_resultados/A2_relatorio.md) |
| A3 · Aplicação em fala | [Código](A3_comparar_fir_iir.py), [análise](A3_resultados/A3_analise.md), gravação de entrada, áudios filtrados, gráficos e medidas |
| Apresentação | [PowerPoint de oito slides](AP2_Tema_9_FIR_IIR.pptx) com roteiro nas notas do apresentador |

## Execução

Usei Python porque permite escrever a equação de diferenças e conferir
os resultados com gráficos no mesmo ambiente. NumPy faz as operações
numéricas e Matplotlib gera os gráficos.

Abra o terminal na pasta extraída e execute:

```bash
python -m pip install -r requirements.txt
python A2_banco_testes.py
python A3_comparar_fir_iir.py
```

O A1 pode ser aberto em Jupyter ou VS Code. Os gráficos e áudios já
gerados estão nas pastas de resultados. Para usar outra gravação,
o A3 aceita `--entrada caminho.wav` com PCM de 16 bits, 22.050 Hz e
duração mínima de dois segundos.

## Material de aula

- *Instruções AP2 · Módulo 5: Funções*, seções 3 a 7.
- *Apostila de Transformada Z · Da Teoria à Prática*, capítulos 2 a 8.
- `transformada_z.py`, diagrama de polos, zeros e círculo unitário.
- *Filtros Passa-Banda e Rejeita-Banda*, Doutornet Tecnologia, vol. 3, 2026.
- *Processamento de Sinais Aplicado a Jogos de Computador*, material de aula sobre Transformada Z, estabilidade e análise espectral.
