# Data Challenge 07 — Site Productive Capacity Classification

Project C531: 9 permanent plots, 43 dominant-height (Hd) measurements, ages 2.6–6.4 years,
Hd 11.1–26.3 m. Five methods build site-index curves from the same data, one folder each.

```
desafio-07-classificacao-sitio/
├── data/Sitio.xls
├── comum/00_funcoes_desafio07.R        # shared settings + functions (sourced by every script)
├── 01_curva_guia/            01_curva_guia.R            + outputs/
├── 02_diferenca_algebrica/   02_diferenca_algebrica.R   + outputs/
├── 03_atribuicao_preliminar/ 03_atribuicao_preliminar.R + outputs/   (needs output of 01)
├── 04_hammer/                04_hammer.R                + outputs/
├── 05_predicao_parametros/   05_predicao_parametros.R   + outputs/
├── 06_comparacao_metodos/    06_comparacao_metodos.R    + outputs/   (extra: report figure)
├── report/relatorio_desafio07.md
└── README.md
```

## Settings used by all five methods (defined once, in `comum/00_funcoes_desafio07.R`)

| Setting | Value | Why |
|---|---|---|
| Index age **Ii** | **4 years** | Inside the observed range (2.6–6.4 yr) and close to the middle of the data, so S is an interpolation, not an extrapolation. |
| Site-index classes **S** | **15, 17, 19, 21 m** (Hd at 4 yr) | The 9 plots have S between about 15.2 and 20.8 m at 4 yr. The PDF's example classes (14, 18, 22, 26) would put 26 m far above every plot, so the classes were centred on the data. Step of 2 m. |
| Base growth model | **Chapman-Richards**, Hd = b0·(1 − e^(−b1·t))^b2 | Used in Parts 1, 2, 3, 5 (and, with a class covariate, Part 4) so the methods are comparable. Its asymptote b0 has a biological meaning. |
| Curve ages | 2.5 to 6.5 years, every 0.25 | Observed range. |

**Why Chapman-Richards was kept (Part 1 comparison).** On the pooled data Chapman-Richards,
Backman and Schumacher are statistically equivalent (Syx 1.812, 1.814 and 1.809 m). The script
keeps Chapman-Richards unless another model has a Syx more than 2% lower, which none does.

## The five methods in one line each

1. **Guide curve** — one curve to all pooled points; S = Hd·g(Ii)/g(t). Anamorphic.
2. **Algebraic difference** — eliminate the site parameter between two ages of the same plot and fit
   the resulting equation to pairs of measurements. Two versions are fitted: isolating b0 (anamorphic)
   and isolating b2 (polymorphic). The polymorphic one fits better and is used.
3. **Preliminary assignment** — S from the guide curve, then S enters the model as a known covariate
   (Hd = f(t, S)); three ways of entering S are compared by AIC; plots are re-classified at the end.
4. **Hammer** — age classes, then productivity classes by rank inside each age class, a model fitted
   to the class means, with the productivity class as covariate.
5. **Parameter prediction** — fit each plot separately, regress the parameters on S, substitute back
   to get one equation Hd = f(t, S).

## Decisions you should know about

- **Plot IDs.** Plot numbers repeat across stands, so a plot is `stand-plot` (e.g. `2-3`).
- **Classifying a plot (methods 1 and 2).** Every measurement is projected to the index age and gives
  its own S. Both the value from the measurement closest to Ii and the mean of all the plot's
  measurements are saved; the class uses the mean.
- **Algebraic difference pairs.** All pairs inside each plot are used (82 pairs). Pairs from one plot
  share measurements, so the standard errors of that fit are optimistic; the estimates themselves are
  fine.
- **Hammer.** Five age classes (one per measurement round: 2.5–3.0, 3.5–4.0, 4.5–5.0, 5.5–6.0,
  6.0–6.5 yr) and three productivity classes (rank thirds) in each. The PDF's last example class
  (5.5–6.5) was split in two because it would have mixed two different measurement rounds.
- **Parameter prediction.** Each plot has 4–5 points and the model has 3 parameters, so a free fit per
  plot is unstable (e.g. b0 from 20.7 to 41.0, b1 from 0.11 to 0.50). b1 is fixed at the pooled value
  and b0, b2 are fitted per plot. The free fits are kept in `parametros_por_parcela.csv`
  (columns `*_livre`) to show the instability.
- **Fit statistics are not comparable across methods.** They are computed on different responses
  (all 43 points; 82 pairs; 15 class means; 43 points with S known). Compare the curves and the
  classifications, not the R² values.

## How to run

For each script edit the `CONFIG` block (3–4 file paths), then run, in this order:
`01` → `02` → `03` (needs 01's output) → `04` → `05` → `06`. Packages: `readxl` only.

## Important notes

- The R scripts were **not executed** while building this (no R available). Every number was computed
  with an independent Python implementation of the same algorithms; the files named
  `*_PYTHON_REFERENCE.*` in each `outputs/` folder are those numbers, to compare with your R run. Small
  differences in the 4th–5th decimal of nonlinear fits are normal.
- The `feixe_*.png` and `comparacao_metodos.png` files in the zip were drawn with Python. Running the R
  scripts overwrites them with the R versions (same file names).
- Chapter 10 of the textbook was not available. The formulations follow the logic in the challenge PDF
  and the standard Chapman-Richards difference equations; check notation against the book.
