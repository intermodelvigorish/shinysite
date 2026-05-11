# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

A Quarto website documenting the **InterModel Vigorish (IMV)** — a statistical metric for comparing binary prediction models (Domingue et al., 2025). The site hosts interactive R Shiny demonstrations and computational examples using the `imv` R package (available on CRAN).

## Build commands

```bash
quarto render          # Render all .qmd files → docs/
quarto render foo.qmd  # Render a single page
quarto preview         # Live preview with hot reload
```

Output goes to `docs/` (GitHub Pages source). The `.quarto/` cache directory and `*.quarto_ipynb` files are gitignored. CLAUDE.md is also gitignored (local-only).

## Architecture

- **Source**: `.qmd` files in the root — each is a standalone Quarto document compiled to HTML
- **Output**: `docs/` — committed rendered HTML (site is served from this directory)
- **Config**: `_quarto.yml` — site structure, navbar, theme (cosmo), and the `shinylive` filter
- **Interactivity**: Embedded Shiny apps use the `shinylive` Quarto extension (`{shinylive-r}` chunks with `standalone: true`). These run entirely client-side via WebAssembly (no server needed).

## Content structure

| Section | Files | Purpose |
|---|---|---|
| Home | `index.qmd` | IMV concept, weighted-coin analogy, interactive tutorial |
| Simulation Examples | `logreg.qmd`, `oracle.qmd`, `overfit.qmd`, `23pl.qmd`, `collapse.qmd`, `implied.qmd` | Shiny-powered simulations illustrating IMV behavior |
| Computational Examples | `example_glm.qmd`, `example_glmer.qmd`, `example_irt_uni.qmd`, `example_irt_multi.qmd` | Real-data walkthroughs using `glm`, `glmer`, and IRT models |

## Key conventions

- Each page is self-contained. Shiny app logic (including IMV helper functions) is inlined directly in `{shinylive-r}` chunks — there are no shared R source files.
- The three core IMV helper functions (`get_coin_weight`, `avg_ll`, `imv_coins`) are duplicated in every shinylive chunk that needs them. When updating the IMV math, update all pages.
- The `imv` package is on CRAN; computational example chunks use `install.packages("imv")` (commented out).
- Computational example pages (`example_*.qmd`) use plain `r` chunks (not shinylive) — they are static walkthroughs, not interactive apps.
- Chart.js is loaded from a CDN (`cdnjs.cloudflare.com`) inside shinylive chunks for canvas-based visualizations; it is not available as a package dependency.
