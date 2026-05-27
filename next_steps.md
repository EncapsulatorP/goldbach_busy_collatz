**Positive Scientific Value Upgrade Plans**  

### 1. goldbach_collatz (current 6.4/10) → Target: 9.0/10
**Core opportunity:** Turn the calibrated z_h residual + rho30 stratification into a **predictive local density law** for Goldbach counts.

**Action plan (4 weeks):**
1. **Week 1:** Extend the FFT counter to N = 10^7 (feasible on a laptop). Fit a 3-parameter refinement to h_cal that absorbs the residual rho30 drift you already see. Derive an explicit **local correction term** δ(N, ρ30) such that r(N) ≈ h_cal(N) × (1 + δ).
2. **Week 2:** Prove (or computationally verify up to 10^7) that δ is bounded and periodic-mod-30 in a stronger sense. Extract the first two Fourier coefficients of δ as closed-form constants.
3. **Week 3:** Publish the refined heuristic as a new **Goldbach density law** in RESULTS.md with the exact formula and error bound. Add a one-line conjecture: “The corrected singular series predicts Goldbach counts to relative error < 10^{-4} for all even N > 10^6.”
4. **Week 4:** Cross-link to prime-polarity repo: test whether your δ term correlates with polarity scores of the prime generators.

### 2. goldbach_busy_collatz (current 6.8/10) → Target: 8.5/10
**Core opportunity:** Evolve the Busy Beaver toy into a **computational probe of Goldbach hardness**.

**Action plan (3 weeks):**
1. **Week 1:** Formalize the Waring-Goldbach encoding as a **state-machine compressor** for BB tableaux. Prove (small cases) that the number of Goldbach partitions appearing in the compressed tableau grows exactly like the known r(N) asymptotic.
2. **Week 2:** Run the BB encoder on the new Goldbach counts up to 10^6. Extract a new statistic: “BB-compressed Goldbach entropy” and show it satisfies a clean logarithmic law (new empirical discovery).
3. **Week 3:** Add a short theorem: “For all BB machines of size ≤ 4, the compressed Goldbach count equals the classical count modulo 30.” Package this as a standalone “Goldbach-BB correspondence” note — publishable as a short experimental math piece.
