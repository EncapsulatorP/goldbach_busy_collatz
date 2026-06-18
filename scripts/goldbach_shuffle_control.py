import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Nmax = 500000

# --- sieve of Eratosthenes ---
sieve = np.ones(Nmax + 1, dtype=bool)
sieve[:2] = False
for i in range(2, int(Nmax**0.5) + 1):
    if sieve[i]:
        sieve[i * i :: i] = False
primes = np.nonzero(sieve)[0]
small_primes = primes[primes <= int(Nmax**0.5) + 2]

C2 = 0.6601618158468695  # twin-prime constant

def odd_prime_factors(n):
    facs = []
    m = n
    while m % 2 == 0:
        m //= 2
    for p in small_primes:
        if p == 2:
            continue
        if p * p > m:
            break
        if m % p == 0:
            facs.append(int(p))
            while m % p == 0:
                m //= p
    if m > 1 and m % 2 == 1:
        facs.append(int(m))
    return facs

def singular_series(n):
    S = 2.0 * C2
    for p in odd_prime_factors(n):
        S *= (p - 1) / (p - 2)
    return S

# --- sample even N, compute ORDERED Goldbach representation count r(N) ---
rng = np.random.default_rng(42)
candidates = np.arange(2000, Nmax, 2)        # even N
sample = np.sort(rng.choice(candidates, size=4000, replace=False))

reps = np.empty(len(sample))
S    = np.empty(len(sample))
for k, N in enumerate(sample):
    jmax = np.searchsorted(primes, N - 2, side="right")
    pl = primes[:jmax]                       # primes p with 2 <= p <= N-2
    reps[k] = np.count_nonzero(sieve[N - pl])  # ordered count: r(N) ~ S(N)*N/ln(N)^2
    S[k] = singular_series(int(N))

lnN = np.log(sample.astype(float))
y = reps * lnN**2 / sample                   # detrended residual ~ S(N)
y_shuf = rng.permutation(y)                  # break the N <-> count pairing

r_real = np.corrcoef(S, y)[0, 1]
r_shuf = np.corrcoef(S, y_shuf)[0, 1]
print(f"corr(S, detrended r(N))      = {r_real:.4f}")
print(f"corr(S, shuffled detrended)  = {r_shuf:.4f}")

div3 = (sample % 3 == 0)

# --- figure ---
fig, ax = plt.subplots(1, 2, figsize=(13.5, 6.0))
lo = min(S.min(), y.min()) - 0.1
hi = max(S.max(), y.max()) + 0.1

for a in ax:
    a.plot([lo, hi], [lo, hi], "k--", lw=1.3, zorder=5, label="y = x  (= S(N))")
    a.set_xlim(lo, hi); a.set_ylim(lo, hi)
    a.set_xlabel("singular series  S(N)   (Möbius / inclusion–exclusion built)")

c_no3 = "#c0392b"; c_y3 = "#2471a3"

ax[0].scatter(S[~div3], y[~div3], s=6, alpha=0.45, c=c_no3, label="3 ∤ N")
ax[0].scatter(S[div3],  y[div3],  s=6, alpha=0.45, c=c_y3,  label="3 | N")
ax[0].set_ylabel("detrended r(N) = r(N)·(ln N)² / N")
ax[0].set_title(f"REAL Goldbach counts — residual tracks S(N)\nr = {r_real:.4f}", fontsize=12)
ax[0].legend(loc="upper left", fontsize=9, framealpha=0.9)

ax[1].scatter(S[~div3], y_shuf[~div3], s=6, alpha=0.45, c=c_no3)
ax[1].scatter(S[div3],  y_shuf[div3],  s=6, alpha=0.45, c=c_y3)
ax[1].set_ylabel("shuffled detrended r(N)  (same marginal, pairing destroyed)")
ax[1].set_title(f"SHUFFLED control — diagonal evaporates\nr = {r_shuf:.4f}", fontsize=12)

fig.suptitle("Shuffle control: is the S(N) alignment data, or an artifact of the transform?",
             fontsize=13, y=0.99)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = "/mnt/user-data/outputs/goldbach_shuffle_control.png"
fig.savefig(out, dpi=130)
print("saved", out)
