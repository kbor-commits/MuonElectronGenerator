#!/usr/bin/env python3
"""Plot energy spectra and angular distributions from cosmic_muons CSV output.

Usage: plot_muons.py <run_dir> <plot_dir> [primary_recipe]
"""
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

run_dir, plot_dir = sys.argv[1], sys.argv[2]
primary = sys.argv[3] if len(sys.argv) > 3 else "Reyna"
os.makedirs(plot_dir, exist_ok=True)

# Reference categorical palette (fixed order) + recessive ink
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({
    "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb",
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 10,
    "axes.titlesize": 11, "axes.titleweight": "bold", "legend.frameon": False,
    "lines.linewidth": 2, "figure.dpi": 150,
})

def load(path):
    return np.genfromtxt(path, delimiter=",", comments="#")

def recipe_data(r):
    g = load(os.path.join(run_dir, f"muons_{r}_nt_generated.csv"))
    pfile = os.path.join(run_dir, f"muons_{r}_nt_plane.csv")
    p = load(pfile) if os.path.exists(pfile) else None
    return g, p

PBINS = np.logspace(-2, np.log10(3000), 81)
def spectrum(p, bins=PBINS, norm=None):
    h, e = np.histogram(p, bins=bins)
    w = np.diff(e); c = np.sqrt(e[1:] * e[:-1])
    n = norm if norm else len(p)
    with np.errstate(divide="ignore", invalid="ignore"):
        y = h / w / n; err = np.sqrt(h) / w / n
    m = h > 0
    return c[m], y[m], err[m]

summary = []
gen, plane = recipe_data(primary)
p_g, th_g = gen[:, 1], gen[:, 3]
N = len(gen)

# ---------------------------------------------------------- 1. momentum spectrum
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
ax = axes[0]
x, y, e = spectrum(p_g)
ax.errorbar(x, y, yerr=e, fmt="o", ms=3, color=C[0], label="Generated (MuonElectronGenerator)")
if plane is not None:
    x2, y2, _ = spectrum(plane[:, 1], norm=N)
    ax.step(x2, y2, where="mid", color=C[1], lw=1.4, label="At scoring plane (after Geant4 transport)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Muon momentum p [GeV/c]"); ax.set_ylabel("dN/dp  [1/(GeV/c), per muon]")
ax.set_title(f"Momentum spectrum - {primary} recipe, all zenith angles")
ax.legend(loc="lower left", fontsize=8.5)

ax = axes[1]
ranges = [(0, 30), (30, 60), (60, 75), (75, 90)]
for i, (a, b) in enumerate(ranges):
    sel = (th_g >= a) & (th_g < b)
    if sel.sum() < 50: continue
    x, y, e = spectrum(p_g[sel])
    ax.plot(x, y * x**2.7, color=C[i], label=f"{a}°-{b}°  ({sel.sum()/N:.1%} of muons)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Muon momentum p [GeV/c]"); ax.set_ylabel("p$^{2.7}$ · dN/dp  (normalised per bin)")
ax.set_title("Spectrum vs zenith angle (p$^{2.7}$-weighted)")
ax.legend(title="Zenith angle", fontsize=8.5, title_fontsize=8.5, loc="lower center")
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, f"{primary}_momentum_spectrum.png")); plt.close(fig)

# ---------------------------------------------------------- 2. angular distribution
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
tb = np.linspace(0, 90, 46); tc = 0.5 * (tb[1:] + tb[:-1])
h, _ = np.histogram(th_g, bins=tb)
dth = np.radians(np.diff(tb))
ax = axes[0]
ax.errorbar(tc, h / dth / N, yerr=np.sqrt(h) / dth / N, fmt="o", ms=3.5, color=C[0], label="Generated")
if plane is not None:
    hp, _ = np.histogram(plane[:, 3], bins=tb)
    ax.step(tc, hp / dth / N, where="mid", color=C[1], lw=1.4, label="At scoring plane")
t = np.radians(np.linspace(0, 90, 400))
f = np.cos(t)**3 * np.sin(t); f /= np.trapz(f, t)
ax.plot(np.degrees(t), f, color=INK2, lw=1.2, ls="--", label=r"$\cos^3\theta\,\sin\theta$ (expected)")
ax.set_xlabel("Zenith angle θ [deg]"); ax.set_ylabel("dN/dθ  [1/rad, per muon]")
ax.set_title(f"Zenith-angle distribution - {primary} recipe")
ax.set_xlim(0, 90); ax.legend(fontsize=8.5)

ax = axes[1]
# Flux per unit solid angle through a horizontal surface, normalised to vertical
dOmega = 2 * np.pi * (np.cos(np.radians(tb[:-1])) - np.cos(np.radians(tb[1:])))
dn = h / dOmega
c3 = np.cos(np.radians(tc))**3
ref = np.sum(dn * c3) / np.sum(c3 * c3)  # least-squares scale to cos^3 over all bins
ax.errorbar(tc, dn / ref, yerr=np.sqrt(h) / dOmega / ref, fmt="o", ms=3.5, color=C[0],
            label="Generated, dN/dΩ (horizontal detector)")
ax.plot(np.degrees(t), np.cos(t)**3, color=INK2, lw=1.2, ls="--", label=r"$\cos^3\theta$")
ax.plot(np.degrees(t), np.cos(t)**2, color=C[2], lw=1.2, ls=":", label=r"$\cos^2\theta$ (intensity)")
ax.set_xlabel("Zenith angle θ [deg]"); ax.set_ylabel("dN/dΩ  (relative to vertical)")
ax.set_title("Per-solid-angle distribution")
ax.set_xlim(0, 90); ax.legend(fontsize=8.5)
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, f"{primary}_angular_distribution.png")); plt.close(fig)

# ---------------------------------------------------------- 3. 2D momentum vs angle
fig, ax = plt.subplots(figsize=(6.5, 4.8))
H, xe, ye = np.histogram2d(th_g, np.log10(p_g), bins=[np.linspace(0, 90, 61), np.linspace(-2, 3.5, 56)])
im = ax.pcolormesh(xe, 10**ye, H.T, cmap="Blues", norm=matplotlib.colors.LogNorm(vmin=1))
ax.set_yscale("log"); ax.grid(False)
ax.set_xlabel("Zenith angle θ [deg]"); ax.set_ylabel("Muon momentum p [GeV/c]")
ax.set_title(f"Momentum vs zenith angle - {primary} recipe")
cb = fig.colorbar(im, ax=ax); cb.set_label("Muons per bin")
fig.tight_layout()
fig.savefig(os.path.join(plot_dir, f"{primary}_momentum_vs_angle.png")); plt.close(fig)

# ---------------------------------------------------------- 4. recipe comparison
recipes = [r for r in ["Reyna", "Blanpied", "Chris", "Jokisch", "Haruo"]
           if os.path.exists(os.path.join(run_dir, f"muons_{r}_nt_generated.csv"))]
if len(recipes) > 1:
    labels = {"Reyna": "Reyna (cos³θ sinθ, p∝1/cosθ)", "Blanpied": "Blanpied (mmuon.f fit)",
              "Chris": "Chris (Los Alamos, 2200 m)", "Jokisch": "Jokisch (75° spectrum)",
              "Haruo": "Haruo (Jokisch → 2200 m)"}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for i, r in enumerate(recipes):
        g, _ = recipe_data(r)
        x, y, _ = spectrum(g[:, 1])
        axes[0].plot(x, y, color=C[i], label=labels[r])
        hh, _ = np.histogram(g[:, 3], bins=tb)
        axes[1].plot(tc, hh / dth / len(g), color=C[i], label=labels[r])
    axes[0].set_xscale("log"); axes[0].set_yscale("log")
    axes[0].set_xlabel("Muon momentum p [GeV/c]"); axes[0].set_ylabel("dN/dp  [1/(GeV/c), per muon]")
    axes[0].set_title("Momentum spectrum by recipe")
    axes[1].set_xlabel("Zenith angle θ [deg]"); axes[1].set_ylabel("dN/dθ  [1/rad, per muon]")
    axes[1].set_title("Zenith-angle distribution by recipe"); axes[1].set_xlim(0, 90)
    axes[0].legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(plot_dir, "recipe_comparison.png")); plt.close(fig)

for r in (recipes if len(recipes) > 1 else [primary]):
    g, _ = recipe_data(r)
    summary.append((r, len(g), np.mean(g[:, 1]), np.median(g[:, 1]), np.mean(g[:, 3]),
                    np.mean(g[:, 0] == -13) / max(np.mean(g[:, 0] == 13), 1e-12)))
with open(os.path.join(plot_dir, "summary.txt"), "w") as fo:
    fo.write(f"{'recipe':<10}{'N':>10}{'<p> GeV/c':>12}{'median p':>11}{'<theta> deg':>13}{'mu+/mu-':>9}\n")
    for s in summary:
        fo.write(f"{s[0]:<10}{s[1]:>10d}{s[2]:>12.3f}{s[3]:>11.3f}{s[4]:>13.2f}{s[5]:>9.3f}\n")
    if plane is not None:
        fo.write(f"\n{primary}: {len(plane)}/{N} primaries reached the scoring plane\n")
print(open(os.path.join(plot_dir, "summary.txt")).read())
