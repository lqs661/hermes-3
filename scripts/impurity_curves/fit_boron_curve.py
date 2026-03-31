from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


ROOT = Path(__file__).resolve().parents[2]
DATASET = ROOT / "scripts" / "impurity_curves" / "radas_output" / "boron" / "output" / "boron.nc"
OUTPUT_PLOT = (
    ROOT
    / "scripts"
    / "impurity_curves"
    / "radas_output"
    / "boron"
    / "output"
    / "boron_curve_comparison.png"
)

BORON_FIT_LOWER = 1.5
BORON_FIT_UPPER = 100.0
BORON_FIT_LOW_VALUE = 3.78400258e-32
BORON_FIT_HIGH_VALUE = 7.29127229e-33
BORON_FIT_COEFFS = np.array(
    [
        -6.72013939e01,
        -3.90954847e01,
        1.20065392e02,
        -1.94724211e02,
        1.95093156e02,
        -1.26633792e02,
        5.37183419e01,
        -1.47341369e01,
        2.51383555e00,
        -2.42260929e-01,
        1.00647912e-02,
    ]
)


def poly_eval(log_t: np.ndarray, coeffs: np.ndarray) -> np.ndarray:
    total = np.zeros_like(log_t)
    for i, coeff in enumerate(coeffs):
        total += coeff * np.power(log_t, i)
    return total


def fit_window(te: np.ndarray, lz: np.ndarray, lower: float, upper: float):
    mask = (te >= lower) & (te <= upper)
    te_fit = te[mask]
    lz_fit = lz[mask]
    coeffs = np.polyfit(np.log(te_fit), np.log(lz_fit), deg=10)[::-1]
    fit = np.exp(poly_eval(np.log(te_fit), coeffs))
    max_rel = float(np.max(np.abs(fit / lz_fit - 1.0)))
    rms_log = float(np.sqrt(np.mean((np.log(fit) - np.log(lz_fit)) ** 2)))
    return {
        "lower": lower,
        "upper": upper,
        "coeffs": coeffs,
        "low_value": float(lz[te <= lower][-1]),
        "high_value": float(lz[te >= upper][0]),
        "max_rel": max_rel,
        "rms_log": rms_log,
    }


def fixed_fraction_boron_curve(te: np.ndarray) -> np.ndarray:
    te = np.asarray(te, dtype=float)
    result = np.empty_like(te)

    lower_mask = te < BORON_FIT_LOWER
    middle_mask = (te >= BORON_FIT_LOWER) & (te <= BORON_FIT_UPPER)
    upper_mask = te > BORON_FIT_UPPER

    result[lower_mask] = BORON_FIT_LOW_VALUE
    result[upper_mask] = BORON_FIT_HIGH_VALUE
    result[middle_mask] = np.exp(poly_eval(np.log(te[middle_mask]), BORON_FIT_COEFFS))

    return result


def save_comparison_plot(te: np.ndarray, lz_radas: np.ndarray, lz_fit: np.ndarray) -> None:
    relative_error = np.abs(lz_fit / lz_radas - 1.0)

    fig, (ax_curve, ax_error) = plt.subplots(
        2,
        1,
        figsize=(7, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 2]},
    )

    ax_curve.loglog(te, lz_radas, label="RADAS equilibrium Lz @ ne~1e20", linewidth=2)
    ax_curve.loglog(te, lz_fit, "--", label="Hermes fixed_fraction_boron fit", linewidth=2)
    ax_curve.axvline(BORON_FIT_LOWER, color="grey", linestyle=":", linewidth=1)
    ax_curve.axvline(BORON_FIT_UPPER, color="grey", linestyle=":", linewidth=1)
    ax_curve.set_ylabel(r"$L_z\ [W\,m^3]$")
    ax_curve.set_title("Boron cooling curve comparison")
    ax_curve.grid(True, which="both", alpha=0.3)
    ax_curve.legend()

    ax_error.semilogx(te, relative_error, color="tab:red", linewidth=2)
    ax_error.axvline(BORON_FIT_LOWER, color="grey", linestyle=":", linewidth=1)
    ax_error.axvline(BORON_FIT_UPPER, color="grey", linestyle=":", linewidth=1)
    ax_error.set_xlabel(r"$T_e\ [eV]$")
    ax_error.set_ylabel("Relative error")
    ax_error.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUTPUT_PLOT, dpi=200)
    plt.close(fig)


def main():
    ds = xr.open_dataset(DATASET)
    te = ds["dim_electron_temp"].values
    lz = (
        ds["equilibrium_Lz"]
        .sel(dim_electron_density=1e20, method="nearest")
        .isel(dim_ne_tau=0)
        .values
    )

    windows = []
    for lower in [1.0, 1.25, 1.5, 2.0]:
        for upper in [100.0, 200.0, 500.0, 1000.0, 1500.0, 3000.0, 10000.0, 50000.0]:
            mask = (te >= lower) & (te <= upper)
            if np.count_nonzero(mask) < 15:
                continue
            windows.append(fit_window(te, lz, lower, upper))

    windows.sort(key=lambda item: (item["rms_log"], item["max_rel"]))

    print("Top candidate windows:")
    for item in windows[:10]:
        print(
            f"lower={item['lower']:.2f} upper={item['upper']:.2f} "
            f"max_rel={item['max_rel']:.4e} rms_log={item['rms_log']:.4e}"
        )

    best = windows[0]
    print("\nSelected window:")
    print(best["lower"], best["upper"], best["low_value"], best["high_value"])
    print("\nCoefficients:")
    for i, coeff in enumerate(best["coeffs"]):
        suffix = "," if i < len(best["coeffs"]) - 1 else ""
        print(f"{coeff:+.8e}{suffix}")

    peak_index = int(np.nanargmax(lz))
    print("\nPeak:")
    print(f"Te={te[peak_index]:.8f} Lz={lz[peak_index]:.8e}")

    lz_fit = fixed_fraction_boron_curve(te)
    relative_error = np.abs(lz_fit / lz - 1.0)
    fit_mask = (te >= BORON_FIT_LOWER) & (te <= BORON_FIT_UPPER)

    print("\nImplemented fixed_fraction_boron error:")
    print(
        "fit window max_rel="
        f"{np.max(relative_error[fit_mask]):.4e} "
        f"full range max_rel={np.max(relative_error):.4e}"
    )

    save_comparison_plot(te, lz, lz_fit)
    print(f"\nSaved plot: {OUTPUT_PLOT}")


if __name__ == "__main__":
    main()
