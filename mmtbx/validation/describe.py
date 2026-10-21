import statistics as st
import math
import sys


def describe(data: list[float]) -> dict:
    """Return meaningful descriptive stats for a list of floats."""
    if not data:
        raise ValueError("data must not be empty")

    xs = sorted(data)
    n = len(xs)

    def percentile(p: float) -> float:
        """Linear-interpolated percentile (p in 0..100)."""
        k = (n - 1) * p / 100
        lo, hi = math.floor(k), math.ceil(k)
        return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)

    mean = st.fmean(xs)
    stdev = st.stdev(xs) if n > 1 else 0.0
    q1, median, q3 = percentile(25), percentile(50), percentile(75)
    iqr = q3 - q1

    # Outliers via the 1.5 * IQR rule
    lower_fence, upper_fence = q1 - 1.5 * iqr, q3 + 1.5 * iqr
    outliers = [x for x in xs if x < lower_fence or x > upper_fence]

    # Skewness (sample, adjusted Fisher-Pearson) and excess kurtosis
    skew = kurt = None
    if n > 2 and stdev > 0:
        m3 = sum((x - mean) ** 3 for x in xs) / n
        m2 = sum((x - mean) ** 2 for x in xs) / n
        skew = (m3 / m2 ** 1.5) * math.sqrt(n * (n - 1)) / (n - 2)
    if n > 3 and stdev > 0:
        m4 = sum((x - mean) ** 4 for x in xs) / n
        m2 = sum((x - mean) ** 2 for x in xs) / n
        kurt = m4 / m2 ** 2 - 3

    return {
        "count": n,
        "sum": math.fsum(xs),
        "mean": mean,
        "median": median,
        "mode": st.multimode(xs) if len(set(xs)) < n else None,
        "min": xs[0],
        "max": xs[-1],
        "range": xs[-1] - xs[0],
        "variance": st.variance(xs) if n > 1 else 0.0,
        "std_dev": stdev,
        "coef_of_variation": stdev / mean if mean else None,
        "q1": q1,
        "q3": q3,
        "iqr": iqr,
        "p5": percentile(5),
        "p95": percentile(95),
        "skewness": skew,
        "excess_kurtosis": kurt,
        "std_error": stdev / math.sqrt(n),
        "ci95_mean": (mean - 1.96 * stdev / math.sqrt(n),
                      mean + 1.96 * stdev / math.sqrt(n)),
        "outliers": outliers,
    }


def print_stats(stats: dict, out=sys.stdout) -> None:
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"{key:>18}: {value:.4f}", file=out)
        else:
            print(f"{key:>18}: {value}", file=out)

import statistics as st


def plot_histogram(data: list[float], bins: int | str = "auto", title: str = "Histogram") -> None:
    import matplotlib.pyplot as plt
    if not data:
        raise ValueError("data must not be empty")

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.hist(data, bins=bins, edgecolor="black", alpha=0.75)

    # Reference lines for mean and median
    mean, median = st.fmean(data), st.median(data)
    ax.axvline(mean, color="red", linestyle="--", linewidth=1.5, label=f"Mean: {mean:.2f}")
    ax.axvline(median, color="green", linestyle="-.", linewidth=1.5, label=f"Median: {median:.2f}")

    ax.set_title(title)
    ax.set_xlabel("Value")
    ax.set_ylabel("Frequency")
    ax.set_yscale("log")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    sample = [2.3, 4.1, 3.8, 5.0, 4.4, 3.9, 4.2, 12.7, 4.0, 3.6]
    plot_histogram(sample, bins=8, title="Sample Distribution")


if __name__ == "__main__":
    sample = [2.3, 4.1, 3.8, 5.0, 4.4, 3.9, 4.2, 12.7, 4.0, 3.6]
    print_stats(describe(sample))

