
from typing import List


def ascii_histogram(data: List[float ], bins : int = 20, width: int = 50, title : str = ""):
    if not data:
        return ("nema podataka")

    lo, hi = min(data), max(data)
    if lo == hi:
        hi = lo + 1e-9
    bin_width = (hi-lo) / bins
    counts = [0] *bins

    for x in data:
        idx = int((x-lo)/bin_width)
        if idx == bins:
            idx -= 1
        counts[idx] += 1

    max_count = max(counts) if counts else 1
    lines = []
    if title:
        lines.append(title)
        lines.append("-" * max(len(title), width))
    for i, c in enumerate(counts):
        bar_len = int((c/ max_count)*width) if max_count else 0
        bin_start = lo+i*bin_width
        label = f"{bin_start : 8.3f}"
        lines.append(f"{label} | {'#' * bar_len} {c}")
    return "\n".join(lines)

