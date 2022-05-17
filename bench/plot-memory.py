# Program to plot theoretical memory requirements.


# Import modules
# -----------------------------------------------------------------------------
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Path to the data folder, in which the benchmark outputs are stored.  The
# output plots will be saved to this folder.
path = "../data/output/"


# File name of the output plot
filename = "bench-memory.png"


DPI = 300
# -----------------------------------------------------------------------------


# Memory
# -----------------------------------------------------------------------------
N = np.array([10, 25, 50, 75, 100, 250, 500, 750, 1000, 2500, 5000, 7500, \
              10000, 25000, 50000])
s = ((N + 1) * (2 * N + 2) + (N + 1)**2) * 4 / 1024**3
d = ((N + 1) * (2 * N + 2) + (N + 1)**2) * 8 / 1024**3
q = ((N + 1) * (2 * N + 2) + (N + 1)**2) * 16 / 1024**3


# Plot the memory requirements
print("Ploting the memory requirements...")
fig, ax = plt.subplots(figsize=(14.0 / 2.54, 10.0 / 2.54))
ax.loglog(N, s, marker='o', label="single")
ax.loglog(N, d, marker='v', label="double")
ax.loglog(N, q, marker='+', label="quadruple")
ax.set_xlabel("Maximum harmonic degree")
ax.set_ylabel("Peak memory (GBs)")


y_major = matplotlib.ticker.LogLocator(base=10.0, numticks=10)
ax.yaxis.set_major_locator(y_major)
y_minor = matplotlib.ticker.LogLocator(base=10.0,
                                       subs=np.arange(1.0, 10.0)* 0.1,
                                       numticks=10)
ax.yaxis.set_minor_locator(y_minor)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax.grid(visible=True)
ax.legend()

outfile = path + filename
print("Saving the plot to %s..." % outfile)
fig.savefig(outfile, dpi=DPI)
# -----------------------------------------------------------------------------


print("Done.")

