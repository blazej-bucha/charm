# Program to plot outputs from the benchmark program.


# Import modules
# -----------------------------------------------------------------------------
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# -----------------------------------------------------------------------------


# Parse input arguments
# -----------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Plots outputs from '
                                 'the benchmark program.')


group = parser.add_mutually_exclusive_group(required=True)


group.add_argument('--single-precision',
                   '-s',
                   action='store_true',
                   help='Plots single-precision outputs.')


group.add_argument('--double-precision',
                   '-d',
                   action='store_true',
                   help='Plots double-precision outputs.')


group.add_argument('--quad-precision',
                   '-q',
                   action='store_true',
                   help='Plots quadruple-precision outputs.')


# Parse the arguments
args = parser.parse_args()
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
if args.single_precision:
    prec = "f"
elif args.double_precision:
    prec = ""
elif args.quad_precision:
    prec = "q"


# Path to the data folder, in which the benchmark outputs are stored.  The
# output plots will be saved to this folder.
path = "../data/output/"


# Suffix of the name of the file that stores the wall-clock times from SHA
file_sha_t = "sha-time.txt"


# Suffix of the name of the file that stores the wall-clock times from SHS
file_shs_t = "shs-time.txt"


# Suffix of the name of the file that stores the accuracy.
file_acc = "acc.txt"


DPI = 300
# -----------------------------------------------------------------------------


# Wall-clock times
# -----------------------------------------------------------------------------
# Load the SHA times and convert them to minutes
infile = path + "bench" + prec + "-" + file_sha_t
print("Reading input file %s..." % infile)
sha_t = np.loadtxt(infile)
sha_t[:, 1] /= 60.0


# Load the SHS times and convert them to minutes
infile = path + "bench" + prec + "-" + file_shs_t
print("Reading input file %s..." % infile)
shs_t = np.loadtxt(infile)
shs_t[:, 1] /= 60.0


# Plot the wall-clock times
print("Ploting the wall-clock times...")
fig, ax = plt.subplots(figsize=(14.0 / 2.54, 10.0 / 2.54))
ax.loglog(sha_t[:, 0], sha_t[:, 1], marker='o', label="SHA")
ax.loglog(shs_t[:, 0], shs_t[:, 1], marker='v', label="SHS")

ymin_fig, ymax_fig = ax.get_ylim()
xmin_fig, xmax_fig = ax.get_xlim()

xmin = sha_t[:, 0].min()
xmax = sha_t[:, 0].max()

y = 1.0 / 60.0
if y <= ymax_fig:
    ax.text(xmin, y + 0.005, '1 sec')
    ax.hlines(y, xmin, xmax, colors='black', linestyles='dashed',
              linewidth=0.75)
y = 1.0
if y <= ymax_fig:
    ax.text(xmin, y + 0.3, '1 min')
    ax.hlines(y, xmin, xmax, colors='black', linestyles='dashed',
              linewidth=0.75)
y = 60.0
if y <= ymax_fig:
    ax.text(xmin, y + 10, '1 hour')
    ax.hlines(y, xmin, xmax, colors='black', linestyles='dashed',
              linewidth=0.75)


ax.set_ylim(ymin_fig, ymax_fig)
ax.set_xlim(xmin_fig, xmax_fig)


ax.set_xlabel("Maximum harmonic degree")
ax.set_ylabel("Wall-clock time (minutes)")


y_major = matplotlib.ticker.LogLocator(base=10.0, numticks=10)
ax.yaxis.set_major_locator(y_major)
y_minor = matplotlib.ticker.LogLocator(base=10.0,
                                       subs=np.arange(1.0, 10.0)* 0.1,
                                       numticks=10)
ax.yaxis.set_minor_locator(y_minor)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax.grid(visible=True)
ax.legend(loc='lower right')

outfile = path + "bench" + prec + "-time.png"
print("Saving the plot to %s..." % outfile)
fig.savefig(outfile, dpi=DPI)
# -----------------------------------------------------------------------------


# Accuracy
# -----------------------------------------------------------------------------
acc = np.loadtxt(path + "bench" + prec + "-" + file_acc)


# Plot the maximum errors of spherical harmonic analysis
print("Ploting the accuracy...")
fig, ax = plt.subplots(figsize=(14.0 / 2.54, 10.0 / 2.54))
ax.loglog(acc[:, 0], acc[:, 1], marker='o', label="$\epsilon_\max$")
ax.loglog(acc[:, 0], acc[:, 2], marker='v', label="$\epsilon_\mathrm{rms}$")
ax.set_xlabel("Maximum harmonic degree")
ax.legend()
ax.grid(visible=True)

outfile = path + "bench" + prec + "-accuracy.png"
print("Saving the plot to %s..." % outfile)
fig.savefig(outfile, dpi=DPI)
# -----------------------------------------------------------------------------


print("Done.")
