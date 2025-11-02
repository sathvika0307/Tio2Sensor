from ase.io import read
from ase.visualize import view
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load TiO2 nanoparticle
tio2 = read('simulation/tio2_nanoparticle.xyz')

# Load ZnO nanoparticle
zno = read('zno_nanoparticle.xyz')

# Function to plot 3D structure
def plot_structure(atoms, title, filename):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    colors = {'Ti': 'blue', 'O': 'red', 'Zn': 'green'}
    for sym, pos in zip(symbols, positions):
        ax.scatter(pos[0], pos[1], pos[2], c=colors.get(sym, 'black'), s=50)
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    ax.set_title(title)
    # Manual legend
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=10) for c in colors.values()]
    ax.legend(handles, colors.keys())
    plt.savefig(filename)
    plt.close()

# Plot TiO2
plot_structure(tio2, 'TiO2 Nanoparticle (181 atoms)', 'tio2_structure.png')

# Plot ZnO
plot_structure(zno, 'ZnO Nanoparticle (334 atoms)', 'zno_structure.png')

print("Generated structure images: tio2_structure.png, zno_structure.png")
