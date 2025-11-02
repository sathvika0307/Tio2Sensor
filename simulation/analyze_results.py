import numpy as np
import matplotlib.pyplot as plt

# Read log file to extract energies
energies = []
temps = []
with open('log.lammps', 'r') as f:
    for line in f:
        if line.startswith('      ') and len(line.split()) == 6:
            parts = line.split()
            try:
                step = int(parts[0])
                temp = float(parts[1])
                pe = float(parts[2])
                energies.append(pe)
                temps.append(temp)
            except:
                pass

# Plot energy vs step
plt.figure(figsize=(10,5))
plt.subplot(1,2,1)
plt.plot(energies)
plt.xlabel('Step')
plt.ylabel('Potential Energy (eV)')
plt.title('Energy Convergence')

plt.subplot(1,2,2)
plt.plot(temps)
plt.xlabel('Step')
plt.ylabel('Temperature (K)')
plt.title('Temperature Evolution')

plt.tight_layout()
plt.savefig('energy_plot.png')
# plt.show()  # Remove for non-interactive

print(f"Final potential energy: {energies[-1] if energies else 'N/A'} eV")
print(f"Average temperature: {np.mean(temps) if temps else 'N/A'} K")
