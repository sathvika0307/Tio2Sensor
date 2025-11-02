import matplotlib.pyplot as plt

# Data from simulations
materials = ['TiO2 (181 atoms)', 'ZnO (334 atoms)']
binding_energies = [-103.7, -1.83]  # eV, note: TiO2 from report, ZnO simulated

plt.figure(figsize=(8, 6))
bars = plt.bar(materials, binding_energies, color=['blue', 'green'])
plt.ylabel('Binding Energy (eV)')
plt.title('Methane Binding Energy Comparison')
plt.axhline(0, color='black', linewidth=0.5)

# Add value labels
for bar, energy in zip(bars, binding_energies):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1, f'{energy} eV', ha='center', va='bottom')

plt.savefig('binding_energy_comparison.png')
plt.close()

print("Generated binding energy comparison plot: binding_energy_comparison.png")
