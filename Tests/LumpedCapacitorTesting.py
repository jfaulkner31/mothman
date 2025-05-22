from Kernels.LumpedCapacitor import *

mat = LumpedMaterialProperty(name='test prop', function_coeffs=[1983, 0.01, 0.001])
ax = mat.plot_prop(xMin=0, xMax=1200)
plt.title("y = 1983 + 0.01*X + 0.001*X^2", fontsize=15)

plt.savefig('Tests/LumpedCapacitorTestingResult.png')

print("Successfully ran lumped capacitor test - check output for validity!")
