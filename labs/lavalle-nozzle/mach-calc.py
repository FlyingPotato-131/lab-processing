import math

k = 7/5
P0 = 1e5

Arel = (1700 / 1560)**2

M = 1
for i in range(10000):
	M = math.sqrt((k+1)/(k-1) * (Arel * M)**(2*(k-1)/(k+1)) - 2/(k-1))

print(f"M1 = {M}") #glass

Prel1 = 0.918

Mexp = M
for i in range(10000):
	Mexp = math.sqrt((k-1)/2/k + (k+1)/2/k * Prel1**(-k+1) * (2/(k+1) * (1/Mexp**2 + (k-1)/2))**(-k))
print(f"M1_exp = {Mexp}")

# Prel = (2/(k+1) / M**2 + (k-1)/(k+1))**(-k/(k-1)) * (2*k/(k+1) * M**2 - (k-1) / (k+1))**(-1 / (k-1))
Prel = (1 + (k-1)/2 * M**2)**(-k/(k-1))
print(f"P_1 = {1 / Prel}")

Arel = (2270 / 1575)**2

for i in range(10000):
	M = math.sqrt((k+1)/(k-1) * (Arel * M)**(2*(k-1)/(k+1)) - 2/(k-1))

print(f"M2 = {M}") #aluminium

Prel1 = 0.624

Mexp = M
for i in range(10000):
	Mexp = math.sqrt((k-1)/2/k + (k+1)/2/k * Prel1**(-k+1) * (2/(k+1) * (1/Mexp**2 + (k-1)/2))**(-k))
print(f"M1_exp = {Mexp}")

# Prel = (2/(k+1) / M**2 + (k-1)/(k+1))**(-k/(k-1)) * (2*k/(k+1) * M**2 - (k-1) / (k+1))**(-1 / (k-1))
Prel = (1 + (k-1)/2 * M**2)**(-k/(k-1))
print(f"P_2 = {1 / Prel}")
