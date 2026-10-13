# CHECK_LAMBDA.py
# NIST values
LAB_NIST = {
    'Hβ_air': 4861.33,
    'Hβ_vacuum': 4862.68,
    'Hα_air': 6562.80,
    'Hα_vacuum': 6564.63,  # corrected
    '[OIII]_4959_air': 4958.91,
    '[OIII]_4959_vacuum': 4960.30,
    '[OIII]_5007_air': 5006.84,
    '[OIII]_5007_vacuum': 5008.24,
}

print("NIST values:")
for k, v in LAB_NIST.items():
    print(f"  {k:20s} = {v:.4f} Å")

print()
print("Air vs vacuum difference:")
for line in ['Hβ', 'Hα', '[OIII]_4959', '[OIII]_5007']:
    d = LAB_NIST[f'{line}_vacuum'] - LAB_NIST[f'{line}_air']
    print(f"  {line:12s}  Δ = {d:+.4f} Å")

# What the pipeline used
print()
print("Values used in the pipeline:")
LAB_PIPE = {'Hb': 4862.68, 'OIII4959': 4960.30, 'OIII5007': 5008.24, 'Ha': 6562.80}
for k, v in LAB_PIPE.items():
    print(f"  {k:12s} = {v:.4f} Å")

print()
print("Is Hα in air or vacuum?")
print(f"  Pipeline uses 6562.80 (air)")
print(f"  Vacuum would be 6564.63")
print(f"  Difference: {6564.63 - 6562.80:+.3f} Å")
print()
print("If the stack is in vacuum and the fit is in air, Hα would have offset +1.83 Å")
print("The measured Hα offset is +1.8858 Å")
print("It matches!")