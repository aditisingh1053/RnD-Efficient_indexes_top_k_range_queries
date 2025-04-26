import numpy as np
import matplotlib.pyplot as plt
import subprocess

# --------- Generate Input ----------
N = 100  # Number of points
k = 10
delta = 2.0
d = 2

mean = [0, 0]
cov = [[5, 1], [1, 5]]
points = np.random.multivariate_normal(mean, cov, N)
weights = np.random.uniform(1, 10, N)

with open('input.txt', 'w') as f:
    f.write(f"{N} {k} {delta} {d}\n")
    for i in range(N):
        f.write(f"{weights[i]} {points[i,0]} {points[i,1]}\n")

print("✅ Input data written to input.txt")

# --------- Compile and Run C++ Code ----------
# Compile the C++ code
print("🛠️ Compiling C++ code...")
compile_process = subprocess.run(["g++", "-o", "greedy_cdt", "CDT_basic.cpp"], capture_output=True, text=True)

if compile_process.returncode != 0:
    print("❌ Compilation failed:")
    print(compile_process.stderr)
    exit(1)
else:
    print("✅ Compilation successful!")

# Now Run with input redirection
print("🚀 Running the C++ executable with input.txt...")
with open('input.txt', 'r') as infile:
    run_process = subprocess.run(["./greedy_cdt"], stdin=infile, capture_output=True, text=True)

if run_process.returncode != 0:
    print("❌ Execution failed:")
    print(run_process.stderr)
    exit(1)

# Save output to output.txt
with open('output.txt', 'w') as outfile:
    outfile.write(run_process.stdout)

print("✅ Output saved to output.txt")

# --------- Plotting ---------
selected_points = []
with open('output.txt', 'r') as f:
    lines = f.readlines()
    for line in lines:
        if line.strip() == "Selected points:":
            continue
        parts = line.strip().split("->")[0]
        x, y = parts.strip("() ").split(",")
        selected_points.append((float(x), float(y)))

# Plot
plt.figure(figsize=(10,10))

# Scatter all points, color by weight
sc = plt.scatter(points[:,0], points[:,1], c=weights, cmap='hot', s=50, alpha=0.7)
plt.colorbar(sc, label='Weight')  # Show colorbar

# Circle selected points
plt.scatter(selected_points[:,0], selected_points[:,1], facecolors='none', edgecolors='blue', s=300, linewidths=2, label='Selected points')

plt.legend()
plt.title("2D Gaussian Data with Weighted Heatmap and Selected Points Circled")
plt.grid(True)
plt.show()
