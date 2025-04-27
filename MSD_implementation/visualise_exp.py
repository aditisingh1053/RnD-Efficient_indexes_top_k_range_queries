import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

# Function to generate random points and weights
def generate_poisson_points(n, d, lambda_poisson=10):
    points = np.random.poisson(lambda_poisson, (n, d))
    weights = np.random.rand(n)
    return points, weights

def generate_gaussian_points(n, d, mean=0, std=1):
    points = np.random.normal(mean, std, (n, d))
    weights = np.random.rand(n)
    return points, weights

def generate_uniform_points(n, d, low=0.0, high=1.0):
    points = np.random.uniform(low, high, (n, d))
    weights = np.random.rand(n)
    return points, weights

def generate_exponential_points(n, d, scale=1.0):
    points = np.random.exponential(scale, (n, d))
    weights = np.random.rand(n)
    return points, weights

def generate_beta_points(n, d, a=2.0, b=5.0):
    points = np.random.beta(a, b, (n, d))
    weights = np.random.rand(n)
    return points, weights

# Function to visualize points and highlight selected ones
def visualize_points(points, weights, selected_indices, distribution_name,force):
    plt.figure(figsize=(8, 6))

    scatter = plt.scatter(points[:, 0], points[:, 1], c=weights, cmap='viridis', s=50, edgecolor='none')
    plt.colorbar(scatter, label='Weight')

    selected_points = points[selected_indices]
    plt.scatter(
        selected_points[:, 0], selected_points[:, 1],
        facecolors='none', edgecolors='red', s=200, linewidths=2, label='Selected Points'
    )

    plt.xlabel('Dimension 1')
    plt.ylabel('Dimension 2')
    plt.title(f'{distribution_name} Distribution, λ = {force}')
    plt.legend()
    os.makedirs("plots_exp", exist_ok=True)
    plt.savefig(f"plots_exp/{distribution_name}_{force}.png", dpi=300)
    plt.close()

# Function to call the C++ script and get the selected indices
def get_selected_points_from_cpp(n, d, k, lambda_param, points, weights):
    input_data = f"{n} {d} {k} {lambda_param}\n"
    for i in range(n):
        point_str = ' '.join(map(str, points[i])) + f" {weights[i]}\n"
        input_data += point_str

    cpp_process = subprocess.Popen(
        ["./msd_exp"], 
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    stdout, stderr = cpp_process.communicate(input=input_data)
    
    if stderr:
        print("Error:", stderr)
        return []

    selected_points = []
    for line in stdout.strip().split("\n"):
        numbers = list(map(float, line.strip().split()))
        coords = numbers[:-1]  # Drop the weight
        selected_points.append(coords)

    selected_points = np.array(selected_points)

    selected_indices = []
    for sel in selected_points:
        for idx, point in enumerate(points):
            if np.allclose(point, sel, atol=1e-6):
                selected_indices.append(idx)
                break

    print(f"----- {len(selected_indices)} points selected -----")
    return selected_indices


n = 200
d = 2
k = 25

# List of generators and their names
generators = [
    ("Poisson", generate_poisson_points),
    ("Gaussian", generate_gaussian_points),
    ("Uniform", generate_uniform_points),
    ("Exponential", generate_exponential_points),
    ("Beta", generate_beta_points),
]

l=[0,50]
for force in l:
    for dist_name, generator in generators:
        print(f"\nGenerating and visualizing {dist_name} distribution...")

        # Generate points
        if dist_name == "Poisson":
            points, weights = generator(n, d, lambda_poisson=10)
        elif dist_name == "Gaussian":
            points, weights = generator(n, d, mean=0, std=3)
        elif dist_name == "Uniform":
            points, weights = generator(n, d, low=-10, high=10)
        elif dist_name == "Exponential":
            points, weights = generator(n, d, scale=2.0)
        elif dist_name == "Beta":
            points, weights = generator(n, d, a=2.0, b=5.0)

        # Get selected points from C++ program
        selected_indices = get_selected_points_from_cpp(n, d, k, force, points, weights)

        # Visualize
        if selected_indices:
            visualize_points(points, weights, selected_indices, dist_name, force)
        else:
            print(f"No selected points for {dist_name}.")