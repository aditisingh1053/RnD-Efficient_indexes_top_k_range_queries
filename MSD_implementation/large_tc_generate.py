import random

n = 50000         
d = 20    
k = 500            
lambda_val = 15.0

input_filename = "./Testcases/input/input2"
output_filename = "./Testcases/output/output2"

input_lines = []
input_lines.append(f"{n} {d} {k} {lambda_val}")

best_points = []
inferior_points = []


for i in range(k):
    coords = [f"{random.uniform(-1000, 1000):.6f}" for _ in range(d-1)]
    weight = f"{random.uniform(40, 50):.6f}"
    point_str = " ".join(coords + [weight])
    best_points.append(point_str)


for i in range(n - k):
    coords = [f"{random.uniform(-1, 1):.6f}" for _ in range(d-1)]
    weight = f"{random.uniform(0, 1):.6f}"
    point_str = " ".join(coords + [weight])
    inferior_points.append(point_str)


all_points = best_points + inferior_points
random.shuffle(all_points)

input_lines.extend(all_points)

with open(input_filename, "w") as f:
    f.write("\n".join(input_lines))

with open(output_filename, "w") as f:
    f.write("\n".join(best_points))

print(f"Generated input file: {input_filename}")
print(f"Generated output file: {output_filename}")
