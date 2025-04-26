#!/bin/bash

INPUT_DIR="Testcases_exp/input"
EXPECTED_DIR="Testcases_exp/output"

# Check if the C++ filename argument is provided
if [ -z "$1" ]; then
    echo "Error: Please provide the C++ file as an argument (e.g., MMD_explicit.cpp)"
    exit 1
fi

CPP_FILE="$1"
EXECUTABLE="${CPP_FILE%.cpp}"

g++ -o "$EXECUTABLE" "$CPP_FILE"
EXECUTABLE="./$EXECUTABLE"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m'

total_score=0
max_score=0
start_time=$(date +%s.%N)

for input_file in "$INPUT_DIR"/input*; do
    filename=$(basename -- "$input_file")
    testcase_num="${filename##input}"
    testcase_num="${testcase_num%%.*}"
    expected_file="$EXPECTED_DIR/output$testcase_num"

    if [[ ! -f "$expected_file" ]]; then
        echo -e "${YELLOW}Expected output file $expected_file not found. Skipping testcase $testcase_num${NC}"
        continue
    fi

    lambda=$(head -n 1 "$input_file" | awk '{print $4}')
    if [[ -z "$lambda" ]]; then
        echo -e "${YELLOW}Could not extract lambda from $input_file. Skipping testcase $testcase_num${NC}"
        continue
    fi

    temp_output="temp_actual.txt"
    timeout 5s $EXECUTABLE < "$input_file" > "$temp_output"
    exit_status=$?
    if [ $exit_status -eq 124 ]; then
        echo -e "${YELLOW}Testcase $testcase_num timed out.${NC}"
        continue
    fi

    result=$(python3 - <<EOF "$input_file" "$temp_output" "$lambda" "$expected_file"
import sys
import numpy as np

def read_points(filename, skip_first=False):
    points = []
    with open(filename) as f:
        if skip_first:
            f.readline()
        for lineno, line in enumerate(f, 1):
            if not line.strip():
                continue
            tokens = list(map(float, line.strip().split()))
            if len(tokens) < 2:
                raise ValueError(f"Line {lineno} in {filename} is malformed: {line.strip()}")
            coords = tokens[:-1]
            weight = tokens[-1]
            points.append((np.array(coords), weight))
    return points

def compute_closest_pair_distance(points):
    """
    Computes the closest pair distance in a set of points.
    """
    min_distance = float('inf')
    n = len(points)
    for i in range(n):
        for j in range(i + 1, n):
            dist = np.linalg.norm(points[i][0] - points[j][0])
            min_distance = min(min_distance, dist)
    return min_distance

def compute_f(points, lam):
    """
    Computes f(T) = mu(T) + lambda * min(p in T) w(p)
    where mu(T) is the closest pair distance of points in T.
    """
    if not points:
        return 0.0
    
    # Closest pair distance
    closest_pair_distance = compute_closest_pair_distance(points)
    
    # Minimum weight
    min_weight = min(p[1] for p in points)
    
    return closest_pair_distance + lam * min_weight

try:
    input_points = read_points(sys.argv[1], skip_first=True)
    actual_points = read_points(sys.argv[2], skip_first=False)
    expected_points = read_points(sys.argv[4], skip_first=False)
    
    lam = float(sys.argv[3])
    expected_f_val = compute_f(expected_points, lam)
    actual_f_val = compute_f(actual_points, lam)
    
    if max(expected_f_val, actual_f_val) == 0:
        ratio = 1.0
    else:
        ratio = min(expected_f_val, actual_f_val) / max(expected_f_val, actual_f_val)
    
    status = "PASS" if ratio >= 0.5 else "FAIL"
    print(f"{status} {expected_f_val:.12f} {actual_f_val:.12f} {ratio:.6f}")
except Exception as e:
    print(f"ERROR: {e}")
    sys.exit(1)
EOF
)

    status=$(echo "$result" | awk '{print $1}')
    expected_val=$(echo "$result" | awk '{print $2}')
    actual_val=$(echo "$result" | awk '{print $3}')
    ratio=$(echo "$result" | awk '{print $4}')
    max_score=$((max_score + 1))

    # Add ratio to total_score safely
    total_score=$(echo "scale=6; $total_score + $ratio" | bc -l)

    if [[ "$status" == "PASS" ]]; then
        echo -e "${GREEN}Testcase $testcase_num passed${NC}"
        echo -e "Expected f(T): $expected_val"
        echo -e "Actual f(T):   $actual_val"
    elif [[ "$status" == "FAIL" ]]; then
        echo -e "${RED}Testcase $testcase_num failed${NC}"
        echo -e "Expected f(T): $expected_val"
        echo -e "Actual f(T):   $actual_val"
        echo "=================================="
    else
        echo -e "${YELLOW}Testcase $testcase_num encountered an error${NC}"
        echo "$result"
        echo "=================================="
    fi
done

end_time=$(date +%s.%N)
runtime=$(echo "scale=6; $end_time - $start_time" | bc -l)

printf "Final Score: ${GREEN}%.4f${NC}/$max_score\n" "$total_score"
printf "Total Runtime: ${GREEN}%.6f${NC} seconds\n" "$runtime"

rm -f temp_actual.txt "$EXECUTABLE"
