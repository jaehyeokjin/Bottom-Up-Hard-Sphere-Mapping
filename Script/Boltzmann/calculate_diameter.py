import math
from scipy.integrate import simps
import sys

def find_diameter(x, v, threshold):
    min_diff = float('inf')
    min_r = None
    for i in range(len(v)):
        diff = abs(v[i] - threshold)**2
        if diff < min_diff:
            min_diff = diff
            min_r = x[i]
    return min_r


def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <temperature>")
        sys.exit(1)
    
    tempval = sys.argv[1]
    temp = float(tempval)
    filename = "../"+tempval + ".out"

    r_end = 20.00
    n_r = int((r_end) / 0.02) + 1
    
    f = open(filename, 'r')

    x = []
    y = []
    r = [0.0 + 0.02 * i for i in range(0, n_r)]
    v = []
    count = 0
    ref_val = 0.0

    for line in f:
        line_e = line.split()
        if count > -1:
            if count == 0:
                ref_val = float(line_e[2])
                ref_x = float(line_e[1])
                ref_ind = int(ref_x / 0.02)
            y.append(float(line_e[2]))
            x.append(float(line_e[1]))
        count += 1

    for i in range(0, len(r)):
        if i < ref_ind:
            v.append(ref_val)
        else:
            v.append(y[i - ref_ind])
    
    f.close()

    # Calculate Boltzmann hard sphere diameter (v(r) = k_B T)
    boltzmann_threshold = 0.596979866/300.0*temp
    boltzmann_diameter = find_diameter(x, v, boltzmann_threshold)
    #if boltzmann_diameter:
    #    print(f"Boltzmann hard sphere diameter at v(r) = k_B T: {boltzmann_diameter:.6f}")
    #else:
    #    print("Boltzmann hard sphere diameter not found within the given range.")
    
    # Calculate Stillinger diameter (v(r) = 0.5)
    stillinger_threshold = 0.5
    stillinger_diameter = find_diameter(x, v, stillinger_threshold)
    #if stillinger_diameter:
    #    print(f"Stillinger diameter at v(r) = 0.5: {stillinger_diameter:.6f}")
    #else:
    #    print("Stillinger diameter not found within the given range.")
    print(boltzmann_diameter, stillinger_diameter)

if __name__ == "__main__":
    main()
