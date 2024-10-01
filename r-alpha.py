from polygon_tally import *

# importing data
with open("temp-check.txt", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        num_sides = float(row[0])
        polygon_radius = float(row[1])
        x = float(row[2])
        y = float(row[3])

# testing

# parameters
r = (x**2 + y**2)**(1/2) 
theta = np.arctan2(y, x)

# creating the class
k1 = KBasis(0, 0, num_sides, polygon_radius)

r_alpha = k1.variable_radius_xy(x, y)
rho = r / r_alpha

print("\n----- Python Output (Known to be Correct) -----")
print("r: ", r)
print("theta: ", theta)
print("variable radius: ", r_alpha)
print("rho: ", rho)
