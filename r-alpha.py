from polygon_tally import *

# importing data
with open("temp-check.txt", 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        num_sides = float(row[0])
        polygon_radius = float(row[1])
        x = float(row[2])
        y = float(row[3])

"""
# testing
num_sides = 5
polygon_radius = 10
x = -1
y = -1
"""

# parameters
def test():
    r = (x**2 + y**2)**(1/2) 
    theta = np.arctan2(y, x) + 2 * np.pi * (y < 0)

    # creating the class
    k1 = KBasis(0, 0, num_sides, polygon_radius)

    r_alpha = k1.variable_radius_xy(x, y)
    rho = r / r_alpha

    u_alpha = k1.u_alpha_theta(theta)

    """
    print("\n----- Python Input -----")
    print("Number of Sides: ", num_sides)
    print("Polygon Radius: ", polygon_radius)
    print("x: ", x)
    print("y: ", y)
    """

    print("\n----- Python Output (Correct) -----")
    print("r: ", r)
    print("theta: ", theta)
    print("variable radius: ", r_alpha)
    print("rho: ", rho)


test()
