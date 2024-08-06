from glip_glorp import *
from square_integration import *

# checking that z00 is 1 for any rho,phi
def test_zernike00():
    a = random.random()
    b = random.random()

    assert ana_zernike(np.array([a]), np.array([2*pi*b]), 0, 0)[0] == 1


### analytical integration tests
def ana_knm(n,m):        
    integrated_k = nquad(k_nm, [[-side_length,side_length],[-side_length,side_length]], args=(n,m))[0]

    return integrated_k

# checking K_00 should be equal to the area
def test_ana_k00():
    integrated_k = ana_knm(0,0)
    area = (2*side_length)**2
    
    assert round(integrated_k) == round(area)

# checking K_11 should be equal to 0
def test_ana_k11():
    integrated_k = ana_knm(1,1)

    assert round(integrated_k) == 0

# checking K_1n1 should be equal to 0
def test_ana_k1n1():
    integrated_k = ana_knm(1,-1)

    assert round(integrated_k) == 0

### numerical integration tests
def num_knm(n,m): 
    rho,phi = gen_rhophi(1000)
    r,theta = gen_rtheta(rho,phi,4,R0)
    x,y = gen_xy(rho,phi,4,R0)
    
    # to integrate
    K = k_m(n,m,r,theta,4,R0)
    
    dr,dtheta = gen_polar_differentials(r,theta)
    dmu = r * dr * dtheta

    return np.sum(K * dmu)

# checking that k00 should be equal to the area of the square
def test_num_k00():
    integrated_k = num_knm(0,0)
    area = (2*R0/2**(1/2))**2 
    assert round(integrated_k) == round(area)

# checking that k11 should be 0
def test_num_k11():
    integrated_k = num_knm(1,1)
    assert round(integrated_k) == 0

# checking that k1,-1 should be 0
def test_num_k1n1():
    integrated_k = num_knm(1,-1)
    assert round(integrated_k) == 0

### Integrating K numerically and analytically
''' Comment this out to remove numerical integration checks
# checking numerical integrated k_nm values
SSk_arr = []

for n in range(51):
    for m in np.arange(-n,n+1,2):
        row = []
        row.append(n)
        row.append(m)

        integrated_k = num_knm(n,m)
        row.append(integrated_k)
        
        SSk_arr.append(row)
        print("n: ", n, "m: ", m, "integral: ", integrated_k)

with open('numerical_k_integration.txt', 'w') as file:
    writer = csv.writer(file)
    writer.writerows(SSk_arr)
'''

''' Comment this out to remove analytical integration checks
# checking the analytical integrated k_nm values

for n in range(51):
    for m in np.arange(-n,n+1,2):
        integrated_k = ana_knm(n,m)
        print("n:",n, " m:",m, " integral:", integrated_k)

        b = str(n) + ',' + str(m) + ',' + str(integrated_k) + '\n'

        f = open('analytical_k_integration.txt', 'a')
        f.write(b)
        f.close()
'''

### checking ck values
def test_ck_constant_f():
    ck = ck_nm(0,0,lambda x,y : 1)
    
    assert round(ck) == 1

#''' Comment this out to remove calculation of ck both numerically and analytically
rho,phi = gen_rhophi(2000)
r,theta = gen_rtheta(rho,phi,p,R0)

for n in range(51):
    for m in np.arange(-n,n+1,2):
        ana_ck = ck_nm(n,m)
        num_ck = ck(n,m,r,theta,p,R0)
        print("n:",n, " m:",m, "| ana:",ana_ck, " num:",num_ck)

        row = str(n) + ',' + str(m) + ',' + str(ana_ck) + ',' + str(num_ck) + '\n'
        
        f = open("ck_num_vs_ana.txt", 'a')
        f.write(row)
        f.close()
#'''
