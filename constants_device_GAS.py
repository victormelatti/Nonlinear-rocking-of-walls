import numpy as np
#geometrical parameters of the wall:


HALF_BASE = 0.50     # meters for free vibrations
#HALF_BASE = 0.45     # meters for seismic acceleration SD
#HALF_BASE = 0.43     # meters for seismic acceleration LD


HALF_HEIGHT = 6.25    # meters for free vibrations
#HALF_HEIGHT = 6.5    # meters for seismic acceleration SD
#HALF_HEIGHT = 7    # meters for seismic acceleration DL

DEPTH = 6.85   # meters

#physical parameters of the wall:
DENSITY = 2.03    # KN/m^3  MASS DENSITY
INTERFACE_STIFFNESS = 6e6  # KPa/m equivalent to 6 MPa/mm
COMPRESSIVE_STRENGTH = 3.6e3  # Pa #KN/m^2 equivalent to 3.2 MPa

#T_STOP = 7.990   #sec
T_STOP = 9.990   #sec  for the paper I am using almost 10 sec of recordings because I am not comparing the results to Abaqus.
                        #for abaqus comparison use 7.990 sec

T_INC = 0.005    #sec
TEST_STOP = int(T_STOP/T_INC)-1
GRAVITY = 9.81    #m/s^2

#0.248 degree correspond to 5cm of displacment of the top anchor
THETA0 = 0.0 * np.pi / 180  # initial angular displacement in radians
OMEGA0 = 0.0 # initial angular velocity

#parameters for the anchor if TAS
#ANCHOR_LENGTH = 0.32   #m
#DISTANCE_FROM_TOP = 0.8
#EPSILON_YIELDING = 0.002
#EPSILON_ULTIMATE = 0.016
#SIGMA_YIELDING = 235e3*1.7   #KN/m^2
#ANCHOR_DIAMETER = 0.016    #m

#parameters for the anchor if GAS
EMBEDMENT_LENGTH = HALF_BASE * 2
ANCHOR_DIAMETER = 0.014    #m
DISTANCE_FROM_TOP = 0.8
HOLE_DIAMETER = 0.07   #m
ELASTIC_MODULUS_STEEL = 210e6    #KN/m^2
STRAIN_MAX_FORCE = 0.005
STRAIN_ULTIMATE_DISPLACEMENT = 0.01


#GROUT PROPERTIES
ELASTIC_MODULUS_CEMENT = 36.283e6    #KN/m^2
FI_j = 0.5
GROUT_COMPRESSIVE_STRENGTH = 50 #MPa # this has to be in MPa do not convert to KN/m^2


#BUCKLING PARAMETER
BETA = 0.5

#PARAMETERS FOR FORCE BASED DESIGN
A_G = 0.24
S = 1.15
CF = 1
E_STAR = 1
Q = 2

#DEVICE PARAMETERS
DEVICE_SLIDING_CAPACITY = 0.03
NUMBER_OF_ANCHORS = int(input('Enter Number of anchors (es. 1) : '))
ACTIVATION_COEFFICIENT = float(input('Enter sliding force as fraction of yielding force (es. 0.8) : '))







