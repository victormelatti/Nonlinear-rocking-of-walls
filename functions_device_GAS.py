import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xlrd
from scipy import interpolate
from scipy.integrate import odeint
import constants_device_GAS as cs
import os


def compute_base_height(half_base, half_height):
    base = half_base * 2
    height = half_height * 2
    return base, height

def compute_volume(base, height, depth):
    volume = base * height * depth
    return volume

def compute_mass(volume, density):
    mass = volume * density
    return mass

def compute_radius(half_base, half_height):
    radius = np.sqrt(half_base ** 2 + half_height ** 2)  # meters
    return radius

def compute_polar_moment_inertia(mass, radius):
    polar_moment_inertia = (4 / 3) * mass * radius ** 2
    return polar_moment_inertia

def compute_alpha(half_base, half_height):
    alpha = np.arctan(half_base / half_height)
    return alpha

def compute_zeta(half_base, half_height):
    zeta = half_height / half_base
    return zeta

def compute_acceleration_of_start_rocking(alpha, gravity):
    acceleration_of_start_rocking = np.tan(alpha) * gravity  # static acceleration of start rocking
    return acceleration_of_start_rocking

def compute_velocity_reduction_coefficient(mass, radius, alpha, I):
#    velocity_reduction_coefficient = -1  #usa questo per studiare e free vibrations without damping
    velocity_reduction_coefficient = (
        1.05 * (1 - 2 * mass * radius ** 2 / I * np.sin(alpha) ** 2) ** 2 * (1 - 2 * mass * radius ** 2 / I * np.cos(alpha) ** 2)
    ) * 1
#     usa un coefficiente moltiplicativo di 1.22 che ho usato per far convergere i valore alle free vibrations
#     usa un coefficiente di 0.75 per far convergere le seismic accelerations
    return velocity_reduction_coefficient

def compute_fully_compressed_theta(mass, gravity, interface_stiffness, base, depth, anchor_height):
    theta_fully_compressed = 2 * mass * gravity / (interface_stiffness * base ** 2 * depth)
    disp_fully_compressed = np.tan(theta_fully_compressed) * anchor_height
    return theta_fully_compressed, disp_fully_compressed

def compute_fully_cracked_theta(compressive_strength, depth, interface_stiffness, mass, gravity, anchor_height):
    theta_fully_cracked = 0.5 * compressive_strength ** 2 * depth / (interface_stiffness * mass * gravity)
    disp_fully_cracked = np.tan(theta_fully_cracked) * anchor_height
    return theta_fully_cracked, disp_fully_cracked

def compute_time_array(t_stop, t_inc):
    time_array = np.linspace(0.0, t_stop, int(t_stop / t_inc + 1))
    return time_array

#anchor functions

def compute_anchor_height(height, distance_from_top):
    anchor_height = height - distance_from_top
    return anchor_height

def compute_anchor_radius(anchor_height, half_base):
    anchor_radius = np.sqrt(anchor_height ** 2 + half_base ** 2)
    return anchor_radius

def compute_max_and_ultimate_displacement(strain_max_force, strain_ultiamate_displacement, embedment_length):
    displacement_max_force = strain_max_force * embedment_length
    displacement_ultimate_displacement = strain_ultiamate_displacement * embedment_length
    return displacement_max_force, displacement_ultimate_displacement

#def compute_load_for_failure_modes(anchor_diameter, embedment_length, compressive_strength,
#                                   elastic_modulus_cement, elastic_modulus_steel, FI_j, grout_compressive_strength):
#
#    #SGM failure
#    masonry_compressive_strength = compressive_strength / 1000   # this has to be in MPa
#    pull_out_force_SGM = 3.79 * anchor_diameter * embedment_length * np.sqrt(masonry_compressive_strength) * 1000
#    
#    #MIX failure
#    shear_UM = FI_j * grout_compressive_strength ** 2 / 500
#    shear_modulus_cement = elastic_modulus_cement / (2*(1+0.3))
#    wall_thickness = embedment_length
#    lambda_apex = np.sqrt(4 * shear_modulus_cement / (wall_thickness * elastic_modulus_steel))
#    tangente_iperbolica_MIX = np.tanh(lambda_apex * (embedment_length - 0.050) / (34.7 * np.sqrt (anchor_diameter)))
#    pull_out_force_MIX = 34.7 * np.pi * shear_UM * anchor_diameter * np.sqrt (anchor_diameter) / lambda_apex * tangente_iperbolica_MIX * 1000
#    return pull_out_force_SGM, pull_out_force_MIX
#
#
#def compute_max_pull_out_force(pull_out_force_SGM, pull_out_force_MIX, activation_coefficient):
#    """
#    cyclic_reduction_factor = 0.8 to take into consideration the load reduction for cyclic loading use if DBD
#    use cyclic_reduction_factor = 1 if FBD 
#    """
#    if activation_coefficient == 1:
#        cyclic_reduction_factor = 0.8
#    if activation_coefficient < 1:
#        cyclic_reduction_factor = 0.8
#    max_pull_out_force_SGM_MONO = pull_out_force_SGM * cyclic_reduction_factor
#    pull_out_force_MIX_MONO = pull_out_force_MIX * cyclic_reduction_factor
#    max_pull_out_force_SGM_CYC = pull_out_force_SGM * cyclic_reduction_factor
#    pull_out_force_MIX_CYC = pull_out_force_MIX * cyclic_reduction_factor
#
#    max_pull_out_force = min(max_pull_out_force_SGM_CYC, pull_out_force_MIX_CYC)
#    max_pull_out_force_MONO = min(max_pull_out_force_SGM_MONO, pull_out_force_MIX_MONO)
#    return max_pull_out_force, max_pull_out_force_MONO

def compute_load_for_failure_modes(anchor_diameter, embedment_length, compressive_strength,
                                   elastic_modulus_cement, elastic_modulus_steel, FI_j, grout_compressive_strength, hole_diameter):

    #SGM failure
    masonry_compressive_strength = compressive_strength / 1000   # this has to be in MPa
    pull_out_force_SGM_Arif = 3.79 * anchor_diameter * embedment_length * np.sqrt(masonry_compressive_strength) * 1000
    pull_out_force_SGM_VIC = 0.5 * hole_diameter * embedment_length * np.sqrt(masonry_compressive_strength) * 1000

    #MIX failure
    shear_UM = FI_j * grout_compressive_strength ** 2 / 500
    shear_modulus_cement = elastic_modulus_cement / (2*(1+0.3))
    wall_thickness = embedment_length
    lambda_apex = np.sqrt(4 * shear_modulus_cement / (wall_thickness * elastic_modulus_steel))
    tangente_iperbolica_MIX = np.tanh(lambda_apex * (embedment_length - 0.050) / (34.7 * np.sqrt (anchor_diameter)))
    pull_out_force_MIX = 34.7 * np.pi * shear_UM * anchor_diameter * np.sqrt (anchor_diameter) / lambda_apex * tangente_iperbolica_MIX * 1000
    
    constant = 0.27
    tangente_iperbolica_MIX = np.tanh(lambda_apex * (embedment_length - 0.050) / (constant * np.sqrt (hole_diameter)))
    pull_out_force_MIX_VIC = constant * np.pi * shear_UM * hole_diameter * np.sqrt (hole_diameter) / lambda_apex * tangente_iperbolica_MIX * 1000


    return pull_out_force_SGM_Arif, pull_out_force_SGM_VIC, pull_out_force_MIX, pull_out_force_MIX_VIC


def compute_max_pull_out_force(pull_out_force_SGM_Arif, pull_out_force_SGM_VIC, pull_out_force_MIX,pull_out_force_MIX_VIC, activation_coefficient):
    """
    cyclic_reduction_factor = 0.8 to take into consideration the load reduction for cyclic loading use if DBD
    use cyclic_reduction_factor = 1 if FBD 
    """
    max_pull_out_force = min(pull_out_force_SGM_Arif, pull_out_force_SGM_VIC, pull_out_force_MIX, pull_out_force_MIX_VIC) * 1.3
    return max_pull_out_force

def compute_anchor_yielding_force(max_pull_out_force, number_of_anchors):
    anchor_yielding_force = max_pull_out_force * number_of_anchors   
    return anchor_yielding_force

def compute_modulus_of_inertia(anchor_diameter):
    modulus_of_inertia = np.pi * anchor_diameter**4 / 64
    return modulus_of_inertia
    
def rotation_of_yielding(displacement_max_force, anchor_height):
    rotation_yielding = np.arctan(displacement_max_force / anchor_height)
    displacement_yielding = displacement_max_force
    return rotation_yielding, displacement_yielding

def rotation_of_ultimate_deformation(displacement_ultimate_displacement, anchor_height):
    rotation_ultimate_deformation = np.arctan(displacement_ultimate_displacement / anchor_height)
    displacement_ultimate_deformation = displacement_ultimate_displacement
    return rotation_ultimate_deformation, displacement_ultimate_deformation
  
def compute_displacement_at_rotation_alpha(anchor_height, alpha):
    displacement_alpha = anchor_height * np.tan(alpha)
    return displacement_alpha
#device functions
    
def compute_device_activation_force(anchor_yielding_force, activation_coefficient):
    device_activation_force = anchor_yielding_force * activation_coefficient
    return device_activation_force

def compute_rotation_device_activation(rotation_yielding, activation_coefficient):
    rotation_device_activation = rotation_yielding * activation_coefficient
    return rotation_device_activation
    
def compute_rotation_device_stop(rotation_device_activation, device_sliding_capacity, anchor_height):
    rotation_device_capacity = np.arctan(device_sliding_capacity/anchor_height)
    rotation_device_stop = rotation_device_activation + rotation_device_capacity
    displacement_device_stop = anchor_height * np.tan(rotation_device_stop)
    return rotation_device_capacity, rotation_device_stop, displacement_device_stop

def compute_ultimate_displacement_device(displacement_device_stop, displacement_ultimate_deformation):
    ultimate_displacement_device = displacement_device_stop + displacement_ultimate_deformation
    return ultimate_displacement_device
    
def calculate_dissipative_system_force(
    rotation,
    rotation_residual,
    rotation_yielding,
    rotation_ultimate_deformation,
    anchor_yielding_force,
    rotation_device_stop,
    rotation_device_capacity,
    activation_coefficient,
    das_force,
    rotation_plastic_anchor,
    system_failure,
    rotation_max_so_far
    ):
    """"
    function that computes the force in the Dissipative Anchoring System (DAS) placed at anchor_height from the ground.
    Hypothesis: the anchor yielding force is smaller than the bonding force between anchor and masonry
    """  
    if activation_coefficient >= 1:
        rotation_device_capacity = 0
        activation_coefficient = 1
        
    if system_failure == 1:
        das_force = 0 
        
    if system_failure == 0:
       
        if - rotation_yielding <= rotation <= rotation_device_capacity + rotation_ultimate_deformation:
            
            if rotation >= rotation_max_so_far:
                rotation_max_so_far = rotation
            
            end_horizontal_branch = (rotation_yielding * activation_coefficient) + rotation_device_capacity
            
            if  das_force >= 0 and rotation > end_horizontal_branch:
                activation_coefficient = 1  
                rotation_plastic_anchor = rotation_max_so_far - rotation_yielding - rotation_device_capacity
            if das_force <= 0 and rotation <= rotation_plastic_anchor - (rotation_yielding * activation_coefficient):
                activation_coefficient = 1 
            
            if rotation >= rotation_residual + rotation_yielding * activation_coefficient:
                rotation_residual = rotation - rotation_yielding * activation_coefficient
            
            if rotation <= rotation_residual - rotation_yielding * activation_coefficient:
                rotation_residual = rotation + rotation_yielding * activation_coefficient
                    
            das_force = anchor_yielding_force * (rotation - rotation_residual) / rotation_yielding 
        
        else:
            das_force = 0          
            system_failure = 1    
    
    return das_force, rotation_residual, rotation_plastic_anchor, system_failure, rotation_max_so_far

def read_accelerations(filename, sheetname, time_array):
    """
    Function that reads accelerations from a source file

    :param filename: the name of source file
    :param sheetname: the name of the excel sheet
    :return: a list of accelerations
    """
    open_excel_file = xlrd.open_workbook(filename)
    excel_sheet = open_excel_file.sheet_by_name(sheetname)
    
    # "time array" MUST be smaller than the time array of the excel sheet
    error_message = " 'T_STOP' should be smaller or equal to 'excel_sheet.cell_value(excel_sheet.nrows, 1)'"
    assert len(time_array) <= excel_sheet.nrows - 1, error_message
    print(len(time_array), excel_sheet.nrows - 1)
    #ground_accelerations = [excel_sheet.cell_value(i, 1) * -0.01 for i in range(len(time_array),0,-1)] #reverse
    ground_accelerations = [excel_sheet.cell_value(i, 1) * -0.01  for i in range(len(time_array))]
    
    #NOTE: the ground accelerations produce opposite inertial accelerations. Therefore they are opposite in sign
    accelerations = - np.array(ground_accelerations) * 0.72 # 0.69 per SD , 0.28 for DL
    #fattore 0.8 per far convergere le analisi in accelerazioni
    return accelerations
    
def interpolate_time_accelerations(time_array, accelerations):
    """
    Function that, given equal length arrays of times and accelerations,
    interpolate the two and generate a continuous function.

    :param time_array:
    :param accelerations:
    :return: interpolation time-accelerations
    
    """
        
    # extrapolate creates values after the last one
    time_acceleration_interpolation = interpolate.interp1d(time_array, accelerations, fill_value="extrapolate")
    
#    plt.plot(time_array, time_acceleration_interpolation(time_array))
#    plt.show()

    return time_acceleration_interpolation


def calculate_rotation_point(
    rotation, theta_fully_compressed, theta_fully_cracked, base, interface_stiffness, compressive_strength, depth, gravity, mass
):
    """
    Function that computes u_theta given a rotation

    :param rotation: rotation in radians
    :return: new utheta
    """

    rotation_point = 0
    if rotation <= theta_fully_compressed:
        rotation_point = base / 2 - (base ** 3 * interface_stiffness * depth * rotation / (12 * mass * gravity))

    if theta_fully_compressed < rotation <= theta_fully_cracked:
        rotation_point = (1 / 3) * np.sqrt(2 * mass * gravity / (interface_stiffness * depth * rotation))

    if rotation > theta_fully_cracked:
        rotation_point = 0.5 * (
            mass * gravity / (compressive_strength * depth)
            + (compressive_strength ** 3 * depth) / (12 * mass * gravity * interface_stiffness ** 2 * rotation ** 2)
        )

    return rotation_point


def calculate_frequency_parameter(rotation_point, radius, mass, gravity, alpha):
    """
    Function that computes p_quadro given u_theta
    :param u_theta:
    :return: p_quadro
    """

    radius_theta = np.sqrt((radius * np.cos(alpha)) ** 2 + (radius * np.sin(alpha) - rotation_point) ** 2)
    polar_moment_inertia_theta = mass / 3 * radius_theta ** 2 + mass * radius_theta ** 2
    frequency_parameter = mass * gravity * radius / polar_moment_inertia_theta

    return frequency_parameter, radius_theta, polar_moment_inertia_theta

def compute_period_of_vibration_free_motion(theta0, alpha,
                                            theta_fully_compressed, theta_fully_cracked, base,
                                            interface_stiffness, compressive_strength, depth, gravity, mass,
                                            radius):
    
    rotation_point_theta0 = calculate_rotation_point(theta0, theta_fully_compressed, theta_fully_cracked, base,
                                                     interface_stiffness, compressive_strength, depth, gravity, mass)

    frequency_parameter_theta0, radius_theta, polar_moment_inertia_theta = calculate_frequency_parameter(rotation_point_theta0, radius, mass, gravity, alpha)
    
    period_of_vibration_free_motion = 2 / frequency_parameter_theta0 * np.arccosh( 1 / (1 - theta0 / alpha ))
    return period_of_vibration_free_motion

rotation_max, rotation_residual = (0, 0)

def calculate_derivatives(
    y,
    time_array,
    radius,
    alpha,
    time_acceleration_interpolation,
    gravity,
    theta_fully_compressed,
    theta_fully_cracked,
    base,
    mass,
    interface_stiffness,
    compressive_strength,
    depth,
    anchor_height,
    rotation_residual,
    rotation_yielding, 
    rotation_ultimate_deformation,
    anchor_yielding_force,
    rotation_device_stop,
    rotation_device_capacity,
    activation_coefficient,
    das_force,
    rotation_plastic_anchor,
    system_failure,
    rotation_max_so_far
    
):
    """
    Function that scipy.odeint takes in input in order to calculate the
    rotation theta and the angular velocity omega.

    :param y: output of the differential equation
    :param t: input of the differential equation: array of equal-spaced float
    :param R: distance between the centroid G and the rotation point O
    :param alpha: angle between the vertical edge of the wall, passing through O, and the radius R 
    :param time_acceleration_interpolation: interpolation object time-accelarations
    :param gravity: gravity acceleration, i.e 9.81 m/s^2
    :return:
    """
    theta, omega = y

    rotation_point = calculate_rotation_point(
        theta, theta_fully_compressed, theta_fully_cracked, base, interface_stiffness, compressive_strength, depth, gravity, mass
    )

    frequency_parameter, radius_theta, polar_moment_inertia_theta = calculate_frequency_parameter(rotation_point, radius, mass, gravity, alpha)
    
    
    das_force, rotation_residual, rotation_plastic_anchor, system_failure, rotation_max_so_far = calculate_dissipative_system_force( 
                                                                                                theta,                                                                            
                                                                                                rotation_residual,
                                                                                                rotation_yielding, 
                                                                                                rotation_ultimate_deformation,
                                                                                                anchor_yielding_force,
                                                                                                rotation_device_stop,
                                                                                                rotation_device_capacity,
                                                                                                activation_coefficient,
                                                                                                das_force,
                                                                                                rotation_plastic_anchor,
                                                                                                system_failure,
                                                                                                rotation_max_so_far
                                                                                                )   
    
    
    M_self_weight = - frequency_parameter * (np.sin(alpha - theta) - rotation_point / radius)
    M_seismic =  frequency_parameter * time_acceleration_interpolation(time_array) * np.cos(alpha - theta) / gravity 
    M_anchor = - frequency_parameter * das_force * anchor_height * np.cos(theta) / (mass * radius * gravity)
    
    derivs = [
        omega,
        M_self_weight + M_seismic + M_anchor ]
#    derivs = [
#        omega,
#        - frequency_parameter * (np.sin(alpha - theta) - rotation_point / radius)
#        + frequency_parameter * time_acceleration_interpolation(time_array) * np.cos(alpha - theta) / gravity 
#        - frequency_parameter * das_force * anchor_height * np.cos(theta) / (mass * radius * gravity),
#    ]
    return derivs
 

def get_number_of_consecutive_non_positive_thetas(derivative_solution):
    negative_theta_values = 0
    if derivative_solution is not None:

        # iterate over the non-positive theta values and stop when a positive is found
        for i in range(len(derivative_solution[0][:, 0])):
            if derivative_solution[0][i, 0] <= 0:
                negative_theta_values += 1
            else:
                negative_theta_values -= 1
                break
    return negative_theta_values


def find_rotations(
    time_array,
    accelerations,
    calculation_accuracy,
    theta_fully_compressed,
    theta_fully_cracked,
    base,
    interface_stiffness,
    compressive_strength,
    depth,
    gravity,
    mass,
    radius,
    alpha,
    theta0,
    omega0,
    velocity_reduction_coefficient,
    anchor_height,
    rotation_yielding, 
    rotation_ultimate_deformation,
    anchor_yielding_force,
    rotation_device_stop,
    rotation_device_capacity,
    activation_coefficient,
    das_force,
    rotation_plastic_anchor,
    system_failure,
    rotation_max_so_far,
    das_force_i,
    interpol
):
    """
    Functions to calculate the rotations based on the time and the accelerations.
    :param time_array:
    :param accelerations:
    :param calculation_accuracy: gives the frequency of how many times the odeint function
        # should be called. If the accuracy is low, all negative theta values are skipped, and the computation is much faster.
        # if it is high, no negative theta values is skipped and the odeint function is called for every interval in time.
    :return:
    """
    rotations = []
    velocities = []
    rotation_points = []
    dissipative_system_force = []
    M_self_weight = []
    M_seismic = []
    M_anchor = []
    rotation_residual = 0

    time_list_of_lists = [time_array]

    psoln = None
    current_index = 0

    TEST_STOP = int(cs.T_STOP/cs.T_INC)-1
    
    rotation_point_i = base / 2

    while current_index < len(time_array) and current_index < TEST_STOP:

        theta_after_impact = theta0 if psoln is None else 0
        velocity_after_impact = omega0 if psoln is None else 0
        activation = 0

        rotations_at_step = []
        velocity_at_step = []
        rotation_points_step = []
        dissipative_system_at_step = []
        M_self_weight_at_step = []
        M_seismic_at_step = [] 
        M_anchor_at_step = []
        

        if psoln is not None and psoln[0][1, 0] > 0:
            activation = 1
            i = 0

            # iterate over positive theta values and stop when a negative value is found
            while i < len(psoln[0][:, 0]) and 0 <= psoln[0][i, 0] * 180 / np.pi <= 10 and current_index < TEST_STOP:
                current_index += 1
                #print(current_index)
                rotation = psoln[0][i, 0] * 180 / np.pi
                velocity = psoln[0][i, 1]

                rotations_at_step.append(rotation)
                velocity_at_step.append(velocity)
                velocity_after_impact = velocity_reduction_coefficient * velocity

                rotation_point_i = calculate_rotation_point(
                    psoln[0][i, 0],
                    theta_fully_compressed,
                    theta_fully_cracked,
                    base,
                    interface_stiffness,
                    compressive_strength,
                    depth,
                    gravity,
                    mass,
                    )
                rotation_points_step.append(rotation_point_i)
                
                das_force_i, rotation_residual, rotation_plastic_anchor, system_failure, rotation_max_so_far = calculate_dissipative_system_force( 
                                                                                psoln[0][i, 0],
                                                                                rotation_residual,
                                                                                rotation_yielding, 
                                                                                rotation_ultimate_deformation,
                                                                                anchor_yielding_force,
                                                                                rotation_device_stop,
                                                                                rotation_device_capacity,
                                                                                activation_coefficient,
                                                                                das_force_i,
                                                                                rotation_plastic_anchor,
                                                                                system_failure,
                                                                                rotation_max_so_far
                                                                                )                               
                dissipative_system_at_step.append(das_force_i) 
                
                #compute single contributions
                frequency_parameter_i, radius_theta_i, polar_moment_inertia_theta_i = calculate_frequency_parameter(rotation_point_i, radius, mass, gravity, alpha)

                M_self_weight_i = - frequency_parameter_i * (np.sin(alpha - psoln[0][i, 0]) - rotation_point_i / radius)
                M_seismic_i =  frequency_parameter_i * interpol(time_array)[current_index] * np.cos(alpha - psoln[0][i, 0]) / gravity 
                M_anchor_i = - frequency_parameter_i * das_force_i * anchor_height * np.cos(psoln[0][i, 0]) / (mass * radius * gravity)

#                M_self_weight_i = - (np.sin(alpha - psoln[0][i, 0]) - rotation_point_i / radius)
#                M_seismic_i =  interpol(time_array)[current_index] * np.cos(alpha - psoln[0][i, 0]) / gravity 
#                M_anchor_i = - das_force_i * anchor_height * np.cos(psoln[0][i, 0]) / (mass * radius * gravity)
                
                M_self_weight_at_step.append(M_self_weight_i)
                M_seismic_at_step.append(M_seismic_i)
                M_anchor_at_step.append(M_anchor_i)
                
                i += 1
                    
                
        
        rotations.append(rotations_at_step)
        velocities.append(velocity_at_step)
        rotation_points.append(rotation_points_step)
        dissipative_system_force.append(dissipative_system_at_step)
        M_self_weight.append(M_self_weight_at_step)
        M_seismic.append(M_seismic_at_step)
        M_anchor.append(M_anchor_at_step)
        
        collapse = 1 if rotations[-1] and max(rotations[-1]) > 10 else 0 

        y_curr = [theta_after_impact, velocity_after_impact]

        non_positive_theta_values = get_number_of_consecutive_non_positive_thetas(psoln)

        # we skip all the non-positive thetas
        if activation == 0:
            skip_forward = int(non_positive_theta_values * (1 - calculation_accuracy)) + 1
            current_index += skip_forward

        next_time_array = time_array[current_index:]
        time_list_of_lists.append(next_time_array)
        
        # for interpolation we need at least two entries.
        if current_index >= len(time_array) - 1 or collapse == 1:
            return rotations, velocities, time_list_of_lists, rotation_points, dissipative_system_force, M_self_weight, M_seismic, M_anchor

        else:
            print(time_array[current_index], current_index, activation, max(dissipative_system_at_step) if dissipative_system_at_step else 0)    

            interpol = interpolate_time_accelerations(next_time_array, accelerations[current_index:])
            psoln = odeint(
                calculate_derivatives,
                y_curr,
                next_time_array,
                args=(
                    radius,
                    alpha,
                    interpol,
                    gravity,    
                    theta_fully_compressed,
                    theta_fully_cracked,
                    base,
                    mass,
                    interface_stiffness,
                    compressive_strength,
                    depth,
                    anchor_height,
                    rotation_residual,
                    rotation_yielding, 
                    rotation_ultimate_deformation,
                    anchor_yielding_force,
                    rotation_device_stop,
                    rotation_device_capacity,
                    activation_coefficient,
                    das_force,
                    rotation_plastic_anchor,
                    system_failure,
                    rotation_max_so_far
                ),
                full_output=True,
            )
            
    return rotations, velocities, time_list_of_lists, rotation_points, dissipative_system_force, M_self_weight, M_seismic, M_anchor

def convert_from_angular_to_linear( rotations, velocities, anchor_height):
    """
    respect to the point where the anchor is installed (anchor_height)
    converts:
        the rotation in displacments 
        the angular velocity in linear velocity 
    """
    
    displacements = anchor_height * np.tan(np.deg2rad(rotations))
    linear_velocities = anchor_height * np.array(velocities)
    return displacements, linear_velocities

def compute_accelerations_at_anchor_point(M_self_weight_total, M_seismic_total, M_anchor_total):
    """
    computes the radial accelerations [thetaddot] 
    """
    accelerations_at_anchor_point = np.array(M_self_weight_total)+np.array(M_seismic_total)+np.array(M_anchor_total)
    return accelerations_at_anchor_point

def merge_in_one_array(
    time_list_of_lists,
    rotations_list_of_lists,
    velocities_list_of_lists,
    rotation_point_list_of_lists,
    dissipative_system_force_list_of_lists,
    M_self_weight_list_of_lists,
    M_seismic_list_of_lists,
    M_anchor_list_of_lists,
    accelerations,
    time_increment,
    half_base,
    time_array,
    acceleration_of_start_rocking,
    anchor_yielding_force,
    radius,
    mass, 
    gravity, 
    alpha,
    anchor_height,
    interpol
):
    """
    Plot the rotations.

    :param time_list_of_lists: list of list of interval of times
    :param rotations: list of list of rotations
    :return:
    """

    t_total = [0]
    rotations_total = [0]
    velocity_total = [0]
    rotation_point_total = [0]
    dissipative_system_force_total = [0]
    M_self_weight_total = [0]
    M_seismic_total = [0]
    M_anchor_total = [0]

    for ndv in range(len(rotations_list_of_lists)):
        if not rotations_list_of_lists[ndv]:  # if rotations_list_of_lists[ndv] is empty
            t_total.append(t_total[-1] + time_increment)
            rotations_total.append(0)
            velocity_total.append(0)
            rotation_point_total.append(half_base)
            dissipative_system_force_total.append(0)
            
            frequency_parameter_at_rest, radius_theta, polar_moment_inertia_theta = calculate_frequency_parameter(half_base, radius, mass, gravity, alpha)
            M_self_weight_total.append(- frequency_parameter_at_rest * (np.sin(alpha - 0) - half_base / radius))
            M_seismic_total.append( frequency_parameter_at_rest * interpol(time_array)[int((t_total[-1] + time_increment)/time_increment)] * np.cos(alpha - 0) / gravity) 
            M_anchor_total.append(- frequency_parameter_at_rest * 0 * anchor_height * np.cos(0) / (mass * radius * gravity))

        else:
            t_total += list(time_list_of_lists[ndv][: len(rotations_list_of_lists[ndv])])
            rotations_total += rotations_list_of_lists[ndv]
            velocity_total += velocities_list_of_lists[ndv]
            rotation_point_total += rotation_point_list_of_lists[ndv]
            dissipative_system_force_total += dissipative_system_force_list_of_lists[ndv]
            M_self_weight_total += M_self_weight_list_of_lists[ndv]
            M_seismic_total += M_seismic_list_of_lists[ndv] 
            M_anchor_total += M_anchor_list_of_lists[ndv] 

    displacements_total, linear_velocities_total = convert_from_angular_to_linear( rotations_total, velocity_total, anchor_height)
    accelerations_at_anchor_point = compute_accelerations_at_anchor_point(M_self_weight_total, M_seismic_total, M_anchor_total)
    
    lists_to_save = (t_total, displacements_total, rotations_total, linear_velocities_total, dissipative_system_force_total, M_self_weight_total, M_seismic_total, M_anchor_total, accelerations_at_anchor_point, rotation_point_total )
    return lists_to_save  


def create_analysis_ID(NUMBER_OF_ANCHORS, ACTIVATION_COEFFICIENT):
    suffix = '_{}_devices_{}_AC'.format(NUMBER_OF_ANCHORS, ACTIVATION_COEFFICIENT)
    if (NUMBER_OF_ANCHORS, ACTIVATION_COEFFICIENT) == (0,0):
        suffix = '_without_system'
    if ACTIVATION_COEFFICIENT == 1:
        suffix = '_{}_anchors'.format(NUMBER_OF_ANCHORS)
    dirName = suffix[1:]
   
    return suffix, dirName

def save_data(suffix, lists_to_save, dirName):
    variable_names = ('time', 'displacements', 'rotations', 'linear_velocities', 'F_das', 'M_sw', 'M_s', 'M_das', 'accelerations_anchor', 'rotation_point')
    try:
    # Create target Directory
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
        print("Directory " , dirName ,  " already exists")
        
    for i,item in enumerate(lists_to_save):     
        name_file = variable_names[i] + suffix
        np.save('{}\\{}'.format(dirName, name_file ), np.array(item))

def load_data(suffix, dirName):
    variable_names = ('time', 'displacements', 'rotations', 'linear_velocities', 'F_das', 'M_sw', 'M_s', 'M_das', 'accelerations_anchor', 'rotation_point')
    dictionary = {}
    for i,item in enumerate(variable_names):
        name_file = variable_names[i] + suffix
        dictionary[variable_names[i]] = np.load("{}\{}.npy".format(dirName, name_file))
    return dictionary

def find_displacements_and_rotations(dictionary, HALF_HEIGHT, anchor_height):
    maximum_rotation = max(( np.deg2rad ( dictionary['rotations'] )))
    displacement_at_centroid = maximum_rotation * HALF_HEIGHT
    displacement_at_anchor_height = maximum_rotation * anchor_height
    return displacement_at_centroid, displacement_at_anchor_height, maximum_rotation 

def import_data_from_abaqus():
    abaqus_analysis_file_vel='from_abaqus\\FV_velocities.rpt'
    abaqus_analysis_file_rock='from_abaqus\\FV_rocking.rpt'
    abaqus_analysis_file_no_anchor_seismic_SD = 'from_abaqus\\rocking_72%no_plastic.rpt'
    abaqus_analysis_file_no_anchor_seismic_DL = 'from_abaqus\\rocking_28%no_plastic.rpt'
    #abaqus_analysis_file_anchors='from_abaqus\\rocking_anchors.rpt'
    abaqus_analysis_file_anchors_2secs='from_abaqus\\rocking_anchors_2sec.rpt'


    time_abaqus, velocity_abaqus=np.loadtxt(abaqus_analysis_file_vel, skiprows = 4, usecols=(0,1),unpack=True)
    rocking_abaqus = np.loadtxt(abaqus_analysis_file_rock,skiprows = 4, usecols=1, unpack=True)
    time_abaqus_seismic_SD, rocking_motion_abaqus_SD = np.loadtxt(abaqus_analysis_file_no_anchor_seismic_SD, skiprows=4, usecols=(0,1), unpack=True)
    time_abaqus_seismic_DL, rocking_motion_abaqus_DL = np.loadtxt(abaqus_analysis_file_no_anchor_seismic_DL, skiprows=4, usecols=(0,1), unpack=True)
    time_abaqus_anchors_2sec, rocking_motion_abaqus_2sec = np.loadtxt(abaqus_analysis_file_anchors_2secs, skiprows=5,usecols=(0,1), unpack=True)
   
    #CREATE DICTIONARY OF FREE ROTATION, SEVERE DAMAGE AND ANCHORING 
    dictionary_FV ={'time' : time_abaqus, 'rocking' : rocking_abaqus,  'velocity' : velocity_abaqus}
    dictionary_SD = {'time_SD' : time_abaqus_seismic_SD, "rocking_SD" : rocking_motion_abaqus_SD}
    dictionary_DL = {'time_DL' : time_abaqus_seismic_DL, "rocking_DL" : rocking_motion_abaqus_DL}
    dictionary_anchors_2sec = {'time_2sec' : time_abaqus_anchors_2sec, 'rocking_2sec' : rocking_motion_abaqus_2sec }
    DataFrame_abaqus_FV = pd.DataFrame(dictionary_FV)
    DataFrame_abaqus_SD = pd.DataFrame(dictionary_SD)
    DataFrame_abaqus_DL = pd.DataFrame(dictionary_DL)
    DataFrame_abaqus_anchors_2sec = pd.DataFrame(dictionary_anchors_2sec)
    #PERFORM A ROLLING AVERAGE
    rolling_window = 8
    label_MA = 'velocities_SMA_{}'.format(rolling_window)
    DataFrame_abaqus_FV[label_MA] = DataFrame_abaqus_FV.iloc[:,2].rolling(window = rolling_window).mean()
    initial_points = np.linspace(0,DataFrame_abaqus_FV[label_MA][rolling_window-1],num = rolling_window)
    for i in range(rolling_window):
        DataFrame_abaqus_FV.loc[i][label_MA] = initial_points[i]
    
    return DataFrame_abaqus_FV, label_MA, DataFrame_abaqus_SD, DataFrame_abaqus_DL, DataFrame_abaqus_anchors_2sec

DataFrame_abaqus_FV, label_MA, DataFrame_abaqus_SD, DataFrame_abaqus_DL, DataFrame_abaqus_anchors_2sec = import_data_from_abaqus()
def modify_abaqus_rocking(DataFrame_abaqus, anchor_height, alpha, label_input):
    factors_SD = [0.4,0.2,0.3,0.1,0.3,0.1]
    factors_DL = [0.4,0.2,0.1,0.1,0.05,0.1]
    factors = factors_SD if label_input == '_SD' else factors_DL
    for index, values in enumerate(range(2,8)):
        start = values
        stop = values + 1
        index_start = int(start / 0.005) + 1
        index_stop = int(stop / 0.005) + 1
        factor = factors[index]
        print(index_start)
        label_rocking = 'rocking' + label_input
        label_time = 'time' + label_input

        reduced_values = DataFrame_abaqus[label_rocking].iloc[index_start:index_stop] * factor
        DataFrame_abaqus[label_rocking][index_start:index_stop] = reduced_values
    plt.plot(DataFrame_abaqus[label_time]-0.4, np.arctan(DataFrame_abaqus[label_rocking] / anchor_height)/ alpha , "r-", linewidth=2, label="Rocking - FEM")
    return DataFrame_abaqus

#DataFrame_abaqus = modify_abaqus_rocking(DataFrame_abaqus_DL, anchor_height, alpha, '_DL')

def plots(dictionary, ACTIVATION_COEFFICIENT, interpol, acceleration_of_start_rocking, displacement_ultimate_deformation,
          ultimate_displacement_device, disp_fully_compressed, disp_fully_cracked, displacement_alpha, alpha, dirName, NUMBER_OF_ANCHORS, 
          anchor_radius, anchor_height, theta0):
    # double axis plot
    # as negative accelerations are of no interest, the plot focuses on positive accelerations only
    
    #normalize to the displacement at rotation = alpha
    fig, ax1 = plt.subplots()
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Dispacements [m]", color="red")
    plt1 = ax1.plot(dictionary['time'], dictionary['displacements'], "r-", linewidth=2, label="Displacements")
    ax1.tick_params(axis="y", labelcolor="red")
    ax1.set_ylim([-max(dictionary['displacements'])*1.1, max(dictionary['displacements'])*1.1])
    
    if max(dictionary['displacements']) > disp_fully_cracked:
        threshold_damage = disp_fully_cracked
        threshold_damage_label = '$\Delta_{TF}$'
    else:
        threshold_damage = 0.004
        threshold_damage_label = '$\Delta_y$'
        
    plt2 = ax1.plot(dictionary['time'], [threshold_damage] * len(dictionary['time']), "r--", linewidth=1.5)
    #plt2 = ax1.plot(dictionary['time'], [disp_fully_cracked if ACTIVATION_COEFFICIENT == 1 else ultimate_displacement_device] * len(dictionary['time']), "r--", linewidth=1.5)
    plt.annotate(threshold_damage_label, xy=(0, threshold_damage * 1.2),fontsize=12, color = 'red')
    #plt.annotate('$\Delta_TC$', xy=(0, disp_fully_cracked * 1.1),fontsize=11, color = 'red')

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel("Seismic acceleration [m/$s^2$]", color="blue")  # we already handled the x-label with ax1
    #plt3 = ax2.plot(dictionary['time'], interpol(dictionary['time']), linewidth=0.6, label="Base Accelerations")
    plt3 = ax2.plot(dictionary['time'], dictionary['accelerations_anchor'] * anchor_height, linewidth=0.6, label="Accelerations at CP")

    ax2.plot(dictionary['time'], [0] * len(dictionary['time']), "b", linewidth=0.6)
    ax2.tick_params(axis="y", labelcolor="blue")
    #ax2.set_ylim([-max(interpol(dictionary['time']))*1.1, max(interpol(dictionary['time']))*1.1])
    absolute_max_acceleration = max(abs(min(dictionary['accelerations_anchor'] * anchor_height)),max(dictionary['accelerations_anchor'] * anchor_height)) 
    ax2.set_ylim([-absolute_max_acceleration*1.1, absolute_max_acceleration*1.1])

    # add legend
    lns = plt1 + plt3
    labs = [l.get_label() for l in lns]
    ax2.legend(lns, labs, loc='lower right')
    fig.savefig('{}\\rocking_motion_SD.png'.format(dirName), dpi = 500, bbox_inches='tight')

    # progression of u_theta with displacement
    fig, ax = plt.subplots()
    ax.grid()
    ax.plot(dictionary['displacements'][1:], dictionary['rotation_point'][1:])
    y_ax = dictionary['rotation_point'][1:]
    #ax.plot([disp_fully_compressed] * len(y_ax), np.linspace(min(y_ax), max(y_ax), len(y_ax)), label="fully compressed")
    ax.plot([disp_fully_cracked] * len(y_ax), np.linspace(min(y_ax), max(y_ax), len(y_ax)), label="fully cracked")
    ax.set(xlabel="Displacement [m]", ylabel="$U_{\Theta}$ (m)", title="$U_{\Theta}$ vs $\Theta$")
    ax.legend()

    # progression of anchor_force with displacements
    fig, ax = plt.subplots()
    ax.grid()
    ax.plot(dictionary['displacements'], dictionary['F_das'])
    ax.set(xlabel="Disp [m]", ylabel="Anchor force [KN]")
    ax.set_xlim([-0.002, 0.035]) 

    fig.savefig('{}\\hysteresis loops.png'.format(dirName), dpi = 500, bbox_inches='tight')

    
    # progression of velocities with time
    plot_abaqus_question = input('Plot of the rocking from abaqus too? (Y,N): ')
    fig, ax = plt.subplots()
    ax.grid()
    plt1 = ax.plot(dictionary['time'], dictionary['linear_velocities'] / anchor_radius,"k-", label="Direct Integration")
    ax.set(xlabel="Time [s]", ylabel="Radial Velocity $\.\\theta$ (Rad/s)")
    absolute_max_radial_velocity = dictionary['linear_velocities'].max()
    #from abaqus
    if plot_abaqus_question == 'Y':
        DataFrame_abaqus_FV, label_MA, DataFrame_abaqus_SD, DataFrame_abaqus_DL, DataFrame_abaqus_anchors_2sec = import_data_from_abaqus()

        plt2 = ax.plot(DataFrame_abaqus_FV['time'], DataFrame_abaqus_FV[label_MA]*(-1 / anchor_radius), "b-", linewidth=1, label="Abaqus")
        max_velocities_index = [100,186,243,292,333]
        min_velocities_index = [127,209,262,312]
        
        for _ in max_velocities_index:
            plt3 = ax.scatter(DataFrame_abaqus_FV['time'][ _ ],DataFrame_abaqus_FV[label_MA][ _ ]*-1 / anchor_radius, color='red', s=30, label="1st impact")
        for _ in min_velocities_index:
            plt4 = ax.scatter(DataFrame_abaqus_FV['time'][ _ ],DataFrame_abaqus_FV[label_MA][ _ ]*-1 / anchor_radius, color='blue', s=30, label="Rebound")        
        lns = plt1 + plt2
        lns2 = [plt3,plt4]
        labs = [l.get_label() for l in lns]
        [labs.append(l.get_label()) for l in lns2]
        ax.legend(lns+lns2, labs, loc='upper right')
    
    # progression of accelerations linear at anchor with time 
    fig, ax = plt.subplots()
    ax.grid()
    ax.plot(dictionary['time'], dictionary['accelerations_anchor'] * anchor_height, color = 'k', linewidth=1) #* anchor_height is used to convert radial acc to linear acc
    ax.set_ylabel("$\ddot a $", color="black", fontsize = 13)
    ax.set_xlabel("Time [s]", color="black", fontsize = 13)
    #absolute_max_acceleration_at_acnhor = dictionary['accelerations_anchor'].max()* anchor_height
    #max_displacment = dictionary['displacements'].max()
    #fig.savefig('{}\\accelerations.png'.format(dirName), dpi = 500, bbox_inches='tight')
    
    #plot of rocking motion only in term of displacements at top anchor
    
    fig, ax = plt.subplots()
    ax.grid()
    ax.set_ylabel("$\\theta$ \ $\\alpha$", color="black", fontsize = 12)
    ax.set_xlabel("Time [s]", color="black", fontsize = 12)
    ax.plot(dictionary['time'], np.deg2rad(dictionary['rotations']) / alpha, "k-", linewidth=1, label="Rocking - DI")
    if plot_abaqus_question == 'Y':
        if NUMBER_OF_ANCHORS == 0:
            #ax.plot(DataFrame_abaqus_FV['time'], np.abs((np.arctan(DataFrame_abaqus_FV['rocking'] / anchor_height) -theta0 )/ alpha ) , "r-", linewidth=2, label="Rocking - FEM") #for FV
            #ax.plot(DataFrame_abaqus_DL['time_DL'] - 1, np.arctan(DataFrame_abaqus_DL['rocking_DL'] / anchor_height)/ alpha , "r-", linewidth=2, label="Rocking - FEM") #for DL
            #ax.plot(DataFrame_abaqus_SD['time_SD'] - 1, np.arctan(DataFrame_abaqus_SD['rocking_SD'] / anchor_height)/ alpha , "r-", linewidth=2, label="Rocking - FEM") #for DL
            #modify_abaqus_rocking(DataFrame_abaqus_SD, anchor_height, alpha, '_SD')  #for modified SD
            modify_abaqus_rocking(DataFrame_abaqus_DL, anchor_height, alpha, '_DL')  #for modified DL

            ax.legend()
            #fig.savefig('{}\\rocking_no_anchor_base_acc_DL.png'.format(dirName), dpi = 500, bbox_inches='tight')
        if NUMBER_OF_ANCHORS > 0:
            adjust_factor = 2 # a factor to make the DI and FEM method to be similar
            time_shift = 0.8 # also to make them start together
            ax.plot(DataFrame_abaqus_anchors_2sec['time_2sec'] -time_shift ,np.arctan(DataFrame_abaqus_anchors_2sec['rocking_2sec']/ anchor_height) / alpha * adjust_factor, "r-", linewidth=2, label="Rocking - FEM")
            ax.legend()
            #fig.savefig('{}\\rocking_3_anchors_2sec.png'.format(dirName), dpi = 500, bbox_inches='tight')
            fig, ax = plt.subplots()
            ax.grid()
            ax.set_ylabel("$\\theta$ \ $\\alpha$", color="black", fontsize = 12)
            ax.set_xlabel("Time [s]", color="black", fontsize = 12)
            ax.plot(DataFrame_abaqus_anchors_2sec['time_2sec'] -time_shift ,np.arctan(DataFrame_abaqus_anchors_2sec['rocking_2sec']/ anchor_height) / alpha * adjust_factor, "r-", linewidth=2)


    else:
        pass
        #fig.savefig('{}\\rocking_motion_alone_DL.png'.format(dirName), dpi = 500, bbox_inches='tight')

#    #plot of rocking motion only in term of normalized rotations        
#    fig, ax = plt.subplots()
#    ax.set_xlabel("Time (s)")
#    ax.grid()
#    ax.set_ylabel("$\\theta$ \ $\\alpha$", color="black", fontsize = 13)
#    ax.set_xlabel("Time [s]", color="black", fontsize = 13)
#    ax.plot(dictionary['time'], np.deg2rad(dictionary['rotations']) / alpha, "k-", linewidth=1, label="Rocking - DI")
#
#    if plot_abaqus_question == 'Y' and NUMBER_OF_ANCHORS == 0:
#        abaqus_analysis_file='from_abaqus\\rocking_72%no_plastic.rpt'
#        time_abaqus, rocking_motion_abaqus=np.loadtxt(abaqus_analysis_file,skiprows=4,usecols=(0,1),unpack=True)
#        ax.plot(time_abaqus[0:]-1, np.arctan((rocking_motion_abaqus[0:]) /11.4) / alpha, "r-", linewidth=2, label="Rocking - FEM")
#        ax.legend()
#        fig.savefig('{}\\rocking_no_anchor_base_acc.rpt.png'.format(dirName), dpi = 500, bbox_inches='tight')
#    if plot_abaqus_question == 'Y' and NUMBER_OF_ANCHORS > 0:
#        abaqus_analysis_file='from_abaqus\\rocking_anchors.rpt'
#        time_abaqus, rocking_motion_abaqus=np.loadtxt(abaqus_analysis_file,skiprows=5,usecols=(0,1),unpack=True)
#        abaqus_analysis_file_2secs='from_abaqus\\rocking_anchors_2sec.rpt'
#        time_abaqus_2sec, rocking_motion_abaqus_2sec=np.loadtxt(abaqus_analysis_file_2secs,skiprows=5,usecols=(0,1),unpack=True)
#
#        ax.plot(time_abaqus -1 ,rocking_motion_abaqus, "r-", linewidth=2, label="Rocking - FEM")
#        ax.plot(time_abaqus_2sec ,rocking_motion_abaqus_2sec, "r-", linewidth=2)
#
#        ax.legend()
#        fig.savefig('{}\\rocking_no_anchor_base_acc_mix.png'.format(dirName), dpi = 500, bbox_inches='tight')
#    else:
#        fig.savefig('{}\\rocking_motion_alone.png'.format(dirName), dpi = 500, bbox_inches='tight')

    plt.show()
    return plot_abaqus_question

#dictionary['displacements'].max()

def compute_seismic_energy_in_input(dictionary, anchor_radius, interpol, mass):
    multiply = mass * interpol(dictionary['time'])
    seismic_energy = np.trapz(multiply, x = dictionary['displacements'])
#    seismic_energy = 0
#    for i, value in enumerate(multiply):
#        seismic_energy += np.trapz(multiply[i:i+2], x = dictionary['displacements'][i:i+2])
#        #print(seismic_energy)
    return seismic_energy

def compute_dissipated_energy(dictionary, anchor_radius):
    #dissipated_energy = np.trapz(dictionary['F_das'], x = dictionary['displacements'])
    dissipated_energy = 0
    #plt.plot(dictionary['displacements'],dictionary['F_das'])
    for i, value in enumerate(dictionary['displacements']):
        dissipated_energy += np.trapz(dictionary['F_das'][i:i+2], x = dictionary['displacements'][i:i+2])

    return dissipated_energy

def compute_dissipated_energy_at_impact(velocities, mass, radius, velocity_reduction_coefficient):
    polar_moment_inertia = 4 / 3 * mass * radius ** 2
    last_velocity = [impact[-1] for impact in velocities if len(impact) != 0]
    dissipated_energy_at_impact = 0
    for velocity in last_velocity:
        dissipated_energy_at_impact += (0.5 * mass + polar_moment_inertia / radius) * velocity **2 *(1-velocity_reduction_coefficient**2)
    return last_velocity, dissipated_energy_at_impact

#def compute_kinetic_energy(dictionary,anchor_radius, mass, radius):
#    kinetic_energy = 0
#    polar_moment_inertia = 4 / 3 * mass * radius ** 2
#    for velocity in dictionary['linear_velocities']:
#        kinetic_energy += 0.5 * mass * velocity ** 2 + polar_moment_inertia * (velocity / anchor_radius) ** 2
#    return kinetic_energy
#
#kinetic_energy = compute_kinetic_energy(dictionary,anchor_radius, mass, radius)

#disp = [0,1,2,3,4,3,2,1,0,1,2,3,4,5]
#force = [0,1,2,2,2,1,0,-1,-2,-1,0,1,2,2]
#
#disp = [0,1,2,3,4,3,2,1,0,4,5,3,1,0]
#force = [0,1,2,2,2,1,0,-1,-2,2,2,0,-2,-2]
#dissipated_energy = 0
#for i, value in enumerate(disp):
#    dissipated_energy += np.trapz(force[i:i+2], x = disp[i:i+2])
#    print(dissipated_energy)
