import constants_device_GAS as cs
import functions_device_GAS as fn
import numpy as np

__name__ = "__main__"
if __name__ == "__main__":
    print("This file solves the equation of motion of a rocking block with anchor")

    time_array = fn.compute_time_array(cs.T_STOP, cs.T_INC)

    #accelerations = fn.read_accelerations("acceleration_laquila_checking.xlsx", "sec_4", time_array)
    accelerations = fn.read_accelerations("acceleration_laquila_checking.xlsx", "test_acc_increasing", time_array)

    base, height = fn.compute_base_height(cs.HALF_BASE, cs.HALF_HEIGHT)

    volume = fn.compute_volume(base, height, cs.DEPTH)

    mass = fn.compute_mass(volume, cs.DENSITY)

    radius = fn.compute_radius(cs.HALF_BASE, cs.HALF_HEIGHT)

    polar_moment_inertia = fn.compute_polar_moment_inertia(mass, radius)

    alpha = fn.compute_alpha(cs.HALF_BASE, cs.HALF_HEIGHT)

    zeta = fn.compute_zeta(cs.HALF_BASE, cs.HALF_HEIGHT)

    acceleration_of_start_rocking = fn.compute_acceleration_of_start_rocking(alpha, cs.GRAVITY)

    velocity_reduction_coefficient = fn.compute_velocity_reduction_coefficient(mass, radius, alpha, polar_moment_inertia)

    interpol = fn.interpolate_time_accelerations(time_array, accelerations)

    anchor_height = fn.compute_anchor_height(height, cs.DISTANCE_FROM_TOP)
    
    anchor_radius = fn.compute_anchor_radius(anchor_height, cs.HALF_BASE)
    
    displacement_alpha = fn.compute_displacement_at_rotation_alpha(anchor_height, alpha)
        
    theta_fully_compressed, disp_fully_compressed = fn.compute_fully_compressed_theta(mass, cs.GRAVITY, cs.INTERFACE_STIFFNESS, base, cs.DEPTH, anchor_height)

    theta_fully_cracked, disp_fully_cracked = fn.compute_fully_cracked_theta(cs.COMPRESSIVE_STRENGTH, cs.DEPTH, cs.INTERFACE_STIFFNESS, mass, cs.GRAVITY, anchor_height)
        
    displacement_max_force, displacement_ultimate_displacement = fn.compute_max_and_ultimate_displacement(cs.STRAIN_MAX_FORCE, cs.STRAIN_ULTIMATE_DISPLACEMENT, cs.EMBEDMENT_LENGTH)
    
#    pull_out_force_SGM, pull_out_force_MIX = fn.compute_load_for_failure_modes(cs.ANCHOR_DIAMETER, cs.EMBEDMENT_LENGTH, cs.COMPRESSIVE_STRENGTH, cs.ELASTIC_MODULUS_CEMENT, cs.ELASTIC_MODULUS_STEEL, cs.FI_j, cs.GROUT_COMPRESSIVE_STRENGTH, cs.HOLE_DIAMETER)
#
#    max_pull_out_force, max_pull_out_force_MONO = fn.compute_max_pull_out_force(pull_out_force_SGM, pull_out_force_MIX, cs.ACTIVATION_COEFFICIENT)

    pull_out_force_SGM_ARIF, pull_out_force_SGM_VIC, pull_out_force_MIX, pull_out_force_MIX_VIC = fn.compute_load_for_failure_modes(cs.ANCHOR_DIAMETER, cs.EMBEDMENT_LENGTH, cs.COMPRESSIVE_STRENGTH,
                                                                                             cs.ELASTIC_MODULUS_CEMENT, cs.ELASTIC_MODULUS_STEEL, cs.FI_j, cs.GROUT_COMPRESSIVE_STRENGTH, cs.HOLE_DIAMETER)
    
    max_pull_out_force = fn.compute_max_pull_out_force(pull_out_force_SGM_ARIF, pull_out_force_SGM_VIC, pull_out_force_MIX,pull_out_force_MIX_VIC, cs.ACTIVATION_COEFFICIENT)
   
    anchor_yielding_force = fn.compute_anchor_yielding_force(max_pull_out_force, cs.NUMBER_OF_ANCHORS)
    
    modulus_of_inertia = fn.compute_modulus_of_inertia(cs.ANCHOR_DIAMETER)
        
    rotation_yielding, displacement_yielding = fn.rotation_of_yielding(displacement_max_force, anchor_height)
    
    rotation_ultimate_deformation, displacement_ultimate_deformation = fn.rotation_of_ultimate_deformation(displacement_ultimate_displacement, anchor_height)

    device_activation_force = fn.compute_device_activation_force(anchor_yielding_force, cs.ACTIVATION_COEFFICIENT)
    
    rotation_device_activation  = fn.compute_rotation_device_activation(rotation_yielding, cs.ACTIVATION_COEFFICIENT) 
    
    rotation_device_capacity, rotation_device_stop, displacement_device_stop = fn.compute_rotation_device_stop(rotation_device_activation, cs.DEVICE_SLIDING_CAPACITY, anchor_height)

    ultimate_displacement_device = fn.compute_ultimate_displacement_device(displacement_device_stop, displacement_ultimate_deformation)
    
    period_of_vibration_free_motion = fn.compute_period_of_vibration_free_motion(cs.THETA0, alpha, theta_fully_compressed, theta_fully_cracked, base,
                                            cs.INTERFACE_STIFFNESS, cs.COMPRESSIVE_STRENGTH, cs.DEPTH, cs.GRAVITY, mass, radius)  
    
    rotations, velocities, ts, rotation_point, dissipative_system_force, M_self_weight, M_seismic, M_anchor = fn.find_rotations(
        time_array,
        accelerations,
        0.8,
        theta_fully_compressed,
        theta_fully_cracked,
        base,
        cs.INTERFACE_STIFFNESS,
        cs.COMPRESSIVE_STRENGTH,
        cs.DEPTH,
        cs.GRAVITY,
        mass,
        radius,
        alpha,
        cs.THETA0,
        cs.OMEGA0,
        velocity_reduction_coefficient,
        anchor_height,
        rotation_yielding, 
        rotation_ultimate_deformation,
        anchor_yielding_force,
        rotation_device_stop,
        rotation_device_capacity,
        cs.ACTIVATION_COEFFICIENT,
        0,
        0,
        0,
        0,
        0,
        interpol        
    )

    lists_to_save = fn.merge_in_one_array(
        ts, rotations, velocities, rotation_point, dissipative_system_force, M_self_weight, M_seismic, M_anchor, accelerations, 
        cs.T_INC, cs.HALF_BASE, time_array, acceleration_of_start_rocking, anchor_yielding_force, radius, 
        mass, cs.GRAVITY, alpha, anchor_height, interpol)

    suffix, dirName = fn.create_analysis_ID(cs.NUMBER_OF_ANCHORS, cs.ACTIVATION_COEFFICIENT)

    #dirName = fn.save_data(suffix, lists_to_save)
    fn.save_data(suffix, lists_to_save, dirName)
    
    dictionary = fn.load_data(suffix, dirName)
    
    displacement_at_centroid, displacement_at_anchor_height, maximum_rotation= fn.find_displacements_and_rotations(dictionary, cs.HALF_HEIGHT, anchor_height)
    import functions_device_GAS as fn
       
    fn.plots(dictionary, cs.ACTIVATION_COEFFICIENT, interpol, acceleration_of_start_rocking,
            displacement_ultimate_deformation, ultimate_displacement_device, disp_fully_compressed,
            disp_fully_cracked, displacement_alpha, alpha, dirName, cs.NUMBER_OF_ANCHORS, anchor_radius, anchor_height, cs.THETA0)
    
    seismic_energy_in_input = fn.compute_seismic_energy_in_input(dictionary, anchor_radius, interpol, mass)
    
    dissipated_energy = fn.compute_dissipated_energy(dictionary, anchor_radius)

    last_velocity, dissipated_energy_at_impact = fn.compute_dissipated_energy_at_impact(velocities, mass, radius, velocity_reduction_coefficient)

    #time_hystory_max_disp_array = np.array([])
    #time_hystory_max_disp_array = np.append(time_hystory_max_disp_array, max(dictionary['displacements']))


#plt.plot(dictionary['time'],dictionary['displacements'])
#plt.plot(dictionary['time'], interpol(dictionary['time'])*0.05)

#alpha_0 = fn.compute_alpha_0(cs.A_G, cs.S, cs.CF, cs.E_STAR, cs.Q)
#hinge_indentation_FB = fn.compute_hinge_indentation_FB( mass, cs.GRAVITY, cs.COMPRESSIVE_STRENGTH, cs.DEPTH)
#tie_force_FB = fn.compute_tie_force_FB(mass, cs.GRAVITY, alpha_0, cs.HALF_BASE, cs.HALF_HEIGHT, hinge_indentation_FB, anchor_height)

#DENSITY*2*(HALF_BASE*HALF_HEIGHT)*DEPTH*GRAVITY    

#relevant points
#list_of_parameters_string = ['rotation_yielding',
#                             'rotation_ultimate_deformation',
#                             'rotation_device_activation',
#                             'rotation_device_capacity',
#                             'rotation_device_stop',
#                             'rotation_device_capacity + rotation_ultimate_deformation']
#
#list_of_parameters_float = [rotation_yielding,
#                            rotation_ultimate_deformation,
#                            rotation_device_activation,
#                            rotation_device_capacity,
#                            rotation_device_stop,
#                            rotation_device_capacity + rotation_ultimate_deformation]    
#dict_of_parameters = {}
#for i in range(len(list_of_parameters_string)):
#    dict_of_parameters[list_of_parameters_string[i]] = np.rad2deg(list_of_parameters_float[i])
#dict_of_parameters
#
#print(max(rotation_total))
#print(max(dissipative_system_force_total))
#
#acc = np.array(M_self_weight_total)+np.array(M_seismic_total)+np.array(M_anchor_total)
#
#plt.plot(t_total,acc )
#start = int(3/0.005)
#stop = int(4/0.005)
#plt.plot(t_total[start:stop],rotation_total[start:stop] )
#plt.plot(t_total[start:stop],M_self_weight_total[start:stop] )
#plt.plot(t_total[start:stop],M_seismic_total[start:stop] )
#plt.plot(t_total[start:stop],M_anchor_total[start:stop] )
#plt.scatter(rotation_total[start:stop],M_anchor_total[start:stop] )

