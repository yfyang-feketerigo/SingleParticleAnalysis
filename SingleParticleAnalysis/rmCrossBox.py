from ast import dump
import os
import numpy as np
displacement_path = './wi10/FlowDisplacement/'
velocity_path = './wi10/FlowAveVelocity/'
displacement_prefix = 'FlowDisplacement.'
velocity_prefix = 'FlowAveVelocity.'
displacement_oprefix = 'FlowDisplacement_nonCrossBox.'
velocity_oprefix = 'FlowAveVelocity_nonCrossBox.'
displacement_opath = './wi10/FlowDisplacement_nonCrossBox/'
if not os.path.exists(displacement_path):
    os.makedirs(displacement_path)
velocity_opath = './wi10/FlowAveVelocity_nonCrossBoundary/'
if not os.path.exists(velocity_opath):
    os.makedirs(velocity_opath)

istart = 1
istop = 50
delta_step = 4800

velocity_skip_row = 9
displacement_skip_row = 10

nCrossBox_col = 9
for i in range(istart, istop + 1):
    str_step = str(i * delta_step)
    fname_veclocity = velocity_prefix + str_step
    velocity = np.loadtxt(velocity_path + fname_veclocity,
                          skiprows=velocity_skip_row)
    row_total = np.shape(velocity)[0]
    row_coutner = 0
    row_pointer = 0
    while row_coutner != row_total:
        if velocity[row_pointer, nCrossBox_col] != 0:
            velocity = np.delete(velocity, row_pointer, axis=0)
        else:
            row_pointer += 1
        row_coutner += 1
    ofile_header = 'id type x y z ix iy iz flow_velocity grad_box_change '
    np.savetxt(velocity_opath + velocity_oprefix + str_step,
               velocity,
               fmt='%i %i %.6f %.6f %.6f %i %i %i %.6e %i',
               header=ofile_header)

nCrossBox_col = 9
for i in range(istart, istop + 1):
    str_step = str(i * delta_step)
    fname_displacement = displacement_prefix + str_step
    displacement = np.loadtxt(displacement_path + fname_displacement,
                              skiprows=displacement_skip_row)
    row_total = np.shape(displacement)[0]
    row_coutner = 0
    row_pointer = 0
    while row_coutner != row_total:
        if displacement[row_pointer, nCrossBox_col] != 0:
            displacement = np.delete(displacement, row_pointer, axis=0)
        else:
            row_pointer += 1
        row_coutner += 1
    ofile_header = 'id type x y z ix iy iz flow_displacement grad_box_change '
    np.savetxt(displacement_opath + displacement_oprefix + str_step,
               displacement,
               fmt='%i %i %.6f %.6f %.6f %i %i %i %.6e %i',
               header=ofile_header)