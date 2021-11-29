import statsmodels.api as smapi
import numpy as np
import math
import json
import os
import time
lowess = smapi.nonparametric.lowess
# 文件路径信息
set_displacement_path = './{isample}/SingleParticleAnalysis/wi{wi}/FlowDisplacement/'
set_velocity_path = './{isample}/SingleParticleAnalysis/wi{wi}/FlowAveVelocity/'
set_displacement_opath = './{isample}/SingleParticleAnalysis/wi{wi}/FlowDisplacement_correctBox/'
set_velocity_opath = './{isample}/SingleParticleAnalysis/wi{wi}/FlowAveVelocity_correctBox/'
displacement_prefix = 'FlowDisplacement.'
velocity_prefix = 'FlowAveVelocity.'
velocity_oprefix = 'FlowAveVelocity_correctBox.'
# 样本起点与终点
sample_start = 1
sample_stop = 10
# 时刻起点与终点
istart = 1
istop = 10


def main():
    wi = 10
    tau_alpha = 165000
    gamma_dot = wi / tau_alpha
    JSON_SETTINGS_FNAME = "PaAn_PerMoment.json"
    # with open(JSON_SETTINGS_FNAME,encoding='utf-8') as fjson
    root = load_json(JSON_SETTINGS_FNAME)
    bandlife_savePath = "./"
    bandlife_fname = "bandlife"
    bandlife_file = bandlife_savePath + bandlife_fname
    output_path = root["output_path"]
    output_v_prefix = "FlowAveVelocity."
    start_step = 0
    delta_step = 6600
    stop_step = 660000
    delta_gamma_step = 66000
    if (stop_step - start_step) % delta_step != 0:
        raise Exception(
            'step not match! start_step:{} stop_stp:{} delta_step:{} '.format(
                start_step, stop_step, delta_step))
    num_fragments = int((stop_step - start_step) / delta_step - 1)
    bandlife = np.zeros((num_fragments, 2))
    if not os.path.exists("./logFiles"):
        os.makedirs("logFiles")
    ifrag = 0
    while True:
        flowFile = output_path + "FlowAveVelocity/" + output_v_prefix
        t0_step = ifrag * delta_step + start_step
        t_step = t0_step + delta_gamma_step
        os.system(
            "./SingleParticleFlowPerMoment.exe {} {} flag_special_t0 > ./logFiles/SingleParticleFlowPerMoment.{}_{}_{}.log"
            .format(t0_step, t_step, t0_step, t_step, delta_step))
        flowFile = flowFile + "{}".format(t_step)
        bandlife[ifrag, 0] = t_step
        bandlife[ifrag,
                 1] = compute_smoothed_Sigma_from_file(flowFile, gamma_dot)
        ifrag += 1
        if t_step >= stop_step:
            break

    os.system("rm {}*".format(output_path + "FlowAveVelocity/" +
                              output_v_prefix))
    np.savetxt(bandlife_file, bandlife, delimiter=',')
    return


def load_json(fname):
    with open(fname, 'r') as fjson:
        root = json.load(fjson)
        return root


def compute_smoothed_Sigma_from_file(fname_v,
                                     rate,
                                     lgrad_lo=0,
                                     lgrad_hi=19.31,
                                     data_flow_col=8,
                                     data_grad_col=3,
                                     vmap_ghost_frac=0.5,
                                     vmap_smooth_frac=0.3):
    lgrad = lgrad_hi - lgrad_lo
    vmap_smooth_frac_ghost = vmap_smooth_frac / (1 + vmap_ghost_frac * 2)

    data_v = correctGradBox_v(fname_v, rate, lgrad=lgrad)
    data_v = gen_ghost_layer(data_v, vmap_ghost_frac, lgrad=lgrad, rate=rate)
    smoothed_vmap = lowess(data_v[:, data_flow_col],
                           data_v[:, data_grad_col],
                           frac=vmap_smooth_frac_ghost,
                           it=0)
    smoothed_vmap = rm_ghost_layer(smoothed_vmap, lgrad_lo, lgrad_hi)
    # np.savetxt("test.smoothed.vmap", smoothed_vmap)
    sigma = computeSigma(smoothed_vmap, rate, lgrad=lgrad)
    # print(sigma)
    return sigma


def computeSigma(vmap_2col, rate, lgrad=19.31):
    std_vmap = np.zeros(vmap_2col.shape[0])
    std_vmap[:] = (vmap_2col[:, 0] - lgrad / 2.) * rate
    sum_sqrdiff = 0.
    for i in range(vmap_2col.shape[0]):
        diff = vmap_2col[i, 1] - std_vmap[i]
        sum_sqrdiff += square(diff)
    ave = sum_sqrdiff / vmap_2col.shape[0]
    ave /= np.square(lgrad / 2. * rate)
    ave = math.sqrt(ave)
    return ave


def correctGradBox_v(fname_veclocity,
                     rate,
                     lgrad=19.31,
                     __DATA_GRADBOX_COL=6,
                     __DATA_FLOW_COL=8,
                     __DATA_SKIP_ROW=9):

    # print("Processing velocity")

    velocity = np.loadtxt(fname_veclocity, skiprows=__DATA_SKIP_ROW)

    velocity[:,
             __DATA_FLOW_COL] = velocity[:,
                                         __DATA_FLOW_COL] - velocity[:,
                                                                     __DATA_GRADBOX_COL] * lgrad * rate
    # print("Processing displacement")
    '''
    if fname_displacement != "NULL":
        displacement = np.loadtxt(displacement_path + fname_displacement,
                                  skiprows=__DIS_SKIP_ROW)
        displacement[:,
                     __DIS_GRAD_COL] = displacement[:,
                                                    __DIS_GRAD_COL] - displacement[:,
                                                                                   __BOX_GRAD_COL] * __L_GRAD * gamma_dot
    '''
    return velocity


def gen_ghost_layer(flow_data_origin,
                    ghost_frac,
                    lgrad,
                    rate,
                    __DATA_FLOW_COL=8,
                    __DATA_GRAD_COL=3,
                    __DATA_ID_COL=0,
                    __TOTAL_PARTICLE=4320):
    '''gen ghost layer'''

    __GHOST_LAYER_WIDTH = ghost_frac * lgrad

    _total_row, _total_col = flow_data_origin.shape
    for row in range(0, _total_row):
        irow_grad_pos = flow_data_origin[row, __DATA_GRAD_COL]
        '''copy particles at bottom of grad direction onto top of box'''
        if irow_grad_pos <= __GHOST_LAYER_WIDTH:
            new_row = flow_data_origin[row, :].copy()
            new_row = new_row.reshape(1, _total_col)
            new_row[0, __DATA_ID_COL] += __TOTAL_PARTICLE
            new_row[0, __DATA_GRAD_COL] += lgrad
            new_row[0, __DATA_FLOW_COL] += lgrad * rate
            flow_data_origin = np.vstack([flow_data_origin, new_row])
        '''copy particle at top of grad direction below bottom of box'''
        if irow_grad_pos >= (lgrad - __GHOST_LAYER_WIDTH):
            new_row = flow_data_origin[row, :].copy()
            new_row = new_row.reshape(1, _total_col)
            new_row[0, __DATA_ID_COL] += __TOTAL_PARTICLE
            new_row[0, __DATA_GRAD_COL] -= lgrad
            new_row[0, __DATA_FLOW_COL] -= lgrad * rate
            flow_data_origin = np.vstack([flow_data_origin, new_row])
    return flow_data_origin


def rm_ghost_layer(flow_data_ghost, lgrad_lo, lgrad_hi):
    '''delete ghost layer'''
    VMAP_GRADPOS_COL = 0
    ghost_atoms_list = []
    row_looper = 0
    _total_row, _total_col = flow_data_ghost.shape
    while row_looper < _total_row:
        row_grad_pos = flow_data_ghost[row_looper, VMAP_GRADPOS_COL]
        if row_grad_pos < lgrad_lo or row_grad_pos > lgrad_hi:
            ghost_atoms_list.append(row_looper)
        row_looper += 1
    flow_data_ghost = np.delete(flow_data_ghost, ghost_atoms_list, axis=0)
    return flow_data_ghost


def square(x):
    return x * x


if __name__ == "__main__":
    main()