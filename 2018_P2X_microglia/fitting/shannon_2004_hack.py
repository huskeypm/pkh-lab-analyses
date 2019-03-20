# Gotran generated code for the  "shannon_2004" model
from __future__ import division

def init_state_values(**values):
    """
    Initialize state values
    """
    # Imports
    import numpy as np
    from modelparameters.utils import Range

    # Init values
    # h=0.9867005, j=0.991562, m=0.001405627, Xr=0.008641386, Xs=0.005412034,
    # X_tos=0.004051574, Y_tos=0.9945511, R_tos=0.9946,
    # X_tof=0.004051574, Y_tof=0.9945511, d=7.175662e-06, f=1.000681,
    # fCaB_SL=0.01452605, fCaB_jct1=0.02421991, R=0.8884332,
    # O=8.156628e-07, I=1.024274e-07, Ca_TroponinC=0.008773191,
    # Ca_TroponinC_Ca_Mg=0.1078283, Mg_TroponinC_Ca_Mg=0.01524002,
    # Ca_Calmodulin=0.0002911916, Ca_Myosin=0.001298754,
    # Mg_Myosin=0.1381982, Ca_SRB=0.002143165, Na_jct1_buf=3.539892,
    # Na_SL_buf=0.7720854, Na_jct1=8.80329, Na_SL=8.80733, Nai=8.80853,
    # Ca_Calsequestrin=1.242988, Ca_SLB_SL=0.1110363,
    # Ca_SLB_jct1=0.009566355, Ca_SLHigh_SL=0.07297378,
    # Ca_SLHigh_jct1=0.007347888, Ca_SR=0.5545201, Ca_jct1=0.0001737475,
    # Ca_SL=0.0001031812, Cai=8.597401e-05, V=-85.56885
    init_values = np.array([0.9867005, 0.991562, 0.001405627, 0.008641386,\
        0.005412034, 0.004051574, 0.9945511, 0.9946, 0.004051574, 0.9945511,\
        7.175662e-06, 1.000681, 0.01452605, 0.02421991, 0.8884332,\
        8.156628e-07, 1.024274e-07, 0.008773191, 0.1078283, 0.01524002,\
        0.0002911916, 0.001298754, 0.1381982, 0.002143165, 3.539892,\
        0.7720854, 8.80329, 8.80733, 8.80853, 1.242988, 0.1110363,\
        0.009566355, 0.07297378, 0.007347888, 0.5545201, 0.0001737475,\
        0.0001031812, 8.597401e-05, -85.56885], dtype=np.float_)

    # State indices and limit checker
    state_ind = dict([("h",(0, Range())), ("j",(1, Range())), ("m",(2,\
        Range())), ("Xr",(3, Range())), ("Xs",(4, Range())), ("X_tos",(5,\
        Range())), ("Y_tos",(6, Range())), ("R_tos",(7, Range())),\
        ("X_tof",(8, Range())), ("Y_tof",(9, Range())), ("d",(10, Range())),\
        ("f",(11, Range())), ("fCaB_SL",(12, Range())), ("fCaB_jct1",(13,\
        Range())), ("R",(14, Range())), ("O",(15, Range())), ("I",(16,\
        Range())), ("Ca_TroponinC",(17, Range())), ("Ca_TroponinC_Ca_Mg",(18,\
        Range())), ("Mg_TroponinC_Ca_Mg",(19, Range())),\
        ("Ca_Calmodulin",(20, Range())), ("Ca_Myosin",(21, Range())),\
        ("Mg_Myosin",(22, Range())), ("Ca_SRB",(23, Range())),\
        ("Na_jct1_buf",(24, Range())), ("Na_SL_buf",(25, Range())),\
        ("Na_jct1",(26, Range())), ("Na_SL",(27, Range())), ("Nai",(28,\
        Range())), ("Ca_Calsequestrin",(29, Range())), ("Ca_SLB_SL",(30,\
        Range())), ("Ca_SLB_jct1",(31, Range())), ("Ca_SLHigh_SL",(32,\
        Range())), ("Ca_SLHigh_jct1",(33, Range())), ("Ca_SR",(34, Range())),\
        ("Ca_jct1",(35, Range())), ("Ca_SL",(36, Range())), ("Cai",(37,\
        Range())), ("V",(38, Range()))])

    for state_name, value in values.items():
        if state_name not in state_ind:
            raise ValueError("{0} is not a state.".format(state_name))
        ind, range = state_ind[state_name]
        if value not in range:
            raise ValueError("While setting '{0}' {1}".format(state_name,\
                range.format_not_in(value)))

        # Assign value
        init_values[ind] = value

    return init_values

def init_parameter_values(**values):
    """
    Initialize parameter values
    """
    # Imports
    import numpy as np
    from modelparameters.utils import Range

    # Param values
    # Cao=1.8, Cli=15, Clo=150, Cm=1.381e-10, F=96485, Ki=135, Ko=5.4, Mgi=1,
    # Nao=140, Rgas=8314.3, T=310, cell_length=100, cell_radius=10.25,
    # Fx_Na_SL=0.89, Fx_Na_jct1=0.11, G_INa=16, Fx_NaBk_SL=0.89,
    # Fx_NaBk_jct1=0.11, G_NaBk=0.000297, Fx_NaK_SL=0.89,
    # Fx_NaK_jct1=0.11, H_NaK=4, I_NaK_max=1.90719, Km_Ko=1.5,
    # Km_Nai=11, Fx_Ks_SL=0.89, Fx_Ks_jct1=0.11, pKNa=0.01833,
    # g_Kp=0.001, G_tos=0.06, G_tof=0.02, Fx_Cl_SL=0.89,
    # Fx_Cl_jct1=0.11, G_Cl=0.109625, Kd_ClCa=0.1, G_ClBk=0.009,
    # Fx_ICaL_SL=0.1, Fx_ICaL_jct1=0.9, PCa=0.00054, PK=2.7e-07,
    # PNa=1.5e-08, Q10_CaL=1.8, gamma_Cai=0.341, gamma_Cao=0.341,
    # gamma_Ki=0.75, gamma_Ko=0.75, gamma_Nai=0.75, gamma_Nao=0.75,
    # lccCaInact=1.7, Fx_NCX_SL=0.89, Fx_NCX_jct1=0.11, HNa=3,
    # K_mCai=0.00359, K_mCao=1.3, K_mNai=12.29, K_mNao=87.5,
    # Kd_act=0.000256, Q10_NCX=1.57, V_max_INaCa=9, eta=0.35,
    # ksat=0.27, Fx_SLCaP_SL=0.89, Fx_SLCaP_jct1=0.11, H_ICap=1.6,
    # Km=0.0005, Q10_SLCaP=2.35, V_maxAF=0.0673, Fx_CaBk_SL=0.89,
    # Fx_CaBk_jct1=0.11, G_CaBk=0.0002513, EC50_SR=0.45, HSR=2.5,
    # Max_SR=15, Min_SR=1, kiCa=0.5, kim=0.005, koCa=10, kom=0.06,
    # ks=25, KSRleak=5.348e-06, H_Jpump=1.787, Kmf=0.000246, Kmr=1.7,
    # Q10_SRCaP=2.6, V_max_Jpump=0.0053114, Bmax_Calsequestrin=0.14,
    # Bmax_SLB_SL=0.0374, Bmax_SLB_jct1=0.0046, Bmax_SLHigh_SL=0.0134,
    # Bmax_SLHigh_jct1=0.00165, koff_Calsequestrin=65, koff_SLB=1.3,
    # koff_SLHigh=0.03, kon_Calsequestrin=100, kon_SL=100,
    # Bmax_Calmodulin=0.024, Bmax_Myosin_Ca=0.14, Bmax_Myosin_Mg=0.14,
    # Bmax_SRB=0.0171, Bmax_TroponinC=0.07, Bmax_TroponinC_Ca_Mg_Ca=0.14,
    # Bmax_TroponinC_Ca_Mg_Mg=0.14, koff_Calmodulin=0.238,
    # koff_Myosin_Ca=0.00046, koff_Myosin_Mg=5.7e-05, koff_SRB=0.06,
    # koff_TroponinC=0.0196, koff_TroponinC_Ca_Mg_Ca=3.2e-05,
    # koff_TroponinC_Ca_Mg_Mg=0.00333, kon_Calmodulin=34,
    # kon_Myosin_Ca=13.8, kon_Myosin_Mg=0.0157, kon_SRB=100,
    # kon_TroponinC=32.7, kon_TroponinC_Ca_Mg_Ca=2.37,
    # kon_TroponinC_Ca_Mg_Mg=0.003, Bmax_SL=1.65, Bmax_jct1=7.561,
    # koff=0.001, kon=0.0001, stim_amplitude=9.5, stim_duration=5,
    # stim_period=1000, stim_start=100
    init_values = np.array([1.8, 15, 150, 1.381e-10, 96485, 135, 5.4, 1, 140,\
        8314.3, 310, 100, 10.25, 0.89, 0.11, 16, 0.89, 0.11, 0.000297, 0.89,\
        0.11, 4, 1.90719, 1.5, 11, 0.89, 0.11, 0.01833, 0.001, 0.06, 0.02,\
        0.89, 0.11, 0.109625, 0.1, 0.009, 0.1, 0.9, 0.00054, 2.7e-07,\
        1.5e-08, 1.8, 0.341, 0.341, 0.75, 0.75, 0.75, 0.75, 1.7, 0.89, 0.11,\
        3, 0.00359, 1.3, 12.29, 87.5, 0.000256, 1.57, 9, 0.35, 0.27, 0.89,\
        0.11, 1.6, 0.0005, 2.35, 0.0673, 0.89, 0.11, 0.0002513, 0.45, 2.5,\
        15, 1, 0.5, 0.005, 10, 0.06, 25, 5.348e-06, 1.787, 0.000246, 1.7,\
        2.6, 0.0053114, 0.14, 0.0374, 0.0046, 0.0134, 0.00165, 65, 1.3, 0.03,\
        100, 100, 0.024, 0.14, 0.14, 0.0171, 0.07, 0.14, 0.14, 0.238,\
        0.00046, 5.7e-05, 0.06, 0.0196, 3.2e-05, 0.00333, 34, 13.8, 0.0157,\
        100, 32.7, 2.37, 0.003, 1.65, 7.561, 0.001, 0.0001, 9.5, 5, 1000,\
        100], dtype=np.float_)

    # Parameter indices and limit checker
    param_ind = dict([("Cao", (0, Range())), ("Cli", (1, Range())), ("Clo",\
        (2, Range())), ("Cm", (3, Range())), ("F", (4, Range())), ("Ki", (5,\
        Range())), ("Ko", (6, Range())), ("Mgi", (7, Range())), ("Nao", (8,\
        Range())), ("Rgas", (9, Range())), ("T", (10, Range())),\
        ("cell_length", (11, Range())), ("cell_radius", (12, Range())),\
        ("Fx_Na_SL", (13, Range())), ("Fx_Na_jct1", (14, Range())), ("G_INa",\
        (15, Range())), ("Fx_NaBk_SL", (16, Range())), ("Fx_NaBk_jct1", (17,\
        Range())), ("G_NaBk", (18, Range())), ("Fx_NaK_SL", (19, Range())),\
        ("Fx_NaK_jct1", (20, Range())), ("H_NaK", (21, Range())),\
        ("I_NaK_max", (22, Range())), ("Km_Ko", (23, Range())), ("Km_Nai",\
        (24, Range())), ("Fx_Ks_SL", (25, Range())), ("Fx_Ks_jct1", (26,\
        Range())), ("pKNa", (27, Range())), ("g_Kp", (28, Range())),\
        ("G_tos", (29, Range())), ("G_tof", (30, Range())), ("Fx_Cl_SL", (31,\
        Range())), ("Fx_Cl_jct1", (32, Range())), ("G_Cl", (33, Range())),\
        ("Kd_ClCa", (34, Range())), ("G_ClBk", (35, Range())), ("Fx_ICaL_SL",\
        (36, Range())), ("Fx_ICaL_jct1", (37, Range())), ("PCa", (38,\
        Range())), ("PK", (39, Range())), ("PNa", (40, Range())), ("Q10_CaL",\
        (41, Range())), ("gamma_Cai", (42, Range())), ("gamma_Cao", (43,\
        Range())), ("gamma_Ki", (44, Range())), ("gamma_Ko", (45, Range())),\
        ("gamma_Nai", (46, Range())), ("gamma_Nao", (47, Range())),\
        ("lccCaInact", (48, Range())), ("Fx_NCX_SL", (49, Range())),\
        ("Fx_NCX_jct1", (50, Range())), ("HNa", (51, Range())), ("K_mCai",\
        (52, Range())), ("K_mCao", (53, Range())), ("K_mNai", (54, Range())),\
        ("K_mNao", (55, Range())), ("Kd_act", (56, Range())), ("Q10_NCX",\
        (57, Range())), ("V_max_INaCa", (58, Range())), ("eta", (59,\
        Range())), ("ksat", (60, Range())), ("Fx_SLCaP_SL", (61, Range())),\
        ("Fx_SLCaP_jct1", (62, Range())), ("H_ICap", (63, Range())), ("Km",\
        (64, Range())), ("Q10_SLCaP", (65, Range())), ("V_maxAF", (66,\
        Range())), ("Fx_CaBk_SL", (67, Range())), ("Fx_CaBk_jct1", (68,\
        Range())), ("G_CaBk", (69, Range())), ("EC50_SR", (70, Range())),\
        ("HSR", (71, Range())), ("Max_SR", (72, Range())), ("Min_SR", (73,\
        Range())), ("kiCa", (74, Range())), ("kim", (75, Range())), ("koCa",\
        (76, Range())), ("kom", (77, Range())), ("ks", (78, Range())),\
        ("KSRleak", (79, Range())), ("H_Jpump", (80, Range())), ("Kmf", (81,\
        Range())), ("Kmr", (82, Range())), ("Q10_SRCaP", (83, Range())),\
        ("V_max_Jpump", (84, Range())), ("Bmax_Calsequestrin", (85,\
        Range())), ("Bmax_SLB_SL", (86, Range())), ("Bmax_SLB_jct1", (87,\
        Range())), ("Bmax_SLHigh_SL", (88, Range())), ("Bmax_SLHigh_jct1",\
        (89, Range())), ("koff_Calsequestrin", (90, Range())), ("koff_SLB",\
        (91, Range())), ("koff_SLHigh", (92, Range())), ("kon_Calsequestrin",\
        (93, Range())), ("kon_SL", (94, Range())), ("Bmax_Calmodulin", (95,\
        Range())), ("Bmax_Myosin_Ca", (96, Range())), ("Bmax_Myosin_Mg", (97,\
        Range())), ("Bmax_SRB", (98, Range())), ("Bmax_TroponinC", (99,\
        Range())), ("Bmax_TroponinC_Ca_Mg_Ca", (100, Range())),\
        ("Bmax_TroponinC_Ca_Mg_Mg", (101, Range())), ("koff_Calmodulin",\
        (102, Range())), ("koff_Myosin_Ca", (103, Range())),\
        ("koff_Myosin_Mg", (104, Range())), ("koff_SRB", (105, Range())),\
        ("koff_TroponinC", (106, Range())), ("koff_TroponinC_Ca_Mg_Ca", (107,\
        Range())), ("koff_TroponinC_Ca_Mg_Mg", (108, Range())),\
        ("kon_Calmodulin", (109, Range())), ("kon_Myosin_Ca", (110,\
        Range())), ("kon_Myosin_Mg", (111, Range())), ("kon_SRB", (112,\
        Range())), ("kon_TroponinC", (113, Range())),\
        ("kon_TroponinC_Ca_Mg_Ca", (114, Range())),\
        ("kon_TroponinC_Ca_Mg_Mg", (115, Range())), ("Bmax_SL", (116,\
        Range())), ("Bmax_jct1", (117, Range())), ("koff", (118, Range())),\
        ("kon", (119, Range())), ("stim_amplitude", (120, Range())),\
        ("stim_duration", (121, Range())), ("stim_period", (122, Range())),\
        ("stim_start", (123, Range()))])

    for param_name, value in values.items():
        if param_name not in param_ind:
            raise ValueError("{0} is not a parameter.".format(param_name))
        ind, range = param_ind[param_name]
        if value not in range:
            raise ValueError("While setting '{0}' {1}".format(param_name,\
                range.format_not_in(value)))

        # Assign value
        init_values[ind] = value

    return init_values

def state_indices(*states):
    """
    State indices
    """
    state_inds = dict([("h", 0), ("j", 1), ("m", 2), ("Xr", 3), ("Xs", 4),\
        ("X_tos", 5), ("Y_tos", 6), ("R_tos", 7), ("X_tof", 8), ("Y_tof", 9),\
        ("d", 10), ("f", 11), ("fCaB_SL", 12), ("fCaB_jct1", 13), ("R", 14),\
        ("O", 15), ("I", 16), ("Ca_TroponinC", 17), ("Ca_TroponinC_Ca_Mg",\
        18), ("Mg_TroponinC_Ca_Mg", 19), ("Ca_Calmodulin", 20), ("Ca_Myosin",\
        21), ("Mg_Myosin", 22), ("Ca_SRB", 23), ("Na_jct1_buf", 24),\
        ("Na_SL_buf", 25), ("Na_jct1", 26), ("Na_SL", 27), ("Nai", 28),\
        ("Ca_Calsequestrin", 29), ("Ca_SLB_SL", 30), ("Ca_SLB_jct1", 31),\
        ("Ca_SLHigh_SL", 32), ("Ca_SLHigh_jct1", 33), ("Ca_SR", 34),\
        ("Ca_jct1", 35), ("Ca_SL", 36), ("Cai", 37), ("V", 38)])

    indices = []
    for state in states:
        if state not in state_inds:
            raise ValueError("Unknown state: '{0}'".format(state))
        indices.append(state_inds[state])
    if len(indices)>1:
        return indices
    else:
        return indices[0]

def parameter_indices(*params):
    """
    Parameter indices
    """
    param_inds = dict([("Cao", 0), ("Cli", 1), ("Clo", 2), ("Cm", 3), ("F",\
        4), ("Ki", 5), ("Ko", 6), ("Mgi", 7), ("Nao", 8), ("Rgas", 9), ("T",\
        10), ("cell_length", 11), ("cell_radius", 12), ("Fx_Na_SL", 13),\
        ("Fx_Na_jct1", 14), ("G_INa", 15), ("Fx_NaBk_SL", 16),\
        ("Fx_NaBk_jct1", 17), ("G_NaBk", 18), ("Fx_NaK_SL", 19),\
        ("Fx_NaK_jct1", 20), ("H_NaK", 21), ("I_NaK_max", 22), ("Km_Ko", 23),\
        ("Km_Nai", 24), ("Fx_Ks_SL", 25), ("Fx_Ks_jct1", 26), ("pKNa", 27),\
        ("g_Kp", 28), ("G_tos", 29), ("G_tof", 30), ("Fx_Cl_SL", 31),\
        ("Fx_Cl_jct1", 32), ("G_Cl", 33), ("Kd_ClCa", 34), ("G_ClBk", 35),\
        ("Fx_ICaL_SL", 36), ("Fx_ICaL_jct1", 37), ("PCa", 38), ("PK", 39),\
        ("PNa", 40), ("Q10_CaL", 41), ("gamma_Cai", 42), ("gamma_Cao", 43),\
        ("gamma_Ki", 44), ("gamma_Ko", 45), ("gamma_Nai", 46), ("gamma_Nao",\
        47), ("lccCaInact", 48), ("Fx_NCX_SL", 49), ("Fx_NCX_jct1", 50),\
        ("HNa", 51), ("K_mCai", 52), ("K_mCao", 53), ("K_mNai", 54),\
        ("K_mNao", 55), ("Kd_act", 56), ("Q10_NCX", 57), ("V_max_INaCa", 58),\
        ("eta", 59), ("ksat", 60), ("Fx_SLCaP_SL", 61), ("Fx_SLCaP_jct1",\
        62), ("H_ICap", 63), ("Km", 64), ("Q10_SLCaP", 65), ("V_maxAF", 66),\
        ("Fx_CaBk_SL", 67), ("Fx_CaBk_jct1", 68), ("G_CaBk", 69), ("EC50_SR",\
        70), ("HSR", 71), ("Max_SR", 72), ("Min_SR", 73), ("kiCa", 74),\
        ("kim", 75), ("koCa", 76), ("kom", 77), ("ks", 78), ("KSRleak", 79),\
        ("H_Jpump", 80), ("Kmf", 81), ("Kmr", 82), ("Q10_SRCaP", 83),\
        ("V_max_Jpump", 84), ("Bmax_Calsequestrin", 85), ("Bmax_SLB_SL", 86),\
        ("Bmax_SLB_jct1", 87), ("Bmax_SLHigh_SL", 88), ("Bmax_SLHigh_jct1",\
        89), ("koff_Calsequestrin", 90), ("koff_SLB", 91), ("koff_SLHigh",\
        92), ("kon_Calsequestrin", 93), ("kon_SL", 94), ("Bmax_Calmodulin",\
        95), ("Bmax_Myosin_Ca", 96), ("Bmax_Myosin_Mg", 97), ("Bmax_SRB",\
        98), ("Bmax_TroponinC", 99), ("Bmax_TroponinC_Ca_Mg_Ca", 100),\
        ("Bmax_TroponinC_Ca_Mg_Mg", 101), ("koff_Calmodulin", 102),\
        ("koff_Myosin_Ca", 103), ("koff_Myosin_Mg", 104), ("koff_SRB", 105),\
        ("koff_TroponinC", 106), ("koff_TroponinC_Ca_Mg_Ca", 107),\
        ("koff_TroponinC_Ca_Mg_Mg", 108), ("kon_Calmodulin", 109),\
        ("kon_Myosin_Ca", 110), ("kon_Myosin_Mg", 111), ("kon_SRB", 112),\
        ("kon_TroponinC", 113), ("kon_TroponinC_Ca_Mg_Ca", 114),\
        ("kon_TroponinC_Ca_Mg_Mg", 115), ("Bmax_SL", 116), ("Bmax_jct1",\
        117), ("koff", 118), ("kon", 119), ("stim_amplitude", 120),\
        ("stim_duration", 121), ("stim_period", 122), ("stim_start", 123)])

    indices = []
    for param in params:
        if param not in param_inds:
            raise ValueError("Unknown param: '{0}'".format(param))
        indices.append(param_inds[param])
    if len(indices)>1:
        return indices
    else:
        return indices[0]

def monitor_indices(*monitored):
    """
    Monitor indices
    """
    monitor_inds = dict([("Vol_Cell", 0), ("Vol_SR", 1), ("Vol_SL", 2),\
        ("Vol_jct1", 3), ("Vol_myo", 4), ("openProb", 5), ("i_Na_jct1", 6),\
        ("i_Na_SL", 7), ("i_Na", 8), ("alpha_h", 9), ("beta_h", 10),\
        ("alpha_j", 11), ("beta_j", 12), ("alpha_m", 13), ("beta_m", 14),\
        ("i_Nab_jct1", 15), ("i_Nab_SL", 16), ("i_Nab", 17), ("sigma", 18),\
        ("f_NaK", 19), ("i_NaK_jct1", 20), ("i_NaK_SL", 21), ("i_NaK", 22),\
        ("G_IKr", 23), ("i_Kr", 24), ("Xr_infinity", 25), ("tau_Xr", 26),\
        ("Rr", 27), ("pCa_jct1", 28), ("pCa_SL", 29), ("G_Ks_jct1", 30),\
        ("G_Ks_SL", 31), ("E_Ks", 32), ("i_Ks_jct1", 33), ("i_Ks_SL", 34),\
        ("i_Ks", 35), ("Xs_infinity", 36), ("tau_Xs", 37), ("i_Kp", 38),\
        ("i_tos", 39), ("X_tos_infinity", 40), ("tau_X_tos", 41),\
        ("Y_tos_infinity", 42), ("tau_Y_tos", 43), ("R_tos_infinity", 44),\
        ("tau_R_tos", 45), ("i_tof", 46), ("X_tof_infinity", 47),\
        ("tau_X_tof", 48), ("Y_tof_infinity", 49), ("tau_Y_tof", 50),\
        ("i_Cl_Ca", 51), ("i_Clb", 52), ("Q_CaL", 53), ("temp", 54),\
        ("i_CaL_Ca_jct1", 55), ("i_CaL_Na_jct1", 56), ("i_CaL_Ca_SL", 57),\
        ("i_CaL_Na_SL", 58), ("i_CaL_K", 59), ("i_CaL", 60), ("fCa_SL", 61),\
        ("fCa_jct1", 62), ("d_infinity", 63), ("tau_d", 64), ("f_infinity",\
        65), ("tau_f", 66), ("temp_jct1", 67), ("temp_SL", 68), ("Q_NCX",\
        69), ("Ka_SL", 70), ("Ka_jct1", 71), ("i_NaCa_jct1", 72),\
        ("i_NaCa_SL", 73), ("i_NaCa", 74), ("Q_SLCaP", 75), ("i_Cap_jct1",\
        76), ("i_Cap_SL", 77), ("i_Cap", 78), ("i_Cab_jct1", 79),\
        ("i_Cab_SL", 80), ("i_Cab", 81), ("kCaSR", 82), ("koSRCa", 83),\
        ("kiSRCa", 84), ("RI", 85), ("j_rel_SR", 86), ("j_leak_SR", 87),\
        ("Q_SRCaP", 88), ("j_pump_SR", 89), ("dCalsequestrin", 90),\
        ("dCa_SLB_SL", 91), ("dCa_SLB_jct1", 92), ("dCa_SLHigh_SL", 93),\
        ("dCa_SLHigh_jct1", 94), ("dCa_jct1_tot_bound", 95),\
        ("dCa_SL_tot_bound", 96), ("i_Ca_jct1_tot", 97), ("i_Ca_SL_tot", 98),\
        ("dCa_TroponinC", 99), ("dCa_TroponinC_Ca_Mg", 100),\
        ("dMg_TroponinC_Ca_Mg", 101), ("dCa_Calmodulin", 102), ("dCa_Myosin",\
        103), ("dMg_Myosin", 104), ("dCa_SRB", 105),\
        ("dCa_cytosol_tot_bound", 106), ("dNa_jct1_buf", 107), ("dNa_SL_buf",\
        108), ("i_Stim", 109), ("E_Na_jct1", 110), ("E_Na_SL", 111),\
        ("E_Ca_jct1", 112), ("E_Ca_SL", 113), ("E_K", 114), ("E_Cl", 115),\
        ("G_K1", 116), ("i_K1", 117), ("alpha_K1", 118), ("beta_K1", 119),\
        ("K1_infinity", 120), ("J_Na_jct1_SL", 121), ("J_Na_SL_myo", 122),\
        ("J_Ca_jct1_SL", 123), ("J_Ca_SL_myo", 124), ("dh_dt", 125),\
        ("dj_dt", 126), ("dm_dt", 127), ("dXr_dt", 128), ("dXs_dt", 129),\
        ("dX_tos_dt", 130), ("dY_tos_dt", 131), ("dR_tos_dt", 132),\
        ("dX_tof_dt", 133), ("dY_tof_dt", 134), ("dd_dt", 135), ("df_dt",\
        136), ("dfCaB_SL_dt", 137), ("dfCaB_jct1_dt", 138), ("dR_dt", 139),\
        ("dO_dt", 140), ("dI_dt", 141), ("dCa_TroponinC_dt", 142),\
        ("dCa_TroponinC_Ca_Mg_dt", 143), ("dMg_TroponinC_Ca_Mg_dt", 144),\
        ("dCa_Calmodulin_dt", 145), ("dCa_Myosin_dt", 146), ("dMg_Myosin_dt",\
        147), ("dCa_SRB_dt", 148), ("dNa_jct1_buf_dt", 149),\
        ("dNa_SL_buf_dt", 150), ("dNa_jct1_dt", 151), ("dNa_SL_dt", 152),\
        ("dNai_dt", 153), ("dCa_Calsequestrin_dt", 154), ("dCa_SLB_SL_dt",\
        155), ("dCa_SLB_jct1_dt", 156), ("dCa_SLHigh_SL_dt", 157),\
        ("dCa_SLHigh_jct1_dt", 158), ("dCa_SR_dt", 159), ("dCa_jct1_dt",\
        160), ("dCa_SL_dt", 161), ("dCai_dt", 162), ("dV_dt", 163)])

    indices = []
    for monitor in monitored:
        if monitor not in monitor_inds:
            raise ValueError("Unknown monitored: '{0}'".format(monitor))
        indices.append(monitor_inds[monitor])
    if len(indices)>1:
        return indices
    else:
        return indices[0]

def rhs(states, t, parameters, values=None):
    """
    Compute the right hand side of the shannon_2004 ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 39)
    h, j, m, Xr, Xs, X_tos, Y_tos, R_tos, X_tof, Y_tof, d, f, fCaB_SL,\
        fCaB_jct1, R, O, I, Ca_TroponinC, Ca_TroponinC_Ca_Mg,\
        Mg_TroponinC_Ca_Mg, Ca_Calmodulin, Ca_Myosin, Mg_Myosin, Ca_SRB,\
        Na_jct1_buf, Na_SL_buf, Na_jct1, Na_SL, Nai, Ca_Calsequestrin,\
        Ca_SLB_SL, Ca_SLB_jct1, Ca_SLHigh_SL, Ca_SLHigh_jct1, Ca_SR, Ca_jct1,\
        Ca_SL, Cai, V = states

    # Assign parameters
    assert(len(parameters) == 124)
    Cao, Cli, Clo, Cm, F, Ki, Ko, Mgi, Nao, Rgas, T, cell_length,\
        cell_radius, Fx_Na_SL, Fx_Na_jct1, G_INa, Fx_NaBk_SL, Fx_NaBk_jct1,\
        G_NaBk, Fx_NaK_SL, Fx_NaK_jct1, H_NaK, I_NaK_max, Km_Ko, Km_Nai,\
        Fx_Ks_SL, Fx_Ks_jct1, pKNa, g_Kp, G_tos, G_tof, Fx_Cl_SL, Fx_Cl_jct1,\
        G_Cl, Kd_ClCa, G_ClBk, lccCaInact, Fx_NCX_SL, Fx_NCX_jct1, HNa,\
        K_mCai, K_mCao, K_mNai, K_mNao, Kd_act, Q10_NCX, V_max_INaCa, eta,\
        ksat, Fx_SLCaP_SL, Fx_SLCaP_jct1, H_ICap, Km, Q10_SLCaP, V_maxAF,\
        Fx_CaBk_SL, Fx_CaBk_jct1, G_CaBk, EC50_SR, HSR, Max_SR, Min_SR, kiCa,\
        kim, koCa, kom, ks, KSRleak, H_Jpump, Kmf, Kmr, Q10_SRCaP,\
        V_max_Jpump, Bmax_Calmodulin, Bmax_Myosin_Ca, Bmax_Myosin_Mg,\
        Bmax_SRB, Bmax_TroponinC, Bmax_TroponinC_Ca_Mg_Ca,\
        Bmax_TroponinC_Ca_Mg_Mg, koff_Calmodulin, koff_Myosin_Ca,\
        koff_Myosin_Mg, koff_SRB, koff_TroponinC, koff_TroponinC_Ca_Mg_Ca,\
        koff_TroponinC_Ca_Mg_Mg, kon_Calmodulin, kon_Myosin_Ca,\
        kon_Myosin_Mg, kon_SRB, kon_TroponinC, kon_TroponinC_Ca_Mg_Ca,\
        kon_TroponinC_Ca_Mg_Mg, Fx_ICaL_SL, Fx_ICaL_jct1, PCa, PK, PNa,\
        Q10_CaL, gamma_Cai, gamma_Cao, gamma_Ki, gamma_Ko, gamma_Nai,\
        gamma_Nao, Bmax_SL, Bmax_jct1, koff, kon, Bmax_Calsequestrin,\
        Bmax_SLB_SL, Bmax_SLB_jct1, Bmax_SLHigh_SL, Bmax_SLHigh_jct1,\
        koff_Calsequestrin, koff_SLB, koff_SLHigh, kon_Calsequestrin, kon_SL,\
        stim_amplitude, stim_duration, stim_period, stim_start = parameters
    # PKH 
    Cao, Cli, Clo, Cm, F, Ki, Ko, Mgi, Nao, Rgas, T, cell_length, cell_radius, Fx_Na_SL, Fx_Na_jct1, G_INa, Fx_NaBk_SL, Fx_NaBk_jct1, G_NaBk, Fx_NaK_SL, Fx_NaK_jct1, H_NaK, I_NaK_max, Km_Ko, Km_Nai, Fx_Ks_SL, Fx_Ks_jct1, pKNa, g_Kp, G_tos, G_tof, Fx_Cl_SL, Fx_Cl_jct1, G_Cl, Kd_ClCa, G_ClBk, Fx_ICaL_SL, Fx_ICaL_jct1, PCa, PK, PNa, Q10_CaL, gamma_Cai, gamma_Cao, gamma_Ki, gamma_Ko, gamma_Nai, gamma_Nao, lccCaInact, Fx_NCX_SL, Fx_NCX_jct1, HNa, K_mCai, K_mCao, K_mNai, K_mNao, Kd_act, Q10_NCX, V_max_INaCa, eta, ksat, Fx_SLCaP_SL, Fx_SLCaP_jct1, H_ICap, Km, Q10_SLCaP, V_maxAF, Fx_CaBk_SL, Fx_CaBk_jct1, G_CaBk, EC50_SR, HSR, Max_SR, Min_SR, kiCa, kim, koCa, kom, ks, KSRleak, H_Jpump, Kmf, Kmr, Q10_SRCaP, V_max_Jpump, Bmax_Calsequestrin, Bmax_SLB_SL, Bmax_SLB_jct1, Bmax_SLHigh_SL, Bmax_SLHigh_jct1, koff_Calsequestrin, koff_SLB, koff_SLHigh, kon_Calsequestrin, kon_SL, Bmax_Calmodulin, Bmax_Myosin_Ca, Bmax_Myosin_Mg, Bmax_SRB, Bmax_TroponinC, Bmax_TroponinC_Ca_Mg_Ca, Bmax_TroponinC_Ca_Mg_Mg, koff_Calmodulin, koff_Myosin_Ca, koff_Myosin_Mg, koff_SRB, koff_TroponinC, koff_TroponinC_Ca_Mg_Ca, koff_TroponinC_Ca_Mg_Mg, kon_Calmodulin, kon_Myosin_Ca, kon_Myosin_Mg, kon_SRB, kon_TroponinC, kon_TroponinC_Ca_Mg_Ca, kon_TroponinC_Ca_Mg_Mg, Bmax_SL, Bmax_jct1, koff, kon, stim_amplitude, stim_duration, stim_period, stim_start = parameters

    # Init return args
    if values is None:
        values = np.zeros((39,), dtype=np.float_)
    else:
        assert isinstance(values, np.ndarray) and values.shape == (39,)

    # Expressions for the Model parameters component
    Vol_Cell = 3.141592654e-15*cell_length*(cell_radius*cell_radius)
    Vol_SR = 0.035*Vol_Cell
    Vol_SL = 0.02*Vol_Cell
    Vol_jct1 = 0.000539*Vol_Cell
    Vol_myo = 0.65*Vol_Cell

    # Expressions for the Reversal potentials component
    E_Na_jct1 = Rgas*T*math.log(Nao/Na_jct1)/F
    E_Na_SL = Rgas*T*math.log(Nao/Na_SL)/F
    E_Ca_jct1 = Rgas*T*math.log(Cao/Ca_jct1)/(2*F)
    E_Ca_SL = Rgas*T*math.log(Cao/Ca_SL)/(2*F)
    E_K = Rgas*T*math.log(Ko/Ki)/F
    E_Cl = Rgas*T*math.log(Cli/Clo)/F

    # Expressions for the INa component
    openProb = (m*m*m)*h*j
    i_Na_jct1 = Fx_Na_jct1*G_INa*(V - E_Na_jct1)*openProb
    i_Na_SL = Fx_Na_SL*G_INa*(V - E_Na_SL)*openProb
    i_Na = i_Na_jct1 + i_Na_SL

    # Expressions for the h gate component
    alpha_h = (1.04951082543e-06*math.exp(-0.147058823529*V) if V < -40 else 0)
    beta_h = (3.56*math.exp(0.079*V) + 310000.0*math.exp(0.35*V) if V < -40 else\
        1.0/(0.13 + 0.0497581410839*math.exp(-0.0900900900901*V)))
    values[0] = -beta_h*h + (1 - h)*alpha_h

    # Expressions for the j gate component
    alpha_j = ((37.78 + V)*(-127140.0*math.exp(0.2444*V) -\
        3.474e-05*math.exp(-0.04391*V))/(1 + 50262745826.0*math.exp(0.311*V))\
        if V < -40 else 0)
    beta_j = (0.1212*math.exp(-0.01052*V)/(1 +\
        0.0039608683399*math.exp(-0.1378*V)) if V < -40 else\
        0.3*math.exp(-2.535e-07*V)/(1 + 0.0407622039784*math.exp(-0.1*V)))
    values[1] = -beta_j*j + (1 - j)*alpha_j

    # Expressions for the m gate component
    alpha_m = (15.0816 + 0.32*V)/(1 - 0.0089778037307*math.exp(-0.1*V))
    beta_m = 0.08*math.exp(-V/11.)
    values[2] = (1 - m)*alpha_m - beta_m*m

    # Expressions for the INab component
    i_Nab_jct1 = Fx_NaBk_jct1*G_NaBk*(V - E_Na_jct1)
    i_Nab_SL = Fx_NaBk_SL*G_NaBk*(V - E_Na_SL)
    i_Nab = i_Nab_jct1 + i_Nab_SL

    # Expressions for the INaK component
    sigma = -1/7. + math.exp(0.0148588410104*Nao)/7.
    f_NaK = 1.0/(1 + 0.0365*math.exp(-F*V/(Rgas*T))*sigma +\
        0.1245*math.exp(-0.1*F*V/(Rgas*T)))
    i_NaK_jct1 = Fx_NaK_jct1*I_NaK_max*Ko*f_NaK/((1 +\
        math.pow(Km_Nai/Na_jct1, H_NaK))*(Km_Ko + Ko))
    i_NaK_SL = Fx_NaK_SL*I_NaK_max*Ko*f_NaK/((1 + math.pow(Km_Nai/Na_SL,\
        H_NaK))*(Km_Ko + Ko))
    i_NaK = i_NaK_SL + i_NaK_jct1

    # Expressions for the Xr gate component
    Xr_infinity = 1.0/(1 + 0.00127263380134*math.exp(-0.133333333333*V))
    tau_Xr = 1.0/((0.0061 + 0.00061*V)/(-1 + 4.26311451517*math.exp(0.145*V))\
        + (0.00966 + 0.00138*V)/(1 - 0.422739131746*math.exp(-0.123*V)))
    values[3] = (Xr_infinity - Xr)/tau_Xr

    # Expressions for the Rr gate component
    Rr = 1.0/(1 + 4.36323731689*math.exp(0.0446428571429*V))

    # Expressions for the IKs component
    pCa_jct1 = 3 - math.log(Ca_jct1)
    pCa_SL = 3 - math.log(Ca_SL)
    G_Ks_jct1 = 0.00399 + 0.0133/(1 +\
        6.14421235333e-06*math.exp(1.66666666667*pCa_jct1))
    G_Ks_SL = 0.00399 + 0.0133/(1 +\
        6.14421235333e-06*math.exp(1.66666666667*pCa_SL))
    E_Ks = Rgas*T*math.log((Nao*pKNa + Ko)/(pKNa*Nai + Ki))/F
    i_Ks_jct1 = Fx_Ks_jct1*(Xs*Xs)*(-E_Ks + V)*G_Ks_jct1
    i_Ks_SL = Fx_Ks_SL*(Xs*Xs)*(-E_Ks + V)*G_Ks_SL
    i_Ks = i_Ks_jct1 + i_Ks_SL

    # Expressions for the Xs gate component
    Xs_infinity = 1.0/(1 + 1.0939777431*math.exp(-0.059880239521*V))
    tau_Xs = 1.0/((0.002157 + 7.19e-05*V)/(1 -\
        0.0117959385198*math.exp(-0.148*V)) + (0.00393 + 0.000131*V)/(-1 +\
        7.85381970442*math.exp(0.0687*V)))
    values[4] = (-Xs + Xs_infinity)/tau_Xs

    # Expressions for the IKp component
    i_Kp = g_Kp*(V - E_K)/(1 + 1786.47556538*math.exp(-0.167224080268*V))

    # Expressions for the Itos component
    i_tos = G_tos*(0.5*R_tos + Y_tos)*(V - E_K)*X_tos

    # Expressions for the X_gate component
    X_tos_infinity = 1.0/(1 + math.exp(-1/5. - V/15.))
    tau_X_tos = 0.5 + 9/(1 + math.exp(1/5. + V/15.))
    values[5] = (-X_tos + X_tos_infinity)/tau_X_tos

    # Expressions for the Y_gate component
    Y_tos_infinity = 1.0/(1 + 28.5027336438*math.exp(V/10.))
    tau_Y_tos = 30 + 3000/(1 + math.exp(6 + V/10.))
    values[6] = (Y_tos_infinity - Y_tos)/tau_Y_tos

    # Expressions for the R_gate component
    R_tos_infinity = 1.0/(1 + 28.5027336438*math.exp(V/10.))
    tau_R_tos = 220 + 2800.0/(1 + math.exp(6 + V/10.))
    values[7] = (-R_tos + R_tos_infinity)/tau_R_tos

    # Expressions for the Itof component
    i_tof = G_tof*(V - E_K)*X_tof*Y_tof

    # Expressions for the Itof X gate component
    X_tof_infinity = 1.0/(1 + math.exp(-1/5. - V/15.))
    tau_X_tof = 1.5 + 3.5*math.exp(-(V*V)/900.)
    values[8] = (-X_tof + X_tof_infinity)/tau_X_tof

    # Expressions for the Itof Y gate component
    Y_tof_infinity = 1.0/(1 + 28.5027336438*math.exp(V/10.))
    tau_Y_tof = 20 + 20/(1 + 28.5027336438*math.exp(V/10.))
    values[9] = (Y_tof_infinity - Y_tof)/tau_Y_tof

    # Expressions for the K1 gate component
    alpha_K1 = 1.02/(1 + 7.35454251046e-07*math.exp(-0.2385*E_K + 0.2385*V))
    beta_K1 = (0.762624006506*math.exp(-0.08032*E_K + 0.08032*V) +\
        1.15340563519e-16*math.exp(-0.06175*E_K + 0.06175*V))/(1 +\
        0.0867722941577*math.exp(0.5143*E_K - 0.5143*V))
    K1_infinity = alpha_K1/(beta_K1 + alpha_K1)

    # Expressions for the ICl Ca component
    i_Cl_Ca = G_Cl*(Fx_Cl_SL/(1 + Kd_ClCa/Ca_SL) + Fx_Cl_jct1/(1 +\
        Kd_ClCa/Ca_jct1))*(V - E_Cl)

    # Expressions for the IClb component
    i_Clb = G_ClBk*(V - E_Cl)

    # Expressions for the d gate component
    d_infinity = 1.0/(1 + 0.0892185174093*math.exp(-V/6.))
    tau_d = (1 - 0.0892185174093*math.exp(-V/6.))*d_infinity/(0.5075 + 0.035*V)
    values[10] = (d_infinity - d)/tau_d

    # Expressions for the f gate component
    f_infinity = 0.6/(1 + math.exp(5/2. - V/20.)) + 1.0/(1 +\
        16964.681259*math.exp(0.277777777778*V))
    tau_f = 1.0/(0.02 + 0.0197*math.exp(-((0.48865 + 0.0337*V)*(0.48865 +\
        0.0337*V))))
    values[11] = (f_infinity - f)/tau_f

    # Expressions for the FCa gate component
    fCa_SL = 1 - fCaB_SL
    fCa_jct1 = 1 - fCaB_jct1
    values[12] = -0.0119*fCaB_SL + lccCaInact*(1 - fCaB_SL)*Ca_SL
    values[13] = lccCaInact*(1 - fCaB_jct1)*Ca_jct1 - 0.0119*fCaB_jct1

    # Expressions for the INaCa component
    temp_jct1 = (-math.pow(Nao, HNa)*Ca_jct1*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)) + Cao*math.pow(Na_jct1,\
        HNa)*math.exp(F*eta*V/(Rgas*T)))/(1 + ksat*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)))
    temp_SL = (-math.pow(Nao, HNa)*Ca_SL*math.exp(F*(-1 + eta)*V/(Rgas*T)) +\
        Cao*math.pow(Na_SL, HNa)*math.exp(F*eta*V/(Rgas*T)))/(1 +\
        ksat*math.exp(F*(-1 + eta)*V/(Rgas*T)))
    Q_NCX = math.pow(Q10_NCX, -31 + T/10.)
    Ka_SL = 1.0/(1 + (Kd_act*Kd_act*Kd_act)/(Ca_SL*Ca_SL*Ca_SL))
    Ka_jct1 = 1.0/(1 + (Kd_act*Kd_act*Kd_act)/(Ca_jct1*Ca_jct1*Ca_jct1))
    i_NaCa_jct1 =\
        Fx_NCX_jct1*V_max_INaCa*Ka_jct1*Q_NCX*temp_jct1/(K_mCai*math.pow(Nao,\
        HNa)*(1 + math.pow(Na_jct1/K_mNai, HNa)) + K_mCao*math.pow(Na_jct1,\
        HNa) + Cao*math.pow(Na_jct1, HNa) + math.pow(Nao, HNa)*Ca_jct1 +\
        math.pow(K_mNao, HNa)*(1 + Ca_jct1/K_mCai)*Ca_jct1)
    i_NaCa_SL = Fx_NCX_SL*V_max_INaCa*Ka_SL*Q_NCX*temp_SL/(math.pow(K_mNao,\
        HNa)*(1 + Ca_SL/K_mCai)*Ca_SL + math.pow(Nao, HNa)*Ca_SL +\
        Cao*math.pow(Na_SL, HNa) + K_mCao*math.pow(Na_SL, HNa) +\
        K_mCai*math.pow(Nao, HNa)*(1 + math.pow(Na_SL/K_mNai, HNa)))
    i_NaCa = i_NaCa_SL + i_NaCa_jct1

    # Expressions for the ICap component
    Q_SLCaP = math.pow(Q10_SLCaP, -31 + T/10.)
    i_Cap_jct1 = Fx_SLCaP_jct1*V_maxAF*Q_SLCaP/(1 + math.pow(Km/Ca_jct1,\
        H_ICap))
    i_Cap_SL = Fx_SLCaP_SL*V_maxAF*Q_SLCaP/(1 + math.pow(Km/Ca_SL, H_ICap))
    i_Cap = i_Cap_SL + i_Cap_jct1

    # Expressions for the ICab component
    i_Cab_jct1 = Fx_CaBk_jct1*G_CaBk*(-E_Ca_jct1 + V)
    i_Cab_SL = Fx_CaBk_SL*G_CaBk*(V - E_Ca_SL)
    i_Cab = i_Cab_jct1 + i_Cab_SL

    # Expressions for the Jrel SR component
    kCaSR = Max_SR - (Max_SR - Min_SR)/(1 + math.pow(EC50_SR/Ca_SR, HSR))
    koSRCa = koCa/kCaSR
    kiSRCa = kiCa*kCaSR
    RI = 1 - I - O - R
    values[14] = -(Ca_jct1*Ca_jct1)*R*koSRCa + kom*O - Ca_jct1*R*kiSRCa +\
        kim*RI
    values[15] = -Ca_jct1*O*kiSRCa + kim*I + (Ca_jct1*Ca_jct1)*R*koSRCa - kom*O
    values[16] = -kim*I + Ca_jct1*O*kiSRCa - kom*I +\
        (Ca_jct1*Ca_jct1)*RI*koSRCa
    j_rel_SR = ks*(Ca_SR - Ca_jct1)*O
    #print "%f %f %f" % (j_rel_SR,(Ca_SR - Ca_jct1),O)

    # Expressions for the Jleak SR component
    j_leak_SR = KSRleak*(Ca_SR - Ca_jct1)

    # Expressions for the Jpump SR component
    Q_SRCaP = math.pow(Q10_SRCaP, -31 + T/10.)
    j_pump_SR = V_max_Jpump*(-math.pow(Ca_SR/Kmr, H_Jpump) +\
        math.pow(Cai/Kmf, H_Jpump))*Q_SRCaP/(1 + math.pow(Ca_SR/Kmr, H_Jpump)\
        + math.pow(Cai/Kmf, H_Jpump))

    #print Cai,Ca_SR,T,j_pump_SR , V_max_Jpump,(-math.pow(Ca_SR/Kmr, H_Jpump) ,   
    #    math.pow(Cai/Kmf, H_Jpump)),Q_SRCaP/(1 + math.pow(Ca_SR/Kmr, H_Jpump)\
    #    + math.pow(Cai/Kmf, H_Jpump))

    # Expressions for the Ion diffusion component
    J_Na_jct1_SL = -1.8313e-14*Na_SL + 1.8313e-14*Na_jct1
    J_Na_SL_myo = 1.6386e-12*Na_SL - 1.6386e-12*Nai
    J_Ca_jct1_SL = 8.2413e-13*Ca_jct1 - 8.2413e-13*Ca_SL
    J_Ca_SL_myo = -3.7243e-12*Cai + 3.7243e-12*Ca_SL

    # Expressions for the Cytosolic component
    dCa_TroponinC = kon_TroponinC*(-Ca_TroponinC + Bmax_TroponinC)*Cai -\
        koff_TroponinC*Ca_TroponinC
    dCa_TroponinC_Ca_Mg = -koff_TroponinC_Ca_Mg_Ca*Ca_TroponinC_Ca_Mg +\
        kon_TroponinC_Ca_Mg_Ca*(-Ca_TroponinC_Ca_Mg - Mg_TroponinC_Ca_Mg +\
        Bmax_TroponinC_Ca_Mg_Ca)*Cai
    dMg_TroponinC_Ca_Mg = Mgi*kon_TroponinC_Ca_Mg_Mg*(-Ca_TroponinC_Ca_Mg +\
        Bmax_TroponinC_Ca_Mg_Mg - Mg_TroponinC_Ca_Mg) -\
        koff_TroponinC_Ca_Mg_Mg*Mg_TroponinC_Ca_Mg
    dCa_Calmodulin = -koff_Calmodulin*Ca_Calmodulin +\
        kon_Calmodulin*(Bmax_Calmodulin - Ca_Calmodulin)*Cai
    dCa_Myosin = kon_Myosin_Ca*(Bmax_Myosin_Ca - Mg_Myosin - Ca_Myosin)*Cai -\
        koff_Myosin_Ca*Ca_Myosin
    dMg_Myosin = Mgi*kon_Myosin_Mg*(Bmax_Myosin_Mg - Mg_Myosin - Ca_Myosin) -\
        koff_Myosin_Mg*Mg_Myosin
    dCa_SRB = kon_SRB*(-Ca_SRB + Bmax_SRB)*Cai - koff_SRB*Ca_SRB
    dCa_cytosol_tot_bound = dCa_SRB + dCa_Calmodulin + dMg_Myosin +\
        dMg_TroponinC_Ca_Mg + dCa_TroponinC + dCa_TroponinC_Ca_Mg +\
        dCa_Myosin
    values[17] = dCa_TroponinC
    values[18] = dCa_TroponinC_Ca_Mg
    values[19] = dMg_TroponinC_Ca_Mg
    values[20] = dCa_Calmodulin
    values[21] = dCa_Myosin
    values[22] = dMg_Myosin
    values[23] = dCa_SRB

    # Expressions for the IKr component
    G_IKr = 0.0129099444874*math.sqrt(Ko)
    i_Kr = (V - E_K)*G_IKr*Rr*Xr

    # Expressions for the IK1 component
    G_K1 = 0.387298334621*math.sqrt(Ko)
    i_K1 = (V - E_K)*G_K1*K1_infinity

    # Expressions for the ICaL component
    Q_CaL = math.pow(Q10_CaL, -31 + T/10.)
    temp = 0.45*(F*F)*Q_CaL*V*d*f/(Rgas*T)
    i_CaL_Ca_jct1 =\
        4*Fx_ICaL_jct1*PCa*(gamma_Cai*Ca_jct1*math.exp(2*F*V/(Rgas*T)) -\
        Cao*gamma_Cao)*fCa_jct1*temp/(-1 + math.exp(2*F*V/(Rgas*T)))
    i_CaL_Na_jct1 = Fx_ICaL_jct1*PNa*(-Nao*gamma_Nao +\
        gamma_Nai*Na_jct1*math.exp(F*V/(Rgas*T)))*fCa_jct1*temp/(-1 +\
        math.exp(F*V/(Rgas*T)))
    i_CaL_Ca_SL = 4*Fx_ICaL_SL*PCa*(gamma_Cai*Ca_SL*math.exp(2*F*V/(Rgas*T))\
        - Cao*gamma_Cao)*fCa_SL*temp/(-1 + math.exp(2*F*V/(Rgas*T)))
    i_CaL_Na_SL = Fx_ICaL_SL*PNa*(-Nao*gamma_Nao +\
        gamma_Nai*Na_SL*math.exp(F*V/(Rgas*T)))*fCa_SL*temp/(-1 +\
        math.exp(F*V/(Rgas*T)))
    i_CaL_K = PK*(Fx_ICaL_SL*fCa_SL + Fx_ICaL_jct1*fCa_jct1)*(-Ko*gamma_Ko +\
        Ki*gamma_Ki*math.exp(F*V/(Rgas*T)))*temp/(-1 +\
        math.exp(F*V/(Rgas*T)))
    i_CaL = i_CaL_Ca_jct1 + i_CaL_Ca_SL + i_CaL_Na_jct1 + i_CaL_K + i_CaL_Na_SL

    # Expressions for the Na buffer component
    dNa_jct1_buf = kon*(Bmax_jct1 - Na_jct1_buf)*Na_jct1 - koff*Na_jct1_buf
    dNa_SL_buf = -koff*Na_SL_buf + kon*(-Na_SL_buf + Bmax_SL)*Na_SL
    values[24] = dNa_jct1_buf
    values[25] = dNa_SL_buf
    values[26] = -dNa_jct1_buf - Cm*(3*i_NaCa_jct1 + 3*i_NaK_jct1 + i_Na_jct1 +\
        i_CaL_Na_jct1 + i_Nab_jct1)/(F*Vol_jct1) - J_Na_jct1_SL/Vol_jct1
    values[27] = -dNa_SL_buf + (-J_Na_SL_myo + J_Na_jct1_SL)/Vol_SL -\
        Cm*(3*i_NaCa_SL + i_Na_SL + i_CaL_Na_SL + i_Nab_SL +\
        3*i_NaK_SL)/(F*Vol_SL)
    values[28] = J_Na_SL_myo/Vol_myo

    # Expressions for the Ca buffer component
    dCalsequestrin = -koff_Calsequestrin*Ca_Calsequestrin +\
        kon_Calsequestrin*(Bmax_Calsequestrin*Vol_myo/Vol_SR -\
        Ca_Calsequestrin)*Ca_SR
    values[29] = dCalsequestrin
    dCa_SLB_SL = kon_SL*(-Ca_SLB_SL + Bmax_SLB_SL*Vol_myo/Vol_SL)*Ca_SL -\
        koff_SLB*Ca_SLB_SL
    dCa_SLB_jct1 = -koff_SLB*Ca_SLB_jct1 +\
        kon_SL*(0.1*Bmax_SLB_jct1*Vol_myo/Vol_jct1 - Ca_SLB_jct1)*Ca_jct1
    dCa_SLHigh_SL = -koff_SLHigh*Ca_SLHigh_SL +\
        kon_SL*(Bmax_SLHigh_SL*Vol_myo/Vol_SL - Ca_SLHigh_SL)*Ca_SL
    dCa_SLHigh_jct1 = -koff_SLHigh*Ca_SLHigh_jct1 +\
        kon_SL*(0.1*Bmax_SLHigh_jct1*Vol_myo/Vol_jct1 -\
        Ca_SLHigh_jct1)*Ca_jct1
    values[30] = dCa_SLB_SL
    values[31] = dCa_SLB_jct1
    values[32] = dCa_SLHigh_SL
    values[33] = dCa_SLHigh_jct1
    dCa_jct1_tot_bound = dCa_SLHigh_jct1 + dCa_SLB_jct1
    dCa_SL_tot_bound = dCa_SLB_SL + dCa_SLHigh_SL
    i_Ca_jct1_tot = i_CaL_Ca_jct1 - 2*i_NaCa_jct1 + i_Cap_jct1 + i_Cab_jct1
    i_Ca_SL_tot = i_Cap_SL + i_CaL_Ca_SL - 2*i_NaCa_SL + i_Cab_SL
    values[34] = -dCalsequestrin - Vol_myo*j_leak_SR/Vol_SR - j_rel_SR +\
        j_pump_SR
    values[35] = Vol_myo*j_leak_SR/Vol_jct1 + Vol_SR*j_rel_SR/Vol_jct1 -\
        J_Ca_jct1_SL/Vol_jct1 - Cm*i_Ca_jct1_tot/(2*F*Vol_jct1) -\
        dCa_jct1_tot_bound
    values[36] = (J_Ca_jct1_SL - J_Ca_SL_myo)/Vol_SL - dCa_SL_tot_bound -\
        Cm*i_Ca_SL_tot/(2*F*Vol_SL)
    values[37] = -dCa_cytosol_tot_bound + J_Ca_SL_myo/Vol_myo -\
        Vol_SR*j_pump_SR/Vol_myo

    # Expressions for the Cell component
    i_Stim = (-stim_amplitude if -stim_period*math.floor(t/stim_period) + t\
        <= stim_start + stim_duration and\
        -stim_period*math.floor(t/stim_period) + t >= stim_start else 0)
    values[38] = -i_CaL - i_Clb - i_tos - i_tof - i_Stim - i_Nab - i_Ks -\
        i_Cab - i_NaCa - i_Na - i_Cap - i_Cl_Ca - i_NaK - i_Kr - i_K1 - i_Kp

    # Return results
    return values

def monitor(states, t, parameters, monitored=None):
    """
    Computes monitored expressions of the shannon_2004 ODE
    """
    # Imports
    import numpy as np
    import math

    # Assign states
    assert(len(states) == 39)
    h, j, m, Xr, Xs, X_tos, Y_tos, R_tos, X_tof, Y_tof, d, f, fCaB_SL,\
        fCaB_jct1, R, O, I, Ca_TroponinC, Ca_TroponinC_Ca_Mg,\
        Mg_TroponinC_Ca_Mg, Ca_Calmodulin, Ca_Myosin, Mg_Myosin, Ca_SRB,\
        Na_jct1_buf, Na_SL_buf, Na_jct1, Na_SL, Nai, Ca_Calsequestrin,\
        Ca_SLB_SL, Ca_SLB_jct1, Ca_SLHigh_SL, Ca_SLHigh_jct1, Ca_SR, Ca_jct1,\
        Ca_SL, Cai, V = states

    # Assign parameters
    assert(len(parameters) == 124)
    Cao, Cli, Clo, Cm, F, Ki, Ko, Mgi, Nao, Rgas, T, cell_length,\
        cell_radius, Fx_Na_SL, Fx_Na_jct1, G_INa, Fx_NaBk_SL, Fx_NaBk_jct1,\
        G_NaBk, Fx_NaK_SL, Fx_NaK_jct1, H_NaK, I_NaK_max, Km_Ko, Km_Nai,\
        Fx_Ks_SL, Fx_Ks_jct1, pKNa, g_Kp, G_tos, G_tof, Fx_Cl_SL, Fx_Cl_jct1,\
        G_Cl, Kd_ClCa, G_ClBk, lccCaInact, Fx_NCX_SL, Fx_NCX_jct1, HNa,\
        K_mCai, K_mCao, K_mNai, K_mNao, Kd_act, Q10_NCX, V_max_INaCa, eta,\
        ksat, Fx_SLCaP_SL, Fx_SLCaP_jct1, H_ICap, Km, Q10_SLCaP, V_maxAF,\
        Fx_CaBk_SL, Fx_CaBk_jct1, G_CaBk, EC50_SR, HSR, Max_SR, Min_SR, kiCa,\
        kim, koCa, kom, ks, KSRleak, H_Jpump, Kmf, Kmr, Q10_SRCaP,\
        V_max_Jpump, Bmax_Calmodulin, Bmax_Myosin_Ca, Bmax_Myosin_Mg,\
        Bmax_SRB, Bmax_TroponinC, Bmax_TroponinC_Ca_Mg_Ca,\
        Bmax_TroponinC_Ca_Mg_Mg, koff_Calmodulin, koff_Myosin_Ca,\
        koff_Myosin_Mg, koff_SRB, koff_TroponinC, koff_TroponinC_Ca_Mg_Ca,\
        koff_TroponinC_Ca_Mg_Mg, kon_Calmodulin, kon_Myosin_Ca,\
        kon_Myosin_Mg, kon_SRB, kon_TroponinC, kon_TroponinC_Ca_Mg_Ca,\
        kon_TroponinC_Ca_Mg_Mg, Fx_ICaL_SL, Fx_ICaL_jct1, PCa, PK, PNa,\
        Q10_CaL, gamma_Cai, gamma_Cao, gamma_Ki, gamma_Ko, gamma_Nai,\
        gamma_Nao, Bmax_SL, Bmax_jct1, koff, kon, Bmax_Calsequestrin,\
        Bmax_SLB_SL, Bmax_SLB_jct1, Bmax_SLHigh_SL, Bmax_SLHigh_jct1,\
        koff_Calsequestrin, koff_SLB, koff_SLHigh, kon_Calsequestrin, kon_SL,\
        stim_amplitude, stim_duration, stim_period, stim_start = parameters
    # PKH 
    Cao, Cli, Clo, Cm, F, Ki, Ko, Mgi, Nao, Rgas, T, cell_length, cell_radius, Fx_Na_SL, Fx_Na_jct1, G_INa, Fx_NaBk_SL, Fx_NaBk_jct1, G_NaBk, Fx_NaK_SL, Fx_NaK_jct1, H_NaK, I_NaK_max, Km_Ko, Km_Nai, Fx_Ks_SL, Fx_Ks_jct1, pKNa, g_Kp, G_tos, G_tof, Fx_Cl_SL, Fx_Cl_jct1, G_Cl, Kd_ClCa, G_ClBk, Fx_ICaL_SL, Fx_ICaL_jct1, PCa, PK, PNa, Q10_CaL, gamma_Cai, gamma_Cao, gamma_Ki, gamma_Ko, gamma_Nai, gamma_Nao, lccCaInact, Fx_NCX_SL, Fx_NCX_jct1, HNa, K_mCai, K_mCao, K_mNai, K_mNao, Kd_act, Q10_NCX, V_max_INaCa, eta, ksat, Fx_SLCaP_SL, Fx_SLCaP_jct1, H_ICap, Km, Q10_SLCaP, V_maxAF, Fx_CaBk_SL, Fx_CaBk_jct1, G_CaBk, EC50_SR, HSR, Max_SR, Min_SR, kiCa, kim, koCa, kom, ks, KSRleak, H_Jpump, Kmf, Kmr, Q10_SRCaP, V_max_Jpump, Bmax_Calsequestrin, Bmax_SLB_SL, Bmax_SLB_jct1, Bmax_SLHigh_SL, Bmax_SLHigh_jct1, koff_Calsequestrin, koff_SLB, koff_SLHigh, kon_Calsequestrin, kon_SL, Bmax_Calmodulin, Bmax_Myosin_Ca, Bmax_Myosin_Mg, Bmax_SRB, Bmax_TroponinC, Bmax_TroponinC_Ca_Mg_Ca, Bmax_TroponinC_Ca_Mg_Mg, koff_Calmodulin, koff_Myosin_Ca, koff_Myosin_Mg, koff_SRB, koff_TroponinC, koff_TroponinC_Ca_Mg_Ca, koff_TroponinC_Ca_Mg_Mg, kon_Calmodulin, kon_Myosin_Ca, kon_Myosin_Mg, kon_SRB, kon_TroponinC, kon_TroponinC_Ca_Mg_Ca, kon_TroponinC_Ca_Mg_Mg, Bmax_SL, Bmax_jct1, koff, kon, stim_amplitude, stim_duration, stim_period, stim_start = parameters

    # Init return args
    if monitored is None:
        monitored = np.zeros((164,), dtype=np.float_)
    else:
        assert isinstance(monitored, np.ndarray) and monitored.shape == (164,)

    # Expressions for the Model parameters component
    monitored[0] = 3.141592654e-15*cell_length*(cell_radius*cell_radius)
    monitored[1] = 0.035*monitored[0]
    monitored[2] = 0.02*monitored[0]
    monitored[3] = 0.000539*monitored[0]
    monitored[4] = 0.65*monitored[0]

    # Expressions for the Reversal potentials component
    monitored[110] = Rgas*T*math.log(Nao/Na_jct1)/F
    monitored[111] = Rgas*T*math.log(Nao/Na_SL)/F
    monitored[112] = Rgas*T*math.log(Cao/Ca_jct1)/(2*F)
    monitored[113] = Rgas*T*math.log(Cao/Ca_SL)/(2*F)
    monitored[114] = Rgas*T*math.log(Ko/Ki)/F
    monitored[115] = Rgas*T*math.log(Cli/Clo)/F

    # Expressions for the INa component
    monitored[5] = (m*m*m)*h*j
    monitored[6] = Fx_Na_jct1*G_INa*(V - monitored[110])*monitored[5]
    monitored[7] = Fx_Na_SL*G_INa*(-monitored[111] + V)*monitored[5]
    monitored[8] = monitored[6] + monitored[7]

    # Expressions for the h gate component
    monitored[9] = (1.04951082543e-06*math.exp(-0.147058823529*V) if V < -40 else\
        0)
    monitored[10] = (3.56*math.exp(0.079*V) + 310000.0*math.exp(0.35*V) if V\
        < -40 else 1.0/(0.13 + 0.0497581410839*math.exp(-0.0900900900901*V)))
    monitored[125] = -h*monitored[10] + (1 - h)*monitored[9]

    # Expressions for the j gate component
    monitored[11] = ((37.78 + V)*(-127140.0*math.exp(0.2444*V) -\
        3.474e-05*math.exp(-0.04391*V))/(1 + 50262745826.0*math.exp(0.311*V))\
        if V < -40 else 0)
    monitored[12] = (0.1212*math.exp(-0.01052*V)/(1 +\
        0.0039608683399*math.exp(-0.1378*V)) if V < -40 else\
        0.3*math.exp(-2.535e-07*V)/(1 + 0.0407622039784*math.exp(-0.1*V)))
    monitored[126] = (1 - j)*monitored[11] - j*monitored[12]

    # Expressions for the m gate component
    monitored[13] = (15.0816 + 0.32*V)/(1 - 0.0089778037307*math.exp(-0.1*V))
    monitored[14] = 0.08*math.exp(-V/11.)
    monitored[127] = (1 - m)*monitored[13] - m*monitored[14]

    # Expressions for the INab component
    monitored[15] = Fx_NaBk_jct1*G_NaBk*(V - monitored[110])
    monitored[16] = Fx_NaBk_SL*G_NaBk*(-monitored[111] + V)
    monitored[17] = monitored[16] + monitored[15]

    # Expressions for the INaK component
    monitored[18] = -1/7. + math.exp(0.0148588410104*Nao)/7.
    monitored[19] = 1.0/(1 + 0.0365*math.exp(-F*V/(Rgas*T))*monitored[18] +\
        0.1245*math.exp(-0.1*F*V/(Rgas*T)))
    monitored[20] = Fx_NaK_jct1*I_NaK_max*Ko*monitored[19]/((1 +\
        math.pow(Km_Nai/Na_jct1, H_NaK))*(Km_Ko + Ko))
    monitored[21] = Fx_NaK_SL*I_NaK_max*Ko*monitored[19]/((1 +\
        math.pow(Km_Nai/Na_SL, H_NaK))*(Km_Ko + Ko))
    monitored[22] = monitored[20] + monitored[21]

    # Expressions for the Xr gate component
    monitored[25] = 1.0/(1 + 0.00127263380134*math.exp(-0.133333333333*V))
    monitored[26] = 1.0/((0.0061 + 0.00061*V)/(-1 +\
        4.26311451517*math.exp(0.145*V)) + (0.00966 + 0.00138*V)/(1 -\
        0.422739131746*math.exp(-0.123*V)))
    monitored[128] = (monitored[25] - Xr)/monitored[26]

    # Expressions for the Rr gate component
    monitored[27] = 1.0/(1 + 4.36323731689*math.exp(0.0446428571429*V))

    # Expressions for the IKs component
    monitored[28] = 3 - math.log(Ca_jct1)
    monitored[29] = 3 - math.log(Ca_SL)
    monitored[30] = 0.00399 + 0.0133/(1 +\
        6.14421235333e-06*math.exp(1.66666666667*monitored[28]))
    monitored[31] = 0.00399 + 0.0133/(1 +\
        6.14421235333e-06*math.exp(1.66666666667*monitored[29]))
    monitored[32] = Rgas*T*math.log((Nao*pKNa + Ko)/(pKNa*Nai + Ki))/F
    monitored[33] = Fx_Ks_jct1*(Xs*Xs)*(V - monitored[32])*monitored[30]
    monitored[34] = Fx_Ks_SL*(Xs*Xs)*(V - monitored[32])*monitored[31]
    monitored[35] = monitored[33] + monitored[34]

    # Expressions for the Xs gate component
    monitored[36] = 1.0/(1 + 1.0939777431*math.exp(-0.059880239521*V))
    monitored[37] = 1.0/((0.002157 + 7.19e-05*V)/(1 -\
        0.0117959385198*math.exp(-0.148*V)) + (0.00393 + 0.000131*V)/(-1 +\
        7.85381970442*math.exp(0.0687*V)))
    monitored[129] = (monitored[36] - Xs)/monitored[37]

    # Expressions for the IKp component
    monitored[38] = g_Kp*(V - monitored[114])/(1 +\
        1786.47556538*math.exp(-0.167224080268*V))

    # Expressions for the Itos component
    monitored[39] = G_tos*(0.5*R_tos + Y_tos)*(V - monitored[114])*X_tos

    # Expressions for the X_gate component
    monitored[40] = 1.0/(1 + math.exp(-1/5. - V/15.))
    monitored[41] = 0.5 + 9/(1 + math.exp(1/5. + V/15.))
    monitored[130] = (-X_tos + monitored[40])/monitored[41]

    # Expressions for the Y_gate component
    monitored[42] = 1.0/(1 + 28.5027336438*math.exp(V/10.))
    monitored[43] = 30 + 3000/(1 + math.exp(6 + V/10.))
    monitored[131] = (-Y_tos + monitored[42])/monitored[43]

    # Expressions for the R_gate component
    monitored[44] = 1.0/(1 + 28.5027336438*math.exp(V/10.))
    monitored[45] = 220 + 2800.0/(1 + math.exp(6 + V/10.))
    monitored[132] = (-R_tos + monitored[44])/monitored[45]

    # Expressions for the Itof component
    monitored[46] = G_tof*(V - monitored[114])*X_tof*Y_tof

    # Expressions for the Itof X gate component
    monitored[47] = 1.0/(1 + math.exp(-1/5. - V/15.))
    monitored[48] = 1.5 + 3.5*math.exp(-(V*V)/900.)
    monitored[133] = (-X_tof + monitored[47])/monitored[48]

    # Expressions for the Itof Y gate component
    monitored[49] = 1.0/(1 + 28.5027336438*math.exp(V/10.))
    monitored[50] = 20 + 20/(1 + 28.5027336438*math.exp(V/10.))
    monitored[134] = (monitored[49] - Y_tof)/monitored[50]

    # Expressions for the K1 gate component
    monitored[118] = 1.02/(1 +\
        7.35454251046e-07*math.exp(-0.2385*monitored[114] + 0.2385*V))
    monitored[119] = (1.15340563519e-16*math.exp(-0.06175*monitored[114] +\
        0.06175*V) + 0.762624006506*math.exp(-0.08032*monitored[114] +\
        0.08032*V))/(1 + 0.0867722941577*math.exp(0.5143*monitored[114] -\
        0.5143*V))
    monitored[120] = monitored[118]/(monitored[118] + monitored[119])

    # Expressions for the ICl Ca component
    monitored[51] = G_Cl*(-monitored[115] + V)*(Fx_Cl_SL/(1 + Kd_ClCa/Ca_SL)\
        + Fx_Cl_jct1/(1 + Kd_ClCa/Ca_jct1))

    # Expressions for the IClb component
    monitored[52] = G_ClBk*(-monitored[115] + V)

    # Expressions for the d gate component
    monitored[63] = 1.0/(1 + 0.0892185174093*math.exp(-V/6.))
    monitored[64] = (1 -\
        0.0892185174093*math.exp(-V/6.))*monitored[63]/(0.5075 + 0.035*V)
    monitored[135] = (monitored[63] - d)/monitored[64]

    # Expressions for the f gate component
    monitored[65] = 0.6/(1 + math.exp(5/2. - V/20.)) + 1.0/(1 +\
        16964.681259*math.exp(0.277777777778*V))
    monitored[66] = 1.0/(0.02 + 0.0197*math.exp(-((0.48865 +\
        0.0337*V)*(0.48865 + 0.0337*V))))
    monitored[136] = (monitored[65] - f)/monitored[66]

    # Expressions for the FCa gate component
    monitored[61] = 1 - fCaB_SL
    monitored[62] = 1 - fCaB_jct1
    monitored[137] = -0.0119*fCaB_SL + lccCaInact*(1 - fCaB_SL)*Ca_SL
    monitored[138] = lccCaInact*(1 - fCaB_jct1)*Ca_jct1 - 0.0119*fCaB_jct1

    # Expressions for the INaCa component
    monitored[67] = (-math.pow(Nao, HNa)*Ca_jct1*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)) + Cao*math.pow(Na_jct1,\
        HNa)*math.exp(F*eta*V/(Rgas*T)))/(1 + ksat*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)))
    monitored[68] = (-math.pow(Nao, HNa)*Ca_SL*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)) + Cao*math.pow(Na_SL,\
        HNa)*math.exp(F*eta*V/(Rgas*T)))/(1 + ksat*math.exp(F*(-1 +\
        eta)*V/(Rgas*T)))
    monitored[69] = math.pow(Q10_NCX, -31 + T/10.)
    monitored[70] = 1.0/(1 + (Kd_act*Kd_act*Kd_act)/(Ca_SL*Ca_SL*Ca_SL))
    monitored[71] = 1.0/(1 + (Kd_act*Kd_act*Kd_act)/(Ca_jct1*Ca_jct1*Ca_jct1))
    monitored[72] =\
        Fx_NCX_jct1*V_max_INaCa*monitored[67]*monitored[69]*monitored[71]/(K_mCai*math.pow(Nao,\
        HNa)*(1 + math.pow(Na_jct1/K_mNai, HNa)) + K_mCao*math.pow(Na_jct1,\
        HNa) + Cao*math.pow(Na_jct1, HNa) + math.pow(Nao, HNa)*Ca_jct1 +\
        math.pow(K_mNao, HNa)*(1 + Ca_jct1/K_mCai)*Ca_jct1)
    monitored[73] =\
        Fx_NCX_SL*V_max_INaCa*monitored[68]*monitored[69]*monitored[70]/(math.pow(K_mNao,\
        HNa)*(1 + Ca_SL/K_mCai)*Ca_SL + math.pow(Nao, HNa)*Ca_SL +\
        Cao*math.pow(Na_SL, HNa) + K_mCao*math.pow(Na_SL, HNa) +\
        K_mCai*math.pow(Nao, HNa)*(1 + math.pow(Na_SL/K_mNai, HNa)))
    monitored[74] = monitored[72] + monitored[73]

    # Expressions for the ICap component
    monitored[75] = math.pow(Q10_SLCaP, -31 + T/10.)
    monitored[76] = Fx_SLCaP_jct1*V_maxAF*monitored[75]/(1 +\
        math.pow(Km/Ca_jct1, H_ICap))
    monitored[77] = Fx_SLCaP_SL*V_maxAF*monitored[75]/(1 + math.pow(Km/Ca_SL,\
        H_ICap))
    monitored[78] = monitored[77] + monitored[76]

    # Expressions for the ICab component
    monitored[79] = Fx_CaBk_jct1*G_CaBk*(V - monitored[112])
    monitored[80] = Fx_CaBk_SL*G_CaBk*(V - monitored[113])
    monitored[81] = monitored[80] + monitored[79]

    # Expressions for the Jrel SR component
    monitored[82] = Max_SR - (Max_SR - Min_SR)/(1 + math.pow(EC50_SR/Ca_SR,\
        HSR))
    monitored[83] = koCa/monitored[82]
    monitored[84] = kiCa*monitored[82]
    monitored[85] = 1 - I - O - R
    monitored[139] = -Ca_jct1*R*monitored[84] -\
        (Ca_jct1*Ca_jct1)*R*monitored[83] + kim*monitored[85] + kom*O
    monitored[140] = kim*I - Ca_jct1*O*monitored[84] +\
        (Ca_jct1*Ca_jct1)*R*monitored[83] - kom*O
    monitored[141] = -kim*I + Ca_jct1*O*monitored[84] - kom*I +\
        (Ca_jct1*Ca_jct1)*monitored[83]*monitored[85]
    monitored[86] = ks*(Ca_SR - Ca_jct1)*O

    # Expressions for the Jleak SR component
    monitored[87] = KSRleak*(Ca_SR - Ca_jct1)

    # Expressions for the Jpump SR component
    monitored[88] = math.pow(Q10_SRCaP, -31 + T/10.)
    monitored[89] = V_max_Jpump*(-math.pow(Ca_SR/Kmr, H_Jpump) +\
        math.pow(Cai/Kmf, H_Jpump))*monitored[88]/(1 + math.pow(Ca_SR/Kmr,\
        H_Jpump) + math.pow(Cai/Kmf, H_Jpump))

    # Expressions for the Ion diffusion component
    monitored[121] = -1.8313e-14*Na_SL + 1.8313e-14*Na_jct1
    monitored[122] = 1.6386e-12*Na_SL - 1.6386e-12*Nai
    monitored[123] = 8.2413e-13*Ca_jct1 - 8.2413e-13*Ca_SL
    monitored[124] = -3.7243e-12*Cai + 3.7243e-12*Ca_SL

    # Expressions for the Cytosolic component
    monitored[99] = kon_TroponinC*(-Ca_TroponinC + Bmax_TroponinC)*Cai -\
        koff_TroponinC*Ca_TroponinC
    monitored[100] = -koff_TroponinC_Ca_Mg_Ca*Ca_TroponinC_Ca_Mg +\
        kon_TroponinC_Ca_Mg_Ca*(-Ca_TroponinC_Ca_Mg - Mg_TroponinC_Ca_Mg +\
        Bmax_TroponinC_Ca_Mg_Ca)*Cai
    monitored[101] = Mgi*kon_TroponinC_Ca_Mg_Mg*(-Ca_TroponinC_Ca_Mg +\
        Bmax_TroponinC_Ca_Mg_Mg - Mg_TroponinC_Ca_Mg) -\
        koff_TroponinC_Ca_Mg_Mg*Mg_TroponinC_Ca_Mg
    monitored[102] = -koff_Calmodulin*Ca_Calmodulin +\
        kon_Calmodulin*(Bmax_Calmodulin - Ca_Calmodulin)*Cai
    monitored[103] = kon_Myosin_Ca*(Bmax_Myosin_Ca - Mg_Myosin -\
        Ca_Myosin)*Cai - koff_Myosin_Ca*Ca_Myosin
    monitored[104] = Mgi*kon_Myosin_Mg*(Bmax_Myosin_Mg - Mg_Myosin -\
        Ca_Myosin) - koff_Myosin_Mg*Mg_Myosin
    monitored[105] = kon_SRB*(-Ca_SRB + Bmax_SRB)*Cai - koff_SRB*Ca_SRB
    monitored[106] = monitored[102] + monitored[100] + monitored[104] +\
        monitored[101] + monitored[103] + monitored[105] + monitored[99]
    monitored[142] = monitored[99]
    monitored[143] = monitored[100]
    monitored[144] = monitored[101]
    monitored[145] = monitored[102]
    monitored[146] = monitored[103]
    monitored[147] = monitored[104]
    monitored[148] = monitored[105]

    # Expressions for the IKr component
    monitored[23] = 0.0129099444874*math.sqrt(Ko)
    monitored[24] = (V - monitored[114])*Xr*monitored[23]*monitored[27]

    # Expressions for the IK1 component
    monitored[116] = 0.387298334621*math.sqrt(Ko)
    monitored[117] = (V - monitored[114])*monitored[116]*monitored[120]

    # Expressions for the ICaL component
    monitored[53] = math.pow(Q10_CaL, -31 + T/10.)
    monitored[54] = 0.45*(F*F)*V*d*f*monitored[53]/(Rgas*T)
    monitored[55] =\
        4*Fx_ICaL_jct1*PCa*(gamma_Cai*Ca_jct1*math.exp(2*F*V/(Rgas*T)) -\
        Cao*gamma_Cao)*monitored[54]*monitored[62]/(-1 +\
        math.exp(2*F*V/(Rgas*T)))
    monitored[56] = Fx_ICaL_jct1*PNa*(-Nao*gamma_Nao +\
        gamma_Nai*Na_jct1*math.exp(F*V/(Rgas*T)))*monitored[54]*monitored[62]/(-1 +\
        math.exp(F*V/(Rgas*T)))
    monitored[57] =\
        4*Fx_ICaL_SL*PCa*(gamma_Cai*Ca_SL*math.exp(2*F*V/(Rgas*T)) -\
        Cao*gamma_Cao)*monitored[54]*monitored[61]/(-1 +\
        math.exp(2*F*V/(Rgas*T)))
    monitored[58] = Fx_ICaL_SL*PNa*(-Nao*gamma_Nao +\
        gamma_Nai*Na_SL*math.exp(F*V/(Rgas*T)))*monitored[54]*monitored[61]/(-1 +\
        math.exp(F*V/(Rgas*T)))
    monitored[59] = PK*(Fx_ICaL_SL*monitored[61] +\
        Fx_ICaL_jct1*monitored[62])*(-Ko*gamma_Ko +\
        Ki*gamma_Ki*math.exp(F*V/(Rgas*T)))*monitored[54]/(-1 +\
        math.exp(F*V/(Rgas*T)))
    monitored[60] = monitored[59] + monitored[58] + monitored[56] +\
        monitored[55] + monitored[57]

    # Expressions for the Na buffer component
    monitored[107] = kon*(Bmax_jct1 - Na_jct1_buf)*Na_jct1 - koff*Na_jct1_buf
    monitored[108] = -koff*Na_SL_buf + kon*(-Na_SL_buf + Bmax_SL)*Na_SL
    monitored[149] = monitored[107]
    monitored[150] = monitored[108]
    monitored[151] = -monitored[107] - monitored[121]/monitored[3] -\
        Cm*(monitored[15] + monitored[56] + monitored[6] + 3*monitored[72] +\
        3*monitored[20])/(F*monitored[3])
    monitored[152] = (monitored[121] - monitored[122])/monitored[2] -\
        Cm*(monitored[16] + monitored[58] + 3*monitored[21] + 3*monitored[73]\
        + monitored[7])/(F*monitored[2]) - monitored[108]
    monitored[153] = monitored[122]/monitored[4]

    # Expressions for the Ca buffer component
    monitored[90] = -koff_Calsequestrin*Ca_Calsequestrin +\
        kon_Calsequestrin*(Bmax_Calsequestrin*monitored[4]/monitored[1] -\
        Ca_Calsequestrin)*Ca_SR
    monitored[154] = monitored[90]
    monitored[91] = kon_SL*(-Ca_SLB_SL +\
        Bmax_SLB_SL*monitored[4]/monitored[2])*Ca_SL - koff_SLB*Ca_SLB_SL
    monitored[92] = kon_SL*(0.1*Bmax_SLB_jct1*monitored[4]/monitored[3] -\
        Ca_SLB_jct1)*Ca_jct1 - koff_SLB*Ca_SLB_jct1
    monitored[93] = -koff_SLHigh*Ca_SLHigh_SL +\
        kon_SL*(Bmax_SLHigh_SL*monitored[4]/monitored[2] -\
        Ca_SLHigh_SL)*Ca_SL
    monitored[94] = -koff_SLHigh*Ca_SLHigh_jct1 + kon_SL*(-Ca_SLHigh_jct1 +\
        0.1*Bmax_SLHigh_jct1*monitored[4]/monitored[3])*Ca_jct1
    monitored[155] = monitored[91]
    monitored[156] = monitored[92]
    monitored[157] = monitored[93]
    monitored[158] = monitored[94]
    monitored[95] = monitored[92] + monitored[94]
    monitored[96] = monitored[93] + monitored[91]
    monitored[97] = -2*monitored[72] + monitored[55] + monitored[76] +\
        monitored[79]
    monitored[98] = monitored[80] + monitored[77] - 2*monitored[73] +\
        monitored[57]
    monitored[159] = -monitored[4]*monitored[87]/monitored[1] - monitored[86]\
        + monitored[89] - monitored[90]
    monitored[160] = -Cm*monitored[97]/(2*F*monitored[3]) +\
        monitored[4]*monitored[87]/monitored[3] - monitored[95] -\
        monitored[123]/monitored[3] + monitored[1]*monitored[86]/monitored[3]
    monitored[161] = -Cm*monitored[98]/(2*F*monitored[2]) - monitored[96] +\
        (-monitored[124] + monitored[123])/monitored[2]
    monitored[162] = -monitored[106] -\
        monitored[1]*monitored[89]/monitored[4] + monitored[124]/monitored[4]

    # Expressions for the Cell component
    monitored[109] = (-stim_amplitude if\
        -stim_period*math.floor(t/stim_period) + t <= stim_start +\
        stim_duration and -stim_period*math.floor(t/stim_period) + t >=\
        stim_start else 0)
    monitored[163] = -monitored[8] - monitored[109] - monitored[74] -\
        monitored[81] - monitored[17] - monitored[35] - monitored[78] -\
        monitored[60] - monitored[52] - monitored[22] - monitored[24] -\
        monitored[46] - monitored[117] - monitored[39] - monitored[38] -\
        monitored[51]

    # Return results
    return monitored
