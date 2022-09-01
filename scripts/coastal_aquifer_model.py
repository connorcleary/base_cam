import numpy as np
import flopy
import os
import pickle
import flopy.utils.binaryfile as bf
from pars import ModelParameters, load_parameters


def build_steady_model(pars):
    '''
        A function to build a coastal aquifer model.

        Input: 
            pars: parameters object

        Outputs:
            swt: groundwater flow and transport model
    '''
    
    model_ws = f"./model_files/{pars.name}"
    if not os.path.exists(model_ws):
        os.makedirs(model_ws)

    modelname = pars.name
    swt = flopy.seawat.Seawat(modelname, model_ws=model_ws, exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")

    delr = pars.Lx/pars.ncol
    delc = pars.Ly/pars.nrow
    delv = pars.Lz/pars.nlay

    top = 0.0
    botm = np.linspace(top - delv, -pars.Lz, pars.nlay)

    ipakcb = 53

    dis = flopy.modflow.ModflowDis(
            swt,
            pars.nlay,
            pars.nrow,
            pars.ncol,
            nper=1,
            itmuni=4,
            delr=delr,
            delc=delc,
            laycbd=0,
            top=top,
            botm=botm,
            perlen=pars.perlen, # COULD BE A POTENTIAL PROBLEM
            nstp=pars.nstp,
            steady=True
        )

    # Variables for the BAS package
    ibound = np.ones((pars.nlay, pars.nrow, pars.ncol), dtype=np.int32)
    ibound[:, :, 0] = -1
    ibound[:, :, -1] = -1

    if not pars.confined:
        ibound[0, :, int(pars.ncol*pars.onshore_proportion):pars.ncol] = -1

    strt = np.zeros((pars.nlay, pars.nrow, pars.ncol))

    bas = flopy.modflow.ModflowBas(swt, ibound, strt=strt)

    laytyp=np.zeros(50)
    laywet=np.zeros(50)
    
    if not pars.confined:
        laytyp[0] = 1
        laywet[0] = 1

    lpf = flopy.modflow.ModflowLpf(swt, hk=0, vka=0, ipakcb=ipakcb, laytyp=laytyp, laywet=laywet)
    pcg = flopy.modflow.ModflowPcg(swt, hclose=1.0e-5, npcond=1, mxiter=500)

    oc_spd = {} 
    for kstp in range(pars.nstp):
            oc_spd[(0, kstp)] = ["save head", "save budget"]

    oc = flopy.modflow.ModflowOc(swt, stress_period_data=oc_spd, compact=True)

    # find inland boundary cells 
    onshore_boundary_cells = []
    offshore_boundary_cells = []
    for k in range(pars.nlay):
        for j in range(pars.nrow):
            onshore_boundary_cells.append([k, j, 0])
            offshore_boundary_cells.append([k, j, pars.ncol-1])

    if not pars.confined:
        for i in range(int(pars.ncol*pars.onshore_proportion), pars.ncol):
            for j in range(pars.nrow):
                offshore_boundary_cells.append([0, j, i])

    # set up constant head stress period data
    chd_data = {}
    chd_sp1 = []
    # wel data
    wel_data = {}
    wel_sp1 = []
    # Set up ssm 
    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    ssm_data = {}
    ssm_sp1 = []

    for cell in onshore_boundary_cells:
        ssm_sp1.append([cell[0], cell[1], cell[2], 0.0, itype["BAS6"]])
        chd_sp1.append([cell[0], cell[1], cell[2], pars.h_b, pars.h_b])

    for cell in offshore_boundary_cells:
        ssm_sp1.append([cell[0], cell[1], cell[2], 35.0, itype["BAS6"]])
        chd_sp1.append([cell[0], cell[1], cell[2], 0.0, 0.0])
        
    ssm_data[0] = ssm_sp1
    chd_data[0] = chd_sp1
    wel_data[0] = wel_sp1 

    chd = flopy.modflow.ModflowChd(swt, stress_period_data=chd_data, ipakcb = ipakcb)

    if type(pars.sconc) == int:
        sconc = pars.sconc*np.ones((pars.nlay, pars.nrow, pars.ncol))
        sconc[:,:,-1] = 35
    else: 
        sconc = pars.sconc

    btn = flopy.mt3d.Mt3dBtn(
        swt,
        nprs=-1,
        prsity=0.35, # can change this
        sconc=sconc, # can change this: to each area having different starti,
        chkmas=False,
        nprobs=10,
        nprmas=10,
        dt0=pars.dt
    )

    adv = flopy.mt3d.Mt3dAdv(swt, 
        mixelm=0,
        dceps=1.0e-5,
        nplane=1,
        npl=16,
        nph=16,
        npmin=4,
        npmax=32,
        dchmoc=1.0e-3,
        nlsink=1,
        npsink=16,
        percel=0.5)
    # sip = flopy.modflow.ModflowSip(swt)
    # lmt = flopy.modflow.ModflowLmt(swt)
    dsp = flopy.mt3d.Mt3dDsp(swt, al=0.4, trpt=0.1, trpv=0.01, dmcoef=1e-9)
    gcg = flopy.mt3d.Mt3dGcg(swt, iter1=500, mxiter=1, isolve=2, cclose=1e-5)
    mxss = int(np.ceil(2*pars.nlay*pars.nrow + pars.nrow*pars.ncol))#*(1-pars.onshore_proportion)+1))
    ssm = flopy.mt3d.Mt3dSsm(swt, stress_period_data=ssm_data, mxss=mxss)

    vdf = flopy.seawat.SeawatVdf(
        swt,
        iwtable=0,
        densemin=0,
        densemax=0,
        denseref=1000.0,
        denseslp=0.7143,
        firstdt=pars.dt,
    )

    fname = r"./MT3D001.UCN"
    if os.path.isfile(fname):
        os.remove(fname)

    swt.write_input()

    return swt 


def run_model(swt):
    """
        A function to run the seawat model

        Inputs: 
            swt: model object
        Outputs:
            None
    """
    swt.write_input()
    success, buff = swt.run_model(silent=False, report=True)
    if not success:
        raise Exception("SEAWAT did not terminate normally.")


def extract_results(name):
    """
        Open model results from binary files

        Inputs:
            name: name of model/realization/scenario
        Outputs:
            head: head matrix [nstp, nlay, nrow, ncol]
            qx: longitudinal flux matrix [nstp, nlay, nrow, ncol]
            qy: transverse flux matrix matrix [nstp, nlay, nrow, ncol]
            qz: vertical flux matrix matrix [nstp, nlay, nrow, ncol]
            concentration: concentration matrix [nstp, nlay, nrow, ncol]
    """
    pars = load_parameters(name)
    name = pars.name
    model_ws = f"./model_files\\{name}"
    nstp = pars.perlen/pars.dt

    # open binary files
    ucnobj = bf.UcnFile(os.path.join(model_ws, "MT3D001.UCN"))
    cbbobj = bf.CellBudgetFile(os.path.join(model_ws, f'{name}.cbc'))
    headobj = bf.HeadFile(os.path.join(model_ws, f'{name}.hds'))

    # get head and concentration data
    concentration = ucnobj.get_alldata()[:]
    head = headobj.get_alldata()[:]
    
    # select every n items
    times = ucnobj.get_times()
    concentration = concentration

    qx = np.zeros_like(concentration)
    qy = np.zeros_like(concentration)
    qz = np.zeros_like(concentration)

    # get fluxes
    for t in range(qx.shape[0]):
        qx[t] = cbbobj.get_data(text="flow right face", totim=times[t])[0]
        if pars.nrow > 1:
            qy[t] = cbbobj.get_data(text="flow front face", totim=times[t])[0]
        qz[t] = cbbobj.get_data(text="flow lower face", totim=times[t])[0]

    save_results(name, concentration, head, qx, qy, qz)
    return concentration, head, qx, qy, qz


def save_results(name, concentration, head, qx, qy, qz):
    """
        Save extracted results to a .npy file

        Inputs:
            name: model name
            concentration, head etc. : numpy arrays of model outputs
        Outputs:
            None
    """
    ws = os.path.join(f'.\\results\\{name}')
    if not os.path.exists(ws):
        os.makedirs(ws)

    with open(os.path.join(ws, f"qx.npy"), 'wb') as f: np.save(f, np.array(qx))
    with open(os.path.join(ws, f"qy.npy"), 'wb') as f: np.save(f, np.array(qy))
    with open(os.path.join(ws, f"qz.npy"), 'wb') as f: np.save(f, np.array(qz))
    with open(os.path.join(ws, f"head.npy"), 'wb') as f: np.save(f, np.array(head))
    with open(os.path.join(ws, f"concentration.npy"), 'wb') as f: np.save(f, np.array(concentration))


def load_results(name):
    """
        Load extracted results from .npy files

        Inputs:
            name: name of the model
        Outputs:
            concentration, head... : numpy matrices of results
    """
    ws = os.path.join(f'.\\results\\{name}')

    with open(os.path.join(ws, f"qx.npy"), 'rb') as f: qx = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"qy.npy"), 'rb') as f: qy = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"qz.npy"), 'rb') as f: qz = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"head.npy"), 'rb') as f: head = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"concentration.npy"), 'rb') as f: concentration = np.load(f, allow_pickle=True)

    return concentration, head, qx, qy, qz, 