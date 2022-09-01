import numpy as np
import flopy
from scipy.io import loadmat
import flopy.utils.binaryfile as bf
import matplotlib.pyplot as plt
from pars import ModelParameters, load_parameters
import os

def load_field(name):
    """
        Load parameter fields from .mat file

        Input:
            name of field
        Output:
            relevant parameter fields
    """
    fields = loadmat(name)
    Kh = fields['model_Kh']
    Kv = fields['model_Kv']
    # others 

    return Kh, Kv


def build_effective_model(pars, hk, vka, direction):
    """
        Build groundwater flow model for calculating effective conductivities
    """
    model_ws = f".\\temp\\{pars.name}"
    if not os.path.exists(model_ws):
        os.makedirs(model_ws)

    gwf = flopy.modflow.Modflow(pars.name, model_ws=model_ws, exe_name=r"C:\Users\ccl124\bin\MF2005.exe")

    delr = pars.Lx/pars.ncol
    delc = pars.Ly/pars.nrow
    delv = pars.Lz/pars.nlay

    top = 0.0 
    botm = np.linspace(top - delv, -pars.Lz, pars.nlay)

    dis = flopy.modflow.ModflowDis(
        gwf, 
        pars.nlay, 
        pars.nrow, 
        pars.ncol, 
        delr=delr, 
        delc=delc, 
        top=top,
        botm=botm,
        itmuni=1,
        perlen= 1,
        tsmult = 1,
        nstp=1,
    )

    ibound = np.ones((pars.nlay, pars.nrow, pars.ncol), dtype=np.int32)
    strt = np.ones((pars.nlay, pars.nrow, pars.ncol), dtype=np.float32)

    if direction.lower()[0] == "h":
        ibound[:, :, 0] = -1
        ibound[:, :, -1] = -1
        strt[:, :, 0] = pars.h_b
        strt[:, :, -1] = 0
    else:
        ibound[0, :, :] = -1
        ibound[-1, :, :] = -1
        strt[0, :, :] = pars.h_b
        strt[-1, :, :] = 0.0

    bas = flopy.modflow.ModflowBas(
        gwf, 
        ibound=ibound, 
        strt=strt
    )

    lpf = flopy.modflow.ModflowLpf(
        gwf,
        hk=hk,
        vka=vka, 
        ipakcb=53,
        laytyp=0
        )

    spd = {(0, 0): ["print head", "print budget", "save head", "save budget"]}
    oc = flopy.modflow.ModflowOc(
        gwf, 
        stress_period_data=spd, 
        compact=True
    )

    pcg = flopy.modflow.ModflowPcg(gwf, mxiter=500)

    return gwf


def run_effective_model(gwf):
    gwf.write_input()
    success, buff = gwf.run_model()
    if not success:
        raise Exception("MODFLOW did not terminate normally.")


def extract_effective_model(modelname):

    model_ws = f".\\temp\\{modelname}"
    if not os.path.exists(model_ws):
        os.makedirs(model_ws)

    cbb = bf.CellBudgetFile(f"{model_ws}\\{modelname}.cbc")
    times = cbb.get_times()
    kstpkper_list = cbb.get_kstpkper()
    frf = cbb.get_data(text="FLOW RIGHT FACE", totim=times[-1])[0]
    flf = cbb.get_data(text="FLOW LOWER FACE", totim=times[-1])[0]
    return frf, flf


def calculate_K_eff(pars, hk, vka):

    for direction in ["horizontal", "vertical"]:
        gwf = build_effective_model(pars, hk, vka, direction)
        run_effective_model(gwf)
        frf, flf = extract_effective_model(pars.name)
        if direction[0] == "v":
            Kv_eff = np.sum(flf[0,:,:])/(pars.h_b/(pars.Lz))/(pars.Lx*pars.Ly)
        else:
            qinflow = np.sum(frf[:,:,0])
            Kh_eff = np.sum(frf[:, :, 0])/(pars.h_b/pars.Lx)/(pars.Lz*pars.Ly)

    return Kv_eff, Kh_eff, qinflow