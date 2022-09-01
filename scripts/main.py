import coastal_aquifer_model as cam
from pars import ModelParameters, load_parameters
import results
import parameter_utils as par_utils
import numpy as np
import post_processing as proc
from multiprocessing import Pool

def run_realization(ppars):
    name,n,pars = ppars 
    swt = cam.build_steady_model(pars)
    Kh, Kv =  par_utils.load_field(f"./fields/pirot3D/80x40x50realization{n}.mat")
    Kh = np.transpose(Kh, (2, 0, 1))[:,:,:]
    Kv = np.transpose(Kv, (2, 0, 1))[:,:,:]
    swt.lpf.hk = Kh
    swt.lpf.vka = Kv
    cam.run_model(swt)
    # concentration, head, qx, qy, qz = cam.extract_results(f"{name}{n}")


def run_ensemble(name=None, realizations=1, **kwargs):
    
    ppars = []
    for n in range(realizations):
        ppars.append([name, n, ModelParameters(f"{name}{n}", **kwargs)])
    
    p = Pool(processes=8)
    p.map(run_realization, ppars)
    # for n in range(realizations):
    #     #pars = ModelParameters(f"{name}{n}", **kwargs)
    #     swt = cam.build_steady_model(pars)
    #     Kh, Kv =  par_utils.load_field(f"./fields/pirot3D/80x40x50realization{n}.mat")
    #     Kh = np.transpose(Kh, (2, 0, 1))[:,:,:]
    #     Kv = np.transpose(Kv, (2, 0, 1))[:,:,:]
    #     swt.lpf.hk = Kh
    #     swt.lpf.vka = Kv
    #     cam.run_model(swt)
    #     concentration, head, qx, qy, qz = cam.extract_results(f"{name}{n}")



def plot_all_rows(name, realizations):
    for i in range(realizations):
        results.plot_results(name+str(i), row=20)

def main():
    # results.plot_effective_K(realizations=30)
    #run_ensemble("test", 1, h_b=0.6, perlen=1e6)
    # plot_all_rows("pirot3d", 1)
    #proc.get_steady_state_time_evolutions("pirot3d", 5)
    #results.compare_time_evolution("pirot3d", 30)
    #results.plot_surface_evolution("pirot3d_multi", 30, isochlors=[0.1, 0.9], nsteps=30)
    cam.extract_results("test0")

if __name__=="__main__":
    main()