from pars import ModelParameters, load_parameters
import numpy as np
import coastal_aquifer_model as cam
import matplotlib.pyplot as plt
import plot_helpers as plth
import os
import post_processing as proc
import parameter_utils as par_utils
from matplotlib.gridspec import GridSpec
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import cm

def plot_results(name, timestep=-1, row=0, return_axs=False, figsize=(12,6),
    cmap="viridis", arrow_c="white", aspect=8, vector_T=10, width=0.002, fmt="%3.2f"):
    """
        Plot the head, concentrations and fluxes

        Inputs:
            name: model to plot
            timestep: timestep to plot, default is -1 (the last one)
            row: row to plot
            return_axs: flag whether to return the axes objects
            figsize: figure dimensions (inches)
            cmap: colormap
            arrow_c: color of arrows
            aspect: vertical exageration 
            vector_T: spaces between vector arrows
            width: arrow width
            fmt: format of contour labels
        Outputs:
            axs: axes objects (optional)
    """

    # load parameters and results
    pars = load_parameters(name)
    concentration, head, qx, qy, qz = cam.load_results(name)

    f, axs = plt.subplots(2, 1, figsize=figsize)

    # set up x and y arrays, to be distance above sea level and distance onshore
    if not pars.confined:
        x = np.linspace(-pars.Lx*pars.offshore_proportion, pars.Lx-pars.Lx*pars.offshore_proportion, pars.ncol)
        y = np.linspace(-pars.sea_level, pars.Lz-pars.sea_level, pars.nlay)
    else:
        x = np.linspace(0, pars.Lx, pars.ncol)
        y = np.linspace(0, pars.Lz, pars.nlay)

    # select relevent slice in time and the alongshore direction, and set values above the water table as nan
    concentration_array = concentration[timestep,:,row,:]
    if not pars.confined:
        head_array = head[timestep,:,row,:] - pars.sea_level*np.ones_like(head[timestep,:,row,:])
        for i in range(pars.nlay):
            for j in range(pars.ncol):
                if concentration_array[i, j] == np.float32(1.e30):
                    head_array[i,j] = np.nan
                    concentration_array[i, j] = np.nan
    else:
        head_array = head[timestep,:,row,:]
    # plot head colormesh
    headcm = axs[0].pcolormesh(x, y, np.flipud(head_array), 
            	                cmap=cmap, vmax=np.nanmax(head_array[:, 1:]), vmin=np.min([0, pars.h_b]))

    # plot head contours
    hc = axs[0].contour(x, y, np.flipud(head_array), colors=arrow_c, 
                        levels=np.linspace(np.min([0, pars.h_b]), np.nanmax(head_array[:, 1:]), 15))

    # label contours
    axs[0].clabel(hc, hc.levels, inline=True, fontsize=10, fmt=fmt)

    # plot concentration colormesh
    conccm = axs[1].pcolormesh(x, y, np.flipud(concentration_array), 
            	                cmap=cmap, vmax=35, vmin=0)

    # plot arrows
    axs[1].quiver(plth.sample_grid_data(x, vector_T, vector_T), plth.sample_grid_data(y, vector_T, vector_T), 
                    plth.sample_grid_data(np.flipud(qx[timestep,:,row,:]), vector_T, vector_T), 
                    -plth.sample_grid_data(np.flipud(qz[timestep,:,row,:]), vector_T, vector_T), 
                    color=arrow_c, width=width)

    axs[0].set_aspect(aspect)
    axs[1].set_aspect(aspect)
    axs[0].set_title("Head")
    axs[1].set_title("Salinity")
    axs[0].set_ylabel("Height above sealevel (m)")
    axs[1].set_ylabel("Height above sealevel (m)")
    axs[1].set_xlabel("Distance onshore (m)")

    f.suptitle(f"Head and salinity distributions for {name}")

    headcb = plt.colorbar(headcm, shrink=1, ax=axs[0])
    conccb = plt.colorbar(conccm, shrink=1, ax=axs[1])
    headcb.ax.set_title('Head (m)', fontsize = 'small')
    conccb.ax.set_title('Salinity (kg/m^3)', fontsize = 'small')
    
    ws = os.path.join(f'.\\figures\\{name}')
    if not os.path.exists(ws):
        os.makedirs(ws)
    plt.savefig(f"{ws}\\head_and_concentration", dpi=300)

    # return axs objects if necessary
    if return_axs: return axs


def plot_effective_K(realizations=30):
    
    hk_array = np.zeros(realizations)
    vk_array = np.zeros(realizations)
    pars = ModelParameters()

    for n in range(realizations):
        Kh, Kv =  par_utils.load_field(f".\\fields\\pirot3D\\80x40x50realization{n}.mat")
        Kh = np.transpose(Kh, (2, 0, 1))[:,:,:]
        Kv = np.transpose(Kv, (2, 0, 1))[:,:,:]
        hk_array[n], vk_array[n], _ = par_utils.calculate_K_eff(pars, Kh, Kv)

    f, ax = plt.subplots()
    ax.scatter(hk_array, vk_array)
    plt.show()


def compare_time_evolution(name, realizations):

    com, toe, mix = proc.load_steady_metric_evolutions(name)
    for n in range(realizations):
        pars = load_parameters(f"{name}{n}")
        
        concentration, head, qx, qy, qz = cam.load_results(f"{name}{n}")
        Kh, Kv =  par_utils.load_field(f".\\fields\\pirot3D\\80x40x50realization{n}.mat")
        Kh = np.transpose(Kh, (2, 0, 1))[:,:,:]

        for row in range(pars.nrow):
            
            f = plt.figure(constrained_layout=True, figsize=(7,7))
            gs = GridSpec(4, 2, figure=f, width_ratios=[1,2])

            axleg = f.add_subplot(gs[0, 0])
            axtime = f.add_subplot(gs[1:, 0])
            axhk = f.add_subplot(gs[0, 1])
            axmax = f.add_subplot(gs[1, 1])
            axmin = f.add_subplot(gs[2, 1])
            axss = f.add_subplot(gs[3, 1])

            axhk.set_aspect(10)
            axmax.set_aspect(10)
            axmin.set_aspect(10)
            axss.set_aspect(10)

            t = 1/365*1e6*np.logspace(1, 3, base=10, num=30).astype(int)
            # load time evolutions
            
            imax = np.atleast_1d(np.argmax(com[n,row,:]))[0]
            imin = np.atleast_1d(np.argmin(com[n,row,imax:]))[0] + imax
            iss = com.shape[-1] -1
            for i in range(com[0, row, imin:].shape[0]):
                if np.average(np.absolute(concentration[i+imin,:,row,:]-concentration[i+imin-1,:,row,:])) < 0.05:
                    iss = i + imin
                    break
            imax = 13

            nmax = int(t[imax]//(1e6/365))
            nmin = int(t[imin]//(1e6/365))
            nss = int(t[iss]//(1e6/365))

            cmax = concentration[nmax,:,:,:]
            cmin = concentration[nmin,:,:,:] 
            css = concentration[nss,:,:,:]
            qxmax = qx[nmax,:,:,:]
            qxmin = qx[nmin,:,:,:] 
            qxss  = qx[nss,:,:,:]
            qzmax  = qz[nmax,:,:,:]
            qzmin  = qz[nmin,:,:,:]
            qzss = qz[nss,:,:,:]

            x = np.linspace(0, 800, 80)
            y = np.linspace(-25, 0, 50)

            cmhk = axhk.pcolormesh(x, y, np.flipud(np.log(Kh[:,row,:])), cmap="coolwarm", vmax=-1, vmin=-13)
            cmmax = axmax.pcolormesh(x, y, np.flipud(cmax[:,row,:]), cmap="viridis", vmax=35, vmin=0)
            cmmin = axmin.pcolormesh(x, y, np.flipud(cmin[:,row,:]), cmap="viridis", vmax=35, vmin=0)
            cmss = axss.pcolormesh(x, y, np.flipud(css[:,row,:]), cmap="viridis", vmax=35, vmin=0)

            axmax.quiver(plth.sample_grid_data(x, 3, 3), plth.sample_grid_data(y, 3, 3), 
                            plth.sample_grid_data(np.flipud(qxmax[:,0,:]), 3,3), -plth.sample_grid_data(np.flipud(qzmax[:,0,:]), 3,3), color="white")
            axmin.quiver(plth.sample_grid_data(x, 3, 3), plth.sample_grid_data(y, 3, 3), 
                            plth.sample_grid_data(np.flipud(qxmin[:,row,:]), 3,3), -plth.sample_grid_data(np.flipud(qzmin[:,row,:]), 3,3), color="white")
            axss.quiver(plth.sample_grid_data(x, 3, 3), plth.sample_grid_data(y, 3, 3), 
                            plth.sample_grid_data(np.flipud(qxss[:,row,:]), 3,3), -plth.sample_grid_data(np.flipud(qzss[:,row,:]), 3,3), color="white")

            
            axtime.invert_yaxis()
            axtime.plot(com[0,row,:], t)
            axtime.axhline(y=t[imax], color = "red", linestyle=(0,(1,2)), label="Time 0")
            axtime.axhline(y=t[imin], color = "blue", linestyle=(1,(1,2)), label="Time 1")
            axtime.axhline(y=t[iss], color = "green", linestyle=(2,(1,2)), label="Time 2")

            axtime.set_yscale('log')

            axhk.set_title("Horizontal hydraulic conductivity")
            axmax.set_title("Time 0")
            axmin.set_title("Time 1")
            axss.set_title("Time 2")
            axtime.set_title("Mixing zone centroid")

            axmin.set_ylabel("Depth (m)")
            axss.set_xlabel("Distance offshore (m)")
            axtime.set_xlabel("Distance offshore (m)")
            axtime.set_ylabel("time (log[years])", labelpad=10, rotation=270)

            h, l = axtime.get_legend_handles_labels()
            axleg.axis("off")
            axleg.legend(h, l, loc="lower center", borderaxespad=0)

            f.suptitle("Saline concentrations at different states")

            cb1 = plt.colorbar(cmhk,ax=axhk, location="right", shrink=0.75, pad=0.05)
            cb2 = plt.colorbar(cmmax,ax=[axmax, axmin, axss], location="right", shrink=0.75/3, pad=0.05)

            cb1.ax.set_title('log10[hk]', fontsize = 'small')
            cb2.ax.set_title('C (kg/m^3)', fontsize = 'small')

            results_location = f'.\\results\\{name}'
            plt.savefig(f"{results_location}\\steady3D_evolution{n}{row}", dpi=300)


def animate_func(num, pars, isochlor_evo, isochlors, Y, Z, ax):
    ax.clear()
    viridis = cm.get_cmap('viridis', 256)
    for i, iso in enumerate(isochlors):
        ax.plot_surface((pars.Lx/pars.ncol)*isochlor_evo[num, i, :, :], Y, Z, alpha=0.5, color=viridis(iso))
    ax.set_box_aspect((2, 1, 1))
    ax.set_xlim(0, 800)
    pass


def plot_surface_evolution(name, realizations, isochlors=[0.1], nsteps=30):

    t = np.logspace(1, 3, base=10, num=nsteps).astype(int)-1
    for n in range(realizations):
        pars = load_parameters(f"{name}{n}")
        concentration, _, _, _, _ = cam.load_results(f"{name}{n}")
        y = np.linspace(0, pars.Ly, pars.nrow)
        z = np.linspace(0, -pars.Lz, pars.nlay)
        Y, Z = np.meshgrid(y, z)
        isochlor_evo= np.zeros((nsteps, len(isochlors), pars.nlay, pars.nrow))
        for i, time in enumerate(t):
            for lay in range(pars.nlay):
                for row in range(pars.nrow):
                        for j, iso in enumerate(isochlors):
                            if np.array(concentration[time, lay, row, :]>iso*35)[-1]:
                                isochlor_evo[i, j, lay, row] = np.argmax(concentration[time, lay, row, :]>iso*35)
                            else:
                                isochlor_evo[i, j, lay, row] = pars.ncol-1

        fig = plt.figure(figsize =(14, 9))
        ax = plt.axes(projection ='3d')
        ani = animation.FuncAnimation(fig, animate_func, fargs = (pars, isochlor_evo, isochlors, Y, Z, ax), interval=100,   
                                   frames=nsteps)
        ax.view_init(30, 225)
        f = f"animate_func_{n}.gif"
        writergif = animation.PillowWriter(fps=nsteps/6)
        ani.save(f, writer=writergif)
        plt.show()

