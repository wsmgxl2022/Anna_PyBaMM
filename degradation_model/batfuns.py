import pybamm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
from scipy import interpolate, integrate 
from scipy.integrate import solve_ivp


def set_rc_params(rcParams):

    rcParams["lines.markersize"] = 5
    rcParams["lines.linewidth"] = 2
    rcParams["xtick.minor.visible"] = True
    rcParams["ytick.minor.visible"] = True
    rcParams["font.size"] = 12
    rcParams["legend.fontsize"] = 10
    rcParams["legend.frameon"] = False
    rcParams["font.family"] = 'serif'
    rcParams['font.serif'] = 'Times New Roman'
    rcParams['mathtext.rm'] = 'serif'
    rcParams['mathtext.it'] = 'serif:italic'
    rcParams['mathtext.bf'] = 'serif:bold'
    rcParams['mathtext.fontset'] = 'custom'
    rcParams["savefig.bbox"] = "tight"
    rcParams["axes.grid"] = True
    rcParams["axes.axisbelow"] = True
    rcParams["grid.linestyle"] = "--"
    rcParams["grid.color"] = (0.8, 0.8, 0.8)
    rcParams["grid.alpha"] = 0.5
    rcParams["grid.linewidth"] = 0.5
    rcParams['figure.dpi'] = 150
    rcParams['savefig.dpi'] = 600
    rcParams['figure.max_open_warning']=False
    
    return rcParams

def nmc_volume_change_mohtat(sto,c_s_max):
    t_change = -1.10/100*(1-sto)
    return t_change

def graphite_volume_change_mohtat(sto,c_s_max):
    stoichpoints = np.array([0,0.12,0.18,0.24,0.50,1])
    thicknesspoints = np.array([0,2.406/100,3.3568/100,4.3668/100,5.583/100,13.0635/100])
    x = [sto]
    t_change = pybamm.Interpolant(stoichpoints, thicknesspoints, x, name=None, interpolator='linear', extrapolate=True, entries_string=None)
    return t_change


def get_parameter_values():
    parameter_values = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Mohtat2020)
    parameter_values.update(
        {
            # mechanical properties
            "Positive electrode Poisson's ratio": 0.3,
            "Positive electrode Young's modulus [Pa]": 375e9,
            "Positive electrode reference concentration for free of deformation [mol.m-3]": 0,
            "Positive electrode partial molar volume [m3.mol-1]": 7.28e-7,
            "Positive electrode volume change": nmc_volume_change_mohtat,
            # Loss of active materials (LAM) model
            "Positive electrode LAM constant exponential term": 2,
            "Positive electrode critical stress [Pa]": 375e6,
            # mechanical properties
            "Negative electrode Poisson's ratio": 0.2,
            "Negative electrode Young's modulus [Pa]": 15e9,
            "Negative electrode reference concentration for free of deformation [mol.m-3]": 0,
            "Negative electrode partial molar volume [m3.mol-1]": 3.1e-6,   
            "Negative electrode volume change": graphite_volume_change_mohtat,
            # Loss of active materials (LAM) model
            "Negative electrode LAM constant exponential term": 2,
            "Negative electrode critical stress [Pa]": 60e6,
            # Other
            "Cell thermal expansion coefficient [m.K-1]": 1.48E-6,
            "Lower voltage cut-off [V]": 3.0,
            # Initializing Particle Concentration
            # "Initial concentration in negative electrode [mol.m-3]": x100*parameter_values["Maximum concentration in negative electrode [mol.m-3]"],
            # "Initial concentration in positive electrode [mol.m-3]": y100*parameter_values["Maximum concentration in positive electrode [mol.m-3]"]
        },
        check_already_exists=False,
    )
    return parameter_values


def split_long_string(title, max_words=None):
    """Get title in a nice format"""
    max_words = max_words or pybamm.settings.max_words_in_line
    words = title.split()
    # Don't split if fits on one line, don't split just for units
    if len(words) <= max_words or words[max_words].startswith("["):
        return title
    else:
        first_line = (" ").join(words[:max_words])
        second_line = (" ").join(words[max_words:])
        return first_line + "\n" + second_line

def cycle_adaptive_simulation(model, parameter_values, experiment,SOC_0=1, save_at_cycles=None,drive_cycle=None):
    experiment_one_cycle = pybamm.Experiment(
        experiment.operating_conditions_cycles[:1],
        termination=experiment.termination_string,
        cccv_handling=experiment.cccv_handling,
        drive_cycles={"DriveCycle": drive_cycle},
    )
    Vmin = 3.0
    Vmax = 4.2
    esoh_model = pybamm.lithium_ion.ElectrodeSOH()
    esoh_sim = pybamm.Simulation(esoh_model, parameter_values=parameter_values)
    param = model.param
    esoh_solver = pybamm.lithium_ion.ElectrodeSOHSolver(parameter_values, param)
    Cn = parameter_values.evaluate(param.n.cap_init)
    Cp = parameter_values.evaluate(param.p.cap_init)
    eps_n = parameter_values["Negative electrode active material volume fraction"]
    eps_p = parameter_values["Positive electrode active material volume fraction"]
    C_over_eps_n = Cn / eps_n
    C_over_eps_p = Cp / eps_p
    c_n_max = parameter_values.evaluate(param.n.prim.c_max)
    c_p_max = parameter_values.evaluate(param.p.prim.c_max)
    n_Li_init = parameter_values.evaluate(param.n_Li_particles_init)
    
    esoh_sol = esoh_sim.solve(
        [0],
        inputs={"V_min": Vmin, "V_max": Vmax, "C_n": Cn, "C_p": Cp, "n_Li": n_Li_init},
        solver=pybamm.AlgebraicSolver(),
    )

    parameter_values.update(
        {
            "Initial concentration in negative electrode [mol.m-3]": esoh_sol[
                "x_100"
            ].data[0]
            * c_n_max,
            "Initial concentration in positive electrode [mol.m-3]": esoh_sol[
                "y_100"
            ].data[0]
            * c_p_max,
            
        }
    )
    parameter_values.update(
        {
            "Negative electrode LAM min stress [Pa]": 0,
            "Negative electrode LAM max stress [Pa]": 0,
            "Positive electrode LAM min stress [Pa]": 0,
            "Positive electrode LAM max stress [Pa]": 0,
        },
    )

    sim_ode = pybamm.Simulation(
        model, experiment=experiment_one_cycle, parameter_values=parameter_values,
        solver=pybamm.CasadiSolver("safe")
    )
    sol0 = sim_ode.solve(initial_soc=SOC_0)
    model = sim_ode.solution.all_models[0]
    cap0 = sol0.summary_variables["Capacity [A.h]"]

    def sol_to_y(sol, loc="end"):
        if loc == "start":
            pos = 0
        elif loc == "end":
            pos = -1
        model = sol.all_models[0]
        n_Li = sol["Total lithium in particles [mol]"].data[pos].flatten()
        Cn = sol["Negative electrode capacity [A.h]"].data[pos].flatten()
        Cp = sol["Positive electrode capacity [A.h]"].data[pos].flatten()
        # y = np.concatenate([n_Li, Cn, Cp])
        y = n_Li
        for var in model.initial_conditions:
            if var.name not in [
                "X-averaged negative particle concentration",
                "X-averaged positive particle concentration",
                "Discharge capacity [A.h]",
                "Porosity times concentration",
            ]:
                value = sol[var.name].data
                if value.ndim == 1:
                    value = value[pos]
                elif value.ndim == 2:
                    value = np.average(value[:, pos])
                elif value.ndim == 3:
                    value = np.average(value[:, :, pos])
                y = np.concatenate([y, value.flatten()])
            elif var.name == "Porosity times concentration":
                for child in var.children:
                    value = sol[child.name].data
                    if value.ndim == 1:
                        value = value[pos]
                    elif value.ndim == 2:
                        value = np.average(value[:, pos])
                    elif value.ndim == 3:
                        value = np.average(value[:, :, pos])
                    y = np.concatenate([y, value.flatten()])
        return y

    def y_to_sol(y, esoh_sim, model):
        n_Li = y[0]
        Cn = C_over_eps_n * y[1]
        Cp = C_over_eps_p * y[2]

        esoh_sol = esoh_sim.solve(
            [0],
            inputs={"V_min": Vmin, "V_max": Vmax, "C_n": Cn, "C_p": Cp, "n_Li": n_Li},
        )
        esoh_sim.built_model.set_initial_conditions_from(esoh_sol)
        ics = {}
        x_100 = esoh_sol["x_100"].data[0]
        y_100 = esoh_sol["y_100"].data[0]
        x_0 = esoh_sol["x_0"].data[0]
        y_0 = esoh_sol["y_0"].data[0]
        start = 1
        for var in model.initial_conditions:
            if var.name == "X-averaged negative particle concentration":
                ics[var.name] = ((x_100-x_0)*SOC_0+x_0) * np.ones((model.variables[var.name].size, 2))
            elif var.name == "X-averaged positive particle concentration":
                ics[var.name] = ((y_100-y_0)*SOC_0+y_0)  * np.ones((model.variables[var.name].size, 2))
            elif var.name == "Discharge capacity [A.h]":
                ics[var.name] = np.zeros(1)
            else:
                if var.name == "Porosity times concentration":
                    for child in var.children:
                        # end = start + model.variables[child.name].size
                        # ics[child.name] = y[start:end, np.newaxis]
                        end = start + 1
                        ics[child.name] = y[start] * np.ones((model.variables[var.name].size, 1))
                        start = end
                else:
                    # end = start + model.variables[var.name].size
                    # ics[var.name] = y[start:end, np.newaxis]
                    end = start + 1
                    ics[var.name] = y[start] * np.ones((model.variables[var.name].size, 1))
                    start = end
        model.set_initial_conditions_from(ics)
        return pybamm.Solution(
            [np.array([0])],
            model.concatenated_initial_conditions.evaluate()[:, np.newaxis],
            model,
            {},
        )

    def dydt(t, y):
        if y[0] < 0 or y[1] < 0 or y[2] < 0:
            return 0 * y

        # print(t)
        # Set up based on current value of y
        y_to_sol(
            y,
            esoh_sim,
            sim_ode.op_conds_to_built_models[
                experiment_one_cycle.operating_conditions[0]["electric"]
            ],
        )

        # Simulate one cycle
        sol = sim_ode.solve()

        dy = sol_to_y(sol) - y

        return dy

    if experiment.termination == {}:
        event = None
    else:

        def capacity_cutoff(t, y):
            sol = y_to_sol(y, esoh_sim, model)
            cap = pybamm.make_cycle_solution([sol], esoh_solver, True)[1]["Capacity [A.h]"]
            return cap / cap0 - experiment_one_cycle.termination["capacity"][0] / 100

        capacity_cutoff.terminal = True

    num_cycles = len(experiment.operating_conditions_cycles)
    if save_at_cycles is None:
        t_eval = np.arange(1, num_cycles + 1)
    elif save_at_cycles == -1:
        t_eval = None
    else:
        t_eval = np.arange(1, num_cycles + 1, save_at_cycles)
    y0 = sol_to_y(sol0, loc="start")
    timer = pybamm.Timer()
    sol = solve_ivp(
        dydt,
        [1, num_cycles],
        y0,
        t_eval=t_eval,
        events=capacity_cutoff,
        first_step=10,
        method="RK23",
        atol=1e-2,
        rtol=1e-2,
    )
    time = timer.time()

    all_sumvars = []
    for idx in range(sol.y.shape[1]):
        fullsol = y_to_sol(sol.y[:, idx], esoh_sim, model)
        sumvars = pybamm.make_cycle_solution([fullsol], esoh_solver, True)[1]
        all_sumvars.append(sumvars)

    all_sumvars_dict = {
        key: np.array([sumvars[key] for sumvars in all_sumvars])
        for key in all_sumvars[0].keys()
    }
    all_sumvars_dict["Cycle number"] = sol.t
    
    all_sumvars_dict["cycles evaluated"] = sol.nfev
    all_sumvars_dict["solution time"] = time
    
    return all_sumvars_dict

def cycle_adaptive_simulation_V2(model, parameter_values, experiment,SOC_0=1, save_at_cycles=None,drive_cycle=None):
    experiment_one_cycle = pybamm.Experiment(
        experiment.operating_conditions_cycles[:1],
        termination=experiment.termination_string,
        cccv_handling=experiment.cccv_handling,
        drive_cycles={"DriveCycle": drive_cycle},
    )
    Vmin = 3.0
    Vmax = 4.2
    esoh_model = pybamm.lithium_ion.ElectrodeSOH()
    esoh_sim = pybamm.Simulation(esoh_model, parameter_values=parameter_values)
    param = model.param
    esoh_solver = pybamm.lithium_ion.ElectrodeSOHSolver(parameter_values, param)
    Cn = parameter_values.evaluate(param.n.cap_init)
    Cp = parameter_values.evaluate(param.p.cap_init)
    eps_n = parameter_values["Negative electrode active material volume fraction"]
    eps_p = parameter_values["Positive electrode active material volume fraction"]
    C_over_eps_n = Cn / eps_n
    C_over_eps_p = Cp / eps_p
    c_n_max = parameter_values.evaluate(param.n.prim.c_max)
    c_p_max = parameter_values.evaluate(param.p.prim.c_max)
    n_Li_init = parameter_values.evaluate(param.n_Li_particles_init)
    
    esoh_sol = esoh_sim.solve(
        [0],
        inputs={"V_min": Vmin, "V_max": Vmax, "C_n": Cn, "C_p": Cp, "n_Li": n_Li_init},
        solver=pybamm.AlgebraicSolver(),
    )

    parameter_values.update(
        {
            "Initial concentration in negative electrode [mol.m-3]": esoh_sol[
                "x_100"
            ].data[0]
            * c_n_max,
            "Initial concentration in positive electrode [mol.m-3]": esoh_sol[
                "y_100"
            ].data[0]
            * c_p_max,
            
        }
    )
    model1 = model
    sim_ode = pybamm.Simulation(
        model, experiment=experiment_one_cycle, parameter_values=parameter_values,
        solver=pybamm.CasadiSolver("safe")
    )
    sol0 = sim_ode.solve(initial_soc=SOC_0)
    model = sim_ode.solution.all_models[0]
    cap0 = sol0.summary_variables["Capacity [A.h]"]

    Omega_n = param.n.Omega
    R_n = parameter_values.evaluate(param.n.prim.R_typ)
    E_n = param.n.E
    nu_n = param.n.nu
    CC_n = parameter_values.evaluate(2*E_n/(1-nu_n)/3)

    Omega_p = param.n.Omega
    R_p = parameter_values.evaluate(param.p.prim.R_typ)
    E_p = param.p.E
    nu_p = param.p.nu
    CC_p = parameter_values.evaluate(2*E_p/(1-nu_p)/3)

    def graphite_volume_change(sto):
        stoichpoints = np.array([0,0.12,0.18,0.24,0.50,1])
        thicknesspoints = np.array([0,2.406/100,3.3568/100,4.3668/100,5.583/100,13.0635/100])
        x = [sto]
        t_change = pybamm.Interpolant(stoichpoints, thicknesspoints, x, name=None, interpolator='linear', extrapolate=True, entries_string=None)
        t_change = np.interp(x,stoichpoints,thicknesspoints)
        return t_change

    def nmc_volume_change(sto):
        t_change = -1.10/100*(1-sto)
        return t_change

    def sigma_hfun(c_s_n,c_s_p):
        Rvec_n = np.linspace(0,R_n,len(c_s_n))
        Rvec_p = np.linspace(0,R_p,len(c_s_p))
        sigma_h_n = []
        sigma_h_p = []
        for nn in range(np.size(c_s_n,1)):
            c_s_n1 = c_s_n[:,nn]
            y_s = np.vectorize(graphite_volume_change)(c_s_n1)
            cube_n = Rvec_n[1:]**3-Rvec_n[:-1]**3
            mul_n = (y_s[1:]+y_s[:-1])/2
            sigma_h_s_n = CC_n*((1/R_n**3)*1/3*np.sum(cube_n*mul_n)-1/3*y_s[-1])
            sigma_h_n.append(sigma_h_s_n)
            c_s_p1 = c_s_p[:,nn]
            x_s = np.vectorize(nmc_volume_change)(c_s_p1)
            cube_p = Rvec_p[1:]**3-Rvec_p[:-1]**3
            mul_p = (x_s[1:]+x_s[:-1])/2
            sigma_h_s_p = CC_p*((1/R_p**3)*1/3*np.sum(cube_p*mul_p)-1/3*x_s[-1])
            sigma_h_p.append(sigma_h_s_p)
        sigma_h_n = np.array(sigma_h_n)
        sigma_h_p = np.array(sigma_h_p)

        return sigma_h_n,sigma_h_p

    def sol_to_y(sol, loc="end"):
        if loc == "start":
            pos = 0
        elif loc == "end":
            pos = -1
        model = sol.all_models[0]
        n_Li = sol["Total lithium in particles [mol]"].data[pos].flatten()
        Cn = sol["Negative electrode capacity [A.h]"].data[pos].flatten()
        Cp = sol["Positive electrode capacity [A.h]"].data[pos].flatten()
        # y = np.concatenate([n_Li, Cn, Cp])
        y = n_Li
        for var in model.initial_conditions:
            if var.name not in [
                "X-averaged negative particle concentration",
                "X-averaged positive particle concentration",
                "Discharge capacity [A.h]",
                "Porosity times concentration",
            ]:
                value = sol[var.name].data
                if value.ndim == 1:
                    value = value[pos]
                elif value.ndim == 2:
                    value = np.average(value[:, pos])
                elif value.ndim == 3:
                    value = np.average(value[:, :, pos])
                y = np.concatenate([y, value.flatten()])
            elif var.name == "Porosity times concentration":
                for child in var.children:
                    value = sol[child.name].data
                    if value.ndim == 1:
                        value = value[pos]
                    elif value.ndim == 2:
                        value = np.average(value[:, pos])
                    elif value.ndim == 3:
                        value = np.average(value[:, :, pos])
                    y = np.concatenate([y, value.flatten()])
        return y

    def y_to_sol(y, esoh_sim, model):
        n_Li = y[0]
        Cn = C_over_eps_n * y[1]
        Cp = C_over_eps_p * y[2]

        esoh_sol = esoh_sim.solve(
            [0],
            inputs={"V_min": Vmin, "V_max": Vmax, "C_n": Cn, "C_p": Cp, "n_Li": n_Li},
        )
        esoh_sim.built_model.set_initial_conditions_from(esoh_sol)
        ics = {}
        x_100 = esoh_sol["x_100"].data[0]
        y_100 = esoh_sol["y_100"].data[0]
        x_0 = esoh_sol["x_0"].data[0]
        y_0 = esoh_sol["y_0"].data[0]
        start = 1
        for var in model.initial_conditions:
            if var.name == "X-averaged negative particle concentration":
                ics[var.name] = ((x_100-x_0)*SOC_0+x_0) * np.ones((model.variables[var.name].size, 2))
            elif var.name == "X-averaged positive particle concentration":
                ics[var.name] = ((y_100-y_0)*SOC_0+y_0)  * np.ones((model.variables[var.name].size, 2))
            elif var.name == "Discharge capacity [A.h]":
                ics[var.name] = np.zeros(1)
            else:
                if var.name == "Porosity times concentration":
                    for child in var.children:
                        # end = start + model.variables[child.name].size
                        # ics[child.name] = y[start:end, np.newaxis]
                        end = start + 1
                        ics[child.name] = y[start] * np.ones((model.variables[var.name].size, 1))
                        start = end
                else:
                    # end = start + model.variables[var.name].size
                    # ics[var.name] = y[start:end, np.newaxis]
                    end = start + 1
                    ics[var.name] = y[start] * np.ones((model.variables[var.name].size, 1))
                    start = end
        model.set_initial_conditions_from(ics)
        return pybamm.Solution(
            [np.array([0])],
            model.concatenated_initial_conditions.evaluate()[:, np.newaxis],
            model,
            {},
        )

    def dydt(t, y):
        if y[0] < 0 or y[1] < 0 or y[2] < 0:
            return 0 * y

        # print(t)
        # Set up based on current value of y
        y_to_sol(
            y,
            esoh_sim,
            sim_ode.op_conds_to_built_models[
                experiment_one_cycle.operating_conditions[0]["electric"]
            ],
        )

        # Simulate one cycle
        sol = sim_ode.solve()
        dy = sol_to_y(sol) - y
        t =  sol["Time [s]"].entries
        t = t/3600
        c_s_n = sol["X-averaged negative particle concentration"].entries
        c_s_p = sol["X-averaged positive particle concentration"].entries
        # sigma_hs_n,sigma_hs_p = sigma_hfun(c_s_n,c_s_p)
        sigma_ts_n = sol["X-averaged negative particle surface tangential stress [Pa]"].entries
        sigma_rs_n = sol["X-averaged negative particle surface radial stress [Pa]"].entries
        sigma_hs_n = (sigma_rs_n+2*sigma_ts_n)/2
        sigma_ts_p = sol["X-averaged positive particle surface tangential stress [Pa]"].entries
        sigma_rs_p = sol["X-averaged positive particle surface radial stress [Pa]"].entries
        sigma_hs_p = (sigma_rs_p+2*sigma_ts_p)/2
        dnli = dy[0]
        beta_LAM_n = param.n.beta_LAM_dimensional
        beta_LAM2_n = param.n.beta_LAM_dimensional2
        m_LAM_n = param.n.m_LAM
        stress_critical_n = param.n.stress_critical_dim
        j_stress_LAM_n = parameter_values.evaluate(-beta_LAM_n*(abs(min(sigma_hs_n)) / stress_critical_n) ** m_LAM_n + beta_LAM2_n*(abs(max(sigma_hs_n)) / stress_critical_n) ** m_LAM_n)
        act_n_loss = j_stress_LAM_n*t[-1]*3600
        # print(f"beta1= {parameter_values.evaluate(param.n.beta_LAM_dimensional)}")
        # print(f"beta2= {parameter_values.evaluate(param.n.beta_LAM_dimensional2)}")
        # print(f"min stress = {abs(min(sigma_hs_n))}")
        # print(f"max stress = {abs(max(sigma_hs_n))}")
        # print(f"j = {j_stress_LAM_n}")
        # print(f"loss = {act_n_loss}")
        C_n_loss = parameter_values.evaluate(act_n_loss*(param.n.L * param.n.prim.c_max * param.F* param.A_cc)/3600)
        dCn = act_n_loss+dy[1]
        c_save_n1 = sol["R-averaged negative particle concentration"].entries
        c_save_n = c_save_n1[1,:]
        dnli += parameter_values.evaluate(3600/param.F)*C_n_loss*np.average(c_save_n)
        beta_LAM_p = param.p.beta_LAM_dimensional
        beta_LAM2_p = param.p.beta_LAM_dimensional2
        m_LAM_p = param.p.m_LAM
        stress_critical_p = param.p.stress_critical_dim
        j_stress_LAM_p = parameter_values.evaluate(-beta_LAM_p*(abs(max(sigma_hs_p)) / stress_critical_p) ** m_LAM_p + beta_LAM2_p*(abs(min(sigma_hs_p)) / stress_critical_p) ** m_LAM_p)
        act_p_loss = j_stress_LAM_p*t[-1]*3600
        C_p_loss = parameter_values.evaluate(act_p_loss*(param.p.L * param.p.prim.c_max * param.F* param.A_cc)/3600)
        dCp = act_p_loss+dy[1]
        c_save_p1 = sol["R-averaged positive particle concentration"].entries
        c_save_p = c_save_p1[1,:]
        dnli += parameter_values.evaluate(3600/param.F)*C_p_loss*np.average(c_save_p)
        dy2 = np.zeros(6)
        dy2[0] = dnli
        dy2[1] = dCn
        dy2[2] = dCp
        dy2[3] = dy[3]
        dy2[4] = dy[4]
        dy2[5] = dy[5]

        # parameter_values.update(
        #     {
        #         "Negative electrode LAM min stress [Pa]": min(sigma_hs_n),
        #         "Negative electrode LAM max stress [Pa]": max(sigma_hs_n),
        #         "Positive electrode LAM min stress [Pa]": min(sigma_hs_p),
        #         "Positive electrode LAM max stress [Pa]": max(sigma_hs_p),
        #     },
        # )

        # sim_ode1 = pybamm.Simulation(
        #     model1, experiment=experiment_one_cycle, parameter_values=parameter_values,
        #     solver=pybamm.CasadiSolver("safe")
        # )
        # sol1 = sim_ode1.solve(initial_soc=SOC_0)
        # # model1 = sim_ode1.solution.all_models[0]

        # y_to_sol(
        #     y,
        #     esoh_sim,
        #     sim_ode1.op_conds_to_built_models[
        #         experiment_one_cycle.operating_conditions[0]["electric"]
        #     ],
        # )

        # # # Simulate one cycle
        # sol = sim_ode1.solve()

        # dy = sol_to_y(sol) - y

        return dy2

    if experiment.termination == {}:
        event = None
    else:

        def capacity_cutoff(t, y):
            sol = y_to_sol(y, esoh_sim, model)
            cap = pybamm.make_cycle_solution([sol], esoh_solver, True)[1]["Capacity [A.h]"]
            return cap / cap0 - experiment_one_cycle.termination["capacity"][0] / 100

        capacity_cutoff.terminal = True

    num_cycles = len(experiment.operating_conditions_cycles)
    if save_at_cycles is None:
        t_eval = np.arange(1, num_cycles + 1)
    elif save_at_cycles == -1:
        t_eval = None
    else:
        t_eval = np.arange(1, num_cycles + 1, save_at_cycles)
    y0 = sol_to_y(sol0, loc="start")
    timer = pybamm.Timer()
    sol = solve_ivp(
        dydt,
        [1, num_cycles],
        y0,
        t_eval=t_eval,
        events=capacity_cutoff,
        first_step=10,
        method="RK23",
        atol=1e-2,
        rtol=1e-2,
    )
    time = timer.time()

    all_sumvars = []
    for idx in range(sol.y.shape[1]):
        fullsol = y_to_sol(sol.y[:, idx], esoh_sim, model)
        sumvars = pybamm.make_cycle_solution([fullsol], esoh_solver, True)[1]
        all_sumvars.append(sumvars)

    all_sumvars_dict = {
        key: np.array([sumvars[key] for sumvars in all_sumvars])
        for key in all_sumvars[0].keys()
    }
    all_sumvars_dict["Cycle number"] = sol.t
    
    all_sumvars_dict["cycles evaluated"] = sol.nfev
    all_sumvars_dict["solution time"] = time
    
    return all_sumvars_dict

def plot(all_sumvars_dict,esoh_data):
    esoh_vars = ["x_0", "y_0", "x_100", "y_100", "C_n", "C_p"]
    # esoh_vars = ["Capacity [A.h]", "Loss of lithium inventory [%]",
    #              "Loss of active material in negative electrode [%]",
    #              "Loss of active material in positive electrode [%]"]
    fig, axes = plt.subplots(3,2,figsize=(7,7))
    for k, name in enumerate(esoh_vars):
        ax = axes.flat[k]
        ax.plot(all_sumvars_dict["Cycle number"],all_sumvars_dict[name],"ro")
        ax.plot(esoh_data["N"],esoh_data[name],"kx")
#         ax.scatter(all_sumvars_dict["Cycle number"],all_sumvars_dict[name],color="r")
    #     ax.plot(long_sol.summary_variables[name],"b-")
        ax.set_title(split_long_string(name))
        if k>3:
            ax.set_xlabel("Cycle number")
    # fig.subplots_adjust(bottom=0.4)
    fig.legend(["Acc Sim"] + ["Reported"], 
           loc="lower center", ncol=1, fontsize=11)
    fig.tight_layout()
    return fig

def plot1(all_sumvars_dict,esoh_data):
    esoh_vars = ["Capacity [A.h]","n_Li"]
    esoh_data["Capacity [A.h]"]=esoh_data["Cap"]
    param = pybamm.LithiumIonParameters()
    esoh_data["n_Li"]= 3600/param.F.value*(esoh_data["y_100"]*esoh_data["C_p"]+esoh_data["x_100"]*esoh_data["C_n"])
    fig, axes = plt.subplots(2,1,figsize=(7,7))
    for k, name in enumerate(esoh_vars):
        ax = axes.flat[k]
        ax.plot(all_sumvars_dict["Cycle number"],all_sumvars_dict[name],"ro")
        ax.plot(esoh_data["N"],esoh_data[name],"kx")
        ax.set_title(split_long_string(name))
        ax.set_xlabel("Cycle number")
    fig.legend(["Acc Sim"] + ["Reported"], 
           loc="upper right", ncol=1, fontsize=11)
    fig.tight_layout()
    return fig

def plotc(all_sumvars_dict,esoh_data):
    esoh_vars = ["x_100", "y_0", "C_n", "C_p", "Capacity [A.h]", "Loss of lithium inventory [%]"]
    fig, axes = plt.subplots(3,2,figsize=(7,7))
    for k, name in enumerate(esoh_vars):
        ax = axes.flat[k]
        ax.plot(all_sumvars_dict["Cycle number"],all_sumvars_dict[name],"ro")
        ax.plot(esoh_data["N"],esoh_data[name],"kx")
        ax.set_title(split_long_string(name))
        # if k ==2 or k==3:
        #     ax.set_ylim([3,6.2])
        if k>3:
            ax.set_xlabel("Cycle number")
    fig.legend(["Sim"] + ["Data"], 
           loc="lower center",bbox_to_anchor=[0.5,-0.05], ncol=1, fontsize=11)
    fig.tight_layout()
    return fig

def plotd(esoh_data):
    esoh_vars = ["x_100", "y_0", "C_n", "C_p", "Capacity [A.h]", "Loss of lithium inventory [%]"]
    fig, axes = plt.subplots(3,2,figsize=(7,7))
    for k, name in enumerate(esoh_vars):
        ax = axes.flat[k]
        ax.plot(esoh_data["N"],esoh_data[name],"kx")
        ax.set_title(split_long_string(name))
        if k ==2 or k==3:
            ax.set_ylim([3,6.2])
        if k>3:
            ax.set_xlabel("Cycle number")
    fig.tight_layout()
    return fig

def plotc2(all_sumvars_dict1,all_sumvars_dict2,esoh_data,leg1="sim1",leg2="sim2"):
    esoh_vars = ["x_100", "y_0", "C_n", "C_p", "Capacity [A.h]", "Loss of lithium inventory [%]"]
    fig, axes = plt.subplots(3,2,figsize=(7,7))
    for k, name in enumerate(esoh_vars):
        ax = axes.flat[k]
        ax.plot(all_sumvars_dict1["Cycle number"],all_sumvars_dict1[name],"b--")
        ax.plot(all_sumvars_dict2["Cycle number"],all_sumvars_dict2[name],"r--")
        ax.plot(esoh_data["N"],esoh_data[name],"kx")
        ax.set_title(split_long_string(name))
        # if k ==2 or k==3:
        #     ax.set_ylim([3,6.2])
        if k>3:
            ax.set_xlabel("Cycle number")
    fig.legend([leg1, leg2 , "Data"], 
           loc="lower center",bbox_to_anchor=[0.5,-0.1], ncol=1, fontsize=11)
    fig.tight_layout()
    return fig

def plotcomp(all_sumvars_dict0,all_sumvars_dict1):
    esoh_vars = ["x_100", "y_0", "C_n", "C_p", "Capacity [A.h]", "Loss of lithium inventory [%]"]
    fig, axes = plt.subplots(3,2,figsize=(7,7))
    for k, name in enumerate(esoh_vars):
        ax = axes.flat[k]
        ax.plot(all_sumvars_dict0["Cycle number"],all_sumvars_dict0[name],"k")
        ax.plot(all_sumvars_dict1["Cycle number"],all_sumvars_dict1[name],"b")
        ax.set_title(split_long_string(name))
        if k ==2 or k==3:
            ax.set_ylim([3,6.2])
        if k>3:
            ax.set_xlabel("Cycle number")
    fig.legend(["Baseline"] + ["Sim"], 
           loc="lower center",bbox_to_anchor=[0.5,-0.05], ncol=1, fontsize=11)
    fig.tight_layout()
    return fig

def plotcomplong(all_sumvars_dict0,all_sumvars_dict1,all_sumvars_dict2):
    esoh_vars = ["x_100", "y_0", "C_n", "C_p", "Capacity [A.h]", "Loss of lithium inventory [%]"]
    fig, axes = plt.subplots(3,2,figsize=(7,7))
    for k, name in enumerate(esoh_vars):
        ax = axes.flat[k]
        ax.plot(all_sumvars_dict0["Cycle number"],all_sumvars_dict0[name],"kx")
        ax.plot(all_sumvars_dict1["Cycle number"],all_sumvars_dict1[name],"bo")
        ax.plot(all_sumvars_dict2["Cycle number"],all_sumvars_dict2[name],"m.")
        ax.set_title(split_long_string(name))
        if k ==2 or k==3:
            ax.set_ylim([3,6.2])
        if k>3:
            ax.set_xlabel("Cycle number")
    fig.legend(["Baseline"] + ["Accl Sim"] + ["Long Sim"], 
            loc="lower center",bbox_to_anchor=[0.5,-0.02], ncol=1, fontsize=11)
    fig.tight_layout()
    return fig


def load_data_calendar(cell,eSOH_DIR,oCV_DIR):
    param = pybamm.LithiumIonParameters() 
    cell_no = f'{cell:02d}'
    dfe=pd.read_csv(eSOH_DIR+"aging_param_cell_"+cell_no+".csv")
    dfo=pd.read_csv(oCV_DIR+"ocv_data_cell_"+cell_no+".csv")
    # if cell_no=='24':
    #     dfe = dfe.drop(dfe.index[0])
    #     dfe = dfe.reset_index(drop=True)
    #     dfe['Time']=dfe['Time']-dfe['Time'][0]
    dfe['N']=dfe['Time']
    dfe["Capacity [A.h]"]=dfe["Cap"]
    dfe["n_Li"]= 3600/param.F.value*(dfe["y_100"]*dfe["C_p"]+dfe["x_100"]*dfe["C_n"])
    dfe["Loss of lithium inventory [%]"]=(1-dfe["n_Li"]/dfe["n_Li"][0])*100
    N =dfe.N.unique()

    # print("Cycle Numbers:")
    # print(*N, sep = ", ") 
    # print(len(N_0))
    # print(len(dfo_0))
    # rev_exp = []
    # irrev_exp = []

    return cell_no,dfe,N,dfo

def init_exp_calendar(cell_no,dfe,param,parameter_values):
    # dfe_0 = dfe[dfe['N']==N[0]]
    C_n_init = dfe['C_n'][0]
    C_p_init = dfe['C_p'][0]
    y_0_init = dfe['y_0'][0] 
    eps_n_data = parameter_values.evaluate(C_n_init*3600/(param.n.L * param.n.prim.c_max * param.F* param.A_cc))
    eps_p_data = parameter_values.evaluate(C_p_init*3600/(param.p.L * param.p.prim.c_max * param.F* param.A_cc))
    # cs_p_init = parameter_values.evaluate(y_0_init* param.c_p_max)
    if cell_no=='22':
        SOC_0 = 1
        Temp = 45
    #     x_init = esoh_sol["x_100"].data[0] 
    #     y_init = esoh_sol["y_100"].data[0] 
    elif cell_no=='23':
        SOC_0 = 0.5
        Temp = 45
    elif cell_no=='24':
        SOC_0 = 0.98
        Temp = -5
    elif cell_no=='25':
        SOC_0 = 0.5
        Temp = -5
        
    return eps_n_data,eps_p_data,SOC_0,Temp#,x_init,y_init


def load_data(cell,eSOH_DIR,oCV_DIR): 
    param = pybamm.LithiumIonParameters() 
    cell_no = f'{cell:02d}'
    dfe=pd.read_csv(eSOH_DIR+"aging_param_cell_"+cell_no+".csv")
    dfe_0=pd.read_csv(eSOH_DIR+"aging_param_cell_"+cell_no+".csv")
    dfo_0=pd.read_csv(oCV_DIR+"ocv_data_cell_"+cell_no+".csv")
    # if cell_no=='13':
    #     dfo_d=dfo_0[dfo_0['N']==dfe['N'].iloc[-5]]
    #     dfo_0=dfo_0.drop(dfo_d.index.values)
    #     dfo_0=dfo_0.reset_index(drop=True)
    #     dfe = dfe.drop(dfe.index[-5])
    #     dfe = dfe.reset_index(drop=True)
    # Remove first RPT
    dfe = dfe.drop(dfe.index[0])
    dfe = dfe.reset_index(drop=True)
    # dfo_d=dfo_0[dfo_0['N']==0]
    # dfo_0=dfo_0.drop(dfo_d.index.values)
    if cell_no=='13':
        dfe = dfe.drop(dfe.index[-1])
        dfe = dfe.reset_index(drop=True)
        dfe_0 = dfe_0.drop(dfe_0.index[-1])
        dfe_0 = dfe_0.reset_index(drop=True)
    dfe['N']=dfe['N']-dfe['N'][0]
    dfe['Ah_th']=dfe['Ah_th']-dfe['Ah_th'][0]
    cycles = np.array(dfe['N'].astype('int'))
    cycles=cycles-1
    cycles[0]=cycles[0]+1
    dfe['N_mod'] = cycles
    N =dfe.N.unique()
    N_0 = dfe_0.N.unique()
    # print("Cycle Numbers:")
    # print(*N, sep = ", ") 
    # print(len(N_0))
    # print(len(dfo_0))
    rev_exp = []
    irrev_exp = []
    dfe["Capacity [A.h]"]=dfe["Cap"]
    dfe["n_Li"]= 3600/param.F.value*(dfe["y_100"]*dfe["C_p"]+dfe["x_100"]*dfe["C_n"])
    # dfe["Total lithium in particles [mol]"]= 3600/param.F.value*(dfe["y_100"]*dfe["C_p"]+dfe["x_100"]*dfe["C_n"])
    dfe["Loss of lithium inventory [%]"]=(1-dfe["n_Li"]/dfe["n_Li"][0])*100
    dfe["Capacity Retention [%]"]=(dfe["Capacity [A.h]"]/dfe["Capacity [A.h]"][0])*100
    dfe_0["Capacity [A.h]"]=dfe_0["Cap"]
    dfe_0["n_Li"]= 3600/param.F.value*(dfe_0["y_100"]*dfe_0["C_p"]+dfe_0["x_100"]*dfe_0["C_n"])
    # dfe_0["Total lithium in particles [mol]"]= 3600/param.F.value*(dfe_0["y_100"]*dfe_0["C_p"]+dfe_0["x_100"]*dfe_0["C_n"])
    dfe_0["Loss of lithium inventory [%]"]=(1-dfe_0["n_Li"]/dfe_0["n_Li"][0])*100
    dfe_0["Capacity Retention [%]"]=(dfe_0["Capacity [A.h]"]/dfe_0["Capacity [A.h]"][0])*100

    for i in range(len(N_0)-1):
        # print(i)
        dfo = dfo_0[dfo_0['N']==N_0[i+1]]
        # print(max(dfo['E'])-min(dfo['E']))
        rev_exp.append(max(dfo['E'])-min(dfo['E']))
    dfe['rev_exp']=rev_exp
    rev_exp = []
    for i in range(len(N_0)):
        # print(i)
        dfo = dfo_0[dfo_0['N']==N_0[i]]
        # print(max(dfo['E'])-min(dfo['E']))
        rev_exp.append(max(dfo['E'])-min(dfo['E']))
    dfe_0['rev_exp']=rev_exp

    dfo_1 = dfo_0[dfo_0['N']==N_0[1]]
    for i in range(len(N_0)-1):
        # print(i)
        dfo = dfo_0[dfo_0['N']==N_0[i+1]]
        # print(max(dfo['E'])-min(dfo['E']))
        irrev_exp.append(min(dfo['E'])-min(dfo_1['E']))
    # if cell == 12:
    #     for ii in range(1,len(irrev_exp)):
    #         irrev_exp[ii]=irrev_exp[ii]+36
    dfe['irrev_exp']=irrev_exp
    irrev_exp = []
    dfo_1 = dfo_0[dfo_0['N']==N_0[0]]
    for i in range(len(N_0)):
        # print(i)
        dfo = dfo_0[dfo_0['N']==N_0[i]]
        # print(max(dfo['E'])-min(dfo['E']))
        irrev_exp.append(min(dfo['E'])-min(dfo_1['E']))
    # if cell == 12:
    #     for ii in range(2,len(irrev_exp)):
    #         irrev_exp[ii]=irrev_exp[ii]+36
    dfe_0['irrev_exp']=irrev_exp
    # print('Irreversible Expansion')
    return cell_no,dfe,dfe_0,dfo_0,N,N_0

def init_exp(cell_no,dfe,spm,parameter_values):
    # dfe_0 = dfe[dfe['N']==N[0]]
    param = spm.param
    C_n_init = dfe['C_n'][0]
    C_p_init = dfe['C_p'][0]
    y_0_init = dfe['y_0'][0] 
    eps_n_data = parameter_values.evaluate(C_n_init*3600/(param.n.L * param.n.prim.c_max * param.F* param.A_cc))
    eps_p_data = parameter_values.evaluate(C_p_init*3600/(param.p.L * param.p.prim.c_max * param.F* param.A_cc))
    # cs_p_init = parameter_values.evaluate(y_0_init* param.c_p_max) 
    if (int(cell_no)-1)//3 ==0:
        c_rate_c = 'C/5'
        c_rate_d = 'C/5'
        dis_set = " until 3V"
    elif (int(cell_no)-1)//3 ==1:
        c_rate_c = '1.5C'
        c_rate_d = '1.5C'
        dis_set = " until 3V"
    elif (int(cell_no)-1)//3 ==2:
        c_rate_c = '2C'
        c_rate_d = '2C'
        dis_set = " until 3V"
    elif (int(cell_no)-1)//3 ==3:
        c_rate_c = 'C/5'
        c_rate_d = '1.5C'
        dis_set = " until 3V"
    elif (int(cell_no)-1)//3 ==4:
        c_rate_c = 'C/5'
        c_rate_d = 'C/5'
        dis_set = " for 150 min"
    elif (int(cell_no)-1)//3 ==5:
        c_rate_c = 'C/5'
        c_rate_d = '1.5C'
        dis_set = " for 20 min"
    elif (int(cell_no)-1)//3 ==6:
        c_rate_c = '1.5C'
        c_rate_d = '1.5C'
        dis_set = " for 20 min"
    if int(cell_no)%3 == 0:
        Temp = 45
    if int(cell_no)%3 == 1:
        Temp = 25
    if int(cell_no)%3 == 2:
        Temp = -5
    SOC_0 = 1
    return eps_n_data,eps_p_data,c_rate_c,c_rate_d,dis_set,Temp,SOC_0

def get_pulse_res(spm,parameter_values,esoh_sol,t_in,SOC):
    param=spm.param
    c_n_max = parameter_values.evaluate(param.n.prim.c_max)
    c_p_max = parameter_values.evaluate(param.p.prim.c_max)
    x_100 = esoh_sol["x_100"].data[0]
    y_100 = esoh_sol["y_100"].data[0]
    x_0 = esoh_sol["x_0"].data[0]
    y_0 = esoh_sol["y_0"].data[0]
    cs_n_0 = (SOC*(x_100-x_0)+x_0)*c_n_max
    cs_p_0 = (SOC*(y_100-y_0)+y_0)*c_p_max
    parameter_values.update(
      {
          "Initial concentration in negative electrode [mol.m-3]": cs_n_0,
          "Initial concentration in positive electrode [mol.m-3]": cs_p_0,        
      }
    )
    sim_pulse = pybamm.Simulation(spm, parameter_values=parameter_values, 
                            solver=pybamm.CasadiSolver(mode="safe", rtol=1e-6, atol=1e-6,dt_max=0.1))
    sol_pulse = sim_pulse.solve(t_eval=t_in)
    I   =  sol_pulse["Current [A]"].entries
    Vt  =  sol_pulse["Terminal voltage [V]"].entries
    idx = np.where(np.diff(np.sign(-I)))[0]
    Rs = abs((Vt[idx+1]-Vt[idx])/(I[idx+1]-I[idx]))[0]
    return Rs

def get_Rs(cyc_no,eSOH,parameter_values,Ns,spm):
  param=spm.param
  model = spm
  Vmin = 3.0
  Vmax = 4.2
  esoh_model = pybamm.lithium_ion.ElectrodeSOH()
  esoh_sim = pybamm.Simulation(esoh_model, parameter_values=parameter_values)
  Cn = eSOH["C_n"][Ns[cyc_no]]
  Cp = eSOH["C_p"][Ns[cyc_no]]
  c_n_max = parameter_values.evaluate(param.n.prim.c_max)
  c_p_max = parameter_values.evaluate(param.p.prim.c_max)
  n_Li_init = eSOH["Total lithium in particles [mol]"][Ns[cyc_no]]
  c_plated_Li = eSOH['X-averaged lithium plating concentration [mol.m-3]'][Ns[cyc_no]]
  eps_n_data = parameter_values.evaluate(Cn*3600/(param.n.L * param.n.prim.c_max * param.F* param.A_cc))
  eps_p_data = parameter_values.evaluate(Cp*3600/(param.p.L * param.p.prim.c_max * param.F* param.A_cc))
  del_sei = eSOH['X-averaged SEI thickness [m]'][Ns[cyc_no]]
  esoh_sol = esoh_sim.solve(
      [0],
      inputs={"V_min": Vmin, "V_max": Vmax, "C_n": Cn, "C_p": Cp, "n_Li": n_Li_init},
      solver=pybamm.AlgebraicSolver(),
  )
  parameter_values.update(
          {
              "Negative electrode active material volume fraction": eps_n_data,
              "Positive electrode active material volume fraction": eps_p_data,
              "Initial inner SEI thickness [m]": 0e-09,
              "Initial outer SEI thickness [m]": del_sei,
              "Initial plated lithium concentration [mol.m-3]": c_plated_Li,
          }
        )
  t_in0 = 1
  t_in1 = 1
  t_inf = t_in0+t_in1
  t_in = np.arange(0,t_inf,0.1)
  I_in = []
  for tt in t_in:
    if tt<t_in0:
        I_in = np.append(I_in,0)
    elif tt>=t_in0 and tt<t_in0+t_in1:
        I_in = np.append(I_in,5)
  timescale = parameter_values.evaluate(spm.timescale)
  current_interpolant = pybamm.Interpolant(
    t_in, -I_in, timescale * pybamm.t
  )
  parameter_values["Current function [A]"] = current_interpolant
  SOC_vals = np.linspace(1,0,11)
  Rs_ch_s = []
  for SOC in SOC_vals[1:10]:
      Rs_t = get_pulse_res(spm,parameter_values,esoh_sol,t_in,SOC)
      Rs_ch_s.append(Rs_t)
  Rs_ch = np.average(Rs_ch_s)
  
  current_interpolant = pybamm.Interpolant(
    t_in, I_in, timescale * pybamm.t
  )
  parameter_values["Current function [A]"] = current_interpolant
  SOC_vals = np.linspace(1,0,11)
  Rs_dh_s = []
  for SOC in SOC_vals[1:10]:
      Rs_t = get_pulse_res(spm,parameter_values,esoh_sol,t_in,SOC)
      Rs_dh_s.append(Rs_t)
  Rs_dh = np.average(Rs_dh_s)
  Rs = (Rs_dh + Rs_ch)/2
  Rs_ave_s = (np.array(Rs_dh_s)+np.array(Rs_ch_s))/2
  return Rs,Rs_ave_s

def load_cycling_data_ch(cell,eSOH_DIR,oCV_DIR,cyc_DIR,cyc_no):
    cell_no,dfe,dfe_0,dfo_0,N,N_0 = load_data(cell,eSOH_DIR,oCV_DIR)
    cycles = np.array(dfe_0['N'].astype('int'))+1
    cycles = cycles[1:]
    # print(cell_no)
    cyc_data_raw1 = pd.read_csv(cyc_DIR+'cycling_data_cell_'+cell_no+'.csv')
    if cell == 1:
        offset = 12
    else:
        offset = 0
    if len(cycles) == cyc_no+1:
        N1 = cyc_data_raw1["Cycle number"].iloc[-1]-offset
    else:
        N1 = cycles[cyc_no]
    # print(N1)
    cyc_data_raw = cyc_data_raw1[ cyc_data_raw1['Cycle number'] == N1 ]
    cyc_data = cyc_data_raw.reset_index(drop=True)
    t_c1 = cyc_data['Time [s]']-cyc_data['Time [s]'][0]
    t_c1 = t_c1.values
    I_c1 = cyc_data['Current [mA]']/1000
    I_c1 = I_c1.values
    V_c1 = cyc_data['Voltage [V]']
    V_c1 = V_c1.values
    E_c1 = cyc_data["Expansion [mu m]"]
    E_c1 = E_c1.values
    idx_I = np.where(np.sign(I_c1[:-1]) != np.sign(I_c1[1:]))[0] 
    idx_I = idx_I[idx_I>50]
    t = t_c1[:idx_I[0]]
    V = V_c1[:idx_I[0]]
    I = I_c1[:idx_I[0]]
    E = E_c1[:idx_I[0]]-E_c1[0]
    # t = t_c1
    # V = V_c1
    # I = I_c1
    Q = integrate.cumtrapz(I,t, initial=0)/3600 #Ah
    
    return t,V,I,Q,E

def cyc_comp_ch(cyc_no,eSOH,t_d,Q_d,V_d,E_d,parameter_values,spm,Ns,c_rate_c,c_rate_d):
    # dfo = dfo_0[dfo_0['N']==N[cyc_no]]
    model = spm
    Vmin = 3.0
    Vmax = 4.2
    esoh_model = pybamm.lithium_ion.ElectrodeSOH()
    esoh_sim = pybamm.Simulation(esoh_model, parameter_values=parameter_values)
    param = model.param
    Cn = eSOH["C_n"][Ns[cyc_no]]
    # print(Cn)
    Cp = eSOH["C_p"][Ns[cyc_no]]
    c_n_max = parameter_values.evaluate(param.n.prim.c_max)
    c_p_max = parameter_values.evaluate(param.p.prim.c_max)
    n_Li_init = eSOH["Total lithium in particles [mol]"][Ns[cyc_no]]
    eps_n_data = parameter_values.evaluate(Cn*3600/(param.n.L * param.n.prim.c_max * param.F* param.A_cc))
    eps_p_data = parameter_values.evaluate(Cp*3600/(param.p.L * param.p.prim.c_max * param.F* param.A_cc))
    del_sei = eSOH['X-averaged SEI thickness [m]'][Ns[cyc_no]]
    c_plated_Li = eSOH['X-averaged lithium plating concentration [mol.m-3]'][Ns[cyc_no]]
    esoh_sol = esoh_sim.solve(
        [0],
        inputs={"V_min": Vmin, "V_max": Vmax, "C_n": Cn, "C_p": Cp, "n_Li": n_Li_init},
        solver=pybamm.AlgebraicSolver(),
    )

    parameter_values.update(
        {
            "Initial concentration in negative electrode [mol.m-3]": esoh_sol[
       "x_100"
            ].data[0]
            * c_n_max,
            "Initial concentration in positive electrode [mol.m-3]": esoh_sol[
                "y_100"
            ].data[0]
            * c_p_max,
            "Negative electrode active material volume fraction": eps_n_data,
            "Positive electrode active material volume fraction": eps_p_data,
            "Initial temperature [K]": 273.15+25,
            "Ambient temperature [K]": 273.15+25,
            "Initial inner SEI thickness [m]": 0e-09,
            # "Initial outer SEI thickness [m]": 5e-09,
            "Initial outer SEI thickness [m]": del_sei,
            "Initial plated lithium concentration [mol.m-3]": c_plated_Li,        
        }
    )
    dis_set = " until 3V"

    if c_rate_d=="C/5":
        timestep = '10 sec'
    else:
        timestep = '1 sec'

    experiment_cyc_comp_ch = pybamm.Experiment(
        [
            "Discharge at "+c_rate_d+dis_set,
            "Rest for 10 sec",
            "Charge at "+c_rate_c+" until 4.2V", 
            "Hold at 4.2V until C/100",
            # "Rest for 10 sec",
            # "Discharge at "+c_rate_d+dis_set,
        ],
        period=timestep,
    )
    sim_exp = pybamm.Simulation(
        model, experiment=experiment_cyc_comp_ch, parameter_values=parameter_values,
        solver=pybamm.CasadiSolver(mode="safe", rtol=1e-6, atol=1e-6,dt_max=0.1),
    )
    sol_exp = sim_exp.solve()
    t_t = sol_exp["Time [s]"].entries
    I_t = sol_exp["Current [A]"].entries
    Q_t = -sol_exp['Discharge capacity [A.h]'].entries
    Vt_t = sol_exp["Terminal voltage [V]"].entries
    exp_t = 30e6*sol_exp["Cell thickness change [m]"].entries
    idx = np.where(np.diff(np.sign(-I_t)))[0]
    I = I_t[idx[0]:]
    t = t_t[idx[0]:]-t_t[idx[0]]
    Q = Q_t[idx[0]:]-Q_t[idx[0]]
    Vt = Vt_t[idx[0]:]
    Exp = exp_t[idx[0]:]-exp_t[idx[0]]

    if max(t)<max(t_d):
        int_V = interpolate.CubicSpline(t_d,V_d,extrapolate=True)
        rmse_V = pybamm.rmse(Vt,int_V(t))
        int_E = interpolate.CubicSpline(t_d,E_d,extrapolate=True)
        rmse_E = pybamm.rmse(Exp,int_E(t))
        # int_VQ = interpolate.CubicSpline(Q_d,V_d,extrapolate=True)
        # rmse_VQ = pybamm.rmse(Vt,int_VQ(Q))
        # int_EQ = interpolate.CubicSpline(Q_d,E_d,extrapolate=True)
        # rmse_EQ = pybamm.rmse(Exp,int_EQ(Q))
    else:
        int_V = interpolate.CubicSpline(t,Vt,extrapolate=True)
        rmse_V = pybamm.rmse(V_d,int_V(t_d))
        int_E = interpolate.CubicSpline(t,Exp,extrapolate=True)
        rmse_E = pybamm.rmse(E_d,int_E(t_d))
        # int_VQ = interpolate.CubicSpline(Q[1:],Vt[1:],extrapolate=True)
        # rmse_VQ = pybamm.rmse(V_d,int_VQ(Q_d))
        # int_EQ = interpolate.CubicSpline(Q[1:],Exp[1:],extrapolate=True)
        # rmse_EQ = pybamm.rmse(E_d,int_EQ(Q_d))
    # rmse_V =0
    # max_V = 0
    rmse_VQ  = 0 ; rmse_EQ = 0
    return t,I,Q,Vt,Exp,sol_exp,rmse_V,rmse_E,rmse_VQ,rmse_EQ

def get_rmse(t_d,V_d,E_d,t,Vt,Exp):
    if max(t)<max(t_d):
        int_V = interpolate.CubicSpline(t_d,V_d,extrapolate=True)
        rmse_V = pybamm.rmse(Vt,int_V(t))
        int_E = interpolate.CubicSpline(t_d,E_d,extrapolate=True)
        rmse_E = pybamm.rmse(Exp,int_E(t))
        # int_VQ = interpolate.CubicSpline(Q_d,V_d,extrapolate=True)
        # rmse_VQ = pybamm.rmse(Vt,int_VQ(Q))
        # int_EQ = interpolate.CubicSpline(Q_d,E_d,extrapolate=True)
        # rmse_EQ = pybamm.rmse(Exp,int_EQ(Q))
    else:
        int_V = interpolate.CubicSpline(t,Vt,extrapolate=True)
        rmse_V = pybamm.rmse(V_d,int_V(t_d))
        int_E = interpolate.CubicSpline(t,Exp,extrapolate=True)
        rmse_E = pybamm.rmse(E_d,int_E(t_d))
        # int_VQ = interpolate.CubicSpline(Q,Vt,extrapolate=True)
        # rmse_VQ = pybamm.rmse(V_d,int_VQ(Q_d))
        # int_EQ = interpolate.CubicSpline(Q,Exp,extrapolate=True)
        # rmse_EQ = pybamm.rmse(E_d,int_EQ(Q_d))
    return rmse_V,rmse_E

def load_cycling_data(cell,eSOH_DIR,oCV_DIR,cyc_DIR,cyc_no):
    cell_no,dfe,dfe_0,dfo_0,N,N_0 = load_data(cell,eSOH_DIR,oCV_DIR)
    cycles = np.array(dfe_0['N'].astype('int'))+1
    cycles = cycles[1:]
    # print(cell_no)
    cyc_data_raw1 = pd.read_csv(cyc_DIR+'cycling_data_cell_'+cell_no+'.csv')
    if len(cycles) == cyc_no+1:
        N1 = cyc_data_raw1["Cycle number"].iloc[-2]-2
    else:
        N1 = cycles[cyc_no]
    # print(N1)
    cyc_data_raw = cyc_data_raw1[ cyc_data_raw1['Cycle number'] == N1 ]
    cyc_data = cyc_data_raw.reset_index(drop=True)
    t_c1 = cyc_data['Time [s]']-cyc_data['Time [s]'][0]
    t_c1 = t_c1.values
    I_c1 = cyc_data['Current [mA]']/1000
    I_c1 = I_c1.values
    V_c1 = cyc_data['Voltage [V]']
    V_c1 = V_c1.values
    E_c1 = cyc_data["Expansion [mu m]"]
    E_c1 = E_c1.values
    idx_I = np.where(np.sign(I_c1[:-1]) != np.sign(I_c1[1:]))[0] 
    idx_I = idx_I[idx_I>50]
    t = t_c1
    V = V_c1
    I = I_c1
    E = E_c1-E_c1[0]
    # t = t_c1
    # V = V_c1
    # I = I_c1
    Q = integrate.cumtrapz(I,t, initial=0)/3600 #Ah
    
    return t,V,I,Q,E

def cyc_comp_dh_ch(cyc_no,eSOH,parameter_values,spm,Ns,c_rate_c,c_rate_d):
    # dfo = dfo_0[dfo_0['N']==N[cyc_no]]
    model = spm
    Vmin = 3.0
    Vmax = 4.2
    esoh_model = pybamm.lithium_ion.ElectrodeSOH()
    esoh_sim = pybamm.Simulation(esoh_model, parameter_values=parameter_values)
    param = model.param
    Cn = eSOH["C_n"][Ns[cyc_no]]
    # print(Cn)
    Cp = eSOH["C_p"][Ns[cyc_no]]
    c_n_max = parameter_values.evaluate(param.n.prim.c_max)
    c_p_max = parameter_values.evaluate(param.p.prim.c_max)
    n_Li_init = eSOH["Total lithium in particles [mol]"][Ns[cyc_no]]
    eps_n_data = parameter_values.evaluate(Cn*3600/(param.n.L * param.n.prim.c_max * param.F* param.A_cc))
    eps_p_data = parameter_values.evaluate(Cp*3600/(param.p.L * param.p.prim.c_max * param.F* param.A_cc))
    del_sei = eSOH['X-averaged SEI thickness [m]'][Ns[cyc_no]]
    c_plated_Li = eSOH['X-averaged lithium plating concentration [mol.m-3]'][Ns[cyc_no]]
    esoh_sol = esoh_sim.solve(
        [0],
        inputs={"V_min": Vmin, "V_max": Vmax, "C_n": Cn, "C_p": Cp, "n_Li": n_Li_init},
        solver=pybamm.AlgebraicSolver(),
    )

    parameter_values.update(
        {
            "Initial concentration in negative electrode [mol.m-3]": esoh_sol[
       "x_100"
            ].data[0]
            * c_n_max,
            "Initial concentration in positive electrode [mol.m-3]": esoh_sol[
                "y_100"
            ].data[0]
            * c_p_max,
            "Negative electrode active material volume fraction": eps_n_data,
            "Positive electrode active material volume fraction": eps_p_data,
            "Initial temperature [K]": 273.15+25,
            "Ambient temperature [K]": 273.15+25,
            "Initial inner SEI thickness [m]": 0e-09,
            # "Initial outer SEI thickness [m]": 5e-09,
            "Initial outer SEI thickness [m]": del_sei,
            "Initial plated lithium concentration [mol.m-3]": c_plated_Li,        
        }
    )
    dis_set = " until 3V"

    if c_rate_d=="C/5":
        timestep = '1 sec'
    else:
        timestep = '1 sec'

    experiment_cyc_comp_ch = pybamm.Experiment(
        [
            "Discharge at "+c_rate_d+dis_set,
            "Rest for 10 sec",
            "Charge at "+c_rate_c+" until 4.2V", 
            "Hold at 4.2V until C/100",
            # "Rest for 10 sec",
            # "Discharge at "+c_rate_d+dis_set,
        ],
        period=timestep,
    )
    sim_exp = pybamm.Simulation(
        model, experiment=experiment_cyc_comp_ch, parameter_values=parameter_values,
        solver=pybamm.CasadiSolver(mode="safe", rtol=1e-6, atol=1e-6,dt_max=0.1),
    )
    sol_exp = sim_exp.solve()
    t_t = sol_exp["Time [s]"].entries
    I_t = sol_exp["Current [A]"].entries
    Q_t = -sol_exp['Discharge capacity [A.h]'].entries
    Vt_t = sol_exp["Terminal voltage [V]"].entries
    exp_t = 30e6*sol_exp["Cell thickness change [m]"].entries
    idx = np.where(np.diff(np.sign(-I_t)))[0]
    I = I_t
    t = t_t-t_t[0]
    Q = Q_t-Q_t[0]
    Vt = Vt_t
    Exp = exp_t-exp_t[0]

    return t,I,Q,Vt,Exp,sol_exp

def cyc_comp_ch_dh(cyc_no,eSOH,t_d,Q_d,V_d,E_d,parameter_values,spm,Ns,c_rate_c,c_rate_d):
    # dfo = dfo_0[dfo_0['N']==N[cyc_no]]
    model = spm
    Vmin = 3.0
    Vmax = 4.2
    esoh_model = pybamm.lithium_ion.ElectrodeSOH()
    esoh_sim = pybamm.Simulation(esoh_model, parameter_values=parameter_values)
    param = model.param
    Cn = eSOH["C_n"][Ns[cyc_no]]
    # print(Cn)
    Cp = eSOH["C_p"][Ns[cyc_no]]
    c_n_max = parameter_values.evaluate(param.n.prim.c_max)
    c_p_max = parameter_values.evaluate(param.p.prim.c_max)
    n_Li_init = eSOH["Total lithium in particles [mol]"][Ns[cyc_no]]
    eps_n_data = parameter_values.evaluate(Cn*3600/(param.n.L * param.n.prim.c_max * param.F* param.A_cc))
    eps_p_data = parameter_values.evaluate(Cp*3600/(param.p.L * param.p.prim.c_max * param.F* param.A_cc))
    del_sei = eSOH['X-averaged SEI thickness [m]'][Ns[cyc_no]]
    c_plated_Li = eSOH['X-averaged lithium plating concentration [mol.m-3]'][Ns[cyc_no]]
    esoh_sol = esoh_sim.solve(
        [0],
        inputs={"V_min": Vmin, "V_max": Vmax, "C_n": Cn, "C_p": Cp, "n_Li": n_Li_init},
        solver=pybamm.AlgebraicSolver(),
    )

    parameter_values.update(
        {
            "Initial concentration in negative electrode [mol.m-3]": esoh_sol[
       "x_100"
            ].data[0]
            * c_n_max,
            "Initial concentration in positive electrode [mol.m-3]": esoh_sol[
                "y_100"
            ].data[0]
            * c_p_max,
            "Negative electrode active material volume fraction": eps_n_data,
            "Positive electrode active material volume fraction": eps_p_data,
            "Initial temperature [K]": 273.15+25,
            "Ambient temperature [K]": 273.15+25,
            "Initial inner SEI thickness [m]": 0e-09,
            # "Initial outer SEI thickness [m]": 5e-09,
            "Initial outer SEI thickness [m]": del_sei,
            "Initial plated lithium concentration [mol.m-3]": c_plated_Li,        
        }
    )
    dis_set = " until 3V"

    if c_rate_d=="C/5":
        timestep = '10 sec'
    else:
        timestep = '1 sec'

    experiment_cyc_comp_ch = pybamm.Experiment(
        [
            "Discharge at "+c_rate_d+dis_set,
            "Rest for 10 sec",
            "Charge at "+c_rate_c+" until 4.2V", 
            "Hold at 4.2V until C/100",
            "Rest for 10 sec",
            "Discharge at "+c_rate_d+dis_set,
        ],
        period=timestep,
    )
    sim_exp = pybamm.Simulation(
        model, experiment=experiment_cyc_comp_ch, parameter_values=parameter_values,
        solver=pybamm.CasadiSolver(mode="safe", rtol=1e-6, atol=1e-6,dt_max=0.1),
    )
    sol_exp = sim_exp.solve()
    t_t = sol_exp["Time [s]"].entries
    I_t = sol_exp["Current [A]"].entries
    Q_t = -sol_exp['Discharge capacity [A.h]'].entries
    Vt_t = sol_exp["Terminal voltage [V]"].entries
    exp_t = 30e6*sol_exp["Cell thickness change [m]"].entries
    idx = np.where(np.diff(np.sign(-I_t)))[0]
    I = I_t[idx[0]:]
    t = t_t[idx[0]:]-t_t[idx[0]]
    Q = Q_t[idx[0]:]-Q_t[idx[0]]
    Vt = Vt_t[idx[0]:]
    Exp = exp_t[idx[0]:]-exp_t[idx[0]]

    if max(t)<max(t_d):
        int_V = interpolate.CubicSpline(t_d,V_d,extrapolate=True)
        rmse_V = pybamm.rmse(Vt,int_V(t))
        int_E = interpolate.CubicSpline(t_d,E_d,extrapolate=True)
        rmse_E = pybamm.rmse(Exp,int_E(t))
        # int_VQ = interpolate.CubicSpline(Q_d,V_d,extrapolate=True)
        # rmse_VQ = pybamm.rmse(Vt,int_VQ(Q))
        # int_EQ = interpolate.CubicSpline(Q_d,E_d,extrapolate=True)
        # rmse_EQ = pybamm.rmse(Exp,int_EQ(Q))
    else:
        int_V = interpolate.CubicSpline(t,Vt,extrapolate=True)
        rmse_V = pybamm.rmse(V_d,int_V(t_d))
        int_E = interpolate.CubicSpline(t,Exp,extrapolate=True)
        rmse_E = pybamm.rmse(E_d,int_E(t_d))
        # int_VQ = interpolate.CubicSpline(Q[1:],Vt[1:],extrapolate=True)
        # rmse_VQ = pybamm.rmse(V_d,int_VQ(Q_d))
        # int_EQ = interpolate.CubicSpline(Q[1:],Exp[1:],extrapolate=True)
        # rmse_EQ = pybamm.rmse(E_d,int_EQ(Q_d))
    # rmse_V =0
    # max_V = 0
    rmse_VQ  = 0 ; rmse_EQ = 0
    return t,I,Q,Vt,Exp,sol_exp,rmse_V,rmse_E,rmse_VQ,rmse_EQ