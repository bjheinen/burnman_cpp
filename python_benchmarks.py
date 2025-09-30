import numpy as np
import timeit
import statistics
import csv
import burnman as bm

def make_mineral_benchmarks(test_mineral):
    return {
        "reset_cache": lambda: (test_mineral.reset(), 0)[1],
        "get_molar_volume": lambda: (test_mineral.reset(), test_mineral.molar_volume)[1],
        "get_density": lambda: (test_mineral.reset(), test_mineral.density)[1],
        "get_molar_internal_energy": lambda: (test_mineral.reset(), test_mineral.molar_internal_energy)[1],
        "get_molar_gibbs": lambda: (test_mineral.reset(), test_mineral.molar_gibbs)[1],
        "get_molar_helmholtz": lambda: (test_mineral.reset(), test_mineral.molar_helmholtz)[1],
        "get_molar_entropy": lambda: (test_mineral.reset(), test_mineral.molar_entropy)[1],
        "get_molar_enthalpy": lambda: (test_mineral.reset(), test_mineral.molar_enthalpy)[1],
        "get_isothermal_bulk_modulus_reuss": lambda: (test_mineral.reset(), test_mineral.isothermal_bulk_modulus_reuss)[1],
        "get_isentropic_bulk_modulus_reuss": lambda: (test_mineral.reset(), test_mineral.isentropic_bulk_modulus_reuss)[1],
        "get_isothermal_compressibility_reuss": lambda: (test_mineral.reset(), test_mineral.isothermal_compressibility_reuss)[1],
        "get_isentropic_compressibility_reuss": lambda: (test_mineral.reset(), test_mineral.isentropic_compressibility_reuss)[1],
        "get_shear_modulus": lambda: (test_mineral.reset(), test_mineral.shear_modulus)[1],
        "get_grueneisen_parameter": lambda: (test_mineral.reset(), test_mineral.grueneisen_parameter)[1],
        "get_thermal_expansivity": lambda: (test_mineral.reset(), test_mineral.thermal_expansivity)[1],
        "get_molar_heat_capacity_v": lambda: (test_mineral.reset(), test_mineral.molar_heat_capacity_v)[1],
        "get_molar_heat_capacity_p": lambda: (test_mineral.reset(), test_mineral.molar_heat_capacity_p)[1],
        "get_isentropic_thermal_gradient": lambda: (test_mineral.reset(), test_mineral.isentropic_thermal_gradient)[1],
        "get_p_wave_velocity": lambda: (test_mineral.reset(), test_mineral.p_wave_velocity)[1],
        "get_bulk_sound_velocity": lambda: (test_mineral.reset(), test_mineral.bulk_sound_velocity)[1],
        "get_shear_wave_velocity": lambda: (test_mineral.reset(), test_mineral.shear_wave_velocity)[1],
    }

def run_benchmarks(test_mineral, repeat=5, number=100):
    results = {}
    for name, func in make_mineral_benchmarks(test_mineral).items():
        print(f"Running benchmark: {name}")
        timings = timeit.repeat(func, number=number, repeat=repeat)
        timings_ns = [t / number * 1.0e9 for t in timings]
        mean_ns = statistics.mean(timings_ns)
        stdev_ns = statistics.stdev(timings_ns) if len(timings_ns) > 1 else 0.0
        results[name] = (mean_ns, stdev_ns)
    return results

def mineral_property_benchmarks():
    # Mineral per property benchmarks
    eos_type = 'mgd3'
    # Make params
    params = {}
    # Set-up common mineral params
    params['T_0'] = 300.0
    params['P_0'] = 0.0
    params['E_0'] = 0.0
    params['F_0'] = 0.0
    params['H_0'] = -1443030.0
    params['V_0'] = 11.24e-6
    params['K_0'] = 161.0e9
    params['Kprime_0'] = 3.8
    params['Kdprime_0'] = -1.6e-11
    params['G_0'] = 131.0e9
    params['Gprime_0'] = 2.1
    params['molar_mass'] = 0.0403
    params['n'] = 2
    params['Debye_0'] = 773.0
    params['grueneisen_0'] = 1.5
    params['q_0'] = 1.5
    params['Cp'] = [149.3, 0.002918, -2983000.0, -799.1]
    params['S_0'] = 62.6
    params['a_0'] = 1.87e-05
    params['eta_s_0'] = 2.565
    params['bel_0'] = 0.00411
    params['gel'] = 1.476

    eos_tags = [
        'mt',
        'vinet',
        'bm2',
        'bm3',
        'mgd2',
        'mgd3',
        'slb2',
        'slb3',
        'slb3-conductive',
        'hp_tmt'
    ]

    eos_map = {
        'mt': 'MT',
        'vinet': 'Vinet',
        'bm2': 'BM2',
        'bm3': 'BM3',
        'mgd2': 'MGD2',
        'mgd3': 'MGD3',
        'slb2': 'SLB2',
        'slb3': 'SLB3',
        'slb3-conductive': 'SLB3Conductive',
        'hp_tmt': 'HPTMT'
    }

    all_results = {}
    for eos in eos_tags:
        print(f"Running benchmarks for {eos}")
        p = params.copy()
        p['equation_of_state'] = eos
        mineral = bm.Mineral(params=p)
        mineral.set_state(55.0e9, 1000.0)
        eos_results = run_benchmarks(mineral)
        key = eos_map[eos]
        for bench_name, (mean_ns, stdev_ns) in eos_results.items():
            if bench_name not in all_results:
                all_results[bench_name] = {}
            all_results[bench_name][key] = mean_ns
            all_results[bench_name][f"{key}_std"] = stdev_ns
    columns = ["Benchmark"]
    for eos in eos_tags:
        key = eos_map[eos]
        columns.extend([key, f"{key}_std"])
    with open("python_mineral_benchmarks.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        for bench_name, vals in all_results.items():
            row = [bench_name] + [vals.get(col, "") for col in columns[1:]]
            writer.writerow(row)

def clear_attrib(obj):
    obj.reset()
    cached_props = [
        "stoichiometric_matrix",
        "stoichiometric_array",
        "reaction_basis",
        "reaction_basis_as_strings",
        "n_reactions",
        "independent_element_indices",
        "dependent_element_indices",
        "reduced_stoichiometric_array",
        "compositional_null_basis",
        "endmember_formulae",
        "endmember_names",
        "endmembers_per_phase",
        "elements",
        "n_endmembers",
        "n_elements"
    ]
    for prop in cached_props:
        obj.__dict__.pop(prop, None)
    phases = getattr(obj, "phases", None)
    if phases:
        for phase in phases:
            clear_attrib(phase)
    endmembers = getattr(obj, "endmembers", None)
    if endmembers:
        for em in endmembers:
            em[0].reset()

def time_func(func, number=5, repeat=10000):
    timings = timeit.repeat(func, number=number, repeat=repeat)
    timings_ns = [t / number * 1e9 for t in timings]
    mean_ns = statistics.mean(timings_ns)
    stdev_ns = statistics.stdev(timings_ns) if len(timings_ns) > 1 else 0.0
    return mean_ns, stdev_ns

def assemblage_benchmarks():
    # MgPv
    p_mgpv = {}
    p_mgpv["name"] = "MgSiO3 perovskite"
    p_mgpv["formula"] = {"Mg": 1.0, "Si": 1.0, "O": 3.0}
    p_mgpv["n"] = 5
    p_mgpv["molar_mass"] = 0.1003887
    p_mgpv["F_0"] = -1368000.0
    p_mgpv["V_0"] = 2.4445e-05
    p_mgpv["K_0"] = 2.51e11
    p_mgpv["Kprime_0"] = 4.1
    p_mgpv["G_0"] = 1.73e11
    p_mgpv["Gprime_0"] = 1.7
    p_mgpv["Debye_0"] = 905.0
    p_mgpv["grueneisen_0"] = 1.57
    p_mgpv["q_0"] = 1.1
    p_mgpv["eta_s_0"] = 2.3
    p_mgpv["equation_of_state"] = 'slb3'
    mgpv = bm.Mineral(params=p_mgpv)

    # FePv
    p_fepv = {}
    p_fepv["name"] = "FeSiO3 perovskite"
    p_fepv["formula"] = {"Fe": 1.0, "Si": 1.0, "O": 3.0}
    p_fepv["n"] = 5
    p_fepv["molar_mass"] = 0.1319287
    p_fepv["F_0"] = -1043000.0
    p_fepv["V_0"] = 2.534e-05
    p_fepv["K_0"] = 2.72e11
    p_fepv["Kprime_0"] = 4.1
    p_fepv["G_0"] = 1.33e11
    p_fepv["Gprime_0"] = 1.4
    p_fepv["Debye_0"] = 871.0
    p_fepv["grueneisen_0"] = 1.57
    p_fepv["q_0"] = 1.1
    p_fepv["eta_s_0"] = 2.3
    p_fepv["equation_of_state"] = 'slb3'
    fepv = bm.Mineral(params=p_fepv)

    # AlPv
    p_alpv = {}
    p_alpv["name"] = "AlAlO3 perovskite"
    p_alpv["formula"] = {"Al": 2.0, "O": 3.0}
    p_alpv["n"] = 5
    p_alpv["molar_mass"] = 0.1019612
    p_alpv["F_0"] = -1533878.0
    p_alpv["V_0"] = 2.494e-05
    p_alpv["K_0"] = 2.58e11
    p_alpv["Kprime_0"] = 4.1
    p_alpv["G_0"] = 1.71e11
    p_alpv["Gprime_0"] = 1.5
    p_alpv["Debye_0"] = 886.0
    p_alpv["grueneisen_0"] = 1.57
    p_alpv["q_0"] = 1.1
    p_alpv["eta_s_0"] = 2.5
    p_alpv["equation_of_state"] = 'slb3'
    alpv = bm.Mineral(params=p_alpv)

    # Periclase
    p_per = {}
    p_per["name"] = "Periclase"
    p_per["formula"] = {"Mg": 1.0, "O": 1.0}
    p_per["n"] = 2
    p_per["molar_mass"] = 0.0403044
    p_per["F_0"] = -569444.6
    p_per["V_0"] = 1.1244e-05
    p_per["K_0"] = 1.613836e11
    p_per["Kprime_0"] = 3.84045
    p_per["G_0"] = 1.309e11
    p_per["Gprime_0"] = 2.1438
    p_per["Debye_0"] = 767.0977
    p_per["grueneisen_0"] = 1.36127
    p_per["q_0"] = 1.7217
    p_per["eta_s_0"] = 2.81765
    p_per["equation_of_state"] = 'slb3'
    per = bm.Mineral(params=p_per)

    # Wuestite
    p_wue = {}
    p_wue["name"] = "Wuestite"
    p_wue["formula"] = {"Fe": 1.0, "O": 1.0}
    p_wue["n"] = 2
    p_wue["molar_mass"] = 0.0718444
    p_wue["F_0"] = -242146.0
    p_wue["V_0"] = 1.2264e-05
    p_wue["K_0"] = 1.794442e11
    p_wue["Kprime_0"] = 4.9376
    p_wue["G_0"] = 59000000000.0
    p_wue["Gprime_0"] = 1.44673
    p_wue["Debye_0"] = 454.1592
    p_wue["grueneisen_0"] = 1.53047
    p_wue["q_0"] = 1.7217
    p_wue["eta_s_0"] = -0.05731
    p_wue["equation_of_state"] = 'slb3'
    pm_wue = [['linear', {'delta_V': 1.0, 'delta_S': 2.0, 'delta_E': 3.0}]]
    wue = bm.Mineral(params=p_wue, property_modifiers=pm_wue)

    # CaPv
    p_capv = {}
    p_capv["name"] = "Ca-perovskite"
    p_capv["formula"] = {"Ca": 1.0, "Si": 1.0, "O": 3.0}
    p_capv["n"] = 5
    p_capv["molar_mass"] = 0.1161617
    p_capv["F_0"] = -1463358.0
    p_capv["V_0"] = 2.745e-05
    p_capv["K_0"] = 2.36e11
    p_capv["Kprime_0"] = 3.9
    p_capv["G_0"] = 1.568315e11
    p_capv["Gprime_0"] = 2.22713
    p_capv["Debye_0"] = 795.779
    p_capv["grueneisen_0"] = 1.88839
    p_capv["q_0"] = 0.89769
    p_capv["eta_s_0"] = 1.28818
    p_capv["equation_of_state"] = 'slb3'
    capv = bm.Mineral(params=p_capv)

    # Bridgmanite Solution
    bdg_endmembers = [
    [mgpv, "[Mg][Si]O3"],
    [fepv, "[Fe][Si]O3"],
    [alpv, "[Al][Al]O3"]]
    bdg_sol = bm.classes.solutionmodel.IdealSolution(endmembers=bdg_endmembers)
    bdg = bm.Solution(solution_model=bdg_sol, molar_fractions=[0.88, 0.07, 0.05])
    bdg.name = "Bridgmanite"

    # Ferropericlase Solution
    fp_endmembers = [
    [per, "[Mg]O"],
    [wue, "[Fe]O"]]
    fp_sol = bm.classes.solutionmodel.SymmetricRegularSolution(endmembers=fp_endmembers, energy_interaction=[[13.0e3]])
    fp = bm.Solution(solution_model=fp_sol, molar_fractions=[0.9, 0.1])
    fp.name = "Ferropericlase"

    # Make assemblage
    # Bdg + Fp + CaPv
    pyr = bm.Composite(
    [bdg, fp, capv],
    fractions=[0.7, 0.2, 0.1],
    fraction_type="molar",
    name="bdg + fp + capv"
    )

    pyr.set_state(50.e9,2000.0)

    benchmarks = {
        "get_formula": "formula",
        "get_endmember_formulae": "endmember_formulae",
        "get_n_endmembers": "n_endmembers",
        "get_endmembers_per_phase": "endmembers_per_phase",
        "get_n_elements": "n_elements",
        "get_n_reactions": "n_reactions",
        "get_elements": "elements",
        "get_independent_element_indices": "independent_element_indices",
        "get_dependent_element_indices": "dependent_element_indices",
        "get_stoichiometric_matrix": "stoichiometric_matrix",
        "get_compositional_null_basis": "compositional_null_basis",
        "get_reaction_basis": "reaction_basis",
        "get_molar_internal_energy": "molar_internal_energy",
        "get_molar_gibbs": "molar_gibbs",
        "get_molar_helmholtz": "molar_helmholtz",
        "get_molar_mass": "molar_mass",
        "get_molar_volume": "molar_volume",
        "get_density": "density",
        "get_molar_entropy": "molar_entropy",
        "get_molar_enthalpy": "molar_enthalpy",
        "get_isothermal_bulk_modulus_reuss": "isothermal_bulk_modulus_reuss",
        "get_isentropic_bulk_modulus_reuss": "isentropic_bulk_modulus_reuss",
        "get_isothermal_compressibility_reuss": "isothermal_compressibility_reuss",
        "get_isentropic_compressibility_reuss": "isentropic_compressibility_reuss",
        "get_shear_modulus": "shear_modulus",
        "get_p_wave_velocity": "p_wave_velocity",
        "get_bulk_sound_velocity": "bulk_sound_velocity",
        "get_shear_wave_velocity": "shear_wave_velocity",
        "get_grueneisen_parameter": "grueneisen_parameter",
        "get_thermal_expansivity": "thermal_expansivity",
        "get_molar_heat_capacity_v": "molar_heat_capacity_v",
        "get_molar_heat_capacity_p": "molar_heat_capacity_p"
    }

    def reset_assemblage_benchmark():
        #clear_attrib(pyr)
        pyr.reset()
        return 0

    def set_state_assemblage_benchmark():
        #clear_attrib(pyr)
        pyr.reset()
        pyr.set_state(50.e9, 2000.0)
        return 0

    results = []
    repeat = 5
    number = 100
    print("Running benchmark: reset_cache")
    reset_mean, reset_std = time_func(reset_assemblage_benchmark, number=number, repeat=repeat)
    results.append(('clear_computed_properties', reset_mean, reset_std))
    print("Running benchmark: set_state")
    set_state_mean, set_state_std = time_func(set_state_assemblage_benchmark, number=number, repeat=repeat)
    results.append(('set_state', set_state_mean, set_state_std))
    for col_name, method_name in benchmarks.items():
        print(f"Running benchmark: {col_name}")
        def wrapped():
            #clear_attrib(pyr)
            pyr.reset()
            return getattr(pyr, method_name)
        mean_ns, stdev_ns = time_func(wrapped, number=number, repeat=repeat)
        results.append((col_name, mean_ns, stdev_ns))
    with open("python_assemblage_benchmarks.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Benchmark", "Mean", "Std_dev"])
        writer.writerows(results)


if __name__ == "__main__":
    print("Running Mineral Property Benchmarks...")
    mineral_property_benchmarks()
    print("Running Assemblage Benchmarks...")
    assemblage_benchmarks()
