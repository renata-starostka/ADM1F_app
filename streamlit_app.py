import streamlit as st

st.title("Poo to Power (Anaerobic Digestion)")
st.write(
    "This app is designed to help wastewater operators of anaerobic digestion systems determine what they should add to their anaerobic digesters. It is a combination of the Anaerobic Digestion Model 1, originally in Matlab, converted to C++ with a petsc solver, and now fully in Python (still using petsc). As you can see, it is still in development. Some future versions will hopefully include:"
)
st.markdown('- Sliders and inputs to model your exact influent parameters, temperature, SRT')
st.markdown('- Default conditions for food waste and sewage sludge digestion')
st.markdown('- Outputs that give you biogas production, composition, and a stability score')
st.markdown('- An option to combine two or more different waste streams in the model')
st.markdown('- Later: A dynamic simulation to show when the biogas increase will occur if you add a new wastestream')
st.markdown('- Later: Calibration tools to adjust parameters to your exact wastewater plant')


import numpy as np
from adm1f import AppCtx

# Thresholds for stability criteria
criteria_min = {'Biogas': 55, 'pH': 6.1, 'Alk': 2000, 'NH3': 0, 'NH4': 0, 'VFA': 0, 'LCFA': 0}
criteria_max = {'Biogas': 100, 'pH': 8.3, 'Alk': 20000, 'NH3': 200, 'NH4': 3250, 'VFA': 5000, 'LCFA': 1400}

def stability(x, min_val, max_val):
    if x < 0.95 * min_val:
        return 0
    elif x < min_val:
        return 0.25
    elif x < 1.05 * min_val:
        return 0.5
    elif x < 1.1 * min_val:
        return 0.75
    elif x <= 0.9 * max_val:
        return 1
    elif x <= 0.95 * max_val:
        return 0.75
    elif x <= max_val:
        return 0.5
    else:
        return 0

def simulate_single_phase_AD(Q, Vliq, t_resx=0):
    """
    Simulate single-phase anaerobic digestion using Python ADM1 (petsc4py).

    Parameters:
    - Q: Influent flow rate (m3/d)
    - Vliq: Reactor volume (m3)
    - t_resx: Optional SRT - HRT difference (d)

    Returns:
    - Dictionary with simulation results and stability scores.
    """
    # Set up ADM1 context
    ctx = AppCtx()

    # Load input files (ensure these files exist in the working dir)
    ctx.LoadInfluent('influent.dat')
    ctx.LoadParams('params_SRT.dat')
    ctx.LoadIC('ic.dat')

    # Update volume, flow rate, and t_resx
    ctx.Vliq = Vliq
    ctx.Q = Q
    ctx.t_resx = t_resx

    # Initialize PETSc solver
    ctx.CreateTS()
    ctx.SetUpTS()

    # Run simulation
    ctx.TSSolve()

    # Extract output values (adjust indices if needed)
    biogas = ctx.indicator[43]
    pH = ctx.indicator[26]
    alk = ctx.indicator[57]
    nh3 = ctx.indicator[58]
    nh4 = ctx.indicator[59]
    vfa = ctx.indicator[54]
    lcfa = ctx.indicator[2]

    # Compute stability scores
    scores = {
        'Biogas': stability(biogas, criteria_min['Biogas'], criteria_max['Biogas']),
        'pH': stability(pH, criteria_min['pH'], criteria_max['pH']),
        'Alk': stability(alk, criteria_min['Alk'], criteria_max['Alk']),
        'NH3': stability(nh3, criteria_min['NH3'], criteria_max['NH3']),
        'NH4': stability(nh4, criteria_min['NH4'], criteria_max['NH4']),
        'VFA': stability(vfa, criteria_min['VFA'], criteria_max['VFA']),
        'LCFA': stability(lcfa, criteria_min['LCFA'], criteria_max['LCFA']),
    }

    stability_score = np.mean(list(scores.values()))

    return {
        'reactor_volume (m3)': Vliq,
        'influent_flow (m3/d)': Q,
        'HRT (d)': Vliq / Q,
        'SRT (d)': Vliq / Q + t_resx,
        'biogas_volume (m3/d)': biogas,
        'composition': {'CH4 (partial pressure)': ctx.indicator[40], 'CO2': ctx.indicator[41]},
        'stability_score': stability_score,
        'parameter_scores': scores
    }

#import simulate_single_phase_AD  # Adjust import as needed

st.header("Single-Phase AD Simulation")

st.subheader("Enter Basic Reactor Information")

# User Inputs
Vliq = st.number_input("Reactor Volume (m³)", min_value=0.1, step=0.1, value=200.0)
Q = st.number_input("Influent Flow Rate (m³/day)", min_value=0.01, step=0.1, value=25.0)
t_resx = st.number_input("SRT-HRT difference (days)", step=1.0, value=0.0)

if st.button("Run Simulation"):
    with st.spinner("Running simulation..."):
        try:
            results = simulate_single_phase_AD(Q=Q, Vliq=Vliq, t_resx=t_resx)

            st.success("Simulation complete!")
            st.subheader("Results:")

            st.metric("Biogas Production", f"{results['biogas_volume (m3/d)']:.2f} m³/day")
            st.metric("CH₄ Composition", f"{results['composition']['CH4 (partial pressure)']:.2f}")
            st.metric("CO₂ Composition", f"{results['composition']['CO2']:.2f}")
            st.metric("HRT", f"{results['HRT (d)']:.1f} days")
            st.metric("SRT", f"{results['SRT (d)']:.1f} days")
            st.metric("Stability Score", f"{results['stability_score']:.2f} / 1.00")

            st.subheader("Individual Parameter Stability Scores:")
            for k, v in results["parameter_scores"].items():
                st.write(f"- {k}: {v:.2f}")
        
        except Exception as e:
            st.error(f"Simulation failed: {e}")
