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
