import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# --- Page Config ---
st.set_page_config(page_title="EIT 3-Level Lambda System", layout="wide")

st.title("Electromagnetically Induced Transparency (EIT)")
st.markdown(r"""
This app simulates the absorption profile of a 3-level $\Lambda$ system. 
A strong **control beam** ($\Omega_c$) creates a transparency window for a weak **probe beam** ($\Omega_p$).
""")

# --- Sidebar Inputs ---
st.sidebar.header("System Parameters")
Gamma = st.sidebar.number_input(r"Decay Rate Γ (MHz)", value=5.0)
Is = st.sidebar.number_input(r"Saturation Intensity Is (mW/cm²)", value=1.0)
gamma21 = st.sidebar.slider(r"Dephasing Rate γ21 (MHz)", 0.0, 1.0, 0.05, step=0.01)
L_alpha0 = st.sidebar.slider("Peak Optical Depth (Lα₀)", 1.0, 20.0, 10.0)

st.sidebar.header("Beam Controls")
Ic = st.sidebar.slider(r"Control Intensity Ic (mW/cm²)", 0.0, 100.0, 15.0)
detuning_c = st.sidebar.slider(r"Control Detuning Δc (MHz)", -10.0, 10.0, 0.0)
scan_range = st.sidebar.slider("Scan Range (MHz)", 10.0, 100.0, 30.0)

# --- Calculation Logic ---
def get_absorption(delta, Gamma, Ic, Is, detuning_c, gamma21, L_alpha0):
    # Rabi frequency: Omega = Gamma * sqrt(I / (2 * Is))
    Omega_c = Gamma * np.sqrt(Ic / (2 * Is))
    delta_p = delta + detuning_c
    
    # Susceptibility formula
    numerator = 1j * (Gamma / 2.0)
    denominator = (Gamma / 2.0 - 1j * delta_p) + (Omega_c**2 / 4.0) / (gamma21 - 1j * delta)
    
    chi = numerator / denominator
    # Normalize absorption
    norm_abs = np.imag(chi) / np.imag((1j * (Gamma/2.0)) / (Gamma/2.0))
    return L_alpha0 * norm_abs

# Generate Data
delta_range = np.linspace(-scan_range/2, scan_range/2, 1000)
absorption = get_absorption(delta_range, Gamma, Ic, Is, detuning_c, gamma21, L_alpha0)

# --- Plotting & Illustration ---
fig = plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2.5], hspace=0.3)

# 1. Diagram Axes
ax_diag = plt.subplot(gs[0])
ax_diag.set_xlim(0, 3)
ax_diag.set_ylim(-0.2, 3)
ax_diag.axis('off')

# Drawing the Lambda System
ax_diag.hlines([0.5, 0.5, 2.5], [0.5, 1.7, 1.1], [1.3, 2.5, 1.9], colors='black', lw=3)
ax_diag.text(0.9, 0.2, r'$|1\rangle$', ha='center', fontsize=12)
ax_diag.text(2.1, 0.2, r'$|2\rangle$', ha='center', fontsize=12)
ax_diag.text(1.5, 2.7, r'$|3\rangle$', ha='center', fontsize=12)

# Transitions
ax_diag.annotate("", xy=(1.4, 2.5), xytext=(0.9, 0.5), arrowprops=dict(arrowstyle="->", color="red", lw=1.5))
ax_diag.annotate("", xy=(1.6, 2.5), xytext=(2.1, 0.5), arrowprops=dict(arrowstyle="->", color="blue", lw=2))
ax_diag.text(0.7, 1.5, r'$\Omega_p, \Delta_p$', color="red", fontsize=11)
ax_diag.text(2.1, 1.5, r'$\Omega_c, \Delta_c$', color="blue", fontsize=11)
ax_diag.set_title(r"Energy Level Diagram ($\Lambda$-System)", fontsize=14)

# 2. Absorption Plot
ax_plot = plt.subplot(gs[1])
ax_plot.plot(delta_range, absorption, color='firebrick', lw=2, label=r'Absorption ($\alpha L$)')

# Labels and Styling
ax_plot.set_xlabel(r'Two-photon Detuning $\delta$ (MHz)', fontsize=12)
ax_plot.set_ylabel(r'Absorption ($\alpha L$)', fontsize=12)
ax_plot.grid(True, alpha=0.3)
ax_plot.axvline(0, color='black', linestyle='--', alpha=0.5)
ax_plot.legend(loc='upper right')

# Display in Streamlit
st.pyplot(fig)

# --- Explanation Table ---
st.subheader("Variable Definitions")
st.table({
    "Variable": ["Gamma (Γ)", "Is", "Ic", "Delta_c (Δc)", "Delta_p (Δp)", "delta (δ)", "gamma21 (γ21)"],
    "Description": [
        "Excited state spontaneous decay rate.",
        "Saturation intensity of the transition.",
        "Intensity of the strong control beam.",
        "Detuning of the control beam from the |2> to |3> transition.",
        "Detuning of the probe beam from the |1> to |3> transition.",
        "Two-photon detuning (Δp - Δc).",
        "Ground state decoherence rate (width of the EIT dip)."
    ]
})
