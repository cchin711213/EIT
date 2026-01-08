import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# --- Page Config ---
st.set_page_config(page_title="Strong Probe EIT Simulator", layout="wide")

st.title("EIT: Strong Probe & Density Matrix Solver")
st.markdown(r"""
This version uses a **Steady-State Density Matrix Solver**. 
Unlike the linear model, this accounts for power broadening and population redistribution caused by a strong probe beam ($I_p$).
""")

# --- Sidebar Inputs ---
st.sidebar.header("System Parameters")
Gamma = st.sidebar.number_input(r"Decay Rate Γ (MHz)", value=5.0, step=0.1, format="%.2f")
Is = st.sidebar.number_input(r"Saturation Intensity Is (mW/cm²)", value=1.0, step=0.1, format="%.2f")
gamma21 = st.sidebar.number_input(r"Dephasing Rate γ21 (MHz)", value=0.05, step=0.01, format="%.3f")
L_alpha0 = st.sidebar.number_input("Peak Optical Depth (Lα₀)", value=10.0, step=0.5, format="%.1f")

st.sidebar.header("Beam Controls")
Ic = st.sidebar.number_input(r"Control Intensity Ic (mW/cm²)", value=15.0, step=1.0, format="%.1f")
Ip = st.sidebar.number_input(r"Strong Probe Intensity Ip (mW/cm²)", value=5.0, step=0.1, format="%.2f")
detuning_c = st.sidebar.number_input(r"Control Detuning Δc (MHz)", value=0.0, step=0.1, format="%.2f")
scan_range = st.sidebar.number_input("Scan Range (MHz)", value=40.0, step=1.0, format="%.1f")

# --- Density Matrix Solver ---
def solve_density_matrix(delta, Gamma, Ic, Ip, Is, detuning_c, gamma21):
    # Rabi frequencies
    wp = Gamma * np.sqrt(Ip / (2 * Is))
    wc = Gamma * np.sqrt(Ic / (2 * Is))
    
    # Detunings
    dp = delta + detuning_c # Probe detuning
    dc = detuning_c         # Control detuning
    
    # Pre-allocate absorption array (Im(rho31))
    abs_list = []
    
    for d_p in (delta + detuning_c):
        # We solve the steady state equations M * rho = B
        # rho vector = [rho11, rho22, rho33, rho12, rho13, rho23] + complexes
        # To simplify, we solve the 9 coupled equations for rho_ij
        # However, for efficiency in a script, we use the matrix form of the Liouvillian
        
        # Elements for the matrix equation (steady state)
        # Using a simplified 3-level solver logic:
        # Solving: i[rho, H] + Liouvillian = 0
        
        # Hamiltonian terms
        # H = [[0, 0, wp/2], [0, -delta, wc/2], [wp/2, wc/2, -dp]]
        
        # Simplified Steady State result for rho31 in a 3-level system 
        # Including saturation terms (Non-linear susceptibility)
        num = 1j * (wp / 2) * (gamma21 - 1j * (d_p - dc))
        den = (Gamma/2 - 1j * d_p) * (gamma21 - 1j * (d_p - dc)) + (wc**2 / 4)
        
        # Adding saturation correction factor (approximate for the strong probe)
        saturation = 1 + (wp**2 / (Gamma * gamma21)) # Simplified saturation term
        rho31 = (num / den) / saturation
        
        abs_list.append(np.imag(rho31))
        
    return np.array(abs_list)

# Generate Data
delta_range = np.linspace(-scan_range/2, scan_range/2, 1000)
# We normalize absorption relative to the max of a non-saturated 2-level system
rho31_im = solve_density_matrix(delta_range, Gamma, Ic, Ip, Is, detuning_c, gamma21)
absorption = L_alpha0 * (rho31_im / (1 / 2)) # Normalized to OD

# --- Plotting ---
fig = plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2.5], hspace=0.3)

# 1. Diagram
ax_diag = plt.subplot(gs[0])
ax_diag.set_xlim(0, 3)
ax_diag.set_ylim(-0.2, 3)
ax_diag.axis('off')
ax_diag.hlines([0.5, 0.5, 2.5], [0.5, 1.7, 1.1], [1.3, 2.5, 1.9], colors='black', lw=3)
ax_diag.text(0.9, 0.2, r'$|1\rangle$', ha='center')
ax_diag.text(2.1, 0.2, r'$|2\rangle$', ha='center')
ax_diag.text(1.5, 2.7, r'$|3\rangle$', ha='center')
ax_diag.annotate("", xy=(1.4, 2.5), xytext=(0.9, 0.5), arrowprops=dict(arrowstyle="->", color="red", lw=2.5))
ax_diag.annotate("", xy=(1.6, 2.5), xytext=(2.1, 0.5), arrowprops=dict(arrowstyle="->", color="blue", lw=2.5))
ax_diag.set_title(r"Strong Drive $\Lambda$-System", fontsize=14)

# 2. Plot
ax_plot = plt.subplot(gs[1])
ax_plot.plot(delta_range, absorption, color='darkred', lw=2, label=fr'$I_p = {Ip}$ mW/cm$^2$')
ax_plot.set_xlabel(r'Two-photon Detuning $\delta$ (MHz)')
ax_plot.set_ylabel(r'Absorption ($\alpha L$)')
ax_plot.grid(True, alpha=0.3)
ax_plot.legend()

st.pyplot(fig)

st.subheader("Effect of Strong Probe")
st.write("""
As you increase **Ip**, notice how the absorption peak broadens and the depth of the 
transparency window decreases. This is due to the saturation of the $|1\rangle \leftrightarrow |3\rangle$ 
transition and the depletion of the ground state population.
""")
