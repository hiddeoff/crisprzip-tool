from nicegui import ui
from crisprzip import *
import plotly.graph_objects as go
import numpy as np
import json
import os
from plotly.subplots import make_subplots
from scipy.optimize import curve_fit
from crisprzip.kinetics import *


def show():
    def check_sequence_input(sequence):
        if len(sequence) != 23:
            raise ValueError(f"Your input of length {len(sequence)} does not "
                             f"follow the required format 5'-target-PAM-3' "
                             f"(23 nts total).")
        if sequence[-2:] != "GG":
            raise ValueError(
                f"Currently, the non-canonical PAM '{sequence[-3:]}' "
                f"is not supported. Please provide a target with a "
                f"canonical 'NGG' PAM.")
        for nt in sequence:
            if nt not in ["A", "C", "G", "T"]:
                raise ValueError(f"Nucleotide '{nt}' could not be recognized. "
                                 f"Please specify A, C, G or T.")

    def process_offtarget_input(inputvalue):
        off_targets = []
        for seq in inputvalue.split(","):
            if seq.strip() != "":
                off_targets.append(seq.strip())
        return off_targets

    def make_stc_list(protospacer, off_targets, parameter_set):
        """Generate SearcherTargetComplexes."""
        bare_protein = load_landscape(parameter_set)
        guided_protein = bare_protein.bind_guide_rna(protospacer=protospacer)

        targets = [protospacer] + off_targets
        protein_sequence_complexes = [guided_protein.probe_sequence(target_seq)
                                      for target_seq in targets]
        return protein_sequence_complexes

    def get_cleavage_rate(stc, binding_rate):
        """Calculate cleavage rate."""
        dt = np.logspace(-2, 6)
        f_clv = stc.get_cleaved_fraction(dt, binding_rate)
        with np.errstate(over='ignore'):  # ignore RuntimeWarning: overflow
            k_eff = np.exp(curve_fit(
                f=lambda t, logk: 1 - np.exp(-np.exp(logk) * t),
                xdata=dt,
                ydata=f_clv,
            )[0][0])
        return k_eff

    def get_all_cleavage_rates(protospacer, off_targets,
                               parameter_set, concentration):
        protein_sequence_complexes = make_stc_list(
            protospacer=protospacer,
            off_targets=[protospacer] + off_targets,
            parameter_set=parameter_set,
        )
        k_on_ref = 1E-2
        binding_rate = k_on_ref * concentration

        k_fit_values = [
            get_cleavage_rate(stc, binding_rate) for
            stc in protein_sequence_complexes
        ]
        return k_fit_values

    # Submit button handler for in vitro cleavage
    def submit_handler_in_vitro_cleavage():
        try:

            # Collect user input
            protospacer = "AGCTAAGCTAAGCTAAGCTAGGG"
            off_targets = ["AGCTAAGCTATGCTAAGCTAGGG",
                           "AGCTAAGCTAAGGTGAGCTAGGG",
                           "AGCTAAGCCCCGCTAAGCTAGGG"]
            parameter_set = "sequence_params"
            concentration = 100.

            # Process user input
            for seq in [protospacer] + off_targets:
                check_sequence_input(seq)

            protein_sequence_complexes = make_stc_list(
                protospacer=protospacer,
                off_targets=off_targets,
                parameter_set=parameter_set,
            )

            k_on_ref = 1E-2
            binding_rate = k_on_ref * concentration

            # Clear previous output content - [IMPORTANT], otherwise plots will stack below on each request
            output_container.clear()

            # FIGURE 0
            k_clv_vals = get_all_cleavage_rates(
                protospacer, off_targets, parameter_set, concentration
            )
            fig0 = go.Figure([go.Bar(x=np.arange(len(off_targets) + 1),
                                     y=k_clv_vals)])


            # FIGURE 1 - Cleaved fraction vs time
            fig1 = go.Figure()
            dt = np.logspace(-1, np.log10(7200))

            for i, stc in enumerate(protein_sequence_complexes):
                f_clv = stc.get_cleaved_fraction(dt, binding_rate)
                fig1.add_trace(go.Scatter(
                    x=dt, y=f_clv,
                    name=('on-target' if i == 0 else f'off-target {i}'),
                    mode='lines')
                )
            fig1.update_layout(
                title='Cleaved fraction vs time',
                xaxis_type="log",
                xaxis_title="time (s)",
                yaxis_title="fraction cleaved (f_clv)",
                height=400,
                margin=dict(l=50, r=50, t=50, b=50)
            )

            # FIGURE 2 - Cleaved rate vs time
            fig2 = go.Figure()
            # Concentration range, currently hardcoded and needs a better solution
            conc_min = .1  # nM
            conc_max = 100  # nM
            dc = np.logspace(np.log10(conc_min), np.log10(conc_max))

            for i, stc in enumerate(protein_sequence_complexes):
                k_fit = [get_cleavage_rate(stc, k_on_ref * concentration) for
                         concentration in dc]
                fig2.add_trace(go.Scatter(
                    x=dc, y=k_fit,
                    name=('on-target' if i == 0 else f'off-target {i}'),
                    mode='lines')
                )
            fig2.update_layout(
                title='Cleaved rate vs time',
                xaxis_type="log",
                yaxis_type="log",
                xaxis_title="Cas9 concentration (nM)",
                yaxis_title="cleavage rate (k_clv)",
                height=400,
                margin=dict(l=50, r=50, t=50, b=50)
            )
            with output_container:
                with ui.column(align_items='center').classes('w-full'):
                    ui.plotly(fig0)
                    with ui.row().classes('w-full no-wrap'):
                        ui.plotly(fig1).classes('w-1/2')
                        ui.plotly(fig2).classes('w-1/2')

        except Exception as e:
            ui.notify(f'Error: {str(e)}', type='negative')

    output_container = ui.element('div').classes('w-full no-wrap')
    submit_handler_in_vitro_cleavage()


@ui.page('/')
def index():
    show()

ui.run()
