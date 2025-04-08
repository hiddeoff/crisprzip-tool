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

    # Submit button handler for in vitro cleavage
    def submit_handler_in_vitro_cleavage():
        try:

            # Collect user input
            protospacer = target_sequence_input.value
            off_targets = process_offtarget_input(off_targets_input.value)
            parameter_set = model_dropdown.value
            concentration = float(rnp_concentration_input.value)

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

            # Cleaved fraction vs time
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
            # Cleaved rate vs time
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
                with ui.row().classes('w-full no-wrap'):
                    ui.plotly(fig1).classes('w-1/2')
                    ui.plotly(fig2).classes('w-1/2')
        except Exception as e:
            ui.notify(f'Error: {str(e)}', type='negative')

    ui.markdown(
        'Predictions of cleavage by active Cas9 in a typical in vitro setting, where RNP is abundant enough to cleave targets independently.').style(
        'font-size: 16px')


    with ui.row().classes('w-full h-full no-wrap'):
        with ui.column().classes('w-2/5 h-full'):

            with ui.grid(columns=3).classes('w-full gap-0'):

                # ROW 1: TARGET SEQUENCE INPUT

                def update_placeholder(value):
                    if value == 'protospacer':
                        target_sequence_input.placeholder = 'Enter DNA sequence (e.g., GACGCATAAAGATGAGACGCTGG)'
                    elif value == 'guide RNA':
                        target_sequence_input.placeholder = 'Enter RNA sequence (e.g., GACGCATAAAGAUCUCGCUGG)'
                    target_sequence_input.update()  # Refresh the input to apply changes

                ui.markdown('**target sequence**').style('font-size: 16px')

                ui.input(
                    placeholder='Enter target sequence (e.g., GACGCATAAAGATGAGACGCTGG)'
                ).classes('w-full')

                ui.radio(['protospacer', 'guide RNA'], value='protospacer',
                         on_change=update_placeholder).classes('m-0 p-0 space-y-0 text-sm')


            # Target seq and off-target seq rows
            with ui.row().classes('w-full h-full no-wrap'):
                ui.markdown('**Target sequence:**').style('font-size: 16px')
                # Input field for target sequence
                target_sequence_input = ui.input(
                    placeholder='Enter target sequence (e.g., GACGCATAAAGATGAGACGCTGG)').classes(
                    'w-full')


                ui.radio(['protospacer', 'guide RNA'], value='protospacer',
                         on_change=update_placeholder)
                # Information balloon for target sequence
                ui.icon('info').tooltip(
                    'Provide the DNA sequence as: 5’-20 nt sequence + PAM-3’, e.g., AGACGCATAAAGATGAGACGCTGG').style(
                    'font-size: 20px')

            # Insert spacer
            ui.markdown('---')

            # # Off-target sequences textarea and upload button
            # with ui.row().classes('w-full h-full no-wrap'):
            #     ui.markdown('**Off-target sequences:**').style('font-size: 16px')
            #     off_target_textarea = ui.textarea(placeholder='Off-target sequence(-es)').props('rows=3').classes('w-1/2')
            #     upload_button = ui.button('Upload').props('icon=upload').classes('w-1/3')
            #     # Information balloon for off-target sequence
            #     ui.icon('info').tooltip('Upload off-target sequences from a local text file.').style('font-size: 20px')

            # Off-target sequences input
            with ui.row().classes('w-full h-full no-wrap'):
                ui.markdown('**Off-target sequences:**').style(
                    'font-size: 16px')
                off_targets_input = ui.textarea(
                    placeholder='Enter up to 5 off-target sequences, comma-separated').props(
                    'rows=3').classes('w-1/2')

                # Function to validate and split input
                def validate_off_targets(text):
                    sequences = [seq.strip() for seq in text.split(',')]
                    if len(sequences) > 5:
                        ui.notify(
                            'Error: Limit exceeded. Only 5 sequences currently allowed.',
                            type='negative')
                        return ','.join(
                            sequences[:5])  # Return only the first 5 sequences
                    return text

                # Attach the validation function to the input
                off_targets_input.on('update:model-value',
                                     lambda e: off_targets_input.set_value(
                                         validate_off_targets(
                                             off_targets_input.value)))

                # Upload button for file-based input (currently disabled), placeholder bnelow
                # upload_button = ui.button('Upload').props('icon=upload').classes('w-1/3')
                ui.icon('info').tooltip(
                    'Specify off-target sequences. You can check to up to 5 off-target sequences, comma-separated.').style(
                    'font-size: 20px')

            # Insert spacer
            ui.markdown('---')

            # Model dropdown with information balloon. Then next to it, RNP conc.
            with ui.row().classes('w-full h-full no-wrap'):
                # Model label
                ui.markdown('**Model:**').style('font-size: 16px')

                # Simple dropdown with model names
                model_dropdown = ui.select(
                    options={
                        'sequence_params': 'sequence-params2 (recommended)',
                        'average_params': 'average-params',
                        'average_params_legacy': 'average-params-legacy'
                    },
                    value='sequence_params'
                ).classes('w-1/4')

                # Not needed
                # # 3. Function to load parameters silently
                # def update_params():
                #     global loaded_params
                #     try:
                #         selected_value = model_dropdown.value
                #         loaded_params = load_model_params(selected_value)
                #         # Optional small notification
                #         ui.notify('Parameters loaded', type='positive', duration=1)
                #     except Exception as exc:
                #         ui.notify(f'Error: {str(exc)}', type='negative')

                # # 4. Connect the event handler
                # model_dropdown.on('update:model-value', lambda _: update_params())
                # Next row has an info icon
                ui.icon('info').tooltip(
                    'Select the model for cleavage predictions. Recommended: sequence-params2.').style(
                    'font-size: 20px')
                # Insert spacer
                ui.markdown('<br>')
                # RNP concentration input
                ui.markdown('**RNP concentration (nM):**').style(
                    'font-size: 16px')
                rnp_concentration_input = ui.input(value='100').props(
                    'type=number').classes('ml-2 w-20')

            # Submit button with event handler
            submit_button = ui.button('Submit',
                                      on_click=submit_handler_in_vitro_cleavage).props(
                'icon=send')

        # Vertical separator between inputs and outputs
        ui.separator().props('vertical')

        with ui.column().classes('w-3/5 h-full no-wrap'):
            ui.markdown('<u>**Output**</u>').style(
                'font-size: 20px; color: gray;')
            output_container = ui.element('div').classes('w-full no-wrap')
