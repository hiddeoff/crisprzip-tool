import os
import json
import re

from nicegui import ui
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.optimize import curve_fit

from crisprzip import *
from crisprzip.kinetics import *


def show():

    def update_placeholder(value):
        if value == 'protospacer':
            target_sequence_input.placeholder = 'GACGCATAAAGATGAGACGCTGG'
        elif value == 'guide RNA':
            target_sequence_input.placeholder = 'GACGCAUAAAGAUGAGACGC'
        target_sequence_input.update()  # Refresh the input to apply changes

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
            concentration = rnp_concentration_input.value

            # Validate concentration input
            concentration_error = concentration_validation(concentration)
            if concentration_error:
                message = "RNP concentration cannot be empty and has to be a valid number. Default value is 100."
                ui.notify(f'Error: {message}', type='negative')
                return  # Stop further processing if there's an error
            else:
                concentration = float(concentration)

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

    def sequence_validation(input, length=None) -> str:
        pattern = r'^[ACGTacgt\,\n\s]*$'
        if not (re.fullmatch(pattern, input)):
            return "Only ACGT nucleotides"
        if length is not None and len(input.strip()) > 0:
            input_length = len(input.strip())
            if input_length < length:
                return f"Too short ({input_length}/{length})"
            elif input_length > length:
                return f"Too long ({input_length}/{length})"

    def concentration_validation(input) -> str| None:
        concentration = str(input)
        if not concentration.strip():
            return "RNP concentration cannot be empty."
        if not re.fullmatch(r"^\d*(\.\d*)?([eE][+-]?\d+)?$", concentration):
            return "Please enter a valid number."
        return None  # Return None if the input is valid

    ui.markdown(
        'Predictions of cleavage by active Cas9 in a typical in vitro setting, '
        'where RNP is abundant enough to cleave targets independently.'
    ).style('font-size: 16px')

    with ui.row().classes('w-full h-full no-wrap'):

        # USER INPUT
        with ui.card().classes('p-4 m-4'):
            with ui.grid(columns=2).style(
                    'grid-template-columns: 300px 100px').classes('gap-0'):
                # TARGET SEQUENCE
                with ui.row(align_items='center').classes('p-0'):
                    ui.markdown('**Target sequence**').classes('p-0 leading-[0.7]')
                    ui.icon('info').tooltip(
                        'Select the model for cleavage predictions. Recommended: sequence-params2.').style(
                        'font-size: 20px')
                ui.element()

                with ui.column().classes('w-full p-0'):
                    target_sequence_input = ui.input(
                        validation=lambda x: sequence_validation(x, 23),
                        # once the update_placeholder() function works, this arg should be omitted
                        placeholder='GACGCATAAAGATGAGACGCTGG'
                    ).classes('w-[280px] font-mono').props('dense')
                    ui.element().classes("h-1")

                with ui.column().classes('w-full p-0'):
                    target_input_select = (
                        ui.select(['protospacer', 'guide RNA'],
                                  value='protospacer',
                                  on_change=update_placeholder)
                        .classes('w-[100px]').props('dense')
                    )
                    update_placeholder(target_input_select.value)

                # OFF-TARGET SEQUENCES
                with ui.row(align_items='center').classes('p-0'):
                    ui.markdown('**Off-target sequences**').classes(
                        'p-0 leading-[0.7]')
                    ui.icon('info').tooltip(
                        'Select the model for cleavage predictions. Recommended: sequence-params2.').style(
                        'font-size: 20px')
                ui.element()

                with ui.column().classes('w-full h-full p-0'):
                    off_targets_input = ui.textarea(
                        placeholder='GACGCATAAAGATGAGACGCTGG,\nGACGCATAAAGATGAGACGCTGG,\n...',
                        validation=lambda x: sequence_validation(x, None)
                    ).props('rows=5 dense').classes('w-[280px] h-2fr font-mono')
                    ui.element().classes("h-1")

                with ui.row(align_items='start').classes('w-full h-full p-0'):
                    with ui.row(align_items='center').classes('w-full h-[52px]'):
                        ui.button('upload').props('outline no-caps')

                # CONCENTRATION
                with ui.element().classes('w-full h-full p-0'):
                    with ui.row(align_items='start').classes('w-[280px] gap-0'):
                        ui.markdown('**RNP concentration**').classes(
                            'w-1/2 leading-[1.7]')
                        rnp_concentration_input = (ui.input(placeholder='100',
                                                            validation=concentration_validation)
                                                   .props('dense suffix="nM"')
                                                   .classes('w-1/2'))
                    ui.element().classes("h-4")

                ui.element()

                # PARAMETER SELECTION
                with ui.row(align_items='center').classes('w-full p-0'):
                    ui.markdown('**Model parameters**').classes('leading-[0.7]')
                    (ui.icon('info')
                     .tooltip('Select the model for cleavage predictions.')
                     .style('font-size: 20px'))
                ui.element()

                with ui.column().classes('w-full h-full p-0'):
                    model_dropdown = ui.select(
                        options={
                            'sequence_params': 'sequence-params2 (recommended)',
                            'average_params': 'average-params',
                            'average_params_legacy': 'average-params-legacy'
                        },
                        value='sequence_params'
                    ).props('dense').classes('w-[280px] p-0 m-0')

                with ui.row(align_items='center').classes('h-full w-full p-0'):
                    (ui.button(icon='search')
                     .classes('h-4 w-4').style('font-size: 12px').props('outline'))

                ui.element().classes("h-6")
                ui.element()

                submit_button = (
                    ui.button('Submit', on_click=submit_handler_in_vitro_cleavage)
                    .props('icon=send').classes('w-[280px]')
                )

        # OUTPUT
        with ui.column().classes('w-full h-full no-wrap'):
            ui.markdown('<u>**Output**</u>').style(
                'font-size: 20px; color: gray;')
            output_container = ui.element('div').classes('flex-grow no-wrap')


@ui.page("/")
def index():
    show()

ui.run()
