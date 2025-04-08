from nicegui import ui
from crisprzip import *
import plotly.graph_objects as go
import numpy as np
import json
import os
from plotly.subplots import make_subplots
from scipy.optimize import curve_fit
from crisprzip.kinetics import *

###############################################################################
# Helper functions: CRISPRzip predictions
###############################################################################

# Not needed
# # Load models (lanscape params)
# def load_model_params(model_name: str) -> dict:
#     """Load parameters from CRISPRzip JSON files"""
#     model_map = {
#         'sequence-params2': 'sequence_params.json',
#         'average-params': 'average_params.json',
#         'average-params-legacy': 'average_params_legacy.json'
#     }
    
#     base_path = os.path.join(os.path.dirname(crisprzip.__file__), 'landscape_params/')
#     with open(os.path.join(base_path, model_map[model_name])) as f:
#         return json.load(f)['param_values']
    
# Function to get cleavage rate
def get_cleavage_rate(stc, binding_rate):
    """Calculate cleavage rate."""
    dt = np.logspace(-2, 6)
    f_clv = stc.get_cleaved_fraction(dt, binding_rate)
    with np.errstate(over='ignore'): # ignore RuntimeWarning: overflow
        k_eff = np.exp(curve_fit(
            f=lambda t, logk: 1 - np.exp(-np.exp(logk) * t),
            xdata=dt,
            ydata=f_clv,
        )[0][0])
    return k_eff

# Submit button handler for in vitro cleavage
def submit_handler_in_vitro_cleavage():
    try:
        parameter_set = model_dropdown.value
        bare_protein = load_landscape(parameter_set)
        
        protospacer = target_sequence_input.value
        off_targets = [seq.strip() for seq in off_targets_input.value.split(',')]
        targets = [protospacer] + off_targets
        
        guided_protein = bare_protein.bind_guide_rna(protospacer=protospacer)
        
        protein_sequence_complexes = []
        for target_seq in targets:
            protein_sequence_complexes += [guided_protein.probe_sequence(target_seq)]
        
        concentration = float(rnp_concentration_input.value)
        k_on_ref = 1E-2
        binding_rate = k_on_ref * concentration

        # Clear previous output content - [IMPORTANT], otherwise plots will stack below on each request
        output_container.clear()
        
        # Targets cleaved over time - needs adjusting and not displayed currently
        # fig_tt = go.Figure()
        # stc0 = protein_sequence_complexes[0]
        # kfit0 = get_cleavage_rate(stc0, k_on_ref)
        # x_values = np.logspace(-2, 6)
        
        # fig_tt.add_trace(go.Scatter(
        #     x=x_values,
        #     y=stc0.get_cleaved_fraction(x_values, k_on_ref),
        #     mode='markers',
        # ))
        
        # fig_tt.add_trace(go.Scatter(
        #     x=x_values,
        #     y=1 - np.exp(-kfit0 * x_values),
        #     mode='lines',
        #     name=f'exp. fit (k={kfit0:.2e} s⁻¹)'
        # ))
        
        # fig_tt.update_layout(
        #     xaxis_type="log",
        #     xaxis_title="time (s)",
        #     yaxis_title="fraction cleaved f_clv",
        #     legend=dict(font=dict(size=10)),
        #     height=300,
        #     margin=dict(l=50, r=50, t=30, b=50)
        # )
        #
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
            k_fit = [get_cleavage_rate(stc, k_on_ref * concentration) for concentration in dc]
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


###############################################################################
# GUI layout with NiceGUI
###############################################################################

# Title

# with ui.header(elevated=True).style('background-color: #D9F0D2'):
#     ui.label('CRISPRzip GUI').style('color: black; font-size: 32px; font-family: Courier New; font-weight: bold;')

with ui.header(elevated=True).style('background-color: #D9F0D2'):
    with ui.row().classes():
        # Image/icon to the left
        ui.image('./img/CRISPRzip_logo_v0.jpg').props('width=45px height=45px')
    
        # Text next to the image
        ui.label('CRISPRzip GUI').style('color: black; font-size: 32px; font-family: Courier New; font-weight: bold;')

# Footer with explicit height and centered content
with ui.footer(elevated=True).classes('py-1 h-10 bg-[#D9F0D2] flex items-center'):
    ui.markdown('CRISPRzip GUI is created with [niceGUI](https://nicegui.io). Licensed under MIT.').style('color: gray; font-size: 10px;')

# Create primary tabs
with ui.tabs().classes('w-full bg-gradient-to-b from-white via-white to-[#6ECFF6] rounded-xl') as tabs:
    one = ui.tab('in vitro', icon='science')
    two = ui.tab('in vivo', icon='pest_control_rodent')

# Create primary tab panels
with ui.tab_panels(tabs, value=one).classes('w-full'):
    # Content for first primary tab (in vitro)
    with ui.tab_panel(one):
        # Create nested tabs inside the in vitro tab panel
        with ui.tabs().classes('justify-start') as in_vitro_tabs:
            vitro_cleavage = ui.tab('cleavage (Cas9)', icon='content_cut').props('no-caps')
            vitro_binding = ui.tab('binding (dCas9)', icon='link').props('no-caps')

        # Create panels for the in vitro nested tabs
        with ui.tab_panels(in_vitro_tabs, value=vitro_cleavage).classes('w-full'):        
###############################################################################            
            # Content for in vitro - cleavage tab
###############################################################################            
            with ui.tab_panel(vitro_cleavage).classes('w-full h-full no-wrap'):
                ui.markdown('Predictions of cleavage by active Cas9 in a typical in vitro setting, where RNP is abundant enough to cleave targets independently.').style('font-size: 16px')
                with ui.row().classes('w-full h-full no-wrap'):
                    with ui.column().classes('w-2/5 h-full'):
                        ui.markdown('<u>**Inputs**</u>').style('font-size: 20px; color: gray;')
                        # Target seq and off-target seq rows
                        with ui.row().classes('w-full h-full no-wrap'):
                            ui.markdown('**Target sequence:**').style('font-size: 16px')
                            # Input field for target sequence
                            target_sequence_input = ui.input(placeholder='Enter target sequence (e.g., GACGCATAAAGATGAGACGCTGG)').classes('w-full')
                            # Radio buttons to select between "protospacer" and "guide RNA"
                            def update_placeholder(value):
                                if value == 'protospacer':
                                    target_sequence_input.placeholder = 'Enter DNA sequence (e.g., GACGCATAAAGATGAGACGCTGG)'
                                elif value == 'guide RNA':
                                    target_sequence_input.placeholder = 'Enter RNA sequence (e.g., GACGCATAAAGAUCUCGCUGG)'
                                target_sequence_input.update()  # Refresh the input to apply changes    
                            ui.radio(['protospacer', 'guide RNA'], value='protospacer', on_change=update_placeholder)
                            # Information balloon for target sequence
                            ui.icon('info').tooltip('Provide the DNA sequence as: 5’-20 nt sequence + PAM-3’, e.g., AGACGCATAAAGATGAGACGCTGG').style('font-size: 20px')
                        
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
                            ui.markdown('**Off-target sequences:**').style('font-size: 16px')
                            off_targets_input = ui.textarea(placeholder='Enter up to 5 off-target sequences, comma-separated').props('rows=3').classes('w-1/2')

                            # Function to validate and split input
                            def validate_off_targets(text):
                                sequences = [seq.strip() for seq in text.split(',')]
                                if len(sequences) > 5:
                                    ui.notify('Error: Limit exceeded. Only 5 sequences currently allowed.', type='negative')
                                    return ','.join(sequences[:5])  # Return only the first 5 sequences
                                return text

                            # Attach the validation function to the input
                            off_targets_input.on('update:model-value', lambda e: off_targets_input.set_value(validate_off_targets(off_targets_input.value)))

                            # Upload button for file-based input (currently disabled), placeholder bnelow
                            # upload_button = ui.button('Upload').props('icon=upload').classes('w-1/3')
                            ui.icon('info').tooltip('Specify off-target sequences. You can check to up to 5 off-target sequences, comma-separated.').style('font-size: 20px')

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
                            ui.icon('info').tooltip('Select the model for cleavage predictions. Recommended: sequence-params2.').style('font-size: 20px')
                            # Insert spacer
                            ui.markdown('<br>')
                            # RNP concentration input
                            ui.markdown('**RNP concentration (nM):**').style('font-size: 16px')
                            rnp_concentration_input = ui.input(value='100').props('type=number').classes('ml-2 w-20')
                        
                        # Submit button with event handler
                        submit_button = ui.button('Submit', on_click=submit_handler_in_vitro_cleavage).props('icon=send')

                    # Vertical separator between inputs and outputs
                    ui.separator().props('vertical')

                    with ui.column().classes('w-3/5 h-full no-wrap'):
                        ui.markdown('<u>**Output**</u>').style('font-size: 20px; color: gray;')
                        output_container = ui.element('div').classes('w-full no-wrap')
                    
###############################################################################    
                # Content for in vitro - binding tab
###############################################################################                
            with ui.tab_panel(vitro_binding):
                ui.label('IN VITRO - BINDING (dCas9) GOES HERE')

###############################################################################
###############################################################################
    # Content for second primary tab (in vivo)
    with ui.tab_panel(two):
        # Create nested tabs inside the in vivo tab panel
        with ui.tabs().classes('justify-start') as in_vivo_tabs:
            vivo_cleavage = ui.tab('cleavage (Cas9)', icon='content_cut').props('no-caps')
            vivo_binding = ui.tab('binding (dCas9)', icon='link').props('no-caps')        
        # Create panels for the in vivo nested tabs
        with ui.tab_panels(in_vivo_tabs, value=vivo_cleavage).classes('w-full'):
###############################################################################              
            # Content for in vivo - cleavage tab
###############################################################################              
            with ui.tab_panel(vivo_cleavage):
                ui.label('IN VIVO - CLEAVAGE (Cas9) PLACEHOLDER')
 
###############################################################################              
            # Content for in vivo - binding tab
###############################################################################
            with ui.tab_panel(vivo_binding):
                ui.label('IN VIVO - BINDING (dCas9) PLACEHOLDER')

###############################################################################
# Run the NiceGUI app
###############################################################################
ui.run()
