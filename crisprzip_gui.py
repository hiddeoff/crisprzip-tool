from nicegui import ui
import crisprzip
import plotly.graph_objects as go
import numpy as np
import json
import os

###############################################################################
# Helper functions: CRISPRzip predictions
###############################################################################

# Load models (lanscape params)
def load_model_params(model_name: str) -> dict:
    """Load parameters from CRISPRzip JSON files"""
    model_map = {
        'sequence-params2': 'sequence_params.json',
        'average-params': 'average_params.json',
        'average-params-legacy': 'average_params_legacy.json'
    }
    
    base_path = os.path.join(os.path.dirname(crisprzip.__file__), 'landscape_params/')
    with open(os.path.join(base_path, model_map[model_name])) as f:
        return json.load(f)['param_values']

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
                            target_sequence_input = ui.input(placeholder='Enter target sequence').classes('w-full')
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

                        # Off-target sequences textarea and upload button
                        with ui.row().classes('w-full h-full no-wrap'):
                            ui.markdown('**Off-target sequences:**').style('font-size: 16px')
                            off_target_textarea = ui.textarea(placeholder='Off-target sequence(-es)').props('rows=3').classes('w-1/2')
                            upload_button = ui.button('Upload').props('icon=upload').classes('w-1/3')
                            # Information balloon for off-target sequence
                            ui.icon('info').tooltip('Upload off-target sequences from a local text file.').style('font-size: 20px')
     
                        # Insert spacer
                        ui.markdown('---')

                        # Model dropdown with information balloon. Then next to it, RNP conc.
                        with ui.row().classes('w-full h-full no-wrap'):
                            # 1. Model label
                            ui.markdown('**Model:**').style('font-size: 16px')

                            # 2. Simple dropdown with model names
                            model_dropdown = ui.select(
                                options={
                                    'sequence-params2': 'sequence-params2 (recommended)',
                                    'average-params': 'average-params', 
                                    'average-params-legacy': 'average-params-legacy'
                                },
                                value='sequence-params2'
                            ).classes('w-1/4')

                            # 3. Function to load parameters silently
                            def update_params():
                                global loaded_params
                                try:
                                    selected_value = model_dropdown.value
                                    loaded_params = load_model_params(selected_value)
                                    # Optional small notification
                                    ui.notify('Parameters loaded', type='positive', duration=1)
                                except Exception as exc:
                                    ui.notify(f'Error: {str(exc)}', type='negative')
                            
                            # 4. Connect the event handler
                            model_dropdown.on('update:model-value', lambda _: update_params())
                            # Next row has an infi icon
                            ui.icon('info').tooltip('Select the model for cleavage predictions. Recommended: sequence-params2.').style('font-size: 20px')
                            # Insert spacer
                            ui.markdown('<br>')
                            # RNP concentration input
                            ui.markdown('**RNP concentration (nM):**').style('font-size: 16px')
                            rnp_concentration_input = ui.input(value='100').props('type=number').classes('ml-2 w-20')
                        
                        # Submit button (one row)
                        submit_button = ui.button('Submit').props('icon=send')

                    # Vertical separator between inputs and outputs
                    ui.separator().props('vertical')

                    with ui.column().classes('w-3/5 h-full no-wrap'):
                        ui.markdown('<u>**Output**</u>').style('font-size: 20px; color: gray;')
                    
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
