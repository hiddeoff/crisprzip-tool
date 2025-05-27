import re
from io import StringIO

from nicegui import ui, events
import pandas as pd

from crisprzip import *
from crisprzip.kinetics import *


initial_input = True  # auto-fills upon load - useful when developing


initial_input = True  # auto-fills upon load - useful when developing


def show_input():
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

    def handle_offtarget_uploads(e: events.UploadEventArguments):
        try:
            off_targets_input.value = ""  # clear the input field
            # read the file content
            content = e.content.read().decode("utf-8")

            # parse the csv (one column comma separated)
            if e.name.lower().endswith('.csv'):
                with StringIO(content) as f:
                    df = pd.read_csv(f)
                    sequences = df.iloc[:, 0].tolist()
                    formatted_sequences = ',\n'.join(sequences)
                    off_targets_input.value = formatted_sequences
            else:
                # otherwise juts use the content then from the off-target box
                off_targets_input.value = content

            ui.notify('File uploaded successfully', type='positive')
        except Exception as e:
            ui.notify(f'Error processing file: {str(e)}', type='negative')

    def process_ontarget_input(inputvalue, inputtype):
        if inputtype == "protospacer":
            return inputvalue
        elif inputtype == "guide RNA":
            return inputvalue.replace("U", "T") + "NGG"
        else:
            raise ValueError(f"Unrecognized input type \'{inputtype}\'.")

    def process_offtarget_input(inputvalue):
        off_targets = []
        for seq in inputvalue.split(","):
            if seq.strip() != "":
                off_targets.append(seq.strip())
        return off_targets

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

    def concentration_validation(input) -> str | None:
        concentration = str(input)
        if not concentration.strip():
            return "RNP concentration cannot be empty."
        if not re.fullmatch(r"^\d*(\.\d*)?([eE][+-]?\d+)?$", concentration):
            return "Please enter a valid number."
        return None  # Return None if the input is valid

    wc1 = 230  # column 1 width
    wc2 = 100  # column 2 width
    fsz = 10  # font size (in pt)
    fsi = 12   # font size for information bubble (in pt)

    # CONTENT
    with ui.grid(columns=2).style(
            f'grid-template-columns: {wc1}px {wc2}px').classes('gap-0'):

        # TARGET SEQUENCE
        with ui.row(align_items='center').classes('p-0'):
            ui.markdown('**Target sequence**').classes(
                'p-0 leading-[0.7]').style(f'font-size: {fsz}pt')
            ui.icon('info').tooltip(
                'Select the model for cleavage predictions. Recommended: sequence-params2.').style(
                f'font-size: {fsi}pt')
        ui.element()

        with ui.column().classes('w-full p-0'):
            target_sequence_input = ui.input(
                validation=lambda x: sequence_validation(x, 23),
                # once the update_placeholder() function works, this arg should be omitted
                placeholder='GACGCATAAAGATGAGACGCTGG',
                value=('GACGCATAAAGATGAGACGCTGG' if initial_input else None),
            ).classes(f'w-[{wc1 - 20}px] font-mono').props('dense').style(
                f'font-size: {fsz}pt')

        with ui.column().classes('w-full p-0'):
            target_input_select = (
                ui.select(['protospacer', 'guide RNA'],
                          value='protospacer',
                          on_change=update_placeholder)
                .classes('w-full').props('dense').style(f'font-size: {fsz}pt')
            )
            update_placeholder(target_input_select.value)

        # OFF-TARGET SEQUENCES
        with ui.row(align_items='center').classes('p-0'):
            ui.markdown('**Off-target sequences**').classes(
                'p-0 leading-[0.7]').style(f'font-size: {fsz}pt')
            ui.icon('info').tooltip(
                'Select the model for cleavage predictions. Recommended: sequence-params2.').style(
                f'font-size: {fsi}pt')
        ui.element()

        with ui.column().classes('w-full h-full p-0'):
            off_targets_input = ui.textarea(
                placeholder='GACGCATAAAGATGAGACGCTGG,\nGACGCATAAAGATGAGACGCTGG,\n...',
                validation=lambda x: sequence_validation(x, None),
                value=(('GACGAACAAAGATGAGACGCTGG,\n' +
                        'GACGCATATATACGAGACGCTGG,\n' +
                        'GACGCATAATTATGAGTCGCTGG,\n' +
                        'GACGCATACCGATGTGTCGCTGG,\n' +
                        'GACGCATAAAGATGGGGCTCTGG,\n')
                       if initial_input else None),
            ).props('rows=5 dense').classes(
                f'w-[{wc1 - 20}px] h-2fr font-mono').style(
                f'font-size: {fsz}pt')

        with ui.row(align_items='start').classes('w-full h-full p-0'):
            with ui.row(align_items='center').classes('w-[52px] h-[52px]'):
                upload_component = ui.upload(
                    on_upload=handle_offtarget_uploads,
                    on_rejected=lambda e: ui.notify('File upload failed',
                                                    type='warning'),
                    auto_upload=True
                ).props("accept=.csv hide-upload-btn").classes(
                    'hidden')  # can add .txt (for example) to .props if you want to uplaod something other than .csv files
                # add custom ui button (instead of the regular upload button)
                ui.button('upload', on_click=lambda: (
                    upload_component.reset(),
                    upload_component.run_method('pickFiles')
                )).props('outline no-caps').style(f'font-size: {fsz}pt')

        with ui.column(align_items='start').classes(
                f'w-[{wc1 - 20}px] p-0 m-0 gap-0'):
            # CONTEXT
            with ui.row(align_items='center').classes('w-full p-0'):
                ui.markdown('**Context**').classes('leading-[0]').style(
                    f'font-size: {fsz}pt')
                (ui.icon('info')
                 .tooltip('Select the application context.')
                 .style(f'font-size: {fsi}pt'))
            context_dropdown = ui.select(
                options={'invitro': 'cell-free (in vitro)',
                         'ecoli': 'E. coli',
                         'mammal': 'mammal', },
                value='invitro',
            ).props('dense').classes(f'w-full p-0 m-0').style(
                f'font-size: {fsz}pt')
            ui.element().classes("h-3")

            # PARAMETER SELECTION
            with ui.row(align_items='center').classes('w-full p-0'):
                ui.markdown('**Landscape parameters**').classes(
                    'leading-[0]').style(f'font-size: {fsz}pt')
                (ui.icon('info')
                 .tooltip('Select the model for cleavage predictions.')
                 .style(f'font-size: {fsi}pt'))
            model_dropdown = ui.select(
                options={
                    'sequence_params': 'sequence (default)',
                    'average_params': 'average',
                    'average_params_legacy': 'average (legacy)'
                },
                value='sequence_params',
            ).props('dense').classes(f'w-full p-0 m-0').style(
                f'font-size: {fsz}pt')
            ui.element().classes("h-6")

        with ui.column(align_items='center'):
            with ui.row(align_items='center').classes('h-full w-full p-0'):
                (ui.button('show').style(f'font-size: {fsz}pt').props(
                    'outline no-caps'))

        submit_button = (
            ui.button('Submit')
            .props('icon=send')
            .classes(f'w-[{wc1 - 20}px]')
            .style(f'font-size: {fsz}pt')
        )

    def get_input_values():
        ontarget = process_ontarget_input(
            target_sequence_input.value,
            target_input_select.value
        )

        input_vals = {
            'on_target': ontarget,
            'off_targets': process_offtarget_input(off_targets_input.value),
            'context': context_dropdown.value,
            'parameter_set': model_dropdown.value
        }
        return input_vals

    # Return the submit button and the function to get input values
    return submit_button, get_input_values, model_dropdown


def get_k_on_off(context):
    if context == 'invitro':
        k_on = 0.1
        k_off = 1.0
    elif context == 'ecoli':
        # E. coli volume:      1 µm³ = 1 fL
        # no. of PAMs:         1E6 (500k / genome, 2 genomes)
        # PAM-binding rate:   60 s⁻¹
        # PAM-unbinding rate: 40 s⁻¹
        k_on = 6.02E23 * 1E-9 * 1E-15 / 1E6 / (1 / 60 + 1 / 40)
        k_off = 40
    elif context == 'mammal':
        # mammal nuclear volume: 500 µm³ = 500 fL
        # no. of PAMs:          13E6 (320 mln / genome, 2 genomes, 2% available)
        # PAM-binding rate:     1.33 s⁻¹
        # PAM-unbinding rate:     40 s⁻¹
        k_on = 6.02E23 * 1E-9 * 500E-15 / 13E6 / (1 / 1.33 + 1 / 40)
        k_off = 40
    else:
        raise ValueError(f"Unknown context '{context}'")
    return k_on, k_off


def make_stc_list(protospacer, off_targets, context, parameter_set):
    """Generate SearcherTargetComplexes."""

    k_on, k_off = get_k_on_off(context)

    if parameter_set == 'sequence_params':
        bare_protein = crisprzip.kinetics.load_landscape(parameter_set)
        bare_protein.internal_rates['k_off'] = k_off
        guided_protein = bare_protein.bind_guide_rna(protospacer=protospacer)

        targets = [protospacer] + off_targets
        protein_sequence_complexes = [guided_protein.probe_sequence(target_seq)
                                      for target_seq in targets]
        return protein_sequence_complexes

    elif (parameter_set == 'average_params') or (
            parameter_set == 'average_params_legacy'):
        protein = crisprzip.kinetics.load_landscape(parameter_set)
        protein.internal_rates['k_off'] = k_off

        targets = [protospacer] + off_targets
        protein_sequence_complexes = []
        for target_seq in targets:
            mm_pattern = (GuideTargetHybrid
                          .from_cas9_offtarget(target_seq, protospacer)
                          .get_mismatch_pattern())
            protein_sequence_complexes += [protein.probe_target(mm_pattern)]
        return protein_sequence_complexes
    else:
        raise ValueError(f"Unrecognized parameter set '{parameter_set}'.")

