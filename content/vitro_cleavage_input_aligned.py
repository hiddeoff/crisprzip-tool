from nicegui import ui
import re

def show():

    show_borders = False
    bc = ' border' if show_borders else ''

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

    def concentration_validation(input) -> str:
        if not re.fullmatch(r"^\d*(\.\d*)?([eE][+-]?\d+)?$", input):
            return "Only numeric input"

    with ui.card().classes('p-4'):
        with ui.grid(columns=2).style('grid-template-columns: 300px 100px').classes('gap-0'):

            # TARGET SEQUENCE
            with ui.row(align_items='center').classes('p-0' + bc):
                ui.markdown('**Target sequence**').classes('p-0 leading-[0.7]')
                ui.icon('info').tooltip(
                    'Select the model for cleavage predictions. Recommended: sequence-params2.').style(
                    'font-size: 20px')
            ui.element()

            with ui.column().classes('w-full p-0' + bc):
                target_sequence_input = ui.input(
                    placeholder='GACGCATAAAGATGAGACGCTGG',
                    validation=lambda x: sequence_validation(x, 23),
                ).classes('w-[280px] font-mono').props('dense')

                ui.element().classes("h-1")

            with ui.column().classes('w-full p-0' + bc):
                target_input_select = ui.select(
                    ['protospacer', 'guide RNA'],
                    value='protospacer').classes('w-[100px]').props('dense')

            # OFF-TARGET SEQUENCES
            with ui.row(align_items='center').classes('p-0' + bc):
                ui.markdown('**Off-target sequences**').classes('p-0 leading-[0.7]')
                ui.icon('info').tooltip(
                    'Select the model for cleavage predictions. Recommended: sequence-params2.').style(
                    'font-size: 20px')
            ui.element()

            with ui.column().classes('w-full h-full p-0' + bc):

                off_targets_input = ui.textarea(
                    placeholder='GACGCATAAAGATGAGACGCTGG,\nGACGCATAAAGATGAGACGCTGG,\n...',
                    validation= lambda x: sequence_validation(x, None)
                ).props('rows=5 dense').classes('w-[280px] h-2fr font-mono')

                ui.element().classes("h-1")

            with ui.row(align_items='start').classes('w-full h-full p-0' + bc):
                with ui.row(align_items='center').classes('w-full h-[52px]'):
                    ui.button('upload').props('outline no-caps')

            # CONCENTRATION
            with ui.element().classes('w-full h-full p-0' + bc):
                with ui.row(align_items='start').classes('w-[280px] gap-0'):
                    ui.markdown('**RNP concentration**').classes('w-1/2 leading-[1.7]')
                    concentration_input = (ui.input(placeholder='100',
                                                    validation=concentration_validation)
                                           .props('dense suffix="nM"')
                                           .classes('w-1/2'))
                ui.element().classes("h-4")

            ui.element()

            # PARAMETER SELECTION
            with ui.row(align_items='center').classes('w-full p-0' + bc):
                ui.markdown('**Model parameters**').classes('leading-[0.7]')
                (ui.icon('info')
                .tooltip('Select the model for cleavage predictions.')
                .style('font-size: 20px'))
            ui.element()

            with ui.column().classes('w-full h-full p-0' + bc):
                parameter_set_input = ui.select(
                    options={
                        'sequence_params': 'sequence-params2 (recommended)',
                        'average_params': 'average-params',
                        'average_params_legacy': 'average-params-legacy'
                    },
                    value='sequence_params'
                ).props('dense').classes('w-[280px] p-0 m-0')

            with ui.row(align_items='center').classes('h-full w-full p-0' + bc):
                ui.button(icon='search').classes('h-4 w-4').style('font-size: 12px').props('outline')

            ui.element().classes("h-6")
            ui.element()

            submit_button = ui.button('Submit').props('icon=send').classes('w-[280px]')

@ui.page("/")
def index():
    show()

ui.run()
