import re
from io import StringIO

from nicegui import ui, events
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt

from crisprzip import *
from crisprzip.kinetics import *

initial_input = True  # auto-fills upon load - useful when developing


def show_input():

    def update_placeholder(e):
        val = e.value
        if val == 'protospacer':
            placeholder = 'GACGCATAAAGATGAGACGCTGG'
        else:
            placeholder = 'GACGCAUAAAGAUGAGACGC'
        target_sequence_input.props(f'placeholder={placeholder}')
        target_sequence_input.update()

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

    def sequence_validation(input, input_type=None) -> str:

        if input_type == "protospacer":
            length = 23
            pattern = r'^[ACGTacgt\,\n\s]*$'
            if not (re.fullmatch(pattern, input)):
                return "Only ACGT nucleotides"
            if length is not None and len(input.strip()) > 0:
                input_length = len(input.strip())
                if input_length < length:
                    return f"Too short: {input_length}/{length}"
                elif input_length > length:
                    return f"Too long: {input_length}/{length}"
                if input[-2:] != "GG":
                    return f"Only canonical PAMs \'NGG\'"

        elif input_type == "guide RNA":
            length = 20
            pattern = r'^[ACGUacgu\,\n\s]*$'
            if not (re.fullmatch(pattern, input)):
                return "Only ACGU nucleotides"
            if length is not None and len(input.strip()) > 0:
                input_length = len(input.strip())
                if input_length < length:
                    return f"Too short ({input_length}/{length})"
                elif input_length > length:
                    return f"Too long ({input_length}/{length})"

        elif input_type == "offtargets":
            off_targets = process_offtarget_input(input)
            for i, seq in enumerate(off_targets):
                output = sequence_validation(seq, "protospacer")
                if output is not None:
                    return f"{output} (target #{i + 1})"

    wc1 = 250  # column 1 width
    wc2 = 100  # column 2 width
    fsz = 10   # font size (in pt)
    fsi = 12   # font size for information bubble (in pt)
    fsb = 10   # font size for text in the information bubble (in pt)

    # CONTENT
    with ui.grid(columns=2).style(
            f'grid-template-columns: {wc1}px {wc2}px').classes('gap-0'):

        # TARGET SEQUENCE
        with ui.row(align_items='center').classes('p-0'):
            ui.markdown('**Target sequence**').classes(
                'p-0 leading-[0.7]').style(f'font-size: {fsz}pt')
            with ui.icon('info').style(f'font-size: {fsi}pt'):
                ui.tooltip(
                    'Sequence of the DNA target/protospacer (5\'-to-3\', '
                    '20 nts + PAM) or the gRNA (5\'-to-3\', 20 nts).'
                ).style(f'font-size: {fsb}pt')
        ui.element()

        with ui.column().classes('w-full p-0'):
            target_sequence_input = (
                ui.input(
                    validation=lambda x: sequence_validation(x, input_type=target_input_select.value),
                    value=('GACGCATAAAGATGAGACGCTGG' if initial_input else None),
                )
                .classes(f'w-[{wc1 - 20}px] font-mono')
                .props('dense')
                .style(f'font-size: {fsz}pt')
            )

        with ui.column().classes('w-full p-0'):
            target_input_select = (
                ui.select(['protospacer', 'guide RNA'],
                          value='protospacer',
                          on_change=update_placeholder)
                .classes('w-full').props('dense').style(f'font-size: {fsz}pt')
            )
            update_placeholder(type('E',(object,),{'value': target_input_select.value})())

        # OFF-TARGET SEQUENCES
        with ui.row(align_items='center').classes('p-0'):
            ui.markdown('**Off-target sequences**').classes(
                'p-0 leading-[0.7]').style(f'font-size: {fsz}pt')
            with ui.icon('info').style(f'font-size: {fsi}pt'):
                ui.tooltip(
                    'Sequences of potential DNA off-targets (5\'-to-3\', '
                    '20 nts + PAM). Optional.'
                ).style(f'font-size: {fsb}pt')
        ui.element()

        with ui.column().classes('w-full h-full p-0'):
            off_targets_input = ui.textarea(
                placeholder='GACGCATAAAGATGAGACGCTGG,\nGACGCATAAAGATGAGACGCTGG,\n...',
                validation=lambda x: sequence_validation(x, input_type="offtargets"),
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
                with ui.icon('info').style(f'font-size: {fsi}pt'):
                    ui.tooltip(
                        'Select the (most similar) application context. This choice '
                        'determines the binding and unbinding kinetics of Cas9 to '
                        'DNA. Click on \'show\' to inspect values.'
                    ).style(f'font-size: {fsb}pt')
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
                with ui.icon('info').style(f'font-size: {fsi}pt'):
                    ui.tooltip(
                        'Select the set of landscape parameters for CRISPRzip. Click on '
                        '\'show\' to inspect values.'
                    ).style(f'font-size: {fsb}pt')
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
                show_button = (
                    ui.button('show')
                    .style(f'font-size: {fsz}pt')
                    .props('outline no-caps')
                )

        submit_button = (
            ui.button('Submit')
            .props('icon=send')
            .classes(f'w-[{wc1 - 20}px]')
            .style(f'font-size: {fsz}pt')
        )

    with ui.dialog().props('full-width') as parameter_popup, ui.row(), ui.card().classes('w-[850px] mx-auto p-3'):
        with ui.column(align_items='end').classes('w-full'):
            (ui.button(icon='close', on_click=parameter_popup.close)
             .props('outline').classes('w-[40px] h-[40px] m-0'))

        ui.add_css('''
            .nicegui-scroll-area .q-scrollarea__content {
                padding: 0;
            }
        ''')

        with ui.scroll_area().classes('w-full h-[350px]'):
            with ui.column(align_items='center').classes('w-full'):
                parameter_plot_container = ui.element()
            ui.markdown('''
                The parameter values that determine the outcomes of a CRISPRzip simulation. In summary,
                
                - The **on-target landscape** $U_b$ describes the free energy of Cas9 target recognition on its perfect target;
                - The **mismatch penalties** $Q_b$ are added to that landscape when a mismatching basepair is encountered;
                - The **binding rate** $k_\mbox{on}^{1\mbox{ nM}}$ describes the rate at which 1 nM of Cas9 binds PAM sites;
                - The **unbinding rate** $k_\mbox{off}$ describes the rate at which PAM-bound Cas9 unbinds into solution (set by _context_ choice);
                - The **forward rate** $k_f$ describes the rate at which the R-loop is extended (set by _context_ choice);
                - The **cleavage rate** $k_\mbox{clv}$ describes the rate at which both strands are cut upon R-loop completion (b=20).
            
                With a sequence-average model ('average' or 'average-legacy'), these parameters fully describe the energy landscape on any
                _mismatch pattern_, e.g. mismatches at position 3 and 16. With sequence-specific model ('sequence'), the on-target landscape and 
                mismatch penalties above are only the protein contribution to the total energy, to which the cost of R-loop formation should still
                be added.
                
                For more details, see [Eslami-Mossalam (2022) Nat Comm](https://www.nature.com/articles/s41467-022-28994-2).
                
            ''', extras=['latex']
            )

    mpl.style.use('seaborn-v0_8')

    def plot_parameter_values():
        plot_size = (800, 250)  # in px
        dpi = plt.rcParams['figure.dpi']  # pixel in inches

        protein = crisprzip.kinetics.load_landscape(model_dropdown.value)
        k_on, k_off = get_k_on_off(context_dropdown.value)

        with (ui.matplotlib(figsize=(plot_size[0] / dpi,
                                     plot_size[1] / dpi))
                .classes(f'w-[{plot_size[0]}px] h-[{plot_size[1]}px]')
                .figure) as fig:

            axs = fig.subplots(1, 3, width_ratios=(3, 3, 2))
            for a in axs:
                a.set_facecolor('#ECF0F1')

            # on-target landscape
            axs[0].plot(
                np.arange(21), np.append(0, protein.on_target_landscape),
            )
            y_extent = np.diff(axs[0].get_ylim())
            axs[0].set_xticks(np.arange(0, 21, 5))
            axs[0].set_xlabel("R-loop size $b$")
            axs[0].set_ylabel("free energy $U_b$ ($k_BT$)")
            axs[0].set_title("(protein) on-target landscape")

            # mismatch penalties
            axs[1].plot(
                np.arange(21), np.append(np.nan, protein.mismatch_penalties),
            )
            axs[1].set_ylim(0, y_extent)
            axs[1].set_xticks(np.append(1, np.arange(5, 21, 5)))
            axs[1].set_xlabel("R-loop size $b$")
            axs[1].set_ylabel("free energy $Q_b$ ($k_BT$)")
            axs[1].set_title("(protein) mismatch penalties")

            # rates
            axs[2].plot(
                np.arange(4),
                [k_on,
                 k_off,
                 protein.internal_rates['k_f'],
                 protein.internal_rates['k_clv']],
                ls='', marker='o', markersize=5
            )
            axs[2].set_yscale('log')
            axs[2].set_xlim(-.5, 3.5)
            axs[2].set_xticks(
                np.arange(4),
                ['$k_{on}^{1\,nM}$',
                 '$k_{off}$',
                 '$k_{f}$',
                 '$k_{clv}$'],
            )
            axs[2].set_ylabel("rates $k$ ($s^{-1}$)")
            axs[2].set_title("rates")

            fig.subplots_adjust(left=0.1, right=.95,
                                bottom=.2, top=.85,
                                wspace=.6)

        return fig

    def show_button_handler():
        parameter_plot_container.clear()
        try:
            with parameter_plot_container:
                plot_parameter_values()
        except Exception as e:
            ui.notify(f'Error: {str(e)}', type='negative')

        parameter_popup.open()

    show_button.on_click(show_button_handler)

    def get_input_values():

        # check for valid target sequence
        err_msg = sequence_validation(input=target_sequence_input.value,
                                      input_type=target_input_select.value)
        if err_msg:
            ui.notify(f"Target sequence error: {err_msg}", type='negative')
            return

        # check for valid off-target sequences
        err_msg = sequence_validation(input=off_targets_input.value,
                                      input_type="offtargets")
        if err_msg:
            ui.notify(f"Off-target sequence error: {err_msg}", type='negative')
            return

        ontarget = process_ontarget_input(
            target_sequence_input.value,
            target_input_select.value
        )
        off_targets = process_offtarget_input(off_targets_input.value)

        input_vals = {
            'on_target': ontarget,
            'off_targets': off_targets,
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

