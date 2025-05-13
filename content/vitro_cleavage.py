import re

from nicegui import ui, events
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from io import StringIO
import pandas as pd

from crisprzip import *
from crisprzip.kinetics import *


initial_input = False  # auto-fills upon load - useful when developing


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
            off_targets_input.value = "" # clear the input field
            # read the file content
            content = e.content.read().decode("utf-8")
            
            # parse the csv (one column comma separated)
            if e.name.lower().endswith('.csv'):
                with StringIO(content) as f:
                    df = pd.read_csv(f)
                    sequences = df.iloc[:, 0].tolist()
                    formatted_sequences = ',\n'.join(sequences)
                    off_targets_input.value= formatted_sequences
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

    def concentration_validation(input) -> str| None:
        concentration = str(input)
        if not concentration.strip():
            return "RNP concentration cannot be empty."
        if not re.fullmatch(r"^\d*(\.\d*)?([eE][+-]?\d+)?$", concentration):
            return "Please enter a valid number."
        return None  # Return None if the input is valid

    wc1 = 230  # column 1 width
    wc2 = 100  # column 2 width
    fsz = 10   # font size (in pt)

    # CONTENT
    with ui.grid(columns=2).style(f'grid-template-columns: {wc1}px {wc2}px').classes('gap-0'):

        # TARGET SEQUENCE
        with ui.row(align_items='center').classes('p-0'):
            ui.markdown('**Target sequence**').classes('p-0 leading-[0.7]').style(f'font-size: {fsz}pt')
            ui.icon('info').tooltip('Select the model for cleavage predictions. Recommended: sequence-params2.').style(f'font-size: {fsz}pt')
        ui.element()

        with ui.column().classes('w-full p-0'):
            target_sequence_input = ui.input(
                validation=lambda x: sequence_validation(x, 23),
                # once the update_placeholder() function works, this arg should be omitted
                placeholder='GACGCATAAAGATGAGACGCTGG',
                value=('GACGCATAAAGATGAGACGCTGG' if initial_input else None),
            ).classes(f'w-[{wc1 - 20}px] font-mono').props('dense').style(f'font-size: {fsz}pt')

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
            ui.markdown('**Off-target sequences**').classes('p-0 leading-[0.7]').style(f'font-size: {fsz}pt')
            ui.icon('info').tooltip('Select the model for cleavage predictions. Recommended: sequence-params2.').style(f'font-size: {fsz}pt')
        ui.element()

        with ui.column().classes('w-full h-full p-0'):
            off_targets_input = ui.textarea(
                placeholder='GACGCATAAAGATGAGACGCTGG,\nGACGCATAAAGATGAGACGCTGG,\n...',
                validation=lambda x: sequence_validation(x, None),
                value=(('GACGAACAAAGATGAGACGCTGG,\n' +
                        'GACGCATATATACGAGACGCTGG,\n' +
                        'GACGCATAATTATGAGTCGCTGG,\n' +
                        'GACGCATACCGATGTGTCGCTGG,\n' +
                        'GACGCATAAAGATGGGGCTCTGG')
                       if initial_input else None) ,
            ).props('rows=5 dense').classes(f'w-[{wc1 - 20}px] h-2fr font-mono').style(f'font-size: {fsz}pt')

        with ui.row(align_items='start').classes('w-full h-full p-0'):
            with ui.row(align_items='center').classes('w-[52px] h-[52px]'):
                upload_component = ui.upload(
                    on_upload=handle_offtarget_uploads,
                    on_rejected=lambda e: ui.notify('File upload failed', type='warning'),
                    auto_upload=True
                ).props("accept=.csv hide-upload-btn").classes('hidden') # can add .txt (for example) to .props if you want to uplaod something other than .csv files
                # add custom ui button (instead of the regular upload button)
                ui.button('Upload', on_click=lambda: (
                    upload_component.reset(),
                    upload_component.run_method('pickFiles')
                )).props('outline no-caps').style(f'font-size: {fsz}pt')
                
        # CONCENTRATION
        with ui.element().classes('w-full h-full p-0'):
            with ui.row(align_items='start').classes(f'w-[{wc1 - 20}px] gap-0'):
                ui.markdown('**RNP concentration**').classes('leading-[1.7]').style(f'font-size: {fsz}pt')
                ui.space()
                rnp_concentration_input = (
                    ui.input(placeholder='100',
                             validation=concentration_validation,
                             value=(100 if initial_input else None))
                    .props('dense suffix="nM"')
                    .classes('w-[70px]')
                    .style(f'font-size: {fsz}pt')
                )

        ui.element()

        # PARAMETER SELECTION
        with ui.row(align_items='center').classes('w-full p-0'):
            ui.markdown('**Model parameters**').classes('leading-[0.7]').style(f'font-size: {fsz}pt')
            (ui.icon('info')
             .tooltip('Select the model for cleavage predictions.')
             .style(f'font-size: {fsz}pt'))
        ui.element()

        with ui.column().classes('w-full h-full p-0'):
            model_dropdown = ui.select(
                options={
                    'sequence_params': 'sequence-params2 (default)',
                    'average_params': 'average-params',
                    'average_params_legacy': 'average-params-legacy'
                },
                value='sequence_params'
            ).props('dense').classes(f'w-[{wc1 - 20}px] p-0 m-0').style(f'font-size: {fsz}pt')

        with ui.row(align_items='center').classes('h-full w-full p-0'):
            (ui.button(icon='search')
             .classes('h-4 w-4').style('font-size: 12px').props('outline'))

        ui.element().classes("h-6")
        ui.element()

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

		# Validate concentration input
        concentration = rnp_concentration_input.value
        concentration_error = concentration_validation(concentration)
        if concentration_error:
            message = "RNP concentration cannot be empty and has to be a valid number. Default value is 100."
            ui.notify(f'Error: {message}', type='negative')
            return  # Stop further processing if there's an error
        else:
            concentration = float(concentration)

        input_vals = {
            'on_target': ontarget,
            'off_targets': process_offtarget_input(off_targets_input.value),
            'concentration': float(rnp_concentration_input.value),
            'parameter_set': model_dropdown.value
        }
        return input_vals

    return submit_button, get_input_values


def show_output(output_container, get_input_values: callable):

    # PROCESSING USER INPUT
    def make_stc_list(protospacer, off_targets, parameter_set):
        """Generate SearcherTargetComplexes."""
        bare_protein = crisprzip.kinetics.load_landscape(parameter_set)
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
            off_targets=off_targets,
            parameter_set=parameter_set,
        )
        k_on_ref = 1E-2
        binding_rate = k_on_ref * concentration

        k_fit_values = [
            get_cleavage_rate(stc, binding_rate) for
            stc in protein_sequence_complexes
        ]
        return k_fit_values

    protospacer, off_targets, concentration, parameter_set = get_input_values().values()

    if len(off_targets) > 25:
        ui.notify("Can't process more than 25 off-targets at once", type='warning')
        return

    protein_sequence_complexes = make_stc_list(
        protospacer=protospacer,
        off_targets=off_targets,
        parameter_set=parameter_set,
    )
    cleavage_rates = get_all_cleavage_rates(
        protospacer=protospacer,
        off_targets=off_targets,
        parameter_set=parameter_set,
        concentration=concentration,
    )

    targets = [protospacer] + off_targets
    values = cleavage_rates
    indices = np.arange(len(targets))

    # VISUALIZATION
    mpl.style.use('seaborn-v0_8')

    output_container.clear()
    with output_container:
        with ui.column(align_items='center').classes('w-full'):
            plot_row = ui.row(align_items='start').classes('gap-0')
        selection_container = ui.element('div').classes('w-full')  # determine width!
        showing_selection = False

    with plot_row:

        # Table
        with ui.column(align_items='center').classes('w-[360px]'):

            def to_sci_html(val):
                scistr = f"{val:.2e}"
                val = scistr[:4]
                pow = str(int(scistr[5:]))
                sci_html = f"{val} &middot 10<sup>{pow}</sup>"
                return sci_html

            ui.add_head_html('''
            <style>
                .ag-cell.monospace-column {
                    font-family: monospace !important;
                    font-size: 13px;
                }
            </style>
            ''')
            default_coldefs = {'suppressMovable': True, 'sortable': False, 'resizable': False}
            column_defs = [
                dict(cd, **default_coldefs) for cd in [
                    {'headerName': '', 'field': 'select', 'width': '40',
                     'checkboxSelection': True,
                     'headerCheckboxSelection': True},
                    {'headerName': '#', 'field': 'index', 'width': '40'},
                    {'headerName': 'sequence', 'field': 'sequence',
                     'width': '220', 'cellClass': 'monospace-column'},
                    {'headerName': 'k_clv (s⁻¹)', 'field': 'k_clv',
                     'width': '90'}
                ]]

            grid = ui.aggrid({
                'columnDefs': column_defs,
                'rowData': [{'index': i, 'sequence': targets[i], 'k_clv': to_sci_html(values[i])}
                            for i in range(len(targets))],
                'rowSelection': 'multiple'
            }, html_columns=[3], auto_size_columns=True)

            async def get_selected_ids():
                selection = await grid.get_selected_rows()
                selected_ids = [r['index'] for r in selection]
                return sorted(selected_ids)

            async def sort_grid_kclv():
                nonlocal grid
                selected_ids = await get_selected_ids()
                grid.options['rowData'] = [
                    {'index': i, 'sequence': targets[i], 'k_clv': to_sci_html(values[i])}
                    for i in reversed(np.argsort(values))
                ]
                for i in range(len(values)):
                    if np.argsort(values)[::-1][i] in selected_ids:
                        grid.run_row_method(i, 'setSelected', True)
                grid.update()

            async def sort_grid_index():
                nonlocal grid
                selected_ids = await get_selected_ids()
                grid.options['rowData'] = [
                    {'index': i, 'sequence': targets[i], 'k_clv': to_sci_html(values[i])}
                    for i in range(len(targets))
                ]
                for i in selected_ids:
                    grid.run_row_method(i, 'setSelected', True)
                grid.update()

            with ui.row(align_items='center').classes('w-full'):
                show_button = ui.button("SHOW").classes('w-[100px]')

                ui.icon('info').tooltip(
                    'Select up to 9 rows with checkboxes or Ctrl/Shift to inspect their'
                    'hybridization landscape, cleavage dynamics and concentration dependence.'
                )
                ui.space()
                sort_button = ui.button().props('no-caps').classes("w-[120px]")
                with sort_button:
                    ui.html("sort by <i>k<sub>clv</sub></i>")
                download_button = ui.button().props('icon=download no-caps outline').classes('w-[1em] h-[1em]')


        ui.element().classes('w-[15px]')

        # Plot
        dpi = plt.rcParams['figure.dpi']  # pixel in inches
        with ui.matplotlib(figsize=(375 / dpi, 250 / dpi)).classes('w-[375px] h-[250px]').figure as fig:
            ax = fig.gca()
            ax.bar(indices, values, width=.8, align='center',
                   color="#5898d4", alpha=.8)
            ax.set_yscale('log')
            ax.set_title("Cleavage rate $k_{clv}$ ($s^{-1}$)")
            ax.set_facecolor('#ECF0F1')
            ax.grid(axis='x')
            ax.set_xticks(indices)
            if len(values) > 16:
                ax.tick_params(axis='both', labelsize='x-small')
            fig.tight_layout()

            async def grid_selection_handler():
                selected_ids = await get_selected_ids()
                if not selected_ids and showing_selection:
                    show_button.set_text("clear")
                else:
                    show_button.set_text("show")
                await highlight_selected_bars()

            async def highlight_selected_bars():
                selected_ids = await get_selected_ids()
                with fig:
                    if selected_ids:
                        for i, bar in enumerate(ax.patches):
                            bar.set_alpha(.9 if i in selected_ids else .4)
                    else:
                        for i, bar in enumerate(ax.patches):
                            bar.set_alpha(.6)

            def sort_plot_kclv():
                sorted_positions = np.argsort(values)[::-1]
                with fig:
                    for xpos, i in enumerate(reversed(np.argsort(values))):
                        bar = ax.patches[i]
                        bar.set_x(xpos - .4)
                    ax.set_xticks(indices, indices[sorted_positions])

            def sort_plot_index():
                with fig:
                    for i in range(len(values)):
                        bar = ax.patches[i]
                        bar.set_x(i - .4)
                    ax.set_xticks(indices, indices)

            grid.on('selectionChanged', grid_selection_handler)

    sorted_kclv = False

    async def handle_sort_click():
        nonlocal sorted_kclv

        if not sorted_kclv:
            await sort_grid_kclv()
            sort_plot_kclv()
            sort_button.clear()
            with sort_button:
                ui.html("sort by index")
            sorted_kclv = True

        else:
            await sort_grid_index()
            sort_plot_index()
            sort_button.clear()
            with sort_button:
                ui.html("sort by <i>k<sub>clv</sub></i>")
            sorted_kclv = False

    sort_button.on_click(handle_sort_click)

    def download_grid():
        df = pd.DataFrame({'sequence': targets, 'k_clv [1/s]': values})
        csv_string = df.to_csv(index=True)
        ui.download.content(csv_string, 'crisprzip_kclv.csv')

    download_button.on_click(download_grid)

    async def handle_show_click():
        nonlocal showing_selection

        selected_ids = await get_selected_ids()
        if selected_ids:
            if len(selected_ids) > 9:
                ui.notify("Select at most 9 targets for inspection.", type='warning')
            else:
                await plot_selection()

        elif not showing_selection:
            ui.notify("No targets selected to show.", type='warning')

        else:
            selection_container.clear()
            showing_selection = False
            show_button.set_text("show")

    async def plot_selection():
        nonlocal showing_selection
        selected_ids = await get_selected_ids()
        showing_selection = True

        try:
            k_on_ref = 1E-2
            binding_rate = k_on_ref * concentration

            # Clear previous output content - [IMPORTANT], otherwise plots will stack below on each request
            selection_container.clear()

            with selection_container:
                with ui.column(align_items='center'):
                    dpi = plt.rcParams['figure.dpi']  # pixel in inches
                    with ui.matplotlib(figsize=(750 / dpi, 600 / dpi)).classes('w-[750px] h-[600px]').figure as fig1:
                        spec0 = fig.add_gridspec(
                            ncols=2, nrows=2,
                            wspace=.5, hspace=.6,
                        )

                        ax0 = fig1.add_subplot(spec0[0, 0])
                        ax0.set_facecolor('#ECF0F1')
                        for i in selected_ids:
                            stc = protein_sequence_complexes[i]
                            landscape = stc._get_off_target_landscape()
                            sol_stab = np.log(
                                k_on_ref * concentration / stc.internal_rates['k_off']
                            )
                            ax0.plot(
                                np.arange(22) - 1,
                                np.concatenate([np.array([sol_stab, 0]), landscape]),
                                label=('target (#0)' if i == 0 else f'off-target #{i}'),
                                zorder=5 - .1 * i,
                            )
                        ax0.set_xlabel('R-loop length $b$')
                        ax0.set_xticks(
                            [-1, 0, 1, 5, 10, 15, 20],
                            ['S', 'P', '1', '5', '10', '15', '20'],
                        )
                        ax0.grid(axis='x')
                        ax0.set_title("R-loop landscape")

                        lims0 = ax0.get_ylim()
                        ax0.vlines([0, 5, 10, 15, 20], -50, 50,
                                   color='white', zorder=0, lw=.8)
                        ax0.set_ylim(*lims0)
                        ax0.set_ylabel(r"free energy $\Delta U_b$ ($k_BT$)")
                        ax0.set_facecolor('#ECF0F1')
                        ax0.legend(bbox_to_anchor=(1.1, 1.05))

                        # FIGURE 1 - Cleaved fraction vs time
                        ax1 = fig1.add_subplot(spec0[1, 0])
                        dt = np.logspace(-1, 6)
                        for i in selected_ids:
                            stc = protein_sequence_complexes[i]
                            f_clv = stc.get_cleaved_fraction(dt, binding_rate)
                            ax1.plot(
                                dt, f_clv,
                                label=('target (#0)' if i == 0 else f'off-target #{i}'),
                                zorder=5 - .1* i,
                            )
                        ax1.set_xscale('log')
                        ax1.set_xlabel('time $t$ (s)')
                        ax1.set_xlim(dt.min(), dt.max())

                        ax1.set_ylim(-.1, 1.1)
                        ax1.set_ylabel("fraction cleaved $f_{clv}$")

                        ax1.set_facecolor('#ECF0F1')
                        ax1.set_title("cleavage vs time")

                        # FIGURE 2 - Cleavage rate vs concentration
                        ax2 = fig1.add_subplot(spec0[1, 1])

                        # Concentration range, currently hardcoded and needs a better solution
                        conc_logmin = -2  # nM
                        conc_logmax = 3  # nM
                        dc = np.logspace(conc_logmin, conc_logmax)

                        for i in selected_ids:
                            stc = protein_sequence_complexes[i]
                            k_fit = [get_cleavage_rate(stc, k_on_ref * c)
                                     for c in dc]
                            ax2.plot(
                                dc, k_fit,
                                label=f'{"on" if i==0 else "off"}-target #{i}',
                                zorder=5 - .1* i
                            )
                            ax2.plot(
                                concentration, get_cleavage_rate(stc, k_on_ref * concentration),
                                marker='o',
                                zorder=5 - .1 * i,
                                color=ax2.lines[-1].get_color()
                            )
                        ax2.set_xscale('log')
                        ax2.set_xlabel('RNP concentration $c$ (nM)')
                        ax2.set_xlim(dc.min(), dc.max())

                        ax2.set_yscale('log')
                        ax2.set_ylabel("cleavage rate $k_{clv}$ ($s^{-1}$)")

                        ax2.set_facecolor('#ECF0F1')
                        ax2.set_title("cleavage vs concentration")

        except Exception as e:
            ui.notify(f'Error: {str(e)}', type='negative')

    show_button.on_click(handle_show_click)



def show_contents():
    with ui.row().classes('w-full h-full no-wrap'):

        # INPUT
        with ui.card().classes('p-4 m-2'):
            submit_button, get_input_values = show_input()

        # OUTPUT
        output_container = ui.column().classes('w-full h-full no-wrap m-2')
        submit_button.on_click(
            lambda: show_output(output_container, get_input_values)
        )


@ui.page("/")
def index():
    show_contents()


ui.run()
