import re
from io import StringIO

import numpy as np
from nicegui import ui, events
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from crisprzip import *
from crisprzip.kinetics import *
from .input import show_input, get_k_on_off, make_stc_list

def get_effective_stab(stc):
    landscape = stc.off_target_landscape
    boltzmann = np.exp(-landscape)
    eff_stab  = np.sum(landscape * boltzmann) / np.sum(boltzmann)
    return eff_stab


def get_all_effective_stabs(protospacer, off_targets,
                            context, parameter_set):
    protein_sequence_complexes = make_stc_list(
        protospacer=protospacer,
        off_targets=off_targets,
        context=context,
        parameter_set=parameter_set,
    )
    u_eff_values = [get_effective_stab(stc) for stc in
                    protein_sequence_complexes]
    return u_eff_values


def get_binding_const(stc, k_on_ref):
    """Calculate cleavage rate."""
    dc = np.logspace(-2, 6)
    f_bnd = stc.get_bound_fraction(time=3600., binding_rate=dc * k_on_ref)
    with np.errstate(over='ignore'):  # ignore RuntimeWarning: overflow
        kd = np.exp(curve_fit(
            f=lambda c, logkd: c / (np.exp(logkd) + c),
            xdata=dc,
            ydata=f_bnd,
        )[0][0])
    return kd


def get_all_binding_const(protospacer, off_targets,
                          context, parameter_set):
    protein_sequence_complexes = make_stc_list(
        protospacer=protospacer,
        off_targets=off_targets,
        context=context,
        parameter_set=parameter_set,
    )
    k_on, k_off = get_k_on_off(context)
    k_fit_values = [
        get_binding_const(stc, k_on) for
        stc in protein_sequence_complexes
    ]
    return k_fit_values


def show_output(output_container, get_input_values: callable):
    protospacer, off_targets, context, parameter_set = get_input_values().values()

    if len(off_targets) > 250:
        ui.notify(
            f"Can't process {len(off_targets)} off-targets at once! (max. 250)",
            type='warning')
        return

    protein_sequence_complexes = make_stc_list(
        protospacer=protospacer,
        off_targets=off_targets,
        context=context,
        parameter_set=parameter_set,
    )
    effective_stabs = get_all_effective_stabs(
        protospacer=protospacer,
        off_targets=off_targets,
        context=context,
        parameter_set=parameter_set,
    )

    targets = [protospacer] + off_targets
    values = effective_stabs
    indices = np.arange(len(targets))

    # VISUALIZATION
    mpl.style.use('seaborn-v0_8')

    output_container.clear()
    with output_container:
        with ui.column(align_items='center').classes('w-full'):
            plot_row = ui.row(align_items='start').classes('gap-0')
        selection_container = ui.element('div').classes(
            'w-full')  # determine width!
        showing_selection = False

    with plot_row:

        # Table
        with ui.column(align_items='center').classes('w-[380px]'):

            ui.add_head_html('''
            <style>
                .ag-cell.monospace-column {
                    font-family: monospace !important;
                    font-size: 13px;
                }
            </style>
            ''')
            default_coldefs = {'suppressMovable': True, 'sortable': False,
                               'resizable': False}
            column_defs = [
                dict(cd, **default_coldefs) for cd in [
                    {'headerName': '', 'field': 'select', 'width': '40',
                     'checkboxSelection': True,
                     'headerCheckboxSelection': True},
                    {'headerName': '#', 'field': 'index', 'width': '40'},
                    {'headerName': 'sequence', 'field': 'sequence',
                     'width': '220', 'cellClass': 'monospace-column'},
                    {'headerName': 'ΔU (kBT)', 'field': 'u_eff',
                     'width': '90'}
                ]]

            grid = ui.aggrid({
                'columnDefs': column_defs,
                'rowData': [{'index': i, 'sequence': targets[i],
                             'u_eff': f"{values[i]:.2f}"}
                            for i in range(len(targets))],
                'rowSelection': 'multiple'
            }, html_columns=[3], auto_size_columns=True)

            async def get_selected_ids():
                selection = await grid.get_selected_rows()
                selected_ids = [r['index'] for r in selection]
                return sorted(selected_ids)

            async def sort_grid_ueff():
                nonlocal grid
                selected_ids = await get_selected_ids()
                grid.options['rowData'] = [
                    {'index': i, 'sequence': targets[i],
                     'u_eff': f"{values[i]:.2f}"}
                    for i in np.argsort(values)
                ]
                for i in range(len(values)):
                    if np.argsort(values)[::-1][i] in selected_ids:
                        grid.run_row_method(i, 'setSelected', True)
                grid.update()

            async def sort_grid_index():
                nonlocal grid
                selected_ids = await get_selected_ids()
                grid.options['rowData'] = [
                    {'index': i, 'sequence': targets[i],
                     'u_eff': f"{values[i]:.2f}"}
                    for i in range(len(targets))
                ]
                for i in selected_ids:
                    grid.run_row_method(i, 'setSelected', True)
                grid.update()

            with ui.row(align_items='center').classes('w-full'):
                show_button = ui.button("SHOW").classes('w-[100px]')

                with ui.icon('info').style(f'font-size: 12pt'):
                    ui.tooltip(
                        'Select up to 9 rows with checkboxes or Ctrl/Shift to inspect their '
                        'hybridization landscape, cleavage dynamics and concentration dependence.'
                    ).style(f'font-size: 10pt')
                ui.space()
                sort_button = ui.button().props('no-caps').classes("w-[120px]")
                with sort_button:
                    ui.html("sort by <i>&Delta;U</i>")
                download_button = ui.button().props(
                    'icon=download no-caps outline').classes('w-[1em] h-[1em]')

        ui.element().classes('w-[15px]')

        # Plot
        with ui.column(align_items='center').classes('gap-0 p-0'):
            dpi = plt.rcParams['figure.dpi']  # pixel in inches

            # title
            with ui.matplotlib(figsize=(365 / dpi, 25 / dpi)).classes(
                    'w-[365px] h-[25px]').figure as tfig:
                t_ax = tfig.gca()
                t_ax.text(.5, .5, r"effective stability $-\Delta U_{eff}$ ($k_BT$)",
                          va='center', ha='center', fontsize='large')
                # tfig.tight_layout()
                t_ax.set_axis_off()

            with ui.row().classes('gap-0 p-0'):

                # y axis
                fig0 = ui.matplotlib(figsize=(40 / dpi, 250 / dpi)).classes(
                    'w-[40px] h-[250px]').figure

                ui.add_css('''
                    .nicegui-scroll-area .q-scrollarea__content {
                        padding: 0;
                    }
                    ''')

                # plot contents
                plotwidth = max(325, 25 * len(targets))
                with ui.scroll_area().classes(f'w-[325px] h-[250px]'):
                    with ui.matplotlib(
                            figsize=(plotwidth / dpi, 250 / dpi)).classes(
                            f'w-[{plotwidth}px] h-[250px]').figure as fig:
                        ax = fig.gca()
                        ax.bar(indices, 20 - np.array(values),
                               bottom=-20, width=.8, align='center',
                               color="#5898d4", alpha=.8)
                        ax.set_ylim(
                            -max(values) - .2 * (max(values) - min(values)),
                            -min(values) + .2 * (max(values) - min(values)),
                        )
                        ax.set_facecolor('#ECF0F1')
                        ax.grid(axis='x')
                        ax.set_xticks(indices)
                        x_margin = (min(15,
                                        len(targets)) - 1) / 15  # margin of 0-1
                        ax.set_xlim(-.5 - x_margin,
                                    indices[-1] + .5 + x_margin)
                        ax.get_yaxis().set_ticklabels([])
                        fig.subplots_adjust(left=0., right=1., bottom=.15,
                                            top=.95)

                with fig0:
                    ax0 = fig0.gca()
                    ax0.set_facecolor('None')
                    ax0.set_yticks(ax.get_yticks())
                    ax0.set_ylim(*ax.get_ylim())
                    ax0.set_xlim(0, 0)
                    ax0.get_xaxis().set_visible(False)
                    fig0.subplots_adjust(left=.99, right=1., bottom=.15,
                                         top=.95)

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

            def sort_plot_ueff():
                sorted_positions = np.argsort(values)
                with fig:
                    for xpos, i in enumerate(np.argsort(values)):
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

    sorted_ueff = False

    async def handle_sort_click():
        nonlocal sorted_ueff

        if not sorted_ueff:
            await sort_grid_ueff()
            sort_plot_ueff()
            sort_button.clear()
            with sort_button:
                ui.html("sort by index")
            sorted_ueff = True

        else:
            await sort_grid_index()
            sort_plot_index()
            sort_button.clear()
            with sort_button:
                ui.html("sort by <i>&Delta;U</i>")
            sorted_ueff = False

    sort_button.on_click(handle_sort_click)

    def download_grid():
        df = pd.DataFrame({'sequence': targets, 'k_clv [1/s]': values})
        csv_string = df.to_csv(index=True)
        ui.download.content(csv_string, 'crisprzip_u_eff.csv')

    download_button.on_click(download_grid)

    async def handle_show_click():
        nonlocal showing_selection

        selected_ids = await get_selected_ids()
        if selected_ids:
            if len(selected_ids) > 9:
                ui.notify("Select at most 9 targets for inspection.",
                          type='warning')
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
            k_on, k_off = get_k_on_off(context)
            concentration = 100
            binding_rate = k_on * concentration

            # Clear previous output content - [IMPORTANT], otherwise plots will stack below on each request
            selection_container.clear()

            with selection_container:
                with ui.column(align_items='center'):
                    dpi = plt.rcParams['figure.dpi']  # pixel in inches
                    with ui.matplotlib(figsize=(750 / dpi, 600 / dpi)).classes(
                            'w-[750px] h-[600px]').figure as fig1:
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
                                k_on * concentration / stc.internal_rates['k_off']
                            )
                            ax0.plot(
                                np.arange(22) - 1,
                                np.concatenate(
                                    [np.array([sol_stab, 0]), landscape]),
                                label=(
                                    'target (#0)' if i == 0 else f'off-target #{i}'),
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

                        # FIGURE 1 - Bound fraction vs time
                        ax1 = fig1.add_subplot(spec0[1, 0])
                        dt = np.logspace(-1, 6)
                        for i in selected_ids:
                            stc = protein_sequence_complexes[i]
                            f_bnd = stc.get_bound_fraction(dt, binding_rate,
                                                           pam_inclusion=0)
                            ax1.plot(
                                dt, f_bnd,
                                label=(
                                    'target (#0)' if i == 0 else f'off-target #{i}'),
                                zorder=5 - .1 * i,
                            )
                        ax1.set_xscale('log')
                        ax1.set_xlabel('time $t$ (s)')
                        ax1.set_xlim(dt.min(), dt.max())

                        ax1.set_ylim(-.1, 1.1)
                        ax1.set_ylabel("fraction bound $f_{bnd}$")

                        ax1.set_facecolor('#ECF0F1')
                        ax1.set_title("binding vs time (100 nM)")

                        # FIGURE 2 - Bound fraction vs concentration
                        ax2 = fig1.add_subplot(spec0[1, 1])

                        # Concentration range, currently hardcoded and needs a better solution
                        conc_logmin = -2  # nM
                        conc_logmax = 3  # nM
                        dc = np.logspace(conc_logmin, conc_logmax)

                        for i in selected_ids:
                            stc = protein_sequence_complexes[i]
                            f_bnd = stc.get_bound_fraction(3600, k_on * dc, pam_inclusion=0)
                            ax2.plot(
                                dc, f_bnd,
                                label=f'{"on" if i == 0 else "off"}-target #{i}',
                                zorder=5 - .1 * i
                            )
                            ax2.plot(
                                dc[np.searchsorted(dc, concentration)],
                                f_bnd[np.searchsorted(dc, concentration)],
                                marker='o',
                                zorder=5 - .1 * i,
                                color=ax2.lines[-1].get_color()
                            )
                        ax2.set_xscale('log')
                        ax2.set_xlabel('RNP concentration $c$ (nM)')
                        ax2.set_xlim(dc.min(), dc.max())

                        # ax2.set_yscale('log')
                        ax2.set_ylabel("bound fraction $f_{bnd}$")

                        ax2.set_facecolor('#ECF0F1')
                        ax2.set_title("binding vs concentration (1 hr)")

        except Exception as e:
            ui.notify(f'Error: {str(e)}', type='negative')

    show_button.on_click(handle_show_click)


def show_contents():
    with ui.row().classes('w-full h-full no-wrap'):
        # INPUT
        with ui.card().classes('p-4 m-2'):
            submit_button, get_input_values, model_dropdown = show_input()

        # OUTPUT
        output_container = ui.column().classes('w-full h-full no-wrap m-2')
        submit_button.on_click(
            lambda: show_output(output_container, get_input_values)
        )
