from matplotlib.patches import FancyBboxPatch
from nicegui import ui
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import crisprzip

mpl.style.use('seaborn-v0_8')


# --- Data ---
protospacer = "AGCTAAGCTAAGCTAAGCTAGGG"
off_targets = ["AGCTAAGCTATGCTAAGCTAGGG",
               "AGCTAAGCTAAGGTGAGCTAGGG",
               "AGCTAAGCCCCGCTAAGCTAGGG"]

indices = list(range(4))
targets = [protospacer] + off_targets
parameter_set = "sequence_params"
concentration = 100.


@ui.page('/')
def page():

    # SUPPORT FUNCTIONS + CAS9 CALCULATIONS
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
    values = cleavage_rates


    # --- Layout ---
    with ui.row(align_items='center').classes('w-full gap-0'):

        # AG Grid
        with ui.column(align_items='center').classes('w-1/2'):
            grid = ui.aggrid({
                'columnDefs': [
                    {'headerName': '', 'field': 'select', 'width': '25',
                     'checkboxSelection': True, 'headerCheckboxSelection': True},
                    {'headerName': '#', 'field': 'index', 'width': '30'},
                    {'headerName': 'sequence', 'field': 'sequence', 'width': '150'},
                    {'headerName': 'k_clv', 'field': 'k_clv'},
                ],
                'rowData': [
                    {'index': i, 'sequence': targets[i], 'k_clv': values[i]}
                    for i in range(len(targets))
                ],
                'rowSelection': 'multiple',
            })

        async def get_selected_ids():
            selection = await grid.get_selected_rows()
            selected_ids =  [r['index'] for r in selection]
            return sorted(selected_ids)

        async def print_selected_ids():
            selected_ids = await get_selected_ids()
            ui.notify(f"Selection: {selected_ids}")

        # Plot
        with ui.matplotlib(figsize=(4, 3)).figure as fig:
            ax = fig.gca()
            ax.bar(indices, values, align='center',
                   color="#5898d4", alpha=.8)
            ax.set_yscale('log')
            ax.set_ylabel("$k_{clv}$ ($s^{-1}$)")
            ax.set_facecolor('#ECF0F1')
            ax.grid(axis='x')
            fig.tight_layout()

            async def highlight_selected_bars():
                selected_ids = await get_selected_ids()

                with fig:
                    if selected_ids:
                        for i, bar in enumerate(ax.patches):
                            bar.set_alpha(.9 if i in selected_ids else .4)
                    else:
                        for i, bar in enumerate(ax.patches):
                            bar.set_alpha(.6)
                ui.update(fig)

            grid.on('selectionChanged', highlight_selected_bars)

    show_button = ui.button("Show selected targets")
    output_container = ui.element('div').classes('w-[1200px]')  # determine width!

    async def plot_selection():
        selected_ids = await get_selected_ids()

        try:
            k_on_ref = 1E-2
            binding_rate = k_on_ref * concentration

            # Clear previous output content - [IMPORTANT], otherwise plots will stack below on each request
            output_container.clear()

            with output_container:
                dpi = plt.rcParams['figure.dpi']  # pixel in inches
                with ui.matplotlib(figsize=(1200/dpi, 250/dpi)).classes('border').figure as fig1:
                    spec1 = fig.add_gridspec(
                        ncols=3, nrows=1, right=.8, wspace=.5, bottom=0.2
                    )
                    spec2 = fig.add_gridspec(
                        ncols=1, nrows=1, left=.825,
                    )

                    ax0 = fig1.add_subplot(spec1[0, 0])
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

                    # FIGURE 1 - Cleaved fraction vs time
                    ax1 = fig1.add_subplot(spec1[0, 1])
                    dt = np.logspace(-1, 6)
                    for i in selected_ids:
                        stc = protein_sequence_complexes[i]
                        f_clv = stc.get_cleaved_fraction(dt, binding_rate)
                        ax1.plot(
                            dt, f_clv,
                            label=('on-target (#0)' if i == 0 else f'off-target #{i}'),
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
                    ax2 = fig1.add_subplot(spec1[0, 2])

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
                    # ax2.legend()
                    ax2.set_title("cleavage vs concentration")

                    ax3 = fig1.add_subplot(spec2[0, 0])
                    h, l = ax1.get_legend_handles_labels()
                    ax3.grid('off')
                    ax3.set_axis_off()
                    ax3.set_facecolor('None')
                    ax3.legend(h, l, loc='center left')

        except Exception as e:
            ui.notify(f'Error: {str(e)}', type='negative')

    show_button.on_click(plot_selection)

ui.run()
