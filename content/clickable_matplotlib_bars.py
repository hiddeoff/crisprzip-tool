from matplotlib.patches import FancyBboxPatch
from nicegui import ui
import numpy as np
import matplotlib as mpl

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
values = [5.521E-1, 5.231E-2, 1.889E-2, 5.224E-5]

@ui.page('/')
def page():

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

    ui.button("Get selected ids", on_click=print_selected_ids)


ui.run()
