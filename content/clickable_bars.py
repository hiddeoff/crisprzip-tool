from nicegui import ui
import plotly.graph_objs as go

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

colors = ['blue'] * len(values)

selected_indices = set()

@ui.page('/')
def page():

    # --- Plotly ---
    def create_figure():
        fig = go.Figure()
        fig.add_trace(go.Bar(
            x=indices,
            y=values,
            marker=dict(color=colors)
        ))
        fig.update_layout(
            clickmode='event+select',
            dragmode=False,
            yaxis={'type': 'log'}
        )
        return fig

    def update_plot():
        plot.figure = create_figure()

    # --- Layout ---
    with ui.row(align_items='center').classes('w-full gap-0'):

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
                    {'select': False, 'index': i, 'sequence': targets[i], 'k_clv': values[i]}
                    for i in range(len(targets))
                ],
                'rowSelection': 'multiple',
            })

            # Working event trigger
            grid.on('selectionChanged', lambda event: ui.notify("Selection changed!"))  # works

            async def output_selected_rows():
                ui.notify('You pressed the button!')
                rows = await grid.get_selected_rows()
                ui.notify('Received selected rows')
                if not rows:
                    ui.notify("No row(s) selected!")
                else:
                    ui.notify(
                        f"You selected rows {', '.join([str(r['index']) for r in rows])}")

            ui.button('Select #2', on_click=lambda: grid.run_row_method(2, 'setSelected', True))  # Select 3rd row

            ui.button('Output selected rows', on_click=output_selected_rows)

        plot = ui.plotly(create_figure()).classes('w-1/2')

ui.run()
