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
        title='Click bars or table rows to select/deselect',
        dragmode=False,
        yaxis={'type': 'log'}
    )
    return fig

def update_plot():
    plot.figure = create_figure()

# --- AG Grid ---
def get_ag_rows():
    # Style rows dynamically based on selection
    return [{
        'index': i,
        'sequence': targets[i],
        'k_clv': values[i],
        # '_style': 'background-color: #ffedcc;' if i in selected_indices else ''
    } for i in range(len(targets))]

def update_grid():
    grid.options['rowData'] = get_ag_rows()
    grid.update()

def toggle_selection(index):
    if index in selected_indices:
        selected_indices.remove(index)
        colors[index] = 'blue'
    else:
        selected_indices.add(index)
        colors[index] = 'orange'
    update_plot()
    update_grid()

# --- Layout ---
with ui.row(align_items='center').classes('w-full gap-0'):

    grid = ui.aggrid({
        'columnDefs': [
            {'headerName': '#', 'field': 'index'},
            {'headerName': 'sequence', 'field': 'sequence'},
            {'headerName': 'k_clv', 'field': 'k_clv'},
        ],
        'rowData': get_ag_rows(),
        'rowSelection': 'multiple',
        # 'getRowId': {'function': 'params.data.index.toString()'},
        # 'suppressRowClickSelection': True,
    }).classes('w-1/2').style("grid-template-columns: auto auto auto")

    grid.options['autoSizeStrategy'] = {'type': 'fitCellContents'}

    plot = ui.plotly(create_figure()).classes('w-1/2')

    update_grid()

# --- Interactions ---
plot.on('plotly_click', lambda e: toggle_selection(e.args['points'][0]['pointIndex']))
grid.on('rowClicked', lambda e: toggle_selection(int(e.args['data']['index'])))

# --- Sum Button ---
with ui.column(align_items='center').classes('w-full gap-0'):
    ui.button('Sum Selected').on('click', lambda: ui.notify(
        f"Selected total: {sum(values[i] for i in selected_indices)}"
    ))

# --- Run ---
ui.run()
