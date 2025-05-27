from nicegui import ui
import content.vitro_cleavage
import content.vitro_binding


@ui.page('/')
def index():

    # HEADER
    with ui.header(elevated=True).style('background-color: #F5F9FF'):
        with ui.tabs().style('color: gray') as tabs:
            one = ui.tab('cleavage', icon='content_cut')
            two = ui.tab('binding', icon='link')
        ui.space()
        ui.image('./img/CRISPRzip_logo_no_background.png').props('width=70px height=70px')
        ui.label('CRISPRzip tool').style('color: gray; font-size: 40px; font-family: Helvetica Neue; font-weight: 300;')
        ui.space()
        ui.space()

    # FOOTER
    with ui.footer(elevated=True).classes('py-1 h-10 bg-[#F5FAF4] flex items-center'):
        ui.markdown('CRISPRzip tool is created with [NiceGUI](https://nicegui.io). Licensed under MIT.').style('color: gray; font-size: 10px;')

    with ui.tab_panels(tabs, value=one).classes('w-full'):

        # TAB 1 - CLEAVAGE
        with ui.tab_panel(one):
            content.vitro_cleavage.show_contents()  # content/vitro_cleavage.py

        # TAB 2 - BINDING
        with ui.tab_panel(two):
            content.vitro_binding.show_contents()  # content/vitro_binding.py

ui.run()
