from nicegui import ui
import content.vitro_cleavage
import content.vitro_binding


@ui.page('/')
def index():

    # HEADER
    with ui.header(elevated=True).style('background-color: #D9F0D2'):
        with ui.row().classes():
            # Image/icon to the left
            ui.image('./img/CRISPRzip_logo_v0.jpg').props('width=45px height=45px')

            # Text next to the image
            ui.label('CRISPRzip GUI').style('color: black; font-size: 32px; font-family: Courier New; font-weight: bold;')

    # FOOTER
    with ui.footer(elevated=True).classes('py-1 h-10 bg-[#D9F0D2] flex items-center'):
        ui.markdown('CRISPRzip GUI is created with [niceGUI](https://nicegui.io). Licensed under MIT.').style('color: gray; font-size: 10px;')

    # TABS
    with ui.tabs().classes('w-full bg-gradient-to-b from-white via-white to-[#6ECFF6] rounded-xl') as tabs:
        one = ui.tab('cleavage', icon='content_cut')
        two = ui.tab('binding', icon='link')

    with ui.tab_panels(tabs, value=one).classes('w-full'):

        # TAB 1 - IN VITRO
        with ui.tab_panel(one):
            content.vitro_cleavage.show_contents()  # content/vitro_cleavage.py

        # TAB 2 - IN VIVO
        with ui.tab_panel(two):
            content.vitro_binding.show_contents()  # content/vitro_binding.py

ui.run()
