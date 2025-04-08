from nicegui import ui
import content.vitro_cleavage

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
        one = ui.tab('in vitro', icon='science')
        two = ui.tab('in vivo', icon='pest_control_rodent')

    with ui.tab_panels(tabs, value=one).classes('w-full'):

        # TAB 1 - IN VITRO
        with ui.tab_panel(one):

            with ui.tabs().classes('justify-start') as in_vitro_tabs:
                vitro_cleavage = ui.tab('cleavage (Cas9)', icon='content_cut').props('no-caps')
                vitro_binding = ui.tab('binding (dCas9)', icon='link').props('no-caps')

            with ui.tab_panels(in_vitro_tabs, value=vitro_cleavage).classes('w-full'):
                with ui.tab_panel(vitro_cleavage).classes('w-full h-full no-wrap'):
                    content.vitro_cleavage.show()  # content/vitro_cleavage.py
                with ui.tab_panel(vitro_binding):
                    ui.label('IN VITRO - BINDING (dCas9) GOES HERE')

        # TAB 2 - IN VIVO
        with ui.tab_panel(two):

            with ui.tabs().classes('justify-start') as in_vivo_tabs:
                vivo_cleavage = ui.tab('cleavage (Cas9)', icon='content_cut').props('no-caps')
                vivo_binding = ui.tab('binding (dCas9)', icon='link').props('no-caps')

            with ui.tab_panels(in_vivo_tabs, value=vivo_cleavage).classes('w-full'):
                with ui.tab_panel(vivo_cleavage):
                    ui.label('IN VIVO - CLEAVAGE (Cas9) PLACEHOLDER')
                with ui.tab_panel(vivo_binding):
                    ui.label('IN VIVO - BINDING (dCas9) PLACEHOLDER')


ui.run()
