from nicegui import native,ui
import content.vitro_cleavage
import content.vitro_binding

# macOS packaging support
#from multiprocessing import freeze_support
#freeze_support()


@ui.page('/')
def index():

    # HEADER
    with ui.header(elevated=True).style('background-color: #F5F9FF').classes('p-1'):
        ui.element().classes('w-[25px]')
        with ui.row().classes('items-center'):
            with ui.tabs().style('color: gray') as tabs:
                one = ui.tab('cleavage', icon='content_cut')
                two = ui.tab('binding', icon='link')

        ui.space()
        ui.space()
        ui.space()

        with ui.row(align_items='center').classes('h-[70px] items-center gap-2'):
            ui.image('https://raw.githubusercontent.com/hiddeoff/crisprzip-tool/refs/heads/main/img/CRISPRzip_logo_v0_gradient_nobg.svg').props('width=60px height=60px')
            (ui.label('CRISPRzip tool')
            .style('color: gray; font-size: 40px; font-weight: 100;'
                    'font-family: Helvetica Neue, Roboto, Inter, sans-serif;'))

        ui.space()
        ui.space()
        ui.space()
        ui.space()

        with ui.link(target="https://github.com/hiddeoff/crisprzip-model", new_tab=True).style('textDecoration: none'):
            with ui.column(align_items='center').classes('h-full opacity-60 gap-0 p-0'):
                (ui.image('https://img.icons8.com/glyph-neue/64/github.png')
                 .props('width=40px height=40px'))
                ui.html("CRISPRzip").style('color: black').classes('leading-[1.0]')
                ui.html("on GitHub").style('color: black').classes('leading-[1.0]')
        ui.element().classes('w-[25px]')

    # FOOTER
    with ui.footer(elevated=True).classes('py-1 h-6 bg-[#F5FAF4] flex items-center'):
        (ui.markdown('CRISPRzip tool is created with [NiceGUI](https://nicegui.io). Licensed under MIT.')
         .style('color: gray; font-size: 10px;')
         .classes('leading-[0.0]'))

    with ui.tab_panels(tabs, value=one).classes('w-full'):

        # TAB 1 - CLEAVAGE
        with ui.tab_panel(one):
            content.vitro_cleavage.show_contents()  # content/vitro_cleavage.py

        # TAB 2 - BINDING
        with ui.tab_panel(two):
            content.vitro_binding.show_contents()  # content/vitro_binding.py

ui.run(
    # Uncomment the next two lines if you want to build an executable, or to run in a contained window
    native=True,
    reload=False,
    title='CRISPRzip tool',
    favicon="img/CRISPRzip_logo_v0_gradient_nobg.svg"
)
