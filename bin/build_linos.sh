lib_path="$HOME/Documents/CRISPRzip/crisprzip-tool/crisprzip_venv/lib/python3.12/"

pyinstaller crisprzip_gui.py \
  --name CRISPRzip \
  --windowed \
  --onefile \
  --add-data "$lib_path/site-packages/nicegui:nicegui/static" \
  --add-data "$lib_path/site-packages/latex2mathml:latex2mathml" \
  --collect-all nicegui \
  --collect-all crisprzip \
  --collect-all matplotlib \
  --collect-all numpy \
  --collect-all pandas \
  --collect-all qtpy \
  --hidden-import uvicorn.logging \
  --hidden-import PySide6.QtWebEngineWidgets \
  --exclude-module gi --exclude-module PyGObject --exclude-module gtk
