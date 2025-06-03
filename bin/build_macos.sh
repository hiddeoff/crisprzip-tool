lib_path="/c/Users/hiddeofferhaus/Projects/crisprzip-tool/cziptool_venv/Lib"

pyinstaller crisprzip_gui.py \
  --name CRISPRzip \
  --windowed \
  --onedir \
  --add-data "{$lib_path}/site-packages/nicegui:nicegui/static" \
  --add-data "{$lib_path}/site-packages/latex2mathml:latex2mathml" \
  --collect-all nicegui \
  --collect-all crisprzip \
  --collect-all matplotlib \
  --collect-all numpy \
  --collect-all pandas \
  --hidden-import uvicorn.logging
