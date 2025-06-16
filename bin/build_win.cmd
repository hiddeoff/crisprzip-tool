set "VENV_PATH=%USERPROFILE%\Anaconda3\envs\cziptool_venv"

pyinstaller crisprzip_gui.py ^
    --name CRISPRzip ^
    --onefile ^
    --windowed ^
    --add-data "%VENV_PATH%\Lib\site-packages\nicegui\static;nicegui/static" ^
    --add-data "%VENV_PATH%\Lib\site-packages\latex2mathml\;latex2mathml" ^
    --collect-all nicegui ^
    --collect-all crisprzip ^
    --collect-all matplotlib ^
    --collect-all numpy ^
    --collect-all pandas ^
    --hidden-import uvicorn.logging