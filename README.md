# CRISPRzip tool
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Welcome to the codebase of the **CRISPRzip tool** from the [Depken Lab](https://depkenlab.tudelft.nl/) at TU
Delft.

![A screenshot of the CRISPRzip tool in action](img/Screenshot_border.png)

## About the project
### The CRISPRzip model
CRISPRzip is a physics-based model to study the target 
recognition dynamics of CRISPR-associated nucleases like Cas9
([Eslami-Mossalam, 2022](#references)). Their interactions with target DNA is represented 
as an energy landscape, with which you can simulate binding and cleavage
kinetics. The parameters have been obtained by machine learning on 
high-throughput data. CRISPRzip makes quantitative predictions of on-target 
efficiency and off-target risks of different guide RNAs. With CRISPRzip, we hope 
to contribute to assessing
the risks that come with particular choices in CRISPR application, and as such
contribute to the development of safe gene editing technology.

### The tool
With the CRISPRzip tool/GUI, you can apply CRISRPzip to the sequences of your 
interest without having to do any programming. It's a matter of installing the
right executable for your computer (Windows/Mac/Linux) and launching it.

The CRISPRzip tool is powered by [NiceGUI](https://nicegui.io/) and 
[PyInstaller](https://pyinstaller.org/en/stable/).

### The Python package
If you are familiar with Python and would like to work with the source code:
good news! It's available as a public [GitHub project](https://github.com/hiddeoff/crisprzip)
and as Python package on [PyPi](https://pypi.org/project/crisprzip/). 

### References
Eslami-Mossallam B et al. (2022) *A kinetic model predicts SpCas9 activity,
improves off-target classification, and reveals the physical basis of
targeting fidelity.* Nature Communications.
[10.1038/s41467-022-28994-2](https://doi.org/10.1038/s41467-022-28994-2)

## Usage
You can download the latest version of the CRISPRzip tool from the 
[Releases page](https://github.com/hiddeoff/crisprizp-tool/releases).
Download and launch!

### Available platforms
| Platform      | File                                                                                                        |
|---------------|-------------------------------------------------------------------------------------------------------------|
| Windows       | [crisprzip.exe](https://github.com/hiddeoff/crisprizp-tool/releases/download/latest/crisprzip-tool-win.exe) |
| macOS         | [crisrpzip](https://github.com/hiddeoff/crisprizp-tool/releases/download/latest/crisprzip-tool-macos)       |
| Linux         | [crisrpzip](https://github.com/hiddeoff/crisprizp-tool/releases/download/latest/crisprzip-tool-unix)        |

### Browser NiceGUI application
As an alternative to the downloadable applications, you could clone this
repository and launch CRISPRzip-tool from the terminal. It will open in your
browser. This way to launch the tool is sometimes (a lot) faster than the 
executables. Follow the developer instructions below.

## Developers
To use and develop the GUI, run the following commands:
1. Clone the repository.
```bash
git clone https://github.com/hiddeoff/crisprzip-tool.git
```
2. Navigate to the cloned repository.
3. Create a virtual environment (venv or conda), activate it, and install the required dependencies.
```bash
conda create -n crisprzip_gui python=3.12
conda activate crisprzip_gui
pip install -r requirements.txt
```
4.  Run the GUI. It should launch in your browser.
```bash
python crisprzip_gui.py
```

## Building the executable
If you want to build the executable for the CRISPRzip tool, you will need to have [PyInstaller](https://pyinstaller.org/en/stable/) installed in addition to the requirements in `requirements.txt`.

You can install it with pip in your virtual environment:
```bash
pip install pyinstaller
```
Then, you can build the executable by running the following command from the root of the project directory, depending on your platform:

<details>
  <summary>Windows</summary>

  - **Instructions for Windows**

</details>

<details>
  <summary>macOS</summary>

```bash
pyinstaller crisprzip_gui.py \
  --name CRISPRzip \
  --windowed \
  --onedir \
  --add-data "/path/to/your/venv/lib/python3.12/site-packages/nicegui:nicegui/static" \
  --add-data "/path/to/your/venv/lib/python3.12/site-packages/latex2mathml:latex2mathml" \
  --collect-all nicegui \
  --collect-all crisprzip \
  --collect-all matplotlib \
  --collect-all numpy \
  --collect-all pandas \
  --hidden-import uvicorn.logging
```
- Important: Replace the paths with the correct ones for your system!

</details>

<details>
  <summary>Linux</summary>

  - **Instructions for Linux**

</details>

These commands will create a `dist` folder containing the executable and a `build` folder with the build files. Alternatively, you can run these commands from the `bin` folder, which contains scripts to build the executable for your platform. Remember to update the paths as needed.

## Contributing
If you would like to contribute to this project: that's great! Have a look at our 
[Contributing guidelines](./CONTRIBUTING.md) and our [Code of Conduct](./CODE_OF_CONDUCT.md).

## Acknowledgements
Many thanks to [Elviss Dvinskis](https://github.com/edvinskis) and 
[Raúl Ortiz](https://github.com/rortizmerino) from the [DCC team at TU Delft](https://www.tudelft.nl/en/library/support/library-for-researchers/setting-up-research/dcc)
for their support to get this GUI released!

## Waiver
Technische Universiteit Delft hereby disclaims all copyright interest in the
program “CRISPRzip-tool” written by the Author(s).
Paulien Herder, Dean of Applied Sciences

(c) 2025, Hidde Offerhaus, Delft, The Netherlands.
