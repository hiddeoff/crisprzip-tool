# Contributing

Thanks for considering to contribute to CRISPRzip-tool! We welcome contributions
of all kinds—whether it's fixing a bug, suggesting an enhancement, or
improving the documentation. Below are some guidelines to help you get started.

## How to Contribute

### 1. Reporting Issues
If you encounter a bug or have a feature request, please:
- Check the issue tracker to see if it has already been reported.
- If not, open a new issue and include:
   - A clear description of the problem or feature.
   - Steps to reproduce the issue, if applicable.
   - Relevant details like version numbers or environment.

### 2. Making Code Changes
#### 2.1 Setup

1. Fork the repository and clone your fork locally:
```shell
git clone https://github.com/username/crisprzip-tool.git
```
2. Navigate to the cloned repository.
3. Create a venv or conda env, activate it, and install requirements.
```bash
conda create -n crisprzip_gui python=3.12
conda activate crisprzip_gui
pip install -r requirements.txt
```
4. Create a new branch for your changes:
```shell
git checkout -b feature/your-feature-name
```
5. Run the GUI. It should launch in your browser.
```bash
python crisprzip_gui.py
```

#### 2.2 Making code changes
- Follow PEP 8 for Python code style.
- Write clear and concise commit messages.
- Add or update tests for any new features or bug fixes.

#### 2.3 Submit a draft pull request
Push your changes to your fork:
```shell
git push origin feature/your-feature-name
```
Open a pull request to the main branch on the original repository.
Ensure your pull request includes:
- A detailed description of the changes.
- Reference to related issues (e.g., Closes #123).

Your pull request will be reviewed by a maintainer.
Please be responsive to feedback and make necessary updates.
Once approved, your contribution will be merged.

## Code of Conduct
Please adhere to our [Code of Conduct](./CODE_OF_CONDUCT.md) in all your 
interactions within the project.

## Questions?
Always feel free to reach out by opening an issue or contacting
[hsofferhaus@gmail.com](mailto:hsofferhaus@gmail.com). We appreciate your contributions! 🎉
