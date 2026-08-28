# Pre-course setup

Please complete the setup before the course starts.

You can use one of the following environments:

1. [**Posit Cloud — recommended**](#1-posit-cloud--recommended)
2. [**macOS/Linux with Terminal**](#2-macos-or-linux)
3. [**Windows with WSL**](#3-windows-with-wsl)

You only need to set up **one** of these.

## 1. Posit Cloud — recommended

The easiest option is to use [Posit Cloud](https://posit.cloud/).

1. Create an account or sign in.
2. Open the course workspace/project when provided.
3. Open Jupyter and check that you can create and run a notebook.

For example:

```python
1 + 1
```

You should get:

```text
2
```

No local installation is required.

---

## 2. macOS or Linux

If you prefer to work locally, use your normal Terminal.

Check that basic Unix commands work:

```bash
pwd
ls
```

### Install Miniforge

Install [Miniforge](https://github.com/conda-forge/miniforge) and follow the instructions for your operating system.

After installation, open a new terminal and check:

```bash
conda --version
```

### Install the course environment

Download the course `env.yml` file.

From the directory containing the file, run:

```bash
conda env create -f env.yml
```

Activate the environment:

```bash
conda activate SysBioPhD2026
```

Check that Python and PLINK2 are available:

```bash
python --version
plink2 --version
```

Start Jupyter:

```bash
jupyter lab
```

---

## 3. Windows with WSL

Windows users who want to work locally should use **WSL 2**.

Open **PowerShell as Administrator** and run:

```powershell
wsl --install
```

Restart your computer if prompted.

Then open **Ubuntu** from the Start menu and create your Linux username and password.

Check that basic Unix commands work:

```bash
pwd
ls
```

### Install Miniforge inside WSL

Install [Miniforge](https://github.com/conda-forge/miniforge) **inside WSL**, following the Linux installation instructions.

After installation, open a new WSL terminal and check:

```bash
conda --version
```

### Install the course environment

Download the course `env.yml` file.

From the directory containing it, run:

```bash
conda env create -f env.yml
```

Activate the environment:

```bash
conda activate SysBioPhD2026
```

Check that Python and PLINK2 are available:

```bash
python --version
plink2 --version
```

Start Jupyter:

```bash
jupyter lab
```

JupyterLab should open in your Windows web browser.

---

## Before the course

Make sure your chosen setup works:

- **Posit Cloud:** you can open and run a Jupyter notebook.
- **macOS/Linux:** `conda activate SysBioPhD2026`, `python --version`, `plink2 --version`, and `jupyter lab` work.
- **Windows/WSL:** the same commands work inside Ubuntu/WSL.

You only need to create the Conda environment once. For later sessions, activate it with:

```bash
conda activate SysBioPhD2026
```

If these steps work, you are ready for the course.