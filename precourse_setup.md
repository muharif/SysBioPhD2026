# Pre-course setup

Please complete the setup before the course starts. 

1. [**Posit Cloud — recommended**](#1-posit-cloud)
2. Choose one of these, depending on your operating system. You only need to set up **one** of these.
- [**macOS/Linux with Terminal**](#2-macos-or-linux)
- [**Windows with WSL**](#3-windows-with-wsl)
3. Complete the [pre-course survey](https://forms.gle/bN4FHu9Djp81MRr47)

NOTE: download the [env.yml](https://github.com/muharif/SysBioPhD2026/blob/main/env.yml) file here

## 1. Posit Cloud

Register for the free account at [Posit Cloud](https://posit.cloud/) using your registered email for this course. This will be used as an option to set up the programming environment. The link to join the course workspace will be provided later this week.

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
