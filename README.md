# ActiveModelH: 2D Pseudo-Spectral Solver

A C++ solver for a two-dimensional **Active Model H**: a conserved scalar order parameter $\psi$ (e.g. a density or composition field) coupled to an incompressible Stokes flow that is driven by an active stress. The equations are integrated with a pseudo-spectral method on a periodic grid, using **Intel MKL** for FFTs and **Armadillo** for array algebra. 

---

<p align="center">
  <img src="plot_steady_states/states_combined.png" width="65%"><br>
  <em> Phase diagram of the active model H: density field (φ) and stream function (ψ) along with velocity field lines (black arrow lines). (Left) Laminar flow state, (Middle) Vortex flow state, and (Right) Turbulent flow state.</em>
</p>

The repository also includes a Python script, `txt_to_matrix_extraction.py`. It collects the snapshot files into matrix format `.mat` file for analysis. 

---

## Contents


<div align="center"> 
  
|1. [Building and running the code](#1-building-and-running-the-code)| 2. [Model](#2-model)| 3. [Input and output files](#3-input-and-output-files)|
|----|----|----|
| 4. [Mian results](#4-main-results)| 5. [Theory](#5-theory)| 6. [Numerical method](#6-numerical-method)|
| 7. [Limitations](#7-limitations)| 8. [References](#8-references)||

</div>

## Table of Contents

1. [Building and running the code](#1-building-and-running-the-code)
- [Model](#model)
- [Numerical Method](#numerical-method)
- [Requirements](#requirements)
- [Building](#1-building-and-running-the-code)
- [Usage](#usage)
- [Input File (`in_data`)](#input-file-in_data)
- [Output Files](#output-files)
- [Post-processing: Export to MATLAB](#post-processing-export-to-matlab)
- [Known Issues and Limitations](#known-issues-and-limitations)
- [License](#license)

---

## 1. Building and running the code

To compile the code, it requires Cmake, Armadillo linear algebra libraries and MKL libraries. After downloading both the libraries, the code can be compiled
as the following.

### Set up the MKL environment (path depends on your install)
The exact MKL link line depends on your compiler, platform and threading choice. Intel's [oneMKL Link Line Advisor](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html) will generate the correct flags for you.


```bash
> source /opt/intel/oneapi/setvars.sh        # or: source /opt/intel/mkl/bin/mklvars.sh intel64

```

### Configure and build

```bash
> mkdir -p build
> cd build

> cmake -S ../ -B .

```

 ### Run

To make the executable file (in the same terminal)

```bash
> make
```

Now the executable file **activeH_heun.exe** will be created. To run the program

```bash
> ./activeH_heun <in/out dir>/ [options]
```

The first argument is the simulation directory. It must contain an `in_data` parameter file, and all output is written there. **Include the trailing slash**, since file names are built by simple string concatenation (`dir + "phi.txt"`).

Running without arguments prints the help message.

### Options

| Flag | Description |
|------|-------------|
| `-C <steps>` | Continue from `phi.txt`, `psi.txt` in the simulation directory. `<steps>` is the number of steps used in the previous average. |
| `-v` | Verbose: write the full state  at every print interval |
| `-n` | Disable noise in the $\phi$ equation |
| `-N` | Disable noise in the $\psi$ equation |

### Examples

```bash
# Fresh run with snapshots
> ./activeH_heun runs/test/ -v

# Continue a run (previous run had 1,000,000 steps)
> ./activeH_heun runs/test/ -C 1000000
```

### Stopping a run

Pressing **Ctrl+C** (SIGINT) stops the simulation gracefully: the current step finishes and the final state is still written to disk, so the run can be resumed with `-c`.

---

## 2. Model

### Order-parameter dynamics

The scalar field φ evolves by advection plus a Cahn–Hilliard-type (Model B) relaxation:

$$
\partial_t \phi + \mathbf{v}\cdot\nabla\phi = M\nabla^2\left(a \phi + b\phi^3 \right) - \kappa \nabla^4 \phi
$$

### Flow field

The velocity is obtained from a stream function ψ, which guarantees incompressibility (∇·**v** = 0):

$$
v_x = \partial_y \psi, \qquad v_y = -\partial_x \psi
$$

The stream function solves a biharmonic (Stokes) equation whose source is the active stress, with strength κ':

$$
\eta\ \nabla^4 \psi = \kappa' \left[\ \partial_x\phi \nabla^2(\partial_y\phi) - \partial_y\phi \nabla^2(\partial_x\phi) \right]
$$


### Parameters

| Symbol | Code variable | `in_data` key | Meaning |
|--------|---------------|---------------|---------|
| λ | `l` | `lambda` | Active coefficient λ (currently only read, not used in the dynamics) |
| κ | `k` | `kappa` | Interfacial stiffness (square-gradient coefficient) |
| κ' | `k1` | `kappa1` | Active stress coefficient |
| ζ | `zeta` | `zeta` | Active coefficient ζ (currently only read, not used in the dynamics) |
| m | `m` | `m` | Mobility (currently multiplies only the `a` term) |
| η | `eta` | `eta` | Viscosity |
| T | `tem` | `tem` | Temperature |
| a, b | `a`, `b` | `a`, `b` | Coefficients of the φ² / φ⁴ bulk free energy |
| D | `Diff` | `D` | Noise strength (read but noise is currently disabled) |

---

## Numerical Method

- **Spatial discretisation:** pseudo-spectral on a periodic `Nx × Ny` grid. Fields are stored in Fourier space as real-to-complex transforms of size `Nx × (Ny/2 + 1)`.
- **Derivatives:** computed exactly in Fourier space by multiplication with *i q*.
- **Nonlinear terms** (φ³, advection **v**·∇φ, the active source term) are evaluated in real space and transformed back.
- **Stream function:** obtained by dividing by |q|⁴ in Fourier space, with a small regulariser (ε = 10⁻⁸) to avoid division by zero at q = 0.
- **Time stepping:** explicit forward Euler with time step `dt`.
- **FFTs:** Intel MKL DFTI, called through external C wrapper functions (see [Building](#building)).

> **Grid spacing.** Wave vectors are built as `2π·n / N`, so the lattice spacing is effectively 1 and the box size is `Nx × Ny`. The `Lx` and `Ly` parameters are read but not currently used.

---

## Requirements

| Dependency | Purpose |
|------------|---------|
| C++11 (or later) compiler | e.g. `g++`, `icpx` |
| [Intel oneMKL](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html) | FFTs (DFTI interface) |
| [Armadillo](https://arma.sourceforge.net/) | Matrix/array operations and file I/O |
| [Boost.Iostreams](https://www.boost.org/doc/libs/release/libs/iostreams/) | "Tee" stream that writes the log to both the terminal and `log.txt` |

### External FFT wrapper

The main file declares, but does not define, the following C functions:

```c
void create_descriptor_handles(int Nx, int Ny);
void fft_forward(double *in, MKL_Complex16 *out);
void fft_backward(MKL_Complex16 *in, double *out);
void fft_padded_forward(double *in, MKL_Complex16 *out);
void fft_padded_backward(MKL_Complex16 *in, double *out);
void free_descriptor_handles();
```

These must be provided by a separate C source file (the MKL FFT wrapper) that is compiled and linked together with the main program.

---



## Input File (`in_data`)

A plain-text file with one `key = value` pair per line. Whitespace is ignored and lines without `=` are skipped (so they can be used as comments). **All 21 keys are required**, and any unrecognised key aborts the program.

```ini
# Solver
steps     = 1000000
pinterval = 10000
Nx        = 128
Ny        = 128
Lx        = 128
Ly        = 128
dt        = 0.001

# Model parameters
D      = 0.0
lambda = 0.0
kappa  = 1.0
kappa1 = 0.5
zeta   = 0.0
eta    = 1.0
tem    = 1.0
m      = 1.0
a      = -0.25
b      = 0.25

# Initial conditions (uniform values)
phi0 = 0.0
psi0 = 0.0
vx0  = 0.0
vy0  = 0.0
```

**Notes**

- Setting `lambda` also sets `kappa` to the same value. Put `kappa` **after** `lambda` so that your chosen value of κ is the one that is used.
- `Nx` and `Ny` should be even; multiples of 4 (ideally powers of 2) are recommended, because the padding buffers are sized `3N/2` and `3N/4 + 1`.
- For phase separation, choose `a < 0` and `b > 0` in the usual φ⁴ convention.

---

## Output Files

All files are written to the simulation directory. Fields are saved as plain-text `Nx × Ny` matrices (Armadillo `raw_ascii` format).

| File | When | Contents |
|------|------|----------|
| `log.txt` | Always | Copy of everything printed to the terminal |
| `phi.txt`, `psi.txt`, `vx.txt`, `vy.txt` | End of run | Final state (also the restart files for `-c` / `-C`) |
| `phi_0.txt`, `psi_0.txt`, `vx_0.txt`, `vy_0.txt` | Start, with `-v` | Initial state |
| `phi_<n>.txt`, `psi_<n>.txt`, `vx_<n>.txt`, `vy_<n>.txt` | Every `pinterval` steps, with `-v` | Snapshots |
| `epr_<n>.txt`, `epr_av_<n>.txt` | Every `pinterval` steps, with `-v -a` | Instantaneous and time-averaged EPR density |
| `epr.txt`, `epr_av.txt` | End of run, with `-a` | Final instantaneous and averaged EPR density |
| `vx_av.txt`, `vy_av.txt`, `v_av.txt`, `phi_av.txt` | Created at start | Reserved for scalar time series (mostly unused at the moment) |

> **Snapshot numbering:** `<n>` is the current step **plus a hard-coded offset of 49,000,000** (`to_string(i+49000000)`). Change or remove this offset in the source to suit your run.

### Quick look with Python

```python
import numpy as np
import matplotlib.pyplot as plt

phi = np.loadtxt("runs/test/phi.txt")
plt.imshow(phi.T, origin="lower", cmap="RdBu_r")
plt.colorbar(label=r"$\phi$")
plt.show()
```

---

## Post-processing: Export to MATLAB

`txt_to_matrix_extraction.py` collects the snapshot files written by a verbose (`-v`) run and saves them, together with the run parameters, in a single MATLAB file called `Data_matrix_Active_H.mat`.

### Requirements

- Python 3
- NumPy
- SciPy (used for `scipy.io.savemat`)

```bash
pip install numpy scipy
```

### What it does

1. **Reads the parameters from `log.txt`.** It parses the `key: value` lines that the solver writes at start-up: `Nx`, `Ny`, `Lx`, `Ly`, `dt`, `lambda`, `kappa`, `kappa1`, `m`, `zeta`, `eta`, `a`, `b` and `tem`.
2. **Loads the snapshots.** It loops over the steps `nstart, nstart + nint, …, nend` and reads `phi_<step>.txt` and `psi_<step>.txt`. It reports any file that is missing.
3. **Flattens each snapshot into one row.** Each `Nx × Ny` field becomes a row of length `Nx*Ny` (row-major / C order). The result is a matrix of shape `(number of snapshots) × (Nx*Ny)`.
4. **Writes `Data_matrix_Active_H.mat`** in the current directory.

### Configuration

Edit the variables at the top of the script before running it:

```python
nstart = 0          # first snapshot step
nint   = 500000     # step between snapshots (should equal pinterval)
nend   = 11000000   # last snapshot step

phi_files = 1       # 1 = extract phi snapshots, 0 = skip
psi_files = 1       # 1 = extract psi snapshots, 0 = skip
epr_files = 1       # currently has no effect (see caveats)
```

The step numbers must match the numbers **in the file names**. Remember that the solver adds an offset of 49,000,000 to each snapshot number. For example, with the offset still in place, a run with `pinterval = 500000` produces `phi_49500000.txt`, `phi_50000000.txt` and so on, so `nstart` and `nend` must be set to those numbers.

### Usage

Run the script **from inside the simulation directory**, because it uses relative paths for `log.txt` and the snapshot files:

```bash
cd runs/test/
python txt_to_matrix_extraction.py
```

### Contents of `Data_matrix_Active_H.mat`

| Variable | Contents |
|----------|----------|
| `data_phi` | φ snapshots, one flattened snapshot per row |
| `data_psi` | ψ snapshots, one flattened snapshot per row |
| `data_epr` | Reserved for EPR snapshots (currently all zeros) |
| `Nx`, `Ny`, `dt`, `lambda`, `kappa`, `kappa1`, `zeta`, `eta`, `m`, `a`, `b`, `temp` | Run parameters taken from `log.txt` |

In MATLAB, you can recover snapshot `n` as a 2D field with:

```matlab
load('Data_matrix_Active_H.mat');
phi_n = reshape(data_phi(n, :), Ny, Nx)';   % transpose undoes MATLAB's column-major reshape
imagesc(phi_n); axis image; colorbar;
```

### Caveats

- **`Ny`, `Lx` and `Ly` are saved incorrectly.** The dictionary passed to `savemat` uses the key `'Ny'` three times (`'Ny':Nx`… `'Ny':Lx, 'Ny':Ly`). As a result, `Lx` and `Ly` are never saved, and `Ny` ends up holding the value of `Ly`. Change the keys to `'Lx':Lx, 'Ly':Ly`.
- **EPR snapshots are not extracted.** The `epr_files` flag exists, but there is no loop that reads `epr_<step>.txt`, so `data_epr` is all zeros.
- **Missing snapshots shift the rows.** When a file is missing, its row is not left empty in place. Later snapshots move up one row and the unused rows remain zero at the end. Check the printed warnings before assuming that row `n` corresponds to step `nstart + n*nint`.
- **Most parameters are stored as strings.** Only `Nx` and `Ny` are converted to integers. The other parameters are saved as text and need converting (e.g. with `str2double` in MATLAB).
- **Memory use.** All snapshots are held in memory at once, as three arrays of `(number of snapshots) × Nx × Ny` doubles. This can be large for fine grids or many snapshots.

---

## Known Issues and Limitations

The code is under active development. Before relying on results, be aware of the following:

1. **Noise is disabled.** Thermal noise terms are commented out, and `sample_noise()` is defined but never called. A perfectly uniform initial state (`phi0` constant) therefore **never evolves**. Start from a perturbed configuration by placing a `phi.txt` in the directory and running with `-c`.
2. **Forward Euler, not Heun.** Despite the executable name, the predictor–corrector (Heun) step is commented out, so the scheme is first order. Keep `dt` small, since the κ∇⁴ term is stiff.
3. **`-x` (AXY preset) bug.** `l = 3/16` and `k = 5/16` use integer division and evaluate to **0**. Use `3.0/16` and `5.0/16`.
4. **`-C` does not truly resume the EPR average.** `epr_av.txt` is loaded into a real-space matrix, but the running average is accumulated in `epr_av_ft`, which starts from zero.
5. **Unused parameters.** `lambda`, `zeta`, `Lx`, `Ly` and `D` are read and logged but do not enter the dynamics. The mobility `m` multiplies only the `a` term.
6. **2π approximation.** Wave vectors use `6.28` instead of 2π; replace it with `2*M_PI` for accuracy.
7. **No dealiasing during the run.** The Nyquist row and column are zeroed only at initialisation. The `Convolution` class (3/2-rule padded convolutions) exists but is currently used only to create and free the FFT descriptors.
8. **Scalar time series are not written.** `phi_av` is never updated, and `epr_av_file` is never opened, so the corresponding outputs are empty or constant.

---

## License

*Add your license here (e.g. MIT, GPL-3.0).*

## Citation

*If this code is associated with a publication, add the reference here.*

## Citation

*If this code is associated with a publication, add the reference here.*
