@echo off
REM =============================================================================
REM  setup_scvi_env.bat — Windows / Miniconda
REM
REM  Creates the scvi Conda environment used by Script 01
REM  (01_gating_clustering.R, via reticulate::use_condaenv("scvi")) for scVI
REM  batch correction. No Jupyter kernel is registered -- unlike pyscenic and
REM  oracle, this environment is called directly from R via reticulate, never
REM  opened as a standalone notebook.
REM
REM  USAGE
REM    Open Anaconda Prompt as your NORMAL user (NOT as Administrator).
REM    Run:  setup_scvi_env.bat
REM    Then in R: reticulate::use_condaenv("scvi", required = TRUE)
REM
REM  CRITICAL VERSION PINS:
REM    numpy==2.2.6, torch==2.12.1, jax==0.6.2,
REM    scvi-tools==1.3.3, Python 3.10.20
REM
REM  CPU vs GPU NOTE:
REM    torch and jaxlib below carry no "+cuXXX" version suffix and were
REM    captured from a "pypi_0"-tagged (default PyPI) build in the source
REM    environment -- both signs of a CPU-only install, not CUDA-enabled.
REM    This script replicates that exactly. If GPU acceleration is actually
REM    needed, torch and jax/jaxlib require different install commands
REM    (specific index URLs) not reflected here.
REM =============================================================================

setlocal enableextensions enabledelayedexpansion
set PIP_OPTS=--no-cache-dir --disable-pip-version-check

REM ---- E1. Create clean scvi env ----
echo.
echo [E1] Removing old scvi env + creating fresh one...
call conda deactivate
call conda env remove -n scvi -y 2>nul
call conda create -n scvi python=3.10.20 -y || ( echo FAILED: scvi env create & exit /b 1 )
call conda activate scvi || ( echo FAILED: activate scvi & exit /b 1 )

REM ---- S1. pip / wheel / packaging ----
echo.
echo [S1] pip / wheel / packaging...
pip install %PIP_OPTS% --no-deps "packaging==26.0" || exit /b 1

REM ---- S2. Scientific base ----
echo.
echo [S2] Scientific base...
pip install %PIP_OPTS% --no-deps "numpy==2.2.6" "scipy==1.15.3" "pandas==2.3.3" "scikit-learn==1.7.2" "numba==0.65.1" "llvmlite==0.47.0" "threadpoolctl==3.6.0" "joblib==1.5.3" "six==1.17.0" "python-dateutil==2.9.0.post0" "pytz==2026.2" "tzdata==2026.2" || exit /b 1

REM ---- S3. Async/networking (dask deps) ----
echo.
echo [S3] Async/networking (dask deps)...
pip install %PIP_OPTS% --no-deps "aiohttp==3.14.1" "aiohappyeyeballs==2.6.2" "aiosignal==1.4.0" "async-timeout==5.0.1" "frozenlist==1.8.0" "multidict==6.7.1" "propcache==0.5.2" "yarl==1.24.2" "attrs==26.1.0" "idna==3.18" || exit /b 1

REM ---- S4. Dask ----
echo.
echo [S4] Dask...
pip install %PIP_OPTS% --no-deps "dask==2026.6.0" "cloudpickle==3.1.2" "fsspec==2026.6.0" "locket==1.0.0" "partd==1.4.2" "toolz==1.1.0" "msgpack==1.2.1" "pyyaml==6.0.3" "click==8.4.1" "psutil==7.2.2" || exit /b 1

REM ---- S5. Single-cell data structures ----
echo.
echo [S5] Single-cell data structures...
pip install %PIP_OPTS% --no-deps "h5py==3.16.0" "anndata==0.11.4" "mudata==0.3.8" "natsort==8.4.0" "networkx==3.4.2" "sparse==0.17.0" "xarray==2025.6.1" "docrep==0.3.2" "scverse-misc==0.0.3" "session-info2==0.4.1" "array-api-compat==1.15.0" "humanize==4.15.0" || exit /b 1

REM ---- S6. JAX ecosystem (CPU build -- no +cuda suffix in source env) ----
echo.
echo [S6] JAX ecosystem (CPU build -- no +cuda suffix in source env)...
pip install %PIP_OPTS% --no-deps "jax==0.6.2" "jaxlib==0.6.2" "ml-dtypes==0.5.4" "opt-einsum==3.4.0" "etils==1.13.0" "ml-collections==1.1.0" "orbax-checkpoint==0.11.39" || exit /b 1

REM ---- S7. Flax / Optax / NumPyro / Pyro (scVI's probabilistic backend) ----
echo.
echo [S7] Flax / Optax / NumPyro / Pyro (scVI's probabilistic backend)...
pip install %PIP_OPTS% --no-deps "flax==0.10.7" "optax==0.2.8" "numpyro==0.19.0" "pyro-api==0.1.2" "pyro-ppl==1.9.1" "multipledispatch==1.0.0" || exit /b 1

REM ---- S8. PyTorch + Lightning (CPU build -- no +cuXXX suffix in source env) ----
echo.
echo [S8] PyTorch + Lightning (CPU build -- no +cuXXX suffix in source env)...
pip install %PIP_OPTS% --no-deps "torch==2.12.1" "torchmetrics==1.9.0" "lightning==2.6.5" "pytorch-lightning==2.6.5" "lightning-utilities==0.15.3" "pillow==12.2.0" || exit /b 1

REM ---- S9. TensorBoard / logging ----
echo.
echo [S9] TensorBoard / logging...
pip install %PIP_OPTS% --no-deps "tensorboard==2.20.0" "tensorboard-data-server==0.7.2" "tensorstore==0.1.78" "treescope==0.1.10" "grpcio==1.81.1" "protobuf==7.35.1" "absl-py==2.4.0" "markdown==3.10.2" "werkzeug==3.1.8" || exit /b 1

REM ---- S10. scvi-tools ----
echo.
echo [S10] scvi-tools...
pip install %PIP_OPTS% --no-deps "scvi-tools==1.3.3" || exit /b 1

REM ---- S11. Misc utilities ----
echo.
echo [S11] Misc utilities...
pip install %PIP_OPTS% --no-deps "rich==15.0.0" "pygments==2.20.0" "markdown-it-py==4.2.0" "mdurl==0.1.2" "colorama==0.4.6" "sympy==1.14.0" "mpmath==1.3.0" "filelock==3.29.4" "jinja2==3.1.6" "markupsafe==3.0.3" "typing-extensions==4.15.0" "exceptiongroup==1.3.1" "importlib-metadata==9.0.0" "importlib-resources==7.1.0" "zipp==4.1.0" "simplejson==4.1.1" "aiofiles==25.1.0" "nest-asyncio==1.6.0" "tqdm==4.68.3" || exit /b 1

echo.
echo =============================================================
echo  Setup complete.
echo  scvi env : scvi  (no Jupyter kernel -- used via reticulate from R)
echo  To run: in R, reticulate::use_condaenv("scvi", required = TRUE)
echo =============================================================
endlocal
