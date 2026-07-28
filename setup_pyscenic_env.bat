@echo off
REM =============================================================================
REM  setup_pyscenic_env.bat — Windows / Miniconda
REM
REM  Creates the pyscenic Conda environment for pyscenic/GRNBoost2 GRN inference
REM  (03_pyscenic.ipynb).
REM
REM  USAGE
REM    Open Anaconda Prompt as your NORMAL user (NOT as Administrator).
REM    cd /d to the folder containing this .bat.
REM    Run:  setup_pyscenic_env.bat
REM    Then: conda activate pyscenic -> jupyter lab -> Kernel "Python (pyscenic)"
REM
REM  CRITICAL VERSION PINS:
REM    numpy==1.23.5, pyarrow==11.0.0 (v10 .feather compat), pandas==1.5.3,
REM    Python 3.10
REM =============================================================================

setlocal enableextensions enabledelayedexpansion
set PIP_OPTS=--no-cache-dir --disable-pip-version-check

REM ---- P1. Create clean pyscenic env ----
echo.
echo [P1] Removing old pyscenic env + creating fresh one...
call conda deactivate
call conda env remove -n pyscenic -y 2>nul
call conda create -n pyscenic python=3.10.13 -y || ( echo FAILED: pyscenic env create & exit /b 1 )
call conda activate pyscenic || ( echo FAILED: activate pyscenic & exit /b 1 )

REM ---- P2. Base: numpy + setuptools ----
echo.
echo [P2] Base: numpy + setuptools...
pip install %PIP_OPTS% "numpy==1.23.5" "setuptools<81" "wheel" || exit /b 1

REM ---- P3. Scientific stack ----
echo.
echo [P3] Scientific stack...
pip install %PIP_OPTS% --no-deps "pandas==1.5.3" "scipy==1.10.1" "numexpr==2.8.4" "bottleneck==1.3.7" "python-dateutil>=2.8" "pytz>=2020" "six>=1.5" || exit /b 1

REM ---- P4. Dask ----
echo.
echo [P4] Dask...
pip install %PIP_OPTS% --no-deps "dask==2023.5.0" "distributed==2023.5.0" "cloudpickle>=1.5" "fsspec>=2021.09" "partd>=1.2" "toolz>=0.10" "tornado>=6.0" "msgpack>=1.0" "zict>=2.0" "tblib>=1.6" "pyyaml>=5.3" "click>=7.0" "psutil>=5.0" "locket>=1.0" "urllib3" "importlib-metadata>=4.13.0" "zipp>=3.0" "sortedcontainers>=2.0.5" || exit /b 1

REM ---- P5. aiohttp ----
echo.
echo [P5] aiohttp...
pip install %PIP_OPTS% --no-deps "aiohttp" "aiohappyeyeballs>=2.5.0" "aiosignal>=1.4.0" "async-timeout>=4.0" "frozenlist>=1.1.1" "multidict>=4.5,<7.0" "propcache>=0.2.0" "yarl>=1.17.0,<2.0" "attrs>=22.2.0" "idna>=2.0" || exit /b 1

REM ---- P6. Numba ----
echo.
echo [P6] Numba...
pip install %PIP_OPTS% --no-deps "numba==0.57.1" "llvmlite==0.40.1" || exit /b 1

REM ---- P7. scikit-learn ----
echo.
echo [P7] scikit-learn...
pip install %PIP_OPTS% --no-deps "scikit-learn==1.3.2" "joblib>=1.0" "threadpoolctl>=2.0" || exit /b 1

REM ---- P8. AnnData + Scanpy ----
REM h5py tightened to ==3.10.0 (was >=3.0; number pulled from fix_env.py's
REM own pin, not changed -- fix_env.py itself is now folded into this script).
echo.
echo [P8] AnnData + Scanpy...
pip install %PIP_OPTS% --no-deps "anndata==0.9.2" "scanpy==1.9.3" "loompy==3.0.7" "h5py==3.10.0" "natsort>=7.0" "packaging>=20" "pynndescent>=0.5" "tqdm>=4.0" "umap-learn==0.5.5" "stdlib_list" "session_info" "patsy" "click-plugins" "networkx" "numpy_groupies>=0.9" || exit /b 1

REM ---- P9. pySCENIC (pyarrow==11.0.0 is CRITICAL) ----
echo.
echo [P9] pySCENIC (pyarrow==11.0.0 is CRITICAL for v10 .feather databases)...
pip install %PIP_OPTS% --no-deps "pyscenic==0.12.1" "ctxcore==0.2.0" "arboreto==0.1.6" "frozendict<2.4" "cytoolz" "multiprocessing-on-dill" "dill" "interlap" "boltons" "requests" "pyarrow==11.0.0" "pandas-flavor" "statsmodels==0.14.0" || exit /b 1

REM ---- P10. Plotting ----
echo.
echo [P10] Matplotlib + seaborn...
pip install %PIP_OPTS% --no-deps "matplotlib==3.7.3" "seaborn==0.12.2" "pillow" "kiwisolver" "cycler" "fonttools" "pyparsing" "contourpy" || exit /b 1

REM ---- P11. Jupyter ----
echo.
echo [P11] Jupyter...
pip install %PIP_OPTS% --no-deps "jupyter==1.0.0" "ipykernel==6.25.2" "ipython==8.14.0" "jupyter-client==8.3.1" "jupyter-core==5.3.1" "jupyterlab==4.0.6" "notebook==7.0.4" "nbformat==5.9.2" "nbconvert==7.8.0" "nbclient==0.8.0" "ipywidgets==8.1.1" "jupyter-server==2.7.3" "jupyterlab-server==2.25.0" "jupyter-lsp==2.2.0" "jupyter-events==0.7.0" "jupyter-server-terminals==0.4.4" "notebook-shim==0.2.3" "widgetsnbextension==4.0.9" "jupyterlab-widgets==3.0.9" "jupyterlab-pygments==0.2.2" "jupyter-console==6.6.3" || exit /b 1
pip install %PIP_OPTS% --no-deps "traitlets==5.10.0" "platformdirs>=2.5" "comm>=0.1.1" "debugpy>=1.6.5" "matplotlib-inline>=0.1" "nest-asyncio>=1.4" "pyzmq>=25" "prompt-toolkit==3.0.39" "pygments>=2.4" "decorator>=4.3" "jedi>=0.16" "parso>=0.8" "wcwidth" "executing" "asttokens" "pure-eval" "stack-data" "colorama" "exceptiongroup" "typing_extensions>=4.6" "tomli" "jinja2" "markupsafe" "anyio" "sniffio" "httpcore" "h11" "httpx" "async-lru" "babel" "json5" "jsonschema" "jsonschema-specifications" "referencing" "rpds-py" "fastjsonschema" "argon2-cffi" "argon2-cffi-bindings" "cffi" "pycparser" "send2trash" "terminado" "websocket-client" "prometheus-client" "pywinpty" "python-json-logger" "rfc3339-validator" "rfc3986-validator" "overrides" "defusedxml" "beautifulsoup4" "soupsieve" "bleach" "tinycss2" "webencodings" "mistune" "pandocfilters" "charset-normalizer" "certifi" "pickleshare" "backcall" || exit /b 1

REM ---- P12. Register kernel ----
echo.
echo [P12] Registering Jupyter kernel (pyscenic)...
python -m ipykernel install --user --name pyscenic --display-name "Python (pyscenic)"

echo.
echo =============================================================
echo  Setup complete.
echo  pyscenic kernel : Python (pyscenic)
echo  To run: conda activate pyscenic -^> jupyter lab -^> scenic_nk_analysis.ipynb
echo =============================================================
endlocal
