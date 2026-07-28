@echo off
REM =============================================================================
REM  setup_oracle_env.bat — Windows / Miniconda
REM
REM  Creates the "oracle" Conda environment for CellOracle KO simulation
REM  (04_oracle_ko_simulation.ipynb). This is the consolidated, SELF-CONTAINED
REM  setup: the environment.yml spec and the (previously split) "PART 2" batch
REM  logic are now folded into this single file, mirroring the layout of
REM  setup_pyscenic_env.bat and setup_scvi_env.bat.
REM
REM  USAGE
REM    Open Anaconda Prompt as your NORMAL user (NOT as Administrator).
REM    cd /d to the folder containing this .bat.
REM    Run:  setup_oracle_env.bat
REM    Then: conda activate oracle -> jupyter lab -> Kernel "Python (oracle)"
REM
REM  HARD CONSTRAINTS (celloracle 0.20.0):
REM    anndata    == 0.10.8   -- celloracle requires <=0.10.8
REM    pandas     == 1.5.3    -- celloracle requires <=1.5.3
REM    matplotlib == 3.6.3    -- celloracle requires <3.7
REM    numpy      == 1.23.5   -- celloracle uses deprecated np.bool/int aliases
REM
REM  The conda scientific stack is created from an environment.yml that this
REM  script WRITES ITSELF (base-conda python), so no external .yml is needed.
REM  celloracle + the Linux-only-package stubs are then layered on with pip,
REM  exactly as before. If `conda env create` reports a solver conflict on a
REM  specific package, relax THAT package to a range and note it in the README.
REM =============================================================================

setlocal enableextensions enabledelayedexpansion
set SCRIPT_DIR=%~dp0
set PIP_OPTS=--no-cache-dir --disable-pip-version-check

echo.
echo =====================================================================
echo  oracle environment (CellOracle KO simulation)
echo =====================================================================

call conda deactivate

REM ---- O0. Clear pip cache ----
echo.
echo [O0] Clearing pip cache...
if exist "%LOCALAPPDATA%\pip" rmdir /s /q "%LOCALAPPDATA%\pip" 2>nul

REM ---- O1. Remove old oracle env ----
echo.
echo [O1] Removing old oracle env (if exists)...
call conda env remove -n oracle -y 2>nul

REM ---- O2. Write environment.yml (self-contained, pure batch -- no Python
REM      needed yet, since no env exists at this point in the script) ----
echo.
echo [O2] Writing environment.yml spec...
(
echo name: oracle
echo channels:
echo   - conda-forge
echo   - defaults
echo dependencies:
echo   - python=3.10.13
echo   - numpy=1.23.5
echo   - setuptools=75.6.0
echo   - wheel=0.45.1
echo   - cython=3.0.11
echo   - pandas=1.5.3
echo   - matplotlib=3.6.3
echo   - anndata=0.10.8
echo   - scipy=1.10.1
echo   - scikit-learn=1.3.2
echo   - scikit-misc
echo   - joblib=1.4.2
echo   - numba=0.57.1
echo   - statsmodels=0.14.0
echo   - scanpy=1.9.3
echo   - h5py=3.10.0
echo   - umap-learn=0.5.5
echo   - pynndescent=0.5.11
echo   - natsort=8.4.0
echo   - tqdm=4.66.5
echo   - python-igraph=0.11.6
echo   - leidenalg=0.10.2
echo   - networkx=3.2.1
echo   - seaborn=0.12.2
echo   - packaging=24.1
echo   - python-dateutil=2.9.0
echo   - pytz=2024.1
echo   - jupyterlab=4.0.6
echo   - ipykernel=6.25.2
echo   - ipywidgets=8.1.1
) > "%SCRIPT_DIR%environment.yml"

if not exist "%SCRIPT_DIR%environment.yml" (
    echo FAILED: environment.yml was not created in "%SCRIPT_DIR%" -- check that this folder is writable
    exit /b 1
)

REM ---- O3. Create env from environment.yml ----
echo.
echo [O3] Creating oracle env from environment.yml...
call conda env create -f "%SCRIPT_DIR%environment.yml"
if errorlevel 1 ( echo FAILED: conda env create & exit /b 1 )
call conda activate oracle
if errorlevel 1 ( echo FAILED: conda activate oracle & exit /b 1 )

REM ---- O4. Stubs for Linux-only packages ----
echo.
echo [O4] Creating stubs for Linux-only packages...
echo     (velocyto, pysam, gimmemotifs, pybedtools, genomepy, louvain)
python -c "import base64; exec(base64.b64decode('aW1wb3J0IG9zCgpTVFVCUyA9IHsKICAgICJ2ZWxvY3l0byI6ICAgIHsidmVyc2lvbiI6IjAuMTcuMTciLCJzdWJzIjpbImFuYWx5c2lzIiwiY291bnRlciIsImRpZmZ1c2lvbiIsImVzdGltYXRpb24iLCJleHBybWF0X2RlbHRhbWF0cml4X2xvYWRlciIsImV4dHJhY3RfcmVhZHMiLCJsb2dpYyIsImxvb20iLCJuZWlnaGJvcnMiLCJzZXJpYWxpemF0aW9uIiwidHJhbnNjcmlwdF9tb2RlbCIsImNvbnN0YW50cyIsImNsdXN0ZXJfc3RhdGlzdGljcyJdfSwKICAgICJweXNhbSI6ICAgICAgIHsidmVyc2lvbiI6IjAuMjQuMCIsICJzdWJzIjpbImxpYmNodHNsaWIiLCJjYWxpZ25tZW50ZmlsZSIsImN0YWJpeCIsImN2Y2YiLCJjZmFpZHgiLCJjYWxpZ25lZHNlZ21lbnQiLCJjc2FtZmlsZSIsImNiY2YiXX0sCiAgICAiZ2ltbWVtb3RpZnMiOiB7InZlcnNpb24iOiIwLjE3LjIiLCAic3VicyI6WyJtb3RpZiIsImZhc3RhIiwic2Nhbm5lciIsImNvbmZpZyIsInBsb3QiLCJ1dGlscyIsIm1vYXAiLCJjb21wYXJpc29uIiwicm9jbWV0cmljcyIsImJhY2tncm91bmQiLCJnZW5vbWVfaW5kZXgiLCJwcmVwcm9jZXNzaW5nIiwic3RhdHMiXX0sCiAgICAicHliZWR0b29scyI6ICB7InZlcnNpb24iOiIwLjEyLjAiLCAic3VicyI6WyJiZWR0b29sIiwiY2JlZHRvb2xzIiwiaGVscGVycyIsImZlYXR1cmVmdW5jcyIsImNvbnRyaWIiLCJzY3JpcHRzIiwiZ2Vub21lX3JlZ2lzdHJ5Il19LAogICAgImdlbm9tZXB5IjogICAgeyJ2ZXJzaW9uIjoiMC4xNi40IiwgInN1YnMiOlsicHJvdmlkZXIiLCJnZW5vbWUiLCJhbm5vdGF0aW9uIiwidXRpbHMiLCJzZXFfc3RhdHMiLCJwbHVnaW4iLCJpbnN0YWxsIiwibWFuYWdlIl19LAogICAgImxvdXZhaW4iOiAgICAgeyJ2ZXJzaW9uIjoiMC44LjIiLCAgInN1YnMiOltdfSwKfQoKSU5JVCA9ICgKICAgICJpbXBvcnQgc3lzIGFzIF9zLCB0eXBlcyBhcyBfdFxuIgogICAgIl9fdmVyc2lvbl9fID0gXCJ7dmVyfVwiXG4iCiAgICAiX19hbGxfXyA9IFtdXG4iCiAgICAiX05PTkUgPSBmcm96ZW5zZXQoKGNocig5NSkqMit4K2Nocig5NSkqMiBmb3IgeCBpbiAoY2hyKDEwMikrY2hyKDEwNSkrY2hyKDEwOCkrY2hyKDEwMSksY2hyKDExNSkrY2hyKDExMikrY2hyKDEwMSkrY2hyKDk5KSxjaHIoMTA4KStjaHIoMTExKStjaHIoOTcpK2NocigxMDApK2NocigxMDEpK2NocigxMTQpLGNocigxMTIpK2Nocig5NykrY2hyKDExNikrY2hyKDEwNCksY2hyKDk5KStjaHIoOTcpK2Nocig5OSkrY2hyKDEwNCkrY2hyKDEwMSkrY2hyKDEwMCksY2hyKDk4KStjaHIoMTE3KStjaHIoMTA1KStjaHIoMTA4KStjaHIoMTE2KStjaHIoMTA1KStjaHIoMTEwKStjaHIoMTE1KSkpKVxuIgogICAgImNsYXNzIF9TTSh0eXBlKTpcbiIKICAgICIgICAgZGVmIF9fZ2V0YXR0cl9fKGNscyxuKTogcmV0dXJuIE5vbmUgaWYgbiBpbiBfTk9ORSBlbHNlIGNsc1xuIgogICAgImNsYXNzIF9TdHViKG1ldGFjbGFzcz1fU00pOlxuIgogICAgIiAgICBkZWYgX19pbml0X18oc2VsZiwqYSwqKmspOiBwYXNzXG4iCiAgICAiICAgIGRlZiBfX2NhbGxfXyhzZWxmLCphLCoqayk6IHJldHVybiBzZWxmXG4iCiAgICAiICAgIGRlZiBfX2dldGF0dHJfXyhzZWxmLG4pOiByZXR1cm4gX1N0dWIoKVxuIgogICAgIiAgICBkZWYgX19pdGVyX18oc2VsZik6IHJldHVybiBpdGVyKFtdKVxuIgogICAgIiAgICBkZWYgX19sZW5fXyhzZWxmKTogcmV0dXJuIDBcbiIKICAgICIgICAgZGVmIF9fYm9vbF9fKHNlbGYpOiByZXR1cm4gRmFsc2VcbiIKICAgICIgICAgZGVmIF9fZW50ZXJfXyhzZWxmKTogcmV0dXJuIHNlbGZcbiIKICAgICIgICAgZGVmIF9fZXhpdF9fKHNlbGYsKmEpOiBwYXNzXG4iCiAgICAiZGVmIF9fZ2V0YXR0cl9fKG4pOiByZXR1cm4gX1N0dWJcbiIKICAgICJjbGFzcyBfTShfdC5Nb2R1bGVUeXBlKTpcbiIKICAgICIgICAgZGVmIF9fZ2V0YXR0cl9fKHNlbGYsbik6IHJldHVybiBOb25lIGlmIG4gaW4gX05PTkUgZWxzZSBfU3R1YlxuIgogICAgImZvciBfc3ViIGluIHtzdWJzfTpcbiIKICAgICIgICAgaWYgXCJ7cGtnfS5cIisgX3N1YiBub3QgaW4gX3MubW9kdWxlczpcbiIKICAgICIgICAgICAgIF9tPV9NKFwie3BrZ30uXCIrIF9zdWIpOyBfbS5fX3BhY2thZ2VfXz1cIntwa2d9XCI7IF9tLl9fc3BlY19fPU5vbmU7IF9zLm1vZHVsZXNbXCJ7cGtnfS5cIisgX3N1Yl09X21cbiIKKQoKVE9NTCA9ICgKICAgICJbYnVpbGQtc3lzdGVtXVxuIgogICAgInJlcXVpcmVzID0gW1wic2V0dXB0b29scz49NjFcIl1cbiIKICAgICJidWlsZC1iYWNrZW5kID0gXCJzZXR1cHRvb2xzLmJ1aWxkX21ldGFcIlxuIgogICAgIltwcm9qZWN0XVxuIgogICAgIm5hbWUgPSBcIntwa2d9XCJcbiIKICAgICJ2ZXJzaW9uID0gXCJ7dmVyfVwiXG4iCiAgICAicmVxdWlyZXMtcHl0aG9uID0gXCI+PTMuOFwiXG4iCiAgICAiZGVwZW5kZW5jaWVzID0gW11cbiIKICAgICJbdG9vbC5zZXR1cHRvb2xzLnBhY2thZ2VzLmZpbmRdXG4iCiAgICAid2hlcmUgPSBbXCIuXCJdXG4iCiAgICAiaW5jbHVkZSA9IFtcIntwa2d9KlwiXVxuIgopCgpmb3IgcGtnLCBpbmZvIGluIFNUVUJTLml0ZW1zKCk6CiAgICB2ZXIgPSBpbmZvWyJ2ZXJzaW9uIl0KICAgIHN1YnMgPSByZXByKHR1cGxlKGluZm9bInN1YnMiXSkpCiAgICBkID0gb3MucGF0aC5qb2luKCJfc3R1YnMiLCBwa2csIHBrZykKICAgIG9zLm1ha2VkaXJzKGQsIGV4aXN0X29rPVRydWUpCiAgICBvcGVuKG9zLnBhdGguam9pbigiX3N0dWJzIixwa2csInB5cHJvamVjdC50b21sIiksInciKS53cml0ZShUT01MLmZvcm1hdChwa2c9cGtnLHZlcj12ZXIpKQogICAgb3Blbihvcy5wYXRoLmpvaW4oZCwiX19pbml0X18ucHkiKSwidyIpLndyaXRlKElOSVQuZm9ybWF0KHBrZz1wa2csdmVyPXZlcixzdWJzPXN1YnMpKQogICAgcHJpbnQoIiAgc3R1YjoiLCBwa2cp').decode())"
if errorlevel 1 ( echo FAILED: stub creation & exit /b 1 )
for %%P in (velocyto pysam gimmemotifs pybedtools genomepy louvain) do (
    pip install %PIP_OPTS% .\_stubs\%%P
    if errorlevel 1 ( echo FAILED: stub %%P & exit /b 1 )
)

REM ---- O5. celloracle (--no-deps to protect stubs) ----
echo.
echo [O5] Installing celloracle==0.20.0 (--no-deps)...
pip install %PIP_OPTS% "celloracle==0.20.0" --no-deps
if errorlevel 1 ( echo FAILED: celloracle & exit /b 1 )

REM ---- O6. Remaining deps ----
echo.
echo [O6] Remaining deps: pyarrow, goatools, jupyter...
pip install %PIP_OPTS% "pyarrow==11.0.0" "goatools==1.4.12" "jupyter==1.0.0"
if errorlevel 1 ( echo FAILED: remaining deps & exit /b 1 )

REM ---- O7. Register kernel ----
echo.
echo [O7] Registering Jupyter kernel (oracle)...
python -m ipykernel install --user --name oracle --display-name "Python (oracle)"

REM ---- O8. Cleanup ----
echo.
echo [O8] Cleaning up temp stub directories...
rmdir /s /q _stubs 2>nul

REM ---- O9. Sanity check ----
echo.
echo [O9] Sanity check oracle env...
python -c "import celloracle as co, anndata as ad, scanpy as sc, numpy as np, igraph, leidenalg; print('celloracle:', co.__version__); print('anndata:   ', ad.__version__); print('scanpy:    ', sc.__version__); print('numpy:     ', np.__version__); print('ALL OK')"
if errorlevel 1 ( echo SANITY CHECK FAILED & exit /b 1 )

REM ---- O10. SSL patch (Windows Python 3.10 malformed cert bug) ----
echo.
echo [O10] Patching Windows SSL cert store (Python 3.10 safety fix)...
python -c "import base64,sys,os; sc=os.path.join(sys.prefix,'Lib','site-packages','sitecustomize.py'); open(sc,'a').write(base64.b64decode('CiMgU1NMIFdpbmRvd3MgY2VydC1zdG9yZSBmaXggKHNldHVwX29yYWNsZV9lbnYuYmF0KQojIFB5dGhvbiAzLjEwIG9uIFdpbmRvd3MgY3Jhc2hlcyB3aGVuIGEgbWFsZm9ybWVkIGNlcnRpZmljYXRlIGlzIGluIHRoZSBzdG9yZS4KIyBUaGlzIHdyYXBzIF9sb2FkX3dpbmRvd3Nfc3RvcmVfY2VydHMgdG8gc2lsZW50bHkgc2tpcCBiYWQgY2VydHMuCmltcG9ydCBzc2wgYXMgX3NzbAppZiBoYXNhdHRyKF9zc2wuU1NMQ29udGV4dCwgJ19sb2FkX3dpbmRvd3Nfc3RvcmVfY2VydHMnKToKICAgIF9vcmlnX2x3c2MgPSBfc3NsLlNTTENvbnRleHQuX2xvYWRfd2luZG93c19zdG9yZV9jZXJ0cwogICAgZGVmIF9zYWZlX2x3c2Moc2VsZiwgc3RvcmVuYW1lLCBwdXJwb3NlKToKICAgICAgICB0cnk6CiAgICAgICAgICAgIF9vcmlnX2x3c2Moc2VsZiwgc3RvcmVuYW1lLCBwdXJwb3NlKQogICAgICAgIGV4Y2VwdCBfc3NsLlNTTEVycm9yOgogICAgICAgICAgICBwYXNzCiAgICBfc3NsLlNTTENvbnRleHQuX2xvYWRfd2luZG93c19zdG9yZV9jZXJ0cyA9IF9zYWZlX2x3c2MK').decode()); print('SSL patch written.')"

echo.
echo =============================================================
echo  Setup complete.
echo  oracle kernel : Python (oracle)
echo  To run: conda activate oracle -^> jupyter lab -^> 04_oracle_ko_simulation.ipynb
echo =============================================================
endlocal
