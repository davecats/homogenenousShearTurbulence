"""Export runtime statistics and mean estimates to a self-describing file.

Usage::

    python3 -m post.python export_statistics /path/to/data -o statistics.nc

The input directory contains ``hst.in`` and normally a ``statistics``
subdirectory.  The two runtime files are also accepted directly beside
``hst.in`` for compatibility with older runs. NetCDF4 is written through
xarray's NetCDF4 backend; install ``xarray`` and its ``netCDF4`` backend in
the active virtual environment. JSON remains available with ``--format json``.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np

from ..config import load_config
from ..io import find_run_file, read_runtimedata, read_text_repaired, read_variances


RUNTIME_COLUMNS = (
    ("meanflowx", "mean streamwise flow integral", "solver units"),
    ("meanflowy", "mean spanwise flow integral", "solver units"),
    ("S", "mean shear", "solver units"),
    ("S2", "second shear parameter", "solver units"),
    ("gamma_x", "streamwise shift", "solver units"),
    ("gamma_y", "spanwise shift", "solver units"),
    ("deltat", "time step", "time units"),
    ("cfl_dt", "CFL-limited time step", "time units"),
    ("energy", "turbulent kinetic energy, <u_i u_i>/2", "velocity squared"),
    ("diss", "dissipation", "solver units"),
    ("uv", "Reynolds stress <u v>", "velocity squared"),
    ("vw", "Reynolds stress <v w>", "velocity squared"),
)
VARIANCE_COLUMNS = (
    ("uu", "velocity second moment <u u>", "velocity squared"),
    ("vv", "velocity second moment <v v>", "velocity squared"),
    ("ww", "velocity second moment <w w>", "velocity squared"),
    ("uw", "Reynolds stress <u w>", "velocity squared"),
)


_find_file = find_run_file


def _meanandvar(
    executable: Path, data_file: Path, number_columns: int, skip_rows: int = 0
) -> list[dict[str, float]]:
    """Run meanandvar over all columns and parse ``mean, uncertainty`` output."""
    cpl = shutil.which("cpl")
    if cpl is None:
        raise RuntimeError("cpl command was not found in PATH; load the CPL environment first")
    try:
        completed = subprocess.run(
            [cpl, "run", str(executable), f"1-{number_columns}", str(skip_rows)],
            input=read_text_repaired(data_file),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    except OSError as exc:
        raise RuntimeError(f"could not run cpl/meanandvar: {exc}") from exc
    if completed.returncode:
        raise RuntimeError(
            f"meanandvar failed with status {completed.returncode}: {completed.stderr.strip()}"
        )

    # CPL's WRITE statement separates fields with tabs but does not add a
    # newline between selected columns.  Therefore the result is commonly one
    # line containing ``mean uncertainty`` repeated for every column (with
    # ``./.`` for unavailable values).  Tokenise the complete stdout rather
    # than assuming one result per line.
    values = []
    for token in completed.stdout.split():
        try:
            values.append(float(token))
        except ValueError:
            pass
    rows = [
        {"mean": values[index], "uncertainty": values[index + 1]}
        for index in range(0, len(values) - 1, 2)
    ]
    if len(rows) != number_columns:
        raise RuntimeError(
            f"meanandvar returned {len(rows)} results, expected {number_columns}; "
            f"stdout was: {completed.stdout!r}"
        )
    return rows


# Quantities that vanish identically in a relaminarized (laminar, steady) flow.
TURBULENT_QUANTITIES = ("energy", "diss", "uv", "vw", "uu", "vv", "ww", "uw")


def _zero_turbulent_estimates(columns, estimates):
    """Zero mean and uncertainty of the turbulent quantities (estimates[0] is time)."""
    out = list(estimates)
    for index, (name, _, _) in enumerate(columns, start=1):
        if name in TURBULENT_QUANTITIES:
            out[index] = {"mean": 0.0, "uncertainty": 0.0}
    return out


def _series(data, columns, estimates, time):
    result = {}
    for (name, description, units), estimate in zip(columns, estimates):
        result[name] = {
            "description": description,
            "units": units,
            "time": time.tolist(),
            "values": np.asarray(data[name], dtype=float).tolist(),
            **estimate,
        }
    return result


def _configuration(cfg):
    values = asdict(cfg)
    values["run_dir"] = str(cfg.run_dir.resolve())
    values["deck_path"] = str(cfg.deck_path.resolve())
    values.update(lx=cfg.lx, lz=cfg.lz, nu=cfg.nu)
    return values


def _document(cfg, runtime_file, variance_file, executable, runtime, variances,
              runtime_estimates, variance_estimates, t_init_stat=None, skip_rows=0):
    """Build the JSON-compatible logical representation used by both writers."""
    return {
        "schema": "hst-statistics",
        "schema_version": 1,
        "description": "Self-describing HST time-series statistics and batch-mean estimates.",
        "source": {
            "run_directory": str(cfg.run_dir.resolve()),
            "deck": str(cfg.deck_path.resolve()),
            "runtime_file": str(runtime_file.resolve()),
            "variance_file": str(variance_file.resolve()),
            "meanandvar": str(executable.resolve()),
        },
        "configuration": _configuration(cfg),
        "statistics": {
            "dt_stat": cfg.dt_stat,
            "t_init_stat": t_init_stat,
            "skipped_rows": skip_rows,
        },
        "datasets": {
            "Runtimedata": _series(
                runtime, RUNTIME_COLUMNS, runtime_estimates[1:], runtime["time"]
            ),
            "variances_runtime.dat": _series(
                variances, VARIANCE_COLUMNS, variance_estimates[1:], variances["time"]
            ),
        },
    }


def _write_netcdf(document, output):
    """Write ``document`` as a compressed NetCDF4 file through xarray.

    The two source files normally share exactly one time axis.  A separate
    dimension is used only if a malformed/partial run gives them different
    lengths or timestamps; this preserves all data instead of silently
    truncating one source.
    """
    try:
        import xarray as xr
    except ImportError as exc:
        raise RuntimeError(
            "NetCDF4 output requires xarray; install xarray and its netCDF4 "
            "backend in the active virtual environment, or use --format json"
        ) from exc

    datasets = document["datasets"]
    rt = datasets["Runtimedata"]
    var = datasets["variances_runtime.dat"]
    rt_time = np.asarray(next(iter(rt.values()))["time"], dtype=np.float64)
    var_time = np.asarray(next(iter(var.values()))["time"], dtype=np.float64)
    shared = rt_time.shape == var_time.shape and np.array_equal(rt_time, var_time)

    data_vars = {}
    encoding = {}
    if shared:
        time_name = {"Runtimedata": "time", "variances_runtime.dat": "time"}
        coords = {
            "time": ("time", rt_time, {
                "units": "simulation time",
                "long_name": "shared sample time",
                "source": "Runtimedata and variances_runtime.dat",
            })
        }
    else:
        time_name = {"Runtimedata": "time_runtime", "variances_runtime.dat": "time_variance"}
        coords = {
            "time_runtime": ("time_runtime", rt_time, {
                "units": "simulation time",
                "long_name": "Runtimedata sample time",
            }),
            "time_variance": ("time_variance", var_time, {
                "units": "simulation time",
                "long_name": "variances_runtime.dat sample time",
            }),
        }

    for dataset_name, series in datasets.items():
        coordinate = time_name[dataset_name]
        for name, info in series.items():
            data_vars[name] = (coordinate, np.asarray(info["values"], dtype=np.float64), {
                "units": info["units"],
                "long_name": info["description"],
                "source_dataset": dataset_name,
            })
            data_vars[f"{name}_mean"] = ((), np.float64(info["mean"]), {
                "units": info["units"],
                "long_name": info["description"],
                "source_dataset": dataset_name,
                "statistic": "mean",
            })
            data_vars[f"{name}_uncertainty"] = ((), np.float64(info["uncertainty"]), {
                "units": info["units"],
                "long_name": info["description"],
                "source_dataset": dataset_name,
                "statistic": "standard deviation of the mean",
            })
            encoding[name] = {"zlib": True, "complevel": 4, "dtype": "float64"}

    attrs = {
        "history": "Created by post.python export_statistics",
        "schema": document["schema"],
        "schema_version": document["schema_version"],
        "description": document["description"],
        "source": json.dumps(document["source"], sort_keys=True),
        "configuration": json.dumps(document["configuration"], sort_keys=True),
        "statistics": json.dumps(document["statistics"], sort_keys=True),
        "dataset_names": "Runtimedata, variances_runtime.dat",
    }
    try:
        dataset = xr.Dataset(data_vars=data_vars, coords=coords, attrs=attrs)
        dataset.to_netcdf(
            str(output),
            mode="w",
            format="NETCDF4",
            engine="netcdf4",
            encoding=encoding,
        )
    except (ImportError, ModuleNotFoundError, ValueError) as exc:
        raise RuntimeError(
            "xarray's netCDF4 backend is unavailable; install xarray and "
            "the netCDF4 package in the active virtual environment"
        ) from exc


def export(run_dir: Path, output: Path, executable: Path, output_format="netcdf",
           t_init_stat=None, relaminarized=False) -> None:
    cfg = load_config(run_dir)
    runtime_file = _find_file(cfg.run_dir, "Runtimedata")
    variance_file = _find_file(cfg.run_dir, "variances_runtime.dat")
    runtime = read_runtimedata(runtime_file)
    variances = read_variances(variance_file)

    if t_init_stat is None:
        skip_rows = 0
    else:
        if t_init_stat < 0:
            raise ValueError("--t-init-stat must be non-negative")
        if cfg.dt_stat <= 0:
            raise ValueError("dt_stat must be positive when --t-init-stat is used")
        skip_rows = round(t_init_stat / cfg.dt_stat) + 1

    runtime_estimates = _meanandvar(
        executable, runtime_file, len(RUNTIME_COLUMNS) + 1, skip_rows
    )
    variance_estimates = _meanandvar(
        executable, variance_file, len(VARIANCE_COLUMNS) + 1, skip_rows
    )

    if relaminarized:
        # the flow is laminar: all turbulent statistics are exactly zero (the decay
        # transient and the roundoff-level tail must not enter the averages)
        runtime_estimates = _zero_turbulent_estimates(RUNTIME_COLUMNS, runtime_estimates)
        variance_estimates = _zero_turbulent_estimates(VARIANCE_COLUMNS, variance_estimates)

    # The first column is time. It is retained as the coordinate, not as a
    # physical statistic, hence the [1:] passed to _series.
    document = _document(cfg, runtime_file, variance_file, executable, runtime,
                         variances, runtime_estimates, variance_estimates,
                         t_init_stat, skip_rows)
    if relaminarized:
        document["statistics"]["relaminarized"] = True
    output.parent.mkdir(parents=True, exist_ok=True)
    if output_format == "json":
        output.write_text(json.dumps(document, indent=2) + "\n")
    else:
        _write_netcdf(document, output)


def run(cfg, argv) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output", default=None, help="output path (default: statistics.nc)")
    parser.add_argument("--format", choices=("netcdf4", "json"), default="netcdf4",
                        help="output format (default: netcdf4; JSON is a fallback)")
    parser.add_argument(
        "--meanandvar", type=Path, default=None,
        help="path to the meanandvar executable (default: post/cpl/meanandvar_utils/meanandvar)",
    )
    parser.add_argument(
        "--t-init-stat", type=float, default=None,
        help=("initial transient duration to exclude from statistics; rows skipped "
              "are round(t_init_stat / dt_stat) + 1 (default: no rows skipped)"),
    )
    parser.add_argument(
        "--relaminarized", action="store_true",
        help=("the flow relaminarized: set mean and uncertainty of all turbulent quantities "
              "(energy, diss, uv, vw, uu, vv, ww, uw) to zero and flag the file"),
    )
    args = parser.parse_args(argv)
    executable = args.meanandvar or (
        Path(__file__).resolve().parents[2] / "cpl" / "meanandvar_utils" / "meanandvar"
    )
    output = Path(args.output) if args.output else cfg.run_dir / (
        "statistics.json" if args.format == "json" else "statistics.nc"
    )
    try:
        export(cfg.run_dir, output, executable, args.format, args.t_init_stat, args.relaminarized)
    except (FileNotFoundError, RuntimeError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"wrote {output}")
    return 0