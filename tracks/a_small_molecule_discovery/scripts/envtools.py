"""Locate executables that live in sibling conda/mamba environments.

No dependency beyond the standard library, so it is safe to import from any of the
pipeline environments (pdac-dd, pdac-md, boltz2, ...). Resolution order:

1. an explicit absolute path or a name found on PATH,
2. <env>/bin/<exe> under every conda envs directory found on this machine.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

_ENV_ROOTS = ("miniforge3", "mambaforge", "miniconda3", "anaconda3", "micromamba", ".conda")


def envs_dirs() -> list[Path]:
    dirs: list[Path] = []
    for item in os.environ.get("CONDA_ENVS_PATH", "").split(os.pathsep):
        if item:
            dirs.append(Path(item).expanduser())
    prefix = os.environ.get("CONDA_PREFIX") or os.environ.get("MAMBA_ROOT_PREFIX")
    if prefix:
        p = Path(prefix)
        dirs.extend([p / "envs", p.parent])
    home = Path.home()
    for name in _ENV_ROOTS:
        dirs.append(home / name / "envs")
    dirs.append(home / "micromamba" / "envs")
    seen, out = set(), []
    for d in dirs:
        if d not in seen and d.is_dir():
            seen.add(d)
            out.append(d)
    return out


def env_tool(env: str | list[str] | tuple[str, ...], exe: str) -> str:
    """Return the path of `exe` inside conda env(s) `env`, else `exe` itself."""
    direct = os.path.expanduser(exe)
    if os.path.isabs(direct) and os.path.isfile(direct):
        return direct
    names = [env] if isinstance(env, str) else list(env)
    for name in names:
        for base in envs_dirs():
            cand = base / name / "bin" / exe
            if cand.is_file():
                return str(cand)
    found = shutil.which(exe)
    return found or exe


def env_prefix(env: str | list[str] | tuple[str, ...]) -> str | None:
    names = [env] if isinstance(env, str) else list(env)
    for name in names:
        for base in envs_dirs():
            if (base / name / "bin").is_dir():
                return str(base / name)
    return None
