import os
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

__version__ = "0.6.3"

# The public interface: one function per command (scarlet/cli.py builds its options from their signatures), the
# model loader and the tutorial datasets' downloader (scarlet.datasets).
_lazy_imports = {
    "learn": ("scarlet.api", "learn"),
    "assign": ("scarlet.api", "assign"),
    "transfer": ("scarlet.api", "transfer"),
    "report": ("scarlet.api", "report"),
    "load_model": ("scarlet.model", "load_model"),
    "Model": ("scarlet.model", "Model"),
    "datasets": ("scarlet.datasets", None),
}
__all__ = list(_lazy_imports)

def __getattr__(name):
    if name in _lazy_imports:
        import importlib
        module_name, attr_name = _lazy_imports[name]
        module = importlib.import_module(module_name)
        val = module if attr_name is None else getattr(module, attr_name)
        globals()[name] = val
        return val
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

def __dir__():
    return sorted(list(globals().keys()) + list(_lazy_imports.keys()))
