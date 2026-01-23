import functools
import inspect
import warnings
from typing import Mapping, Callable, Optional, Type

def deprecate_kwargs(
    renames: Mapping[str, str],
    *,
    since: Optional[str] = None,
    remove_in: Optional[str] = None,
    warn_cls: Type[Warning] = FutureWarning,
) -> Callable:
    """
    Decorator to deprecate keyword argument names.

    Parameters
    ----------
    renames : {old_name: new_name}
        Dictionary mapping old argument names to new argument names.
    since : str, optional
        Specify the version since when the deprecation started.
    remove_in : str, optional
        Specify the version when the deprecated argument will be removed.
    warn_cls : Warning subclass
        The warning class to use (default: FutureWarning).
    """
    def decorator(func: Callable) -> Callable:
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            # If deprecated args are used, issue warnings and map to new names
            for old, new in renames.items():
                if old in kwargs:
                    if new in kwargs:
                        raise TypeError(
                            f"{func.__name__}() received both '{old}' (deprecated) "
                            f"and '{new}'. Use only '{new}'."
                        )
                    msg = f"'{old}' is deprecated; use '{new}' instead"
                    if since:
                        msg += f" (since {since})"
                    if remove_in:
                        msg += f"; support will be removed in {remove_in}"
                    warnings.warn(msg, warn_cls, stacklevel=2)
                    kwargs[new] = kwargs.pop(old)

            # Bind arguments to check validity
            sig.bind_partial(*args, **kwargs)
            return func(*args, **kwargs)

        # Preserve the original function signature
        wrapper.__signature__ = sig
        return wrapper
    return decorator

# Alias for easier usage
def rename_kwargs(**old_to_new: str) -> Callable:
    return deprecate_kwargs(old_to_new)
