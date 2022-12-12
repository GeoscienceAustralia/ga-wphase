from functools import lru_cache, wraps
from os import getpid
from typing import Callable, Optional, TypeVar
from typing_extensions import ParamSpec
from uuid import uuid4

class ForkSharedMemoryException(Exception):
    pass

class ForkSharedMemory:
    storage = {}
    key: str
    content: object

    def __init__(self, content):
        self.key = str(uuid4())
        self.content = content
        ForkSharedMemory.storage[self.key] = self

    def __getstate__(self):
        return self.key

    def __setstate__(self, state):
        try:
            stored = ForkSharedMemory.storage[state]
        except KeyError:
            raise ForkSharedMemoryException(
                f"Could not find shared memory key {state}!"
                " Note that this module is only compatible with forking worker pools."
            )
        self.key = state
        self.content = stored.content

T = TypeVar("T")
P = ParamSpec("P")
def per_proc_cache(maxsize: Optional[int] = None):
    """lru_cache variant that clears the cache upon forking."""
    def decorator(getter: Callable[P, T]) -> Callable[P, T]:
        pid = getpid()
        lru = lru_cache(maxsize=maxsize)(getter)

        @wraps(getter)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
            nonlocal pid
            newpid = getpid()
            if newpid != pid:
                # We've forked since the last call, so purge the cache!
                lru.cache_clear()
                pid = newpid
            return lru(*args, **kwargs)

        return wrapper
    return decorator
