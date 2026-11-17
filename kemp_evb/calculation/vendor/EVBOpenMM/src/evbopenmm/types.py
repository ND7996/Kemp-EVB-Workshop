from __future__ import annotations

import os
from typing import Tuple

PathLike = os.PathLike[str] | str
AtomReference = Tuple[str, str, int, str]
CanonicalAtomKey = Tuple[AtomReference, ...]
CanonicalPairKey = Tuple[AtomReference, AtomReference]

__all__ = [
    "PathLike",
    "AtomReference",
    "CanonicalAtomKey",
    "CanonicalPairKey",
]
