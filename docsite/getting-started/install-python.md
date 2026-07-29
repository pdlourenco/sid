# Install — Python

Requires **Python 3.10+**, **NumPy**, and **SciPy**. Matplotlib is optional
(needed only for plotting functions and notebook examples).

## From the repository

```bash
git clone https://github.com/pdlourenco/sid.git
cd sid
pip install -e ./python
```

With plotting support (recommended — needed for the example notebooks):

```bash
pip install -e "./python[plot]"
```

For development (tests + linting + notebook validation):

```bash
pip install -e "./python[dev]"
```

## Quick smoke test

```python
import sid
print(sid.__version__)
```

If that prints a version string, you're set. Head to
[Quick start](quick-start.md) for a runnable example, or jump straight to
the [Python notebooks](../examples/python/index.md).

## Compatibility

| Dependency | Minimum version |
|---|---|
| Python | 3.10 |
| NumPy | 1.22 |
| SciPy | 1.8 |
| Matplotlib | 3.5 (optional) |
