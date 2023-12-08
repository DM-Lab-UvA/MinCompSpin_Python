# Requirements
**C++:**  The code uses the C++11 version of C++.

**Python:** requires `pybind11` to be installed.

# Install

```console
  $ git clone git@github.com:DM-Lab-UvA/MinCompSpin_Python.git
  $ cd MinCompSpin_Python
  $ make
  $ python test.py
```

# Usage

The package can be built as shown in [Install](#install).
You can then import it into python like so:
```python
  import MinCompSpin
  MinCompSpin.main("input/path.bin", n)
```
