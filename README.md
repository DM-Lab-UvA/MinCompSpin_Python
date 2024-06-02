# Requirements
**C++:**  The code uses the C++11 version of C++.

**Python:** The code works with Python 3 and requires `pybind11` to be installed. See [here](https://pybind11.readthedocs.io/en/stable/basics.html) for the requirements for pybind11 depending on your environment. 

# Install MinCompSpin package

```console
  $ git clone git@github.com:DM-Lab-UvA/MinCompSpin_Python.git
  $ cd MinCompSpin_Python
  $ make py_module
```

# Usage

The package can be built as shown in [Install](#install).
You can then import it into python using:
```python
  import MinCompSpin
  MinCompSpin.main("input/path.bin", n)
```

For quick examples on how to use the package check:
 - the Jupyter notebook: `Test_Binding.ipynb`
 - the python script: `test.py`
