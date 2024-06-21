from setuptools import setup, Extension

extension_mod = Extension('MinCompSpin', sources = ['MinCompSpin.cp312-win_amd64.pyd'])

setup(name = 'MinCompSpin', version = '0.1', ext_modules = [extension_mod])
