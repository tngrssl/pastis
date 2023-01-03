from setuptools import setup

setup(name='pastis',
version='0.1.0',
description="",
author="Tangi Roussel",
packages=["pastis"],
install_requires=['pygamma','numpy', 'pandas', 'pymapvbvd', 'xlrd', 'termcolor', 'nibabel', 'IPython', 'brukerapi', 'ismrmrd']
)
