import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='dlotter',
     version='0.0.2.1',
     author="Kasper Hintz",
     author_email="kah@dmi.dk",
     description="Quick and dirty static plots of NWP output",
     long_description=long_description,
     long_description_content_type="text/markdown",
     packages=setuptools.find_packages(),
     setup_requires=[
         'wheel',
         ],
     install_requires=[
         'wheel',
         'dmit',
         'xarray',
         'numpy',
         'cartopy',
         'eccodes',
         'pygrib',
         ],
     url="https://dmidk.github.io/dlotter/",
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ]
 )