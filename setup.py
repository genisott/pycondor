from setuptools import setup, find_packages

setup(name='condor',
    version='2.0',
    description='Python implementation of CONDOR with the BRIM algorithm.',
    url='https://github.com/genisott/pycondor',
    author='Genís Calderer',
    author_email='genis.calderer@gmail.com',
    license='MIT',
    packages=['condor'],
    install_requires=['pandas',
    'numpy',
    'python-igraph'],
    zip_safe=False)
