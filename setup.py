from setuptools import setup, find_packages

setup(
    name='memdiff',
    packages=find_packages(include=[
        'memdiff',
        'memdiff.oseen',
        'memdiff.immersed',
        'memdiff.analysis',
    ]),
    version='0.1.0',
    description='',
    author='Martin Vögele',
    license='MIT',
    install_requires=[
        'argparse',
        'numpy',
        'scipy',
        'matplotlib',
    ],
)
