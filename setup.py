from setuptools import find_packages, setup

with open("README.rst", "r") as fh:
    long_description = fh.read()

# Get admin.
admin = {}
with open("hvec_flow/admin.py") as fp:
    exec(fp.read(), admin)

setup(
    name = 'hvec_flow',
    version = admin['__version__'],
    author = admin['__author__'],
    author_email = admin['__author_email__'],
    description = 'Python package for hydraulic engineers '
                'dealing with flow and turbulence.',
    long_description=long_description,
#    url='https://github.com/pastas/pastas',
#    project_urls={
#        'Source': 'https://github.com/pastas/pastas',
#        'Documentation': 'http://pastas.readthedocs.io/en/latest/',
#        'Tracker': 'https://github.com/pastas/pastas/issues',
#        'Help': 'https://github.com/pastas/pastas/discussions'
#    },
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Other Audience',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering'
    ],
    platforms='Windows',
    install_requires=['numpy>=1.17'],
    packages=find_packages(exclude=[]),
)