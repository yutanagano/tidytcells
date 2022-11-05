from pathlib import Path
from setuptools import setup, find_packages


here = Path(__file__).parent.resolve()


readme = (here/'README.md').read_text(encoding='utf-8')


setup(
    name='tidytcells',
    version='0.0.3',
    description='Standardise TCR/MHC gene names to IMGT nomenclature.',
    long_description=readme,
    long_description_content_type='text/markdown',
    author='Yuta Nagano',
    author_email='yutanagano51@proton.me',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    keywords='immunology, bioinformatics, TCR, MHC, HLA, IMGT',
    package_dir={'': 'src'},
    packages=find_packages('src'),
    package_data={'tidytcells': ['resources/*']},
    python_requires='>=3.10'
)