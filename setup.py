from pathlib import Path
from setuptools import setup, find_packages


HERE = Path(__file__).parent.resolve()
README = (HERE/'README.rst').read_text(encoding='utf-8')
VERSION = (HERE/'VERSION.txt').read_text(encoding='utf-8')


setup(
    name='tidytcells',
    version=VERSION,
    description='Standardise TCR/MHC data.',
    long_description=README,
    long_description_content_type='text/x-rst',
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
    package_data={'tidytcells': ['_resources/*']},
    python_requires='>=3.6',
    extras_require={
        'dev': [
            'odfpy>=1.4.1',
            'pandas>=1.1.5',
            'pytest>=7.0.1',
            'pytest-cov>=4.0.0'
        ]
    }
)