from pathlib import Path
from setuptools import setup, find_packages


HERE = Path(__file__).parent.resolve()
README = (HERE / "README.md").read_text(encoding="utf-8")
VERSION = (HERE / "VERSION.txt").read_text(encoding="utf-8")


setup(
    name="tidytcells",
    version=VERSION,
    description="Standardise TR/MH data",
    long_description=README,
    long_description_content_type="text/markdown",
    author="Yuta Nagano",
    author_email="yutanagano51@proton.me",
    url="https://tidytcells.readthedocs.io",
    download_url="https://github.com/yutanagano/tidytcells",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="immunology, bioinformatics, TCR, TR, MHC, MH, HLA, T cell, IMGT",
    package_dir={"": "src"},
    packages=find_packages("src"),
    package_data={"tidytcells": ["_resources/*"]},
    python_requires=">=3.7",
    extras_require={
        "dev": [
            "beautifulsoup4",
            "build",
            "odfpy>=1.4.1",
            "pandas>=1.1.5",
            "pip",
            "pytest>=7",
            "pytest-cov",
            "requests",
            "setuptools",
            "tox>=4",
            "twine",
            "wheel",
        ],
        "docs": ["sphinx-book-theme"],
    },
)
