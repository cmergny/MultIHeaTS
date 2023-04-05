from setuptools import setup, find_packages

# read the contents of your README file
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="multiheats",
    version="0.401",
    packages=find_packages("src"),
    package_dir={"": "src"},
    readme="README.md",
    description="MultIHeaTS is a Multi-layered Implicit Heat Transfer Solver.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        "wheel",
        "ipython",
        "pytest",
        "numpy",
        "matplotlib",
        "scipy",
        "pandas",
        "tqdm",
        "xlrd",
    ],
)
