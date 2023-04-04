from setuptools import setup, find_packages

setup(
    name="multiheats",
    version="0.2",
    packages=find_packages("src"),
    package_dir={"": "src"},
    description="MultIHeaTS is a Multi-layered Implicit Heat Transfer Solver. It uses finite differences to  simulates and predicts the surface temperature in 1D multi-layered planetary surfaces exposed to solar radiation.",
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
