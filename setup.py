from setuptools import setup, find_packages

setup(
    name="multiheats",
    version="0.1",
    packages=find_packages("src"),
    package_dir={"": "src"},
    install_requires=[
        # list any dependencies your package has here
    ],
)
