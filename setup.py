# Not sure whether this works or not. Modified from atotickov's fork and AI suggestion.
from pathlib import Path
from setuptools import setup, find_packages

setup(
    name="quarTeT",
    version="1.2.5",
    author="Echoring",
    author_email="linyunzhi20@gmail.com",
    description="A telomere-to-telomere toolkit for gap-free genome assembly and centromeric repeat identification",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/aaranyue/quarTeT",
    packages=find_packages(),
    scripts=list(map(str, sorted(Path('./').rglob("*.py")))),
    python_requires='>=3.6',
    # classifiers=[
    #     "Programming Language :: Python :: 3",
    #     "License :: OSI Approved :: MIT License",
    #     "Operating System :: OS Independent",
    # ],
    # install_requires=[
    #     "minimap2>=2.24",
    #     "unimap>=0.1",
    #     "mummer4>=4.0.0rc1",
    #     "trf>=4.09",
    #     "cd-hit>=4.8.1",
    #     "blast>=2.11.0",
    #     "tidk>=0.2.31",
    #     "gnuplot>=5.4",
    #     "r-base>=3.5.0",
    #     "r-rideogram>=0.2.2",
    #     "r-ggplot2>=3.4.4"
    # ],
    # entry_points={
    #     'console_scripts': [
    #         'quartet=quartet:main',
    #     ],
    # },
)
