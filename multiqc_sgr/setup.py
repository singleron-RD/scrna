from setuptools import find_packages, setup

version = "1.0.0"

setup(
    name="multiqc_sgr",
    version=version,
    author="Yiqi Zhou",
    description="MultiQC plugin for Singleron",
    long_description=__doc__,
    keywords="bioinformatics",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["multiqc==1.21"],
    entry_points={
        "multiqc.modules.v1": [
            "scsnp = multiqc_sgr.scsnp:MultiqcModule",
            "scrna = multiqc_sgr.scrna:MultiqcModule",
        ],
        "multiqc.hooks.v1": [
            "before_config = multiqc_sgr:multiqc_sgr_config",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Environment :: Web Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: JavaScript",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
