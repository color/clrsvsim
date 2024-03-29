try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


setup(
    name="clrsvsim",
    version="0.1.1",
    description="Color Genomics Structural Variant Simulator",
    author="Color",
    author_email="dev@color.com",
    url="https://github.com/color/clrsvsim",
    packages=["clrsvsim"],
    install_requires=[
        "cigar>=0.1.3",
        "numpy>=1.10.1",
        "pyfasta>=0.5.2",
        "pysam>=0.10.0",
    ],
    tests_require=[
        # NOTE: `mock` is not actually needed in Python 3.
        # `unittest.mock` can be used instead.
        "mock>=2.0.0",
        "nose>=1.3.7",
    ],
    license="Apache-2.0",
)
