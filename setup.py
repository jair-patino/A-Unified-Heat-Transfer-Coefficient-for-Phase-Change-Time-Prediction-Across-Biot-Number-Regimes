from setuptools import setup, find_packages

setup(
    name="unified-heat-transfer-model",
    version="1.0.0",
    author="Jair Patiño B.",
    author_email="jairpalejandrov@gmail.com",
    description="Unified heat transfer coefficient model for phase change time prediction",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/jair-patino/A-Unified-Heat-Transfer-Coefficient-for-Phase-Change-Time-Prediction-Across-Biot-Number-Regimes",
    packages=find_packages(include=["scripts", "scripts.*", "tests"]),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "pandas>=1.3.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.11.0",
        "tqdm>=4.62.0",
        "pyyaml>=6.0",
        "joblib>=1.1.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.2.5",
            "pytest-cov>=3.0.0",
            "black>=22.3.0",
            "flake8>=4.0.0",
            "sphinx>=4.3.0",
        ],
        "docs": [
            "sphinx>=4.3.0",
            "sphinx-rtd-theme>=1.0.0",
            "nbsphinx>=0.8.7",
        ],
    },
    entry_points={
        "console_scripts": [
            "unified-heat-transfer=scripts.numerical_solver.validation_simulations:main",
        ],
    },
)
