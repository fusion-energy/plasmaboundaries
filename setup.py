import setuptools

setuptools.setup(
    name="plasma_boundaries",
    version="0.0.1",
    author="Remi Delaporte Mathurin",
    author_email="delaporte.mathurin@laposte.net",
    description="Determine plasma flux functions for various plasma parameters and plasma configurations",
    url="https://github.com/RemiTheWarrior/plasma-boundaries",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'matplotlib',
        'scipy',
        'numpy',
        'sympy',
    ],
)
