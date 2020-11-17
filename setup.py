import setuptools

setuptools.setup(
    name="plasmaboundaries",
    version="0.1.5",
    author="Remi Delaporte-Mathurin",
    author_email="rdelaportemathurin@gmail.com",
    description="Determine plasma flux functions for various plasma parameters and plasma configurations",
    url="https://github.com/RemiTheWarrior/plasma-boundaries",
    packages=setuptools.find_packages(),
    package_dir={"plasmaboundaries": "plasmaboundaries"},
    data_files=[("", ["LICENSE.txt"])],
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
    zip_safe=True,
)
