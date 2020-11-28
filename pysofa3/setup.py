import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pysofa3", # Replace with your own username
    version="0.0.1",
    author="LiYangLiu",
    author_email="895479558@qq.com",
    description="A library of functions on astronomy",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/YYjoke/pysofa.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)