import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
      name='acidoseq',
      version='1.3.6',
      description='Studying Acidobacteria reads',
      url='https://github.com/sap218/acidoseq',
      author='Samanthe C Pendle',
      author_email='samanfapendle@outlook.com',
      license='MIT',
      packages=setuptools.find_packages(),
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"],
      zip_safe=False,
      entry_points = {'console_scripts': ['acidoseq=acidoseq.acidoseq:main', 'acidomap=map.map:main']},
      include_package_data = True
)
