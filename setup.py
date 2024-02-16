from setuptools import setup, find_packages

setup(
    name='gcmsnote',
    version='0.1.0',
    packages=find_packages(),
    package_data={
        'gcmsnote': ['data/gcms_lib.csv'],
    },
    include_package_data=True,
    description='A comprehensive toolkit to annotate GCMS metabolomics data.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Boris Minasenko',
    author_email='boris.minasenko@emory.edu',
    url='https://github.com/BM-Boris/gcmsnote',
    install_requires=[
        'numpy',  
        'pandas',  
        'tqdm'
        
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
)
