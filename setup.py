from setuptools import setup, find_packages

setup(name='pyDIAUtils',
      version=0.1,
      author='Aaron Maurais',
      packages=find_packages(),
      package_dir={'pyDIAUtils': 'pyDIAUtils'},
      python_requires='>=3.8',
      install_requires=['matplotlib>=0.1.6',
                        'scikit-learn>=1.3.2',
                        'numpy>=1.26.2',
                        'pandas>=2.1.4']

)

