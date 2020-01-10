from setuptools import setup
from setuptools import find_packages

setup(name='gcncc',
      version='0.1',
      description='Graph Convolutional Network for Clustering and Classification',
      author='Omar Maddouri',
      author_email='omar.maddouri@gmail.com',
      url='https://github.com/omarmaddouri',
      download_url='...',
      license='MIT',
      install_requires=['keras'],
      extras_require={
          'model_saving': ['json', 'h5py'],
      },
      package_data={'gcncc': ['README.md']},
      packages=find_packages())