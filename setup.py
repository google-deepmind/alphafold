# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Install script for setuptools."""

from alphafold import version
from setuptools import find_packages
from setuptools import setup

setup(
    name='alphafold',
    version=version.__version__,
    description=(
        'An implementation of the inference pipeline of AlphaFold v2.0. This is'
        ' a completely new model that was entered as AlphaFold2 in CASP14 and'
        ' published in Nature.'
    ),
    author='DeepMind',
    author_email='alphafold@deepmind.com',
    license='Apache License, Version 2.0',
    url='https://github.com/deepmind/alphafold',
    packages=find_packages(),
    install_requires=[
        'absl-py',
        'biopython',
        'chex',
        'dm-haiku',
        'dm-tree',
        'docker',
        'immutabledict',
        'jax',
        'ml-collections',
        'numpy',
        'pandas',
        'scipy',
        'tensorflow-cpu',
    ],
    tests_require=[
        'matplotlib',  # For notebook_utils_test.
        'mock',
    ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
    ],
)
