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

from setuptools import find_packages
from setuptools import setup

setup(
    name='alphafold',
    version='2.3.1',
    description='An implementation of the inference pipeline of AlphaFold v2.0.'
    'This is a completely new model that was entered as AlphaFold2 in CASP14 '
    'and published in Nature.',
    author='DeepMind',
    author_email='alphafold@deepmind.com',
    license='Apache License, Version 2.0',
    url='https://github.com/deepmind/alphafold',
    packages=find_packages(),
    install_requires=[
        'absl-py==1.0.0',
        'biopython==1.79',
        'chex==0.0.7',
        'dm-haiku==0.0.7',
        'dm-tree==0.1.6',
        'immutabledict==2.0.0',
        'jax==0.3.17',
        'ml-collections==0.1.0',
        'numpy==1.21.6',
        'pandas==1.3.4',
        'protobuf==3.20.1',
        'scipy==1.7.0',
        'tensorflow-cpu==2.9.0'
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
        # Python version <= 3.7 is not recommended due to multithread issue.
        #'Programming Language :: Python :: 3.6',
        #'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
    ],
)
