from setuptools import setup, find_packages

with open("lcvk/__init__.py") as f:
    exec([x for x in f.readlines() if '__version__' in x][0])

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as f:
    requirements = f.read()


setup(
  name='lcvk',
  version=__version__,
  author='Shuzhao Li',
  author_email='shuzhao.li@gmail.com',
  description='Li-circular van Krevelen diagram',
  long_description=long_description,
  long_description_content_type="text/markdown",
  url='https://github.com/shuzhao-li-lab/Li_CVK_diagram',
  license='BSD 3-Clause',
  keywords='bioinformatics metabolic pathway visualization',

  classifiers=[
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Topic :: Software Development :: Libraries :: Python Modules',
  ],

  packages=find_packages(
    include=['*', '']
  ),
  data_files=[  ],
  include_package_data=True,
  zip_safe=True,
  entry_points = {
    },

  python_requires='>=3.7',
  install_requires=requirements,

)
