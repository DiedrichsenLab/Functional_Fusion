# Functional_Fusion
This repository provides the data structures and code to efficiently implement multi-data set project.
It is written and maintained by the Diedrichsen Lab, Western University.


## Installation and dependencies
This project depends on several third party libraries, including:

[numpy](https://numpy.org/) (version>=1.22.2)

nibabel []

[nilearn](https://nilearn.github.io/stable/index.html) (version>=0.9.0), ...

	pip install numpy nilearn ...

[nitools]
    pip install neuroimagingtools

Or you can install the package manually from the original binary source as above links.

Once you clone the functional fusion repository, you need to add it to your PYTHONPATH, so you can import the functionality. Add these lines to your .bash_profile, .bash_rc .zsh_profile file...

```
PYTHONPATH=<abspath_of_repo_parentdir>:${PYTHONPATH}
export PYTHONPATH
```

## Documentation

For detailed documentation see: [https://functional-fusion.readthedocs.io/en/latest/](https://functional-fusion.readthedocs.io/en/latest/).

Source files for the online documentation can be found in the docs folder.
