PSAMM metabolic modeling tools
==============================

Tools related to metabolic modeling, reconstruction, data parsing and
formatting, consistency checking, automatic gap filling, and model simulations.

See [NEWS](NEWS.md) for information on recent changes. The `master` branch
tracks the latest release while the `develop` branch is the latest version in
development. Please apply any pull requests to the `develop` branch when
creating the pull request.

Install
-------

Use `pip` to install (it is recommended to use a
[Virtualenv](https://virtualenv.pypa.io/)):

``` shell
$ pip install git+ssh://git@github.com/zhanglab/model_script.git
```

Documentation
-------------

See the bundled documentation for more information. The documentation can be
generated using [Sphinx](http://sphinx-doc.org/):

``` shell
$ git clone git@github.com:zhanglab/model_script.git
$ pip install -e model_script[docs]
$ cd model_script/docs/
$ make html
```

Then open `model_script/docs/_build/html/index.html`.
