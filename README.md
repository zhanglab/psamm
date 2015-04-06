Metabolic modelling tools
=========================

Tools related to metabolic modelling, reconstruction, data parsing and
formatting.

See [NEWS](NEWS.md) for information on recent changes.

Use `pip` to install (it is recommended to use a
[Virtualenv](https://virtualenv.pypa.io/)):

``` shell
$ pip install git+ssh://git@github.com/zhanglab/model_script.git
```

See the bundled documentation for more information. The documentation can be
generated using [Sphinx](http://sphinx-doc.org/):

``` shell
$ git clone git@github.com:zhanglab/model_script.git
$ pip install -e model_script[docs]
$ cd model_script/docs/
$ make html
```

Then open `model_script/docs/_build/html/index.html`.
