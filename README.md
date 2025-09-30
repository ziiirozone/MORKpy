MORKpy is used to create python bindings of the rust crate [MORK](https://github.com/ziiirozone/MORK) using [pyo3](https://github.com/PyO3/pyo3) and [maturin](https://github.com/PyO3/maturin). 
To actually create a programm that can be exported and used in a python program, one needs to install python, then activate the virtual environnement by using the command :
```
python -m venv .env
```
The program can be compiled using the command :
```
maturin build --release
```
This creates a wheel in [/target/wheels](/target/wheels), a distribution of the program that can be imported in a python program. An example of how to use a wheel is given in [example](/example).
