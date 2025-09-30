#![allow(non_snake_case)]
use pyo3::prelude::*;
use pyo3::types::PyAny;

use MORK::GMORK::GMORK;
use MORK::NDMORK::NDMORK;
use MORK::RK::RK;
use MORK::create_computation_order;
use MORK::SCC;

const ERROR_FRACTION: f64 = 0.001;
const MAX_ITER: u32 = 500;
const MIN_ITER: u32 = 10;

#[pyclass]
pub struct GMORKPy {
    pub nodes: Vec<f64>,
    pub main_weights: Vec<Vec<Vec<f64>>>,
    main_weights_f: Box<dyn Fn(u32) -> Vec<Vec<f64>> + Send + Sync>,
    pub secondary_weights: Vec<Vec<Vec<f64>>>,
    secondary_weights_f: Box<dyn Fn(u32) -> Vec<Vec<f64>> + Send + Sync>,
    pub queue: Vec<SCC>,
    pub s: usize,
    pub stored_length: usize,
    pub cyclic_derivatives: Vec<Vec<bool>>, // [N-1][j]
    pub factorial: Vec<f64>,
    pub h_powers: Vec<f64>,
    pub h: f64,
    pub error_fraction: f64,
    pub min_iter: u32,
    pub max_iter: u32,
}

#[pymethods]
impl GMORKPy {
    #[new]
    pub fn new(
        main_weights_function: Py<PyAny>,
        secondary_weights_function: Py<PyAny>,
        nodes: Vec<f64>,
        maximum_weight_graph: Vec<Vec<bool>>,
    ) -> Self {
        let main_weights_function2 = Box::new(move |N| {
            Python::with_gil(|py| {
                main_weights_function
                    .call1(py, (N,))
                    .expect("Couldn't call weights function")
                    .extract::<Vec<Vec<f64>>>(py)
                    .expect("Couldn't extract weight matrix from weights function")
            })
        });
        let secondary_weights_function2 = Box::new(move |N| {
            Python::with_gil(|py| {
                secondary_weights_function
                    .call1(py, (N,))
                    .expect("Couldn't call weights function")
                    .extract::<Vec<Vec<f64>>>(py)
                    .expect("Couldn't extract weight matrix from weights function")
            })
        });
        let main_weights = vec![main_weights_function2(1)];
        let secondary_weights = vec![secondary_weights_function2(1)];
        let s = main_weights[0].len() - 1;
        let mut cyclic_derivatives: Vec<Vec<bool>> = vec![vec![false; s]];
        let queue = create_computation_order(&maximum_weight_graph);
        for task in queue.iter() {
            if let SCC::Implicit(J, _) = task {
                for &j in J {
                    for &j1 in J {
                        if main_weights[0][j][j1] != 0. {
                            cyclic_derivatives[0][j] = true;
                            break;
                        }
                    }
                }
            }
        }
        GMORKPy {
            s,
            stored_length: 1,
            nodes,
            main_weights,
            main_weights_f: main_weights_function2,
            secondary_weights,
            secondary_weights_f: secondary_weights_function2,
            factorial: vec![1., 1.],
            h: 0.,
            h_powers: vec![1., 0.],
            queue,
            cyclic_derivatives,
            error_fraction: ERROR_FRACTION,
            min_iter: MIN_ITER,
            max_iter: MAX_ITER,
        }
    }
}

impl GMORKPy {
    pub fn approximate(
        &mut self,
        t: f64,
        h: f64,
        f: &dyn Fn(f64, &Vec<Vec<f64>>) -> Vec<f64>,
        y0: &Vec<Vec<f64>>,
    ) -> Vec<Vec<f64>> {
        if h != self.h {
            self.h = h;
            for N in 1..=self.stored_length {
                self.h_powers[N] = h.powi(N as i32)
            }
        }
        // calculate difference threshold for picard iterations
        let mut threshold = y0[0][0].abs();
        for k in 0..y0.len() {
            // verify the length of the method is enough
            if y0[k].len() > self.stored_length {
                GMORK::set_minimum_length(
                    y0[k].len(),
                    self.s,
                    &mut self.stored_length,
                    &mut self.factorial,
                    &mut self.main_weights,
                    &self.main_weights_f,
                    &mut self.secondary_weights,
                    &self.secondary_weights_f,
                    &mut self.h_powers,
                    h,
                    &mut self.cyclic_derivatives,
                    &self.queue,
                );
            }
            for N in 0..y0[k].len() {
                if threshold < y0[k][N].abs() {
                    threshold = y0[k][N].abs()
                }
            }
        }
        threshold *= self.error_fraction;
        GMORK::approximate_GMORK(
            t,
            h,
            f,
            y0,
            threshold,
            self.s,
            &self.queue,
            &self.nodes,
            &self.main_weights,
            &self.secondary_weights,
            &self.h_powers,
            &self.factorial,
            self.min_iter,
            self.max_iter,
            &self.cyclic_derivatives,
        )
    }
}

#[pyclass]
pub struct NDMORKPy {
    pub s: usize,
    pub stored_length: usize,
    pub nodes: Vec<f64>,
    pub weights: Vec<Vec<Vec<f64>>>,
    pub weights_function: Box<dyn Fn(u32) -> Vec<Vec<f64>> + Send + Sync>,
    pub factorial: Vec<f64>,
    pub coefficients: Vec<Vec<f64>>,
    pub h: f64,
    pub h_powers: Vec<f64>,
    pub queue: Vec<SCC>,
    pub cycle_derivative: Vec<Vec<bool>>, // [N-1][j]
    pub error_fraction: f64,
    pub min_iter: u32,
    pub max_iter: u32,
}

#[pymethods]
impl NDMORKPy {
    #[new]
    pub fn new(
        weights_function: Py<PyAny>,
        nodes: Vec<f64>,
        maximum_weight_graph: Vec<Vec<bool>>,
    ) -> Self {
        let weights_function2 = Box::new(move |N| {
            Python::with_gil(|py| {
                weights_function
                    .call1(py, (N,))
                    .expect("Couldn't call weights function")
                    .extract::<Vec<Vec<f64>>>(py)
                    .expect("Couldn't extract weight matrix from weights function")
            })
        });
        let s = nodes.len() - 1;
        let weights = vec![weights_function2(1)];
        let queue = create_computation_order(&maximum_weight_graph);
        let mut cycle_derivative: Vec<Vec<bool>> = vec![vec![false; s]];
        for task in queue.iter() {
            if let SCC::Implicit(J, _) = task {
                for &j in J {
                    for &j1 in J {
                        if weights[0][j][j1] != 0. {
                            cycle_derivative[0][j] = true;
                            break;
                        }
                    }
                }
            }
        }
        NDMORKPy {
            s,
            stored_length: 1,
            nodes,
            weights,
            weights_function: weights_function2,
            factorial: vec![1., 1.],
            coefficients: vec![vec![1.; s + 1]],
            h: 0.,
            h_powers: vec![1., 0.],
            queue,
            cycle_derivative,
            error_fraction: ERROR_FRACTION,
            min_iter: MIN_ITER,
            max_iter: MAX_ITER,
        }
    }

    fn approximate(&mut self, t: f64, h: f64, f: Py<PyAny>, y0: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
        let f_py: &dyn Fn(f64, &Vec<Vec<f64>>) -> Vec<f64> = &|t, y| {
            Python::with_gil(|py| {
                f.call1(py, (t, y))
                    .expect("Couldn't call f function")
                    .extract::<Vec<f64>>(py)
                    .expect("Couldn't extract result of f function as Vec<f64>")
            })
        };
        if h == 0. {
            return y0.clone();
        }
        if h != self.h {
            self.h = h;
            for N in 1..=self.stored_length {
                self.h_powers[N] = h.powi(N as i32)
            }
        }
        // calculate difference threshold for picard iterations
        let mut threshold = y0[0][0].abs();
        for k in 0..y0.len() {
            // verify the length of the method is enough
            // /*
            if y0[k].len() > self.stored_length {
                NDMORK::set_minimum_length(
                    &mut self.h_powers,
                    &self.weights_function,
                    &mut self.weights,
                    h,
                    self.s,
                    &mut self.factorial,
                    &mut self.stored_length,
                    y0[k].len(),
                    &mut self.cycle_derivative,
                    &self.queue,
                    &mut self.coefficients,
                    &self.nodes,
                );
            }
            // */
            for N in 0..y0[k].len() {
                if threshold < y0[k][N].abs() {
                    threshold = y0[k][N].abs()
                }
            }
        }
        threshold *= self.error_fraction;
        NDMORK::approximate_ND(
            t,
            h,
            f_py,
            &y0,
            threshold,
            self.min_iter,
            self.max_iter,
            self.s,
            &self.queue,
            &self.cycle_derivative,
            &self.nodes,
            &self.weights,
            &self.coefficients,
            &self.h_powers,
            &self.factorial,
        )
    }
}

#[pyclass]
pub struct RKPy {
    pub s: usize,
    pub nodes: Vec<f64>,
    pub weights: Vec<Vec<f64>>,
    pub queue: Vec<SCC>,
    pub error_fraction: f64,
    pub min_iter: u32,
    pub max_iter: u32,
}

#[pymethods]
impl RKPy {
    #[new]
    pub fn new(weights: Vec<Vec<f64>>, nodes: Vec<f64>) -> Self {
        let s = weights.len() - 1;
        let weight_graph = (0..=s)
            .map(|j| {
                (0..=s)
                    .map(|j1| {
                        if j1 == s || weights[j][j1] == 0. {
                            false
                        } else {
                            true
                        }
                    })
                    .collect()
            })
            .collect();
        let queue = create_computation_order(&weight_graph);
        let s = nodes.len() - 1;
        RKPy {
            s,
            nodes,
            weights,
            queue,
            error_fraction: ERROR_FRACTION,
            min_iter: MIN_ITER,
            max_iter: MAX_ITER,
        }
    }

    fn approximate(&mut self, t: f64, h: f64, f: Py<PyAny>, y0: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
        if h == 0. {
            return y0.clone();
        }
        let f_py: &dyn Fn(f64, &Vec<Vec<f64>>) -> Vec<f64> = &|t, y| {
            Python::with_gil(|py| {
                f.call1(py, (t, y))
                    .expect("Couldn't call f function")
                    .extract::<Vec<f64>>(py)
                    .expect("Couldn't extract result of f function as Vec<f64>")
            })
        };
        // calculate difference threshold for picard iterations
        let mut threshold = y0[0][0].abs();
        for k in 0..y0.len() {
            for N in 0..y0[k].len() {
                if threshold < y0[k][N].abs() {
                    threshold = y0[k][N].abs()
                }
            }
        }
        threshold *= self.error_fraction;
        RK::approximate_RK(
            t,
            h,
            f_py,
            &y0,
            self.s,
            threshold,
            &self.queue,
            &self.weights,
            &self.nodes,
            self.min_iter,
            self.max_iter,
        )
    }
}

#[pymodule]
fn MORKpy(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<NDMORKPy>()?;
    m.add_class::<RKPy>()?;
    Ok(())
}
