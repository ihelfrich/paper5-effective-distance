//! PyO3 bindings for paper5-core. Exposes the LCP solver to Python as
//! `import paper5_core`.
//!
//! Built with `maturin develop --release --manifest-path crates/paper5-py/Cargo.toml`.

use numpy::{IntoPyArray, PyArray2, PyReadonlyArray1};
use paper5_core::{
    chokepoints::ChokepointState, graph::{Edge, Graph, Mode, Node}, solve_many_to_many,
    solve_pair,
};
use pyo3::prelude::*;

/// Build a `Graph` from flat Python arrays. Arrays are zero-copy-viewed but
/// we allocate owned Vec<T> for the Rust `Graph` since it's held across
/// many queries and should outlive the Python caller's arrays.
#[pyfunction]
fn build_graph(
    node_lon: PyReadonlyArray1<f32>,
    node_lat: PyReadonlyArray1<f32>,
    node_mode: PyReadonlyArray1<u8>,
    edge_from: PyReadonlyArray1<u32>,
    edge_to: PyReadonlyArray1<u32>,
    edge_length_km: PyReadonlyArray1<f32>,
    edge_time_hours: PyReadonlyArray1<f32>,
    edge_mode: PyReadonlyArray1<u8>,
) -> PyResult<GraphHandle> {
    let lon = node_lon.as_slice()?;
    let lat = node_lat.as_slice()?;
    let nmode = node_mode.as_slice()?;
    let nodes: Vec<Node> = (0..lon.len() as u32)
        .map(|i| Node {
            id: i,
            lon: lon[i as usize],
            lat: lat[i as usize],
            mode: mode_from_u8(nmode[i as usize]),
        })
        .collect();

    let ef = edge_from.as_slice()?;
    let et = edge_to.as_slice()?;
    let el = edge_length_km.as_slice()?;
    let eh = edge_time_hours.as_slice()?;
    let em = edge_mode.as_slice()?;
    let edges: Vec<Edge> = (0..ef.len())
        .map(|i| Edge {
            from: ef[i],
            to: et[i],
            length_km: el[i],
            time_hours: eh[i],
            mode: mode_from_u8(em[i]),
        })
        .collect();

    Ok(GraphHandle { inner: Graph::from_edges(nodes, edges) })
}

fn mode_from_u8(x: u8) -> Mode {
    match x {
        0 => Mode::Road,
        1 => Mode::Maritime,
        2 => Mode::Air,
        _ => Mode::Transfer,
    }
}

#[pyclass]
pub struct GraphHandle {
    inner: Graph,
}

#[pymethods]
impl GraphHandle {
    #[getter]
    fn num_nodes(&self) -> usize {
        self.inner.num_nodes()
    }
    #[getter]
    fn num_edges(&self) -> usize {
        self.inner.num_edges()
    }

    /// Single-pair cost. Chokepoint state is specified via kwargs (see
    /// `make_chokepoint_state`).
    fn solve_pair(&self, src: u32, dst: u32, cp: &ChokepointHandle) -> f32 {
        solve_pair(&self.inner, src, dst, &cp.inner)
    }

    /// Many-to-many: returns a (len(sources), len(targets)) 2D numpy array.
    fn solve_many_to_many<'py>(
        &self,
        py: Python<'py>,
        sources: PyReadonlyArray1<u32>,
        targets: PyReadonlyArray1<u32>,
        cp: &ChokepointHandle,
    ) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let s = sources.as_slice()?;
        let t = targets.as_slice()?;
        let flat = solve_many_to_many(&self.inner, s, t, &cp.inner);
        let arr = ndarray::Array2::from_shape_vec((s.len(), t.len()), flat)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e.to_string()))?;
        Ok(arr.into_pyarray_bound(py))
    }
}

#[pyclass]
pub struct ChokepointHandle {
    inner: ChokepointState,
}

#[pymethods]
impl ChokepointHandle {
    #[new]
    #[pyo3(signature = (suez_open=true, panama_capacity_factor=1.0, red_sea_risk=1.0, global_risk=1.0))]
    fn new(
        suez_open: bool,
        panama_capacity_factor: f32,
        red_sea_risk: f32,
        global_risk: f32,
    ) -> Self {
        Self {
            inner: ChokepointState {
                suez_open,
                panama_capacity_factor,
                red_sea_risk,
                global_risk,
            },
        }
    }

    #[staticmethod]
    fn suez_ever_given_2021() -> Self {
        Self { inner: ChokepointState::suez_ever_given_2021() }
    }
    #[staticmethod]
    fn panama_drought_2023() -> Self {
        Self { inner: ChokepointState::panama_drought_2023() }
    }
    #[staticmethod]
    fn houthi_red_sea_2024() -> Self {
        Self { inner: ChokepointState::houthi_red_sea_2024() }
    }
}

#[pymodule]
fn paper5_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", paper5_core::VERSION)?;
    m.add_class::<GraphHandle>()?;
    m.add_class::<ChokepointHandle>()?;
    m.add_function(wrap_pyfunction!(build_graph, m)?)?;
    Ok(())
}
