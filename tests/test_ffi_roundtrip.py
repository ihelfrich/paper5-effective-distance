"""test_ffi_roundtrip.py — verify the Rust↔Python FFI surface.

These tests only run once `paper5_core` is installed via

    maturin develop --release --manifest-path crates/paper5-py/Cargo.toml

On CI and in dev, they are skipped if the extension module is absent so we
don't block pytest on Rust toolchain availability during pure-Python work.
"""
from __future__ import annotations

import numpy as np
import pytest

paper5_core = pytest.importorskip("paper5_core")


def _triangle_graph():
    """Three road nodes, three edges — direct A→C is slower than A→B→C."""
    node_lon = np.array([0.0, 1.0, 1.0], dtype=np.float32)
    node_lat = np.array([0.0, 0.0, 1.0], dtype=np.float32)
    node_mode = np.zeros(3, dtype=np.uint8)  # all road

    edge_from = np.array([0, 1, 0], dtype=np.uint32)
    edge_to = np.array([1, 2, 2], dtype=np.uint32)
    edge_length_km = np.array([100.0, 100.0, 300.0], dtype=np.float32)
    edge_time_hours = np.array([1.0, 1.0, 3.0], dtype=np.float32)
    edge_mode = np.zeros(3, dtype=np.uint8)  # all road

    return paper5_core.build_graph(
        node_lon, node_lat, node_mode,
        edge_from, edge_to, edge_length_km, edge_time_hours, edge_mode,
    )


def test_version_exposed():
    assert paper5_core.__version__.startswith("0.1")


def test_triangle_solve_pair():
    g = _triangle_graph()
    cp = paper5_core.ChokepointHandle()
    # 0->1->2 is cheaper than 0->2 direct (2h vs 3h @ $30/h road rate).
    cost = g.solve_pair(0, 2, cp)
    assert abs(cost - 60.0) < 1e-3


def test_triangle_many_to_many():
    g = _triangle_graph()
    cp = paper5_core.ChokepointHandle()
    sources = np.array([0, 1], dtype=np.uint32)
    targets = np.array([1, 2], dtype=np.uint32)
    m = g.solve_many_to_many(sources, targets, cp)
    assert m.shape == (2, 2)
    # 0->1 = 1h*30 = 30; 0->2 = 60; 1->1 = 0; 1->2 = 30.
    np.testing.assert_allclose(m, [[30.0, 60.0], [0.0, 30.0]], atol=1e-3)


def test_chokepoint_preset_constructors():
    cp_suez = paper5_core.ChokepointHandle.suez_ever_given_2021()
    cp_panama = paper5_core.ChokepointHandle.panama_drought_2023()
    cp_houthi = paper5_core.ChokepointHandle.houthi_red_sea_2024()
    # Smoke-test: they exist and can be passed to solve_pair without error.
    g = _triangle_graph()
    for cp in (cp_suez, cp_panama, cp_houthi):
        _ = g.solve_pair(0, 2, cp)
