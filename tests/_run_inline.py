"""Quick inline smoke test for distance.py — no pytest needed."""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parents[1] / "src"))

import numpy as np
from paper5.distance import great_circle_km, ces_aggregate, compute_d_ovdl, ChokepointState, Agglomeration

PASS = "\033[92mPASS\033[0m"
FAIL = "\033[91mFAIL\033[0m"
errors = 0

def check(name, condition, detail=""):
    global errors
    if condition:
        print(f"  {PASS}  {name}")
    else:
        print(f"  {FAIL}  {name}  {detail}")
        errors += 1

print("=== Paper5 distance.py smoke tests ===")

d = great_circle_km(0, 0, 0, 0)
check("gc same-point = 0", d == 0)

d = great_circle_km(0, 0, 1, 0)
check(f"gc 1-deg equatorial ≈ 111.32 km  (got {d:.2f})", abs(d - 111.32) < 0.3)

d = great_circle_km(-0.1278, 51.5074, -74.006, 40.7128)
check(f"gc LHR→JFK 5500–5650 km  (got {d:.0f})", 5500 < d < 5650)

check("gc symmetry", abs(
    great_circle_km(10, 50, 100, 20) - great_circle_km(100, 20, 10, 50)
) < 1e-9)

costs = np.ones((3,3))*100; w = np.ones(3)/3
for t in (-1, 1, 2):
    r = ces_aggregate(costs, w, w, theta=t)
    check(f"ces uniform θ={t}", abs(r - 100) < 0.01, f"got {r}")

costs = np.array([[100.,200.],[300.,400.]]); w = np.array([0.5,0.5])
expected = 1/(0.25/100+0.25/200+0.25/300+0.25/400)
r = ces_aggregate(costs, w, w, theta=-1)
check(f"ces θ=-1 harmonic (got {r:.2f}, want {expected:.2f})", abs(r-expected) < 1)

r = ces_aggregate(costs, w, w, theta=1)
check(f"ces θ=1 arithmetic (got {r:.2f}, want 250.0)", abs(r - 250.0) < 0.01)

try:
    ces_aggregate(costs, w, w, theta=0)
    check("ces θ=0 raises ValueError", False)
except ValueError:
    check("ces θ=0 raises ValueError", True)

cp = ChokepointState.for_year(2023)
check(f"panama 2023 factor={cp.panama_capacity_factor}", cp.panama_capacity_factor == 1.8)
cp2 = ChokepointState.for_year(2024)
check(f"red sea 2024 risk={cp2.red_sea_risk}", cp2.red_sea_risk == 2.5)
check("normal year factor=1", ChokepointState.for_year(2015).global_risk == 1.0)

ao = [Agglomeration("A",0, 0.0, 0.0, 1.0)]
ad = [Agglomeration("B",0, 10.0, 0.0, 1.0)]
r = compute_d_ovdl(ao, ad, theta=-1)
gc = great_circle_km(0,0,10,0)
check(f"ovdl single-agglom = gc  (r={r:.1f}, gc={gc:.1f})", abs(r-gc) < 1)

ao2 = [Agglomeration("A",0,0.0,0.0,pop=1.0,viirs=10.0),
       Agglomeration("A",1,5.0,0.0,pop=9.0,viirs=1.0)]
ad2 = [Agglomeration("B",0,100.0,0.0,pop=1.0,viirs=1.0)]
d_pop = compute_d_ovdl(ao2, ad2, theta=-1, use_viirs=False)
d_viirs = compute_d_ovdl(ao2, ad2, theta=-1, use_viirs=True)
check("ovdl viirs≠pop when weights differ", abs(d_viirs - d_pop) > 1e-3,
      f"both={d_pop:.2f}")

print()
if errors:
    print(f"=== {errors} FAILURE(S) ===")
    sys.exit(1)
else:
    print(f"=== ALL {14-errors} TESTS PASS ===")
