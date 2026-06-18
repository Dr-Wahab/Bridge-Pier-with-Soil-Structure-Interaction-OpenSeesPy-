"""
===============================================================================
 Single-column RC bridge pier on a large-diameter drilled-shaft foundation
 with soil-structure interaction (SSI) modelled by Beam-on-Nonlinear-Winkler-
 Foundation (BNWF) p-y / t-z / q-z springs.

 Analyses:   (1) gravity   (2) modal   (3) displacement-controlled pushover
             (4) nonlinear time-history (uniform free-field excitation)

 Units: kN, m, s, Mg (tonne), kPa.   g = 9.81 m/s^2
 OpenSeesPy.  2D model: ndm=2, ndf=3 (UX, UY, RZ).  Y is vertical (up).
===============================================================================
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import openseespy.opensees as ops

OUT = "/home/claude"
g = 9.81

# ----------------------------------------------------------------------------
# 1. GLOBAL GEOMETRY AND MATERIAL/SOIL PARAMETERS
# ----------------------------------------------------------------------------
# --- Pier column (above ground) ---
H_col   = 8.0          # column clear height [m]
D_col   = 1.5          # column diameter [m]
n_col   = 8            # column elements
# --- Drilled shaft / pile (below ground) ---
L_pile  = 15.0         # embedded shaft length [m]
D_pile  = 1.8          # shaft diameter [m]
n_pile  = 15           # pile elements (1 m each -> soil node every 1 m)
# --- Concrete (C35) ---
fc      = 35.0e3       # f'c [kPa] (unconfined)
Ec      = 5000.0*np.sqrt(35.0)*1.0e3            # [kPa]
fcc     = 1.30*fc                                # confined strength (Mander-type)
eps_c0  = 0.002; eps_cu = 0.005                  # unconfined strains
eps_cc  = 0.006; eps_ccu = 0.020                 # confined strains
ft      = 0.10*fc                                 # tensile strength
# --- Reinforcing steel (Grade 420) ---
fy      = 420.0e3      # [kPa]
Es      = 200.0e6      # [kPa]
b_hard  = 0.01         # strain-hardening ratio
cover   = 0.05         # clear cover [m]
# longitudinal reinforcement
n_bar_col   = 28;  dbar_col  = 0.032
n_bar_pile  = 32;  dbar_pile = 0.036
Abar_col    = np.pi/4*dbar_col**2
Abar_pile   = np.pi/4*dbar_pile**2
# --- masses / densities ---
rho_c   = 2.4          # concrete mass density [Mg/m^3]
A_col   = np.pi/4*D_col**2
A_pile  = np.pi/4*D_pile**2
m_col_pl  = rho_c*A_col      # mass / length [Mg/m]
m_pile_pl = rho_c*A_pile
M_super = 600.0        # tributary superstructure mass on pier [Mg] (~5.9 MN)

# --- Soil profile (single medium-dense sand layer; water table at surface) ---
phi_deg = 35.0
gamma_e = 10.0         # effective unit weight [kN/m^3]
k_sub   = 16300.0      # initial subgrade modulus for p-y [kN/m^3]
Kc      = 0.8          # lateral earth-pressure coeff. for t-z
delta   = 0.8*phi_deg  # pile-soil interface friction angle [deg]
Nq      = 40.0         # bearing capacity factor for shaft tip

# ----------------------------------------------------------------------------
# 2. SOIL SPRING CAPACITY FUNCTIONS (API RP2A, illustrative calibration)
# ----------------------------------------------------------------------------
def api_sand_py_coeffs(phi_d):
    phi = np.radians(phi_d)
    Ko  = 0.4
    Ka  = np.tan(np.radians(45-phi_d/2))**2
    al  = phi/2.0
    be  = np.radians(45+phi_d/2)
    C1  = (Ko*np.tan(phi)*np.sin(be))/(np.tan(be-phi)*np.cos(al)) \
          + (np.tan(be)**2*np.tan(al))/np.tan(be-phi) \
          + Ko*np.tan(be)*(np.tan(phi)*np.sin(be)-np.tan(al))
    C2  = np.tan(be)/np.tan(be-phi) - Ka
    C3  = Ka*(np.tan(be)**8 - 1.0) + Ko*np.tan(phi)*np.tan(be)**4
    return C1, C2, C3

C1, C2, C3 = api_sand_py_coeffs(phi_deg)

def py_capacity(depth, D, trib):
    """PySimple1 (pult [kN], y50 [m]) for sand at given depth."""
    pus = (C1*depth + C2*D)*gamma_e*depth      # shallow [kN/m]
    pud = C3*D*gamma_e*depth                    # deep    [kN/m]
    pu  = max(min(pus, pud), 1.0)
    A   = 0.9                                    # cyclic/static factor
    pult = A*pu*trib
    y50  = np.clip(0.5*A*pu/(k_sub*max(depth,0.5)), 0.002, 0.04)
    return pult, y50

def tz_capacity(depth, D, trib):
    """TzSimple1 (tult [kN], z50 [m]) skin friction, sand."""
    sv   = gamma_e*depth
    fs   = Kc*sv*np.tan(np.radians(delta))      # unit skin friction [kPa]
    fs   = min(fs, 100.0)
    tult = max(fs*np.pi*D*trib, 1.0)
    return tult, 0.0025

def qz_capacity(depth, D):
    """QzSimple1 (qult [kN], z50 [m]) end bearing, sand."""
    sv   = gamma_e*depth
    qp   = min(Nq*sv, 12.0e3)                    # cap 12 MPa
    qult = qp*(np.pi/4*D**2)
    return qult, 0.05*D

# ----------------------------------------------------------------------------
# 3. MODEL BUILDER  (returns geometry dictionaries for plotting)
# ----------------------------------------------------------------------------
TAG = dict(coreCol=1, coverCol=2, corePile=3, coverPile=4, steel=5,
           secCol=10, secPile=11, transf=1)

def build_model():
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 3)

    # ---- materials: concrete & steel ----
    ops.uniaxialMaterial('Concrete02', TAG['coreCol'],  -fcc, -eps_cc, -0.2*fcc, -eps_ccu, 0.1, ft, ft/eps_c0)
    ops.uniaxialMaterial('Concrete02', TAG['coverCol'], -fc,  -eps_c0, 0.0,      -eps_cu,  0.1, ft, ft/eps_c0)
    ops.uniaxialMaterial('Concrete02', TAG['corePile'], -fcc, -eps_cc, -0.2*fcc, -eps_ccu, 0.1, ft, ft/eps_c0)
    ops.uniaxialMaterial('Concrete02', TAG['coverPile'],-fc,  -eps_c0, 0.0,      -eps_cu,  0.1, ft, ft/eps_c0)
    ops.uniaxialMaterial('Steel02', TAG['steel'], fy, Es, b_hard, 18.0, 0.925, 0.15)

    # ---- fiber sections (circular) ----
    def circ_section(secTag, D, coreMat, coverMat, nbar, Abar):
        r_o  = D/2.0
        r_c  = r_o - cover                # core radius (to bar centre ~ cover)
        ops.section('Fiber', secTag)
        ops.patch('circ', coreMat, 12, 4, 0., 0., 0.0,  r_c, 0., 360.)
        ops.patch('circ', coverMat,12, 2, 0., 0., r_c,  r_o, 0., 360.)
        ops.layer('circ', TAG['steel'], nbar, Abar, 0., 0., r_c-0.01)
    circ_section(TAG['secCol'],  D_col,  TAG['coreCol'],  TAG['coverCol'],  n_bar_col,  Abar_col)
    circ_section(TAG['secPile'], D_pile, TAG['corePile'], TAG['coverPile'], n_bar_pile, Abar_pile)

    ops.geomTransf('PDelta', TAG['transf'])
    ops.beamIntegration('Lobatto', TAG['secCol'],  TAG['secCol'],  5)
    ops.beamIntegration('Lobatto', TAG['secPile'], TAG['secPile'], 5)

    # ---- nodes (top -> tip), one per metre ----
    node_coords = {}
    ys_col  = np.linspace(H_col, 0.0, n_col+1)         # +8 ... 0
    ys_pile = np.linspace(0.0, -L_pile, n_pile+1)[1:]  # -1 ... -15
    ys = np.concatenate([ys_col, ys_pile])
    struct_nodes = []
    for i, y in enumerate(ys, start=1):
        ops.node(i, 0.0, float(y))
        node_coords[i] = (0.0, float(y))
        struct_nodes.append(i)
    top_node  = struct_nodes[0]
    grnd_node = struct_nodes[n_col]      # y = 0
    tip_node  = struct_nodes[-1]

    # ---- beam-column elements ----
    beam_elems = []
    for e in range(len(struct_nodes)-1):
        ni, nj = struct_nodes[e], struct_nodes[e+1]
        integ  = TAG['secCol'] if e < n_col else TAG['secPile']
        eid    = 100+e
        ops.element('forceBeamColumn', eid, ni, nj, TAG['transf'], integ)
        beam_elems.append((eid, ni, nj))

    # ---- lumped mass + gravity nodal loads ----
    node_mass = {}
    for k, nd in enumerate(struct_nodes):
        y = node_coords[nd][1]
        mpl = m_col_pl if y >= 0 else m_pile_pl
        trib = 0.0
        if k > 0:                trib += abs(node_coords[nd][1]-node_coords[struct_nodes[k-1]][1])/2
        if k < len(struct_nodes)-1: trib += abs(node_coords[struct_nodes[k+1]][1]-node_coords[nd][1])/2
        node_mass[nd] = mpl*trib
    node_mass[top_node] += M_super

    # ---- SSI springs: pile node <-> fixed free-field node ----
    spring_elems = []
    soil_nodes   = []
    mat_id = 1000
    ele_id = 2000
    pile_nodes = struct_nodes[n_col+1:]   # y = -1 .. -15
    for nd in pile_nodes:
        depth = abs(node_coords[nd][1])
        is_tip = (nd == tip_node)
        trib = 1.0
        ff = 9000+nd                       # free-field (soil) node, same coords
        ops.node(ff, 0.0, node_coords[nd][1])
        ops.fix(ff, 1, 1, 1)
        soil_nodes.append(ff)
        # lateral p-y
        pult, y50 = py_capacity(depth, D_pile, trib)
        ops.uniaxialMaterial('PySimple1', mat_id, 2, pult, y50, 0.0)
        py_mat = mat_id; mat_id += 1
        # vertical t-z (shaft) or q-z (tip)
        if is_tip:
            qult, z50 = qz_capacity(depth, D_pile)
            ops.uniaxialMaterial('QzSimple1', mat_id, 2, qult, z50, 0.0, 0.0)
        else:
            tult, z50 = tz_capacity(depth, D_pile, trib)
            ops.uniaxialMaterial('TzSimple1', mat_id, 2, tult, z50, 0.0)
        vz_mat = mat_id; mat_id += 1
        ops.element('zeroLength', ele_id, ff, nd, '-mat', py_mat, vz_mat, '-dir', 1, 2)
        spring_elems.append((ele_id, ff, nd)); ele_id += 1

    return dict(node_coords=node_coords, beam_elems=beam_elems,
                spring_elems=spring_elems, soil_nodes=soil_nodes,
                struct_nodes=struct_nodes, node_mass=node_mass,
                top=top_node, ground=grnd_node, tip=tip_node)

# ----------------------------------------------------------------------------
# 4. ANALYSIS HELPERS
# ----------------------------------------------------------------------------
def assign_mass(M):
    for nd, m in M['node_mass'].items():
        ops.mass(nd, m, m, 0.0)

def run_gravity(M):
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for nd, m in M['node_mass'].items():
        ops.load(nd, 0.0, -m*g, 0.0)
    ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Transformation')
    ops.test('NormDispIncr', 1e-8, 30); ops.algorithm('Newton')
    ops.integrator('LoadControl', 0.1); ops.analysis('Static')
    ops.analyze(10)
    ops.loadConst('-time', 0.0)

# ----------------------------------------------------------------------------
# 5. BUILD + GRAVITY + MODAL
# ----------------------------------------------------------------------------
M = build_model()
assign_mass(M)
run_gravity(M)
ops.wipeAnalysis()

nmodes = 3
eig = ops.eigen('-fullGenLapack', nmodes)
periods = [2*np.pi/np.sqrt(w) for w in eig]
print("Modal periods [s]:", [round(t,3) for t in periods])

# Rayleigh damping (5% on modes 1 & 2)
zeta = 0.05
w1, w2 = np.sqrt(eig[0]), np.sqrt(eig[1])
a0 = zeta*2*w1*w2/(w1+w2)
a1 = zeta*2/(w1+w2)
ops.rayleigh(a0, 0.0, a1, 0.0)

# ----------------------------------------------------------------------------
# 6. PUSHOVER  (displacement control at pier top, +X)
# ----------------------------------------------------------------------------
ops.timeSeries('Linear', 2)
ops.pattern('Plain', 2, 2)
ops.load(M['top'], 1.0, 0.0, 0.0)

ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Transformation')
ops.test('NormDispIncr', 1e-6, 60); ops.algorithm('Newton')
dU = 0.002
ops.integrator('DisplacementControl', M['top'], 1, dU)
ops.analysis('Static')

target = 0.30                     # ~3.75% drift of column
nsteps = int(target/dU)
po_disp, po_shear = [0.0], [0.0]
for i in range(nsteps):
    if ops.analyze(1) != 0:
        ops.algorithm('NewtonLineSearch')
        if ops.analyze(1) != 0:
            print("Pushover stopped at step", i); break
        ops.algorithm('Newton')
    ops.reactions()
    V = -sum(ops.nodeReaction(n)[0] for n in M['soil_nodes'])
    po_disp.append(ops.nodeDisp(M['top'], 1))
    po_shear.append(V)

# capture deformed shape at end of pushover
push_def = {nd: (ops.nodeDisp(nd,1), ops.nodeDisp(nd,2)) for nd in M['struct_nodes']}
print(f"Pushover: peak base shear = {max(po_shear):.0f} kN at top disp "
      f"{po_disp[po_shear.index(max(po_shear))]*1000:.0f} mm")

# ----------------------------------------------------------------------------
# 7. SYNTHETIC GROUND MOTION (Kanai-Tajimi filtered noise, Jennings envelope)
# ----------------------------------------------------------------------------
def synthetic_motion(dt=0.01, T=20.0, pga_g=0.40, seed=7):
    rng = np.random.default_rng(seed)
    t = np.arange(0, T, dt); n = len(t)
    white = rng.standard_normal(n)
    wg, zg = 2*np.pi*2.5, 0.6      # Kanai-Tajimi soil filter
    a = np.zeros(n); v = 0.0; x = 0.0
    for i in range(1, n):
        acc = -2*zg*wg*v - wg**2*x + white[i]
        v += acc*dt; x += v*dt; a[i] = wg**2*x + 2*zg*wg*v
    env = np.where(t < 2, (t/2)**2, np.where(t < 12, 1.0, np.exp(-0.25*(t-12))))
    a *= env
    a *= pga_g*g/np.max(np.abs(a))
    return t, a

t_gm, a_gm = synthetic_motion()

# ----------------------------------------------------------------------------
# 8. NONLINEAR TIME-HISTORY  (uniform free-field excitation, +X)
# ----------------------------------------------------------------------------
M = build_model(); assign_mass(M); run_gravity(M)
ops.wipeAnalysis()
eig = ops.eigen('-fullGenLapack', 3)
w1, w2 = np.sqrt(eig[0]), np.sqrt(eig[1])
ops.rayleigh(zeta*2*w1*w2/(w1+w2), 0.0, zeta*2/(w1+w2), 0.0)

dt = t_gm[1]-t_gm[0]
ops.timeSeries('Path', 3, '-dt', dt, '-values', *a_gm.tolist(), '-factor', 1.0)
ops.pattern('UniformExcitation', 3, 1, '-accel', 3)

ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Transformation')
ops.test('NormDispIncr', 1e-6, 50); ops.algorithm('Newton')
ops.integrator('Newmark', 0.5, 0.25); ops.analysis('Transient')

th_t, th_top, th_shear = [], [], []
nT = len(t_gm)
ok = 0
for i in range(nT):
    ok = ops.analyze(1, dt)
    if ok != 0:                                   # fallback strategy
        ops.test('NormDispIncr', 1e-4, 100)
        ops.algorithm('NewtonLineSearch')
        ok = ops.analyze(1, dt)
        ops.test('NormDispIncr', 1e-6, 50); ops.algorithm('Newton')
        if ok != 0:
            print("THA stopped at t=%.2f s" % (i*dt)); break
    ops.reactions()
    V = -sum(ops.nodeReaction(n)[0] for n in M['soil_nodes'])
    th_t.append(ops.getTime()); th_top.append(ops.nodeDisp(M['top'],1)); th_shear.append(V)

print(f"THA: peak top drift = {max(np.abs(th_top))*1000:.0f} mm, "
      f"peak base shear = {max(np.abs(th_shear)):.0f} kN")

# ----------------------------------------------------------------------------
# 9. PUBLICATION-STYLE PLOTS
# ----------------------------------------------------------------------------
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib as mpl

# --- consistent style ------------------------------------------------------
mpl.rcParams.update({
    "font.family": "serif", "font.size": 10, "axes.linewidth": 0.8,
    "axes.edgecolor": "#333333", "axes.labelcolor": "#222222",
    "xtick.color": "#333333", "ytick.color": "#333333",
    "axes.titlesize": 11, "legend.frameon": True, "legend.framealpha": 0.95,
    "legend.edgecolor": "#999999", "legend.fontsize": 8.5,
})
# single, shared colour palette (ColorBrewer-derived, print-safe)
COL   = "#2166ac"   # column
SHAFT = "#b2701f"   # drilled shaft
MASS  = "#1a1a1a"   # superstructure mass
PY    = "#1b7837"   # lateral p-y spring
TZ    = "#762a83"   # axial t-z spring
QZ    = "#d6604d"   # tip q-z spring
SOIL  = "#efe2c0"   # soil domain
GRND  = "#555555"   # ground / fixity
DISP  = "#2166ac"   # response: displacement (always blue)
SHEAR = "#9e2a2b"   # response: base shear (always dark red)
# single line-width rule: structural members scale with diameter
LW_RULE = 3.2                       # lw per metre of diameter
LW_COL   = LW_RULE*D_col            # column line width
LW_SHAFT = LW_RULE*D_pile           # shaft line width
LW_SPR   = 1.3                      # all soil springs (consistent)
LW_RESP  = 1.6                      # all response curves (consistent)

NC = M["node_coords"]

def hspring(ax, x0, x1, y, color, n=5, amp=0.10, lw=LW_SPR):
    """Horizontal zig-zag spring symbol between x0 and x1 at height y."""
    xs = np.linspace(x0, x1, 2*n+2)
    ys = np.full_like(xs, y)
    ys[1:-1:2] += amp; ys[2:-1:2] -= amp
    ax.plot(xs, ys, color=color, lw=lw, solid_joinstyle="round", zorder=3)

def vspring(ax, x, y0, y1, color, n=4, amp=0.10, lw=LW_SPR):
    """Vertical zig-zag spring symbol between y0 and y1 at abscissa x."""
    ys = np.linspace(y0, y1, 2*n+2)
    xs = np.full_like(ys, x)
    xs[1:-1:2] += amp; xs[2:-1:2] -= amp
    ax.plot(xs, ys, color=color, lw=lw, solid_joinstyle="round", zorder=3)

# ============================================================ (a) MODEL ====
fig, ax = plt.subplots(figsize=(7.2, 8.2))

# soil domain + ground surface
ax.axhspan(-L_pile-1.2, 0, xmin=0.0, xmax=1.0, color=SOIL, zorder=0)
ax.axhline(0, color=GRND, lw=1.0, ls=(0, (6, 3)), zorder=1)

# structural members (line width scaled to diameter by a single rule)
for eid, ni, nj in M["beam_elems"]:
    above = NC[ni][1] >= 0 and NC[nj][1] >= 0
    ax.plot([NC[ni][0], NC[nj][0]], [NC[ni][1], NC[nj][1]],
            color=COL if above else SHAFT,
            lw=LW_COL if above else LW_SHAFT,
            solid_capstyle="round", zorder=4)

# soil springs (drawn clear of the geometry; consistent width)
r_pile = D_pile/2
for eid, ff, nd in M["spring_elems"]:
    yy = NC[nd][1]
    is_tip = (nd == M["tip"])
    hspring(ax, -1.7, -r_pile, yy, PY)                       # p-y (lateral)
    ax.plot([-1.85, -1.7], [yy+0.12, yy-0.12], color=GRND, lw=1.0, zorder=3)
    if is_tip:
        vspring(ax, 0.0, yy-1.05, yy, QZ)                    # q-z (tip)
        ax.plot([-0.22, 0.22], [yy-1.12, yy-1.12], color=GRND, lw=1.0)
    else:
        vspring(ax, r_pile+0.45, yy-0.45, yy+0.45, TZ, n=3, amp=0.08)  # t-z

# superstructure mass
ax.plot(0, H_col, marker="s", color=MASS, ms=15, zorder=5)

# fixity ground symbol (left, free-field)
ax.plot([-1.95, -1.95], [-L_pile-1.0, -0.2], color=GRND, lw=1.2, zorder=2)
for yy in np.arange(-L_pile-1.0, 0, 1.0):
    ax.plot([-2.05, -1.95], [yy-0.15, yy+0.10], color=GRND, lw=0.8, zorder=2)

ax.set_xlim(-2.6, 2.4); ax.set_ylim(-L_pile-2.2, H_col+1.6)
ax.set_xlabel("horizontal x [m]"); ax.set_ylabel("elevation y [m]")
ax.set_title("Bridge pier on drilled shaft with BNWF soil springs")
ax.set_aspect("equal"); ax.grid(alpha=0.25, lw=0.5)

# legend: one entry per element type, placed OUTSIDE the geometry
handles = [
    Line2D([0],[0], color=COL,   lw=LW_COL,   label=f"RC column  (\u00d8{D_col:.1f} m, H={H_col:.0f} m)"),
    Line2D([0],[0], color=SHAFT, lw=LW_SHAFT, label=f"Drilled shaft  (\u00d8{D_pile:.1f} m, L={L_pile:.0f} m)"),
    Line2D([0],[0], color=MASS,  lw=0, marker="s", ms=9,
           label=f"Superstructure mass  ({M_super:.0f} Mg)"),
    Line2D([0],[0], color=PY, lw=LW_SPR, label="$p$\u2013$y$ lateral springs"),
    Line2D([0],[0], color=TZ, lw=LW_SPR, label="$t$\u2013$z$ skin-friction springs"),
    Line2D([0],[0], color=QZ, lw=LW_SPR, label="$q$\u2013$z$ tip-bearing spring"),
    Line2D([0],[0], color=GRND, lw=1.0, ls=(0,(6,3)), label="Ground surface"),
    Patch(facecolor=SOIL, edgecolor="none", label=f"Sand ($\\phi$={phi_deg:.0f}\u00b0, $\\gamma'$={gamma_e:.0f} kN/m$^3$)"),
]
ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
          borderaxespad=0., title="Model components")
plt.tight_layout(); plt.savefig(f"{OUT}/fig1_model.png", dpi=160, bbox_inches="tight"); plt.close()

# ===================================================== (b) DEFORMED SHAPE ==
fig, ax = plt.subplots(figsize=(5.0, 8.2))
scale = 4.0
for eid, ni, nj in M["beam_elems"]:
    ax.plot([NC[ni][0], NC[nj][0]], [NC[ni][1], NC[nj][1]],
            color="0.72", lw=2.0, zorder=1)
    xi = NC[ni][0]+scale*push_def[ni][0]; xj = NC[nj][0]+scale*push_def[nj][0]
    yi = NC[ni][1]+scale*push_def[ni][1]; yj = NC[nj][1]+scale*push_def[nj][1]
    above = NC[ni][1] >= 0 and NC[nj][1] >= 0
    ax.plot([xi, xj], [yi, yj], color=COL if above else SHAFT,
            lw=LW_COL if above else LW_SHAFT, solid_capstyle="round", zorder=2)
ax.axhline(0, color=GRND, lw=1.0, ls=(0,(6,3)))
ax.set_xlabel("horizontal x [m]"); ax.set_ylabel("elevation y [m]")
ax.set_title(f"Deformed shape at end of pushover (\u00d7{scale:.0f})")
ax.set_aspect("equal"); ax.grid(alpha=0.25, lw=0.5)
handles = [
    Line2D([0],[0], color="0.72", lw=2.0, label="Undeformed"),
    Line2D([0],[0], color=COL,   lw=4.0, label="Deformed \u2013 column"),
    Line2D([0],[0], color=SHAFT, lw=5.0, label="Deformed \u2013 shaft"),
    Line2D([0],[0], color=GRND,  lw=1.0, ls=(0,(6,3)), label="Ground surface"),
]
ax.legend(handles=handles, loc="lower right", borderaxespad=0.6)
plt.tight_layout(); plt.savefig(f"{OUT}/fig2_deformed.png", dpi=160, bbox_inches="tight"); plt.close()

# ======================================================== (c) PUSHOVER =====
fig, ax = plt.subplots(figsize=(6.6, 4.6))
ax.plot(np.array(po_disp)*1000, np.array(po_shear), color=SHEAR, lw=LW_RESP,
        label="Pushover (with SSI)")
ax.set_xlabel("Top displacement [mm]"); ax.set_ylabel("Base shear [kN]")
ax.set_title("Pushover capacity curve")
ax.grid(alpha=0.25, lw=0.5); ax.legend(loc="lower right")
ax2 = ax.secondary_xaxis("top", functions=(lambda d: d/(H_col*1000)*100,
                                            lambda dr: dr*H_col*1000/100))
ax2.set_xlabel("Column drift ratio [%]")
plt.tight_layout(); plt.savefig(f"{OUT}/fig3_pushover.png", dpi=160, bbox_inches="tight"); plt.close()

# ==================================================== (d) TIME-HISTORY =====
fig, axs = plt.subplots(3, 1, figsize=(8.2, 8.0), sharex=True)
axs[0].plot(t_gm, a_gm/g, color=GRND, lw=0.9, label="Input ground motion")
axs[0].set_ylabel("Accel. [g]"); axs[0].grid(alpha=0.25, lw=0.5)
axs[0].set_title(f"Nonlinear time-history response  (PGA = {max(abs(a_gm))/g:.2f} g)")
axs[0].legend(loc="upper right")
axs[1].plot(th_t, np.array(th_top)*1000, color=DISP, lw=LW_RESP, label="Pier-top displacement")
axs[1].set_ylabel("Disp. [mm]"); axs[1].grid(alpha=0.25, lw=0.5); axs[1].legend(loc="upper right")
axs[2].plot(th_t, np.array(th_shear), color=SHEAR, lw=LW_RESP, label="Base shear")
axs[2].set_ylabel("Shear [kN]"); axs[2].set_xlabel("Time [s]")
axs[2].grid(alpha=0.25, lw=0.5); axs[2].legend(loc="upper right")
plt.tight_layout(); plt.savefig(f"{OUT}/fig4_timehistory.png", dpi=160, bbox_inches="tight"); plt.close()

# ====================================================== (e) HYSTERESIS =====
fig, ax = plt.subplots(figsize=(6.2, 5.0))
ax.plot(np.array(th_top)*1000, np.array(th_shear), color=SHEAR, lw=0.8,
        alpha=0.9, label="Base shear \u2013 top displacement")
ax.set_xlabel("Top displacement [mm]"); ax.set_ylabel("Base shear [kN]")
ax.set_title("Hysteretic response (time-history)")
ax.grid(alpha=0.25, lw=0.5); ax.legend(loc="upper left")
plt.tight_layout(); plt.savefig(f"{OUT}/fig5_hysteresis.png", dpi=160, bbox_inches="tight"); plt.close()

print("Done. Publication-style figures written.")
