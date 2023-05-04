from dolfin import *
import numpy as np
import os

"""
    Code to solve the Poisson-Nernst-Planck equations
    in an axisymmetric geometry. Solver uses a Newton solver
    from the FEniCS library (fenics 2019).
    !!! Solver has been tested and validated for the 
    parameters presented in the associated article. The 
    mesh provided here is for test purposes only!!!
    @ author: Christian Pedersen (May 2023)
"""


### SOLVER PARAMETERS

N = 100

pore_radius = 1.0
pore_length = 5.0
bulk_radius = 10.0
bulk_length = 10.0    # half bulk width

psi_wall = 0.1        # applied pore potential

DT = float(0.005)
dt = Constant(DT)

time = 0
end_time = 2


folder = 'results/'

if not os.path.exists(folder):
    os.makedirs(folder)

### LOAD MESH
mesh = Mesh('mesh/mesh.xml')
X = mesh.coordinates()

r = Expression('x[1]', degree=1)

reservoir_walls  = CompiledSubDomain('on_boundary')
axis_of_symmetry = CompiledSubDomain('x[1] < DOLFIN_EPS && on_boundary')
left_pore        = CompiledSubDomain('x[1] < (pore_radius+DOLFIN_EPS) && x[1] > (pore_radius-DOLFIN_EPS) && x[0] < 0 && on_boundary', pore_radius=pore_radius)
right_pore       = CompiledSubDomain('x[1] < (pore_radius+DOLFIN_EPS) && x[1] > (pore_radius-DOLFIN_EPS) && x[0] > 0 && on_boundary', pore_radius=pore_radius)
pore_ends        = CompiledSubDomain('x[0] < (-pore_end+DOLFIN_EPS) or x[0] > (pore_end-DOLFIN_EPS) && x[1] < pore_radius && on_boundary', pore_end=(bulk_length+pore_length), pore_radius=pore_radius)

sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
sub_domains.set_all(0)

reservoir_walls.mark(sub_domains, 1)
axis_of_symmetry.mark(sub_domains, 5)
pore_ends.mark(sub_domains, 2)
left_pore.mark(sub_domains, 3)
right_pore.mark(sub_domains, 4)

ds=Measure('ds', domain=mesh, subdomain_data=sub_domains)

### CHECK THAT BC POSITION ARE CORRECT
domain = File(folder+'domain.pvd')
domain << sub_domains

### FUNCTIONS AND FUNCTIONSPACES

elms = FiniteElement('CG', mesh.ufl_cell(), 1)

mixedspace = MixedElement([elms, elms, elms])

U, Pp, Pm= FunctionSpace(mesh, elms), FunctionSpace(mesh, elms), FunctionSpace(mesh, elms)
UP = FunctionSpace(mesh, mixedspace)

Q = FunctionSpace(mesh, 'DG', 0)
q = TestFunction(Q)
h = FacetArea(mesh)

up = Function(UP)
up0 = Function(UP)
vq = TestFunction(UP)

### INITIAL CONDITION
pp_init = interpolate(Constant(1), Pp)
pm_init = interpolate(Constant(1), Pm)

FunctionAssigner(UP.sub(1), Pp).assign(up0.sub(1), pp_init)
FunctionAssigner(UP.sub(1), Pp).assign(up.sub(1), pp_init)
FunctionAssigner(UP.sub(2), Pm).assign(up0.sub(2), pm_init)
FunctionAssigner(UP.sub(2), Pm).assign(up.sub(2), pm_init)


def Jrp(u, pp):
   return -(pp.dx(1) + pp*u.dx(1))

def Jrm(u, pm):
   return -(pm.dx(1) - pm*u.dx(1))

def Jzp(u, pp):
   return -(pp.dx(0) + pp*u.dx(0))

def Jzm(u, pm):
   return -(pm.dx(0) - pm*u.dx(0))


u, pp, pm = split(up)
u0, pp0, pm0  = split(up0)
v, qp, qm = split(vq)

n = FacetNormal(mesh)

### BOUNDARY CONDITIONS
bc_left_pore = DirichletBC(UP.sub(0), Constant(-psi_wall), left_pore)
bc_right_pore = DirichletBC(UP.sub(0), Constant(psi_wall), right_pore)

bcs = [bc_left_pore, bc_right_pore]


### VARIATIONAL FORMULATION

eq_phi = inner(r*u.dx(1), v.dx(1))*dx + inner(r*u.dx(0), v.dx(0))*dx - r*u.dx(1)*v*ds(3) - r*u.dx(1)*v*ds(4) - r*N**2*(pp - pm)*v*dx
eq_pp = r*(pp - pp0)*qp*dx - dt*inner(r*Jrp(u, pp), qp.dx(1))*dx - dt*inner(r*Jzp(u, pp), qp.dx(0))*dx
eq_pm = r*(pm - pm0)*qm*dx - dt*inner(r*Jrm(u, pm), qm.dx(1))*dx - dt*inner(r*Jzm(u, pm), qm.dx(0))*dx

eq = eq_phi + eq_pp + eq_pm

u_temp, pp_temp, pm_temp = up.split(deepcopy=True)
u0_temp, pp0_temp, pm0_temp = up0.split(deepcopy=True)


### STORE DATA
u_ = XDMFFile(folder+'psi.xdmf')
pp_ = XDMFFile(folder+'rho_p.xdmf')
pm_ = XDMFFile(folder+'rho_m.xdmf')

u_.parameters['flush_output'] = True
pp_.parameters['flush_output'] = True
pm_.parameters['flush_output'] = True

u_.write_checkpoint(u0_temp, "u", 0.0, XDMFFile.Encoding.HDF5, False)
pp_.write_checkpoint(pp_temp, "pp", 0.0, XDMFFile.Encoding.HDF5, False)
pm_.write_checkpoint(pm_temp, "pm", 0.0, XDMFFile.Encoding.HDF5, False)


frame = 100         # how often to store data
no = frame - 1      # store first time step

### SOLVER LOOP
while time < end_time+DOLFIN_EPS:

    no += 1
    time += DT
    solve(eq == 0, up, bcs)

    if no == frame:

        u_temp, pp_temp, pm_temp = up.split(deepcopy=True)          
        u_.write_checkpoint(u_temp, "u", time, XDMFFile.Encoding.HDF5, True)
        pp_.write_checkpoint(pp_temp, "pp", time, XDMFFile.Encoding.HDF5, True)
        pm_.write_checkpoint(pm_temp, "pm", time, XDMFFile.Encoding.HDF5, True)
        no = 0


    up0.assign(up)

