import pandas as pd
from pyomo.environ import *

# ================= Initialization =========================

def build_model(args):
    """
    args: input arguments
    """
    m = ConcreteModel()

    # ================== Sets ===================

    m.V = Set(initialize = args['V'], doc='Sets of nodes')
    m.Vs = Set(within = m.V, initialize = args['Vs'], doc='subset, nodes of plants')
    m.Vt = Set(within = m.V, initialize = args['Vt'], doc='subset, nodes of transshipment')
    m.Vd = Set(within = m.V, initialize = args['Vd'], doc='subset, nodes of markets')

    m.A = Set(within = m.V*m.V, initialize = args['A'], doc = 'Available links')

    # ================== Parameters ===================
    cinit = {}
    for i in m.V:
        for j in m.V:
            cinit[i,j] = args['C'][i-1][j-1]

    uinit = {}
    for i in m.V:
        for j in m.V:
            uinit[i,j] = args['U'][i-1][j-1]

    def c_init(m,i,j):
        return cinit[i,j]

    def u_init(m,i,j):
        return uinit[i,j]

    m.D = Param(m.Vd, initialize = args['D'], doc = 'annual demand from market')
    m.K = Param(m.Vs, initialize = args['K'], doc = 'annual prod. capacity of plant i')
    m.c = Param(m.V, m.V, initialize = c_init, doc = 'cost of producing, handling and/or shipping between nodes ij')
    m.u = Param(m.V, m.V, initialize = u_init, doc = 'capacity between ij')
    # m.D.pprint()
    # m.c.pprint()

    # ================== Variables ===================
    
    m.x = Var(m.V, m.V, domain = NonNegativeReals, doc = 'quantity shipped from node i to node j')

    # ================== Equations ===================

    def _Eq1(m, i):
        """
        Supply constraint LB
        """
        return 0 <= sum(m.x[i,j] for j in m.V) - sum(m.x[j,i] for j in m.V)
    m.Eq1 = Constraint(m.Vs, rule=_Eq1, doc= 'Supply constraint LB')

    def _Eq2(m, i):
        """
        Supply constraint UB
        """
        return m.K[i] >= sum(m.x[i,j] for j in m.V) - sum(m.x[j,i] for j in m.V)
    m.Eq2 = Constraint(m.Vs, rule=_Eq2, doc= 'Supply constraint UB')

    def _Eq3(m, i):
        """
        Transshipment constraint
        """
        return sum(m.x[i,j] for j in m.V) - sum(m.x[j,i] for j in m.V) == 0 
    m.Eq3 = Constraint(m.Vt, rule=_Eq3, doc= 'Transshipment constraint')

    def _Eq4(m, i):
        """
        Transshipment constraint
        """
        return sum(m.x[i,j] for j in m.V) - sum(m.x[j,i] for j in m.V) == -m.D[i] 
    m.Eq4 = Constraint(m.Vd, rule=_Eq4, doc= 'Market constraint')

    m.cost = Var(domain = NonNegativeReals)

    def _Eq5(m):
        tmp = 0
        for i in m.V:
            for j in m.V:
                tmp += m.c[i,j]*m.x[i,j]
        return tmp == m.cost
    m.Eq5 = Constraint(rule=_Eq5)


    # m.obj = Objective(expr = sum(sum(m.c[i,j]*m.x[i,j] for i in m.V) for j in m.V),\
    #     sense = minimize)
    m.obj = Objective(expr = m.cost, sense = minimize)
    return m

def init_model(m):
    for (i,j) in (m.V*m.V):
        if (i,j) in m.A:
            # print(i,j)
            m.x[i,j].setub(m.u[i,j])
        else:
            m.x[i,j].fix(0)
    
    # m.x.pprint()
    return 

def get_data(tab):
    args = {} 
    args['V'] = [i for i in range(1,8)] 
    args['Vs'] = [i for i in range(1,4)] 
    args['Vt'] = [i for i in range(4,6)] 
    args['Vd'] = [i for i in range(6,8)] 
    args['K'] = {1: 200, 2:300, 3:100}
    args['D'] = {6: 400, 7:180} 
    args['U'] = [[200 for j in args['V']] for i in args['V']]

    raw_C = pd.read_excel('HW5_data.xlsx', index_col=0)
    args['C'] = [[0 for j in args['V']] for i in args['V']]
    args['A'] = []
    def is_NAN(x):
        return x != x
    for j in args['V']:
        for i in args['V']:
            if not is_NAN(raw_C[j][i]):
                args['C'][i-1][j-1] = raw_C[j][i]
                args['A'].append((i,j))
    return args

dat = get_data('HW5_data.xlsx')
m = build_model(dat)
init_model(m)

opt = SolverFactory('gams')
io_options = dict() 

io_options['solver'] = "cplex"

res = opt.solve(m,
    tee=True,
    keepfiles=True,
    symbolic_solver_labels=True,
    #add_options = ['GAMS_MODEL.optfile = 1; option reslim=120; option optcr=0.0;'],
    add_options = ['option reslim=180; option optcr=0.0;'],
    tmpdir='/home/zyuliu/06653',
    io_options=io_options)
# print(value(sum(sum(m.c[i,j]*m.x[i,j] for i in m.V) for j in m.V)))

#Exporting to excel
res_list = []

for i in m.V:
    row = []
    for j in m.V:
        row.append(value(m.x[i,j]))
    res_list.append(row)

df1 = pd.DataFrame(res_list)
with pd.ExcelWriter('HW5_result.xlsx') as writer:  
    df1.to_excel(writer)

print('printing infeasible constraints')

from pyomo.util.infeasible import log_infeasible_constraints
log_infeasible_constraints(m)
