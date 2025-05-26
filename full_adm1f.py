{\rtf1\ansi\ansicpg1252\cocoartf2580
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;\f1\fswiss\fcharset0 ArialMT;}
{\colortbl;\red255\green255\blue255;\red26\green26\blue26;\red255\green255\blue255;}
{\*\expandedcolortbl;;\cssrgb\c13333\c13333\c13333;\cssrgb\c100000\c100000\c100000;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 f
\f1\fs20 \cf2 \cb3 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 rom petsc4py import PETSc\
\pard\pardeftab720\partightenfactor0
\cf2 import numpy as np\
\
MAXLINE = 1000\
\
class AppCtx:\
\'a0 \'a0 def __init__(self):\
\'a0 \'a0 \'a0 \'a0 self.initialconditions = None\
\'a0 \'a0 \'a0 \'a0 self.params = None\
\'a0 \'a0 \'a0 \'a0 self.interface_params = None\
\'a0 \'a0 \'a0 \'a0 self.influent = None\
\'a0 \'a0 \'a0 \'a0 self.adm1_output = None\
\'a0 \'a0 \'a0 \'a0 self.asm1_output = None\
\'a0 \'a0 \'a0 \'a0 self.indicator = None\
\'a0 \'a0 \'a0 \'a0 self.V = [0.0, 0.0]\
\'a0 \'a0 \'a0 \'a0 self.adctx = None \'a0# Placeholder for automatic differentiation\
\'a0 \'a0 \'a0 \'a0 self.rwork = 0.0\
\'a0 \'a0 \'a0 \'a0 self.debug = False\
\'a0 \'a0 \'a0 \'a0 self.Cat_mass = 0.0\
\'a0 \'a0 \'a0 \'a0 self.set_Cat_mass = False\
\'a0 \'a0 \'a0 \'a0 self.t_resx = 0.0\
\'a0 \'a0 \'a0 \'a0 self.set_t_resx = False\
\
def read_vector_from_file(comm, filename, size, debug=False):\
\'a0 \'a0 vec = PETSc.Vec().create(comm=comm)\
\'a0 \'a0 vec.setSizes(size)\
\'a0 \'a0 vec.setFromOptions()\
\
\'a0 \'a0 if comm.getRank() == 0:\
\'a0 \'a0 \'a0 \'a0 print(f"Reading values from file: \{filename\}")\
\'a0 \'a0 values = np.loadtxt(filename, dtype=float, max_rows=size)\
\'a0 \'a0 vec.setValues(range(size), values)\
\'a0 \'a0 vec.assemblyBegin()\
\'a0 \'a0 vec.assemblyEnd()\
\
\'a0 \'a0 if debug:\
\'a0 \'a0 \'a0 \'a0 vec.view()\
\
\'a0 \'a0 return vec\
\
def ReadParams(ctx, filename):\
\'a0 \'a0 num_params = 103\
\'a0 \'a0 ctx.params = read_vector_from_file(PETSc.COMM_WORLD, filename, num_params, ctx.debug)\
\
def ReadInfluent(ctx, filename):\
\'a0 \'a0 num_influent = 28\
\'a0 \'a0 ctx.influent = read_vector_from_file(PETSc.COMM_WORLD, filename, num_influent, ctx.debug)\
\
def ReadInitialConditions(ctx, filename):\
\'a0 \'a0 num_ic = 45\
\'a0 \'a0 ctx.initialconditions = read_vector_from_file(PETSc.COMM_WORLD, filename, num_ic, ctx.debug)\
\
def DigestParToInterfacePar(ctx):\
\'a0 \'a0 num_interface_params = 19\
\'a0 \'a0 ctx.interface_params = PETSc.Vec().create(comm=PETSc.COMM_WORLD)\
\'a0 \'a0 ctx.interface_params.setSizes(num_interface_params)\
\'a0 \'a0 ctx.interface_params.setFromOptions()\
\
\'a0 \'a0 params = ctx.params.getArray(readonly=True)\
\'a0 \'a0 interface_params = np.zeros(num_interface_params)\
\
\'a0 \'a0 interface_params[0] = 0.0\
\'a0 \'a0 interface_params[1] = 0.79\
\'a0 \'a0 interface_params[2] = 0.0\
\'a0 \'a0 interface_params[3] = params[7] * 14\
\'a0 \'a0 interface_params[4] = params[5] * 14\
\'a0 \'a0 interface_params[5] = params[22] * 14\
\'a0 \'a0 interface_params[6] = params[6] * 14\
\'a0 \'a0 interface_params[7] = params[6] * 14\
\
\'a0 \'a0 for i in range(8, 19):\
\'a0 \'a0 \'a0 \'a0 interface_params[i] = params[i + 69]\
\
\'a0 \'a0 ctx.interface_params.setValues(range(num_interface_params), interface_params)\
\'a0 \'a0 ctx.interface_params.assemblyBegin()\
\'a0 \'a0 ctx.interface_params.assemblyEnd()\
\
\'a0 \'a0 if ctx.debug:\
\'a0 \'a0 \'a0 \'a0 print("Interface Params:")\
\'a0 \'a0 \'a0 \'a0 ctx.interface_params.view()\
\
def ProcessADM1(adm1_sol, ctx):\
\'a0 \'a0 x = adm1_sol.getArray(readonly=True)\
\'a0 \'a0 params = ctx.params.getArray(readonly=True)\
\'a0 \'a0 influent = ctx.influent.getArray(readonly=True)\
\
\'a0 \'a0 if ctx.adm1_output is None:\
\'a0 \'a0 \'a0 \'a0 ctx.adm1_output = PETSc.Vec().create(comm=PETSc.COMM_WORLD)\
\'a0 \'a0 \'a0 \'a0 ctx.adm1_output.setSizes(54)\
\'a0 \'a0 \'a0 \'a0 ctx.adm1_output.setFromOptions()\
\'a0 \'a0 \'a0 \'a0 ctx.adm1_output.set(0.0)\
\
\'a0 \'a0 y = ctx.adm1_output.getArray()\
\
\'a0 \'a0 R = params[77]\
\'a0 \'a0 T_base = params[78]\
\'a0 \'a0 T_op = params[79]\
\'a0 \'a0 P_atm = params[93]\
\'a0 \'a0 V_liq = ctx.V[0]\
\'a0 \'a0 kLa = params[94]\
\'a0 \'a0 pK_w_base = params[80]\
\'a0 \'a0 k_P = params[99]\
\'a0 \'a0 factor = (1.0 / T_base - 1.0 / T_op) / (100.0 * R)\
\
\'a0 \'a0 K_H_h2 = 1.0 / 10**(-187.04 / T_op + 5.473) * 55.6 / 1.01325\
\'a0 \'a0 K_H_ch4 = 1.0 / 10**(-675.74 / T_op + 6.880) * 55.6 / 1.01325\
\'a0 \'a0 K_H_co2 = 1.0 / 10**(-1012.40 / T_op + 6.606) * 55.6 / 1.01325\
\'a0 \'a0 K_w = 10**(-pK_w_base) * np.exp(55700.0 * factor)\
\'a0 \'a0 p_gas_h2o = 10**(5.20389 - 1733.926 / (T_op - 39.485))\
\
\'a0 \'a0 for i in range(26):\
\'a0 \'a0 \'a0 \'a0 y[i] = x[i]\
\'a0 \'a0 y[26] = influent[26]\
\'a0 \'a0 y[27] = T_op - 273.15\
\'a0 \'a0 y[28:33] = x[37:42]\
\
\'a0 \'a0 p_gas_h2 = x[32] * R * T_op / 16.0\
\'a0 \'a0 p_gas_ch4 = x[33] * R * T_op / 64.0\
\'a0 \'a0 p_gas_co2 = x[34] * R * T_op\
\'a0 \'a0 P_gas = p_gas_h2 + p_gas_ch4 + p_gas_co2 + p_gas_h2o\
\'a0 \'a0 q_gas = k_P * (P_gas - P_atm)\
\
\'a0 \'a0 y[33] = -np.log10(x[42])\
\'a0 \'a0 y[34] = x[42]\
\'a0 \'a0 y[35] = x[26]\
\'a0 \'a0 y[36] = x[27]\
\'a0 \'a0 y[37] = x[28]\
\'a0 \'a0 y[38] = x[29]\
\'a0 \'a0 y[39] = x[30]\
\'a0 \'a0 y[40] = x[43]\
\'a0 \'a0 y[41] = x[31]\
\'a0 \'a0 y[42] = x[44]\
\'a0 \'a0 y[43] = x[32]\
\'a0 \'a0 y[44] = x[33]\
\'a0 \'a0 y[45] = x[34]\
\'a0 \'a0 y[46] = p_gas_h2\
\'a0 \'a0 y[47] = p_gas_ch4\
\'a0 \'a0 y[48] = p_gas_co2\
\'a0 \'a0 y[49] = P_gas\
\'a0 \'a0 y[50] = q_gas * P_gas / P_atm\
\'a0 \'a0 y[51] = x[7]\
\'a0 \'a0 y[52] = V_liq / influent[26]\
\'a0 \'a0 y[53] = ctx.rwork\
\
\'a0 \'a0 ctx.adm1_output.restoreArray()\
\'a0 \'a0 ctx.params.restoreArrayRead()\
\'a0 \'a0 ctx.influent.restoreArrayRead()\
\'a0 \'a0 adm1_sol.restoreArrayRead()\
\
def ProcessIndicators(ctx):\
\'a0 \'a0 if ctx.indicator is None:\
\'a0 \'a0 \'a0 \'a0 ctx.indicator = PETSc.Vec().create(comm=PETSc.COMM_WORLD)\
\'a0 \'a0 \'a0 \'a0 ctx.indicator.setSizes(67)\
\'a0 \'a0 \'a0 \'a0 ctx.indicator.setFromOptions()\
\'a0 \'a0 \'a0 \'a0 ctx.indicator.set(0.0)\
\
\'a0 \'a0 adm1_out = ctx.adm1_output.getArray(readonly=True)\
\'a0 \'a0 infl = ctx.influent.getArray(readonly=True)\
\'a0 \'a0 y = ctx.indicator.getArray()\
\
\'a0 \'a0 for i in range(min(54, len(adm1_out))):\
\'a0 \'a0 \'a0 \'a0 y[i] = adm1_out[i] * 1000.0 \'a0# mg/L\
\
\'a0 \'a0 y[54] = y[3] + y[4] + y[5] + y[6]\
\'a0 \'a0 y[66] = adm1_out[52]\
\
\'a0 \'a0 ctx.adm1_output.restoreArrayRead()\
\'a0 \'a0 ctx.influent.restoreArrayRead()\
\'a0 \'a0 ctx.indicator.restoreArray()\
\
def PostProcess(ts, ctx):\
\'a0 \'a0 adm1_sol = ts.getSolution()\
\'a0 \'a0 ProcessADM1(adm1_sol, ctx)\
\'a0 \'a0 ProcessIndicators(ctx)\
\
\'a0 \'a0 viewer = PETSc.ViewerASCII().create(PETSc.COMM_WORLD)\
\'a0 \'a0 viewer.setFileName("adm1_output.out")\
\'a0 \'a0 ctx.adm1_output.view(viewer)\
\
\'a0 \'a0 viewer.setFileName("indicator.out")\
\'a0 \'a0 ctx.indicator.view(viewer)\
\
def IFunctionPassive(ts, t, U, Udot, F, ctx):\
\'a0 \'a0 u = U.getArray(readonly=True)\
\'a0 \'a0 udot = Udot.getArray(readonly=True)\
\'a0 \'a0 f = F.getArray()\
\
\'a0 \'a0 for i in range(U.getSize()):\
\'a0 \'a0 \'a0 \'a0 f[i] = udot[i] + u[i] \'a0# dummy ODE: du/dt = -u\
\
\'a0 \'a0 U.restoreArrayRead()\
\'a0 \'a0 Udot.restoreArrayRead()\
\'a0 \'a0 F.restoreArray()\
\
def setup_ts_solver(ctx):\
\'a0 \'a0 ts = PETSc.TS().create(PETSc.COMM_WORLD)\
\'a0 \'a0 ts.setProblemType(PETSc.TS.ProblemType.NONLINEAR)\
\'a0 \'a0 ts.setType('beuler') \'a0# Backward Euler\
\'a0 \'a0 ts.setIFunction(IFunctionPassive, ctx)\
\
\'a0 \'a0 ts.setTime(0.0)\
\'a0 \'a0 ts.setMaxTime(10.0)\
\'a0 \'a0 ts.setTimeStep(0.1)\
\'a0 \'a0 ts.setMaxSteps(1000)\
\
\'a0 \'a0 ts.setSolution(ctx.initialconditions.copy())\
\'a0 \'a0 ts.setFromOptions()\
\
\'a0 \'a0 return ts\
\
def main():\
\'a0 \'a0 OptDB = PETSc.Options()\
\'a0 \'a0 debug_flag = OptDB.getBool("debug", default=False)\
\
\'a0 \'a0 ctx = AppCtx()\
\'a0 \'a0 ctx.debug = debug_flag\
\
\'a0 \'a0 param_file = OptDB.getString("param_file", default="params.txt")\
\'a0 \'a0 influent_file = OptDB.getString("influent_file", default="influent.txt")\
\'a0 \'a0 ic_file = OptDB.getString("ic_file", default="initialconditions.txt")\
\
\'a0 \'a0 ReadParams(ctx, param_file)\
\'a0 \'a0 ReadInfluent(ctx, influent_file)\
\'a0 \'a0 ReadInitialConditions(ctx, ic_file)\
\'a0 \'a0 DigestParToInterfacePar(ctx)\
\
\'a0 \'a0 ts = setup_ts_solver(ctx)\
\'a0 \'a0 ts.solve(ctx.initialconditions)\
\'a0 \'a0 PostProcess(ts, ctx)\
\
if __name__ == "__main__":\
\'a0 \'a0 main()}