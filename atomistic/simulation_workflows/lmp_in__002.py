
# Original npt template (commented out)
# npt = '''variable dt equal 0.003
# variable StepNPT equal 100000
# variable PreStep equal 0
# variable nevery equal 100
# variable out_dump_file_0 string force.0.dump
# # variable out_dump_file string
# # variable mlip_ini string
#
# units metal
# boundary p p p
# atom_style atomic
#
# box tilt large
# read_data data.in
#
#
# pair_style  mlip ${mlip_ini}
# pair_coeff * *
#
# neighbor 2.0 bin
# neigh_modify delay 10 check yes
#
# timestep    ${dt}
# variable Tdamp equal "v_dt * 100"
# variable Pdamp equal "v_dt * 1000"
# velocity        all create ${T} ${random} mom yes rot yes dist gaussian
# thermo 100
#
# fix 100 all npt temp ${T} ${T} ${Tdamp} aniso 0.0 0.0 ${Pdamp}
# velocity all scale ${T}
# velocity all zero linear
#
# run             ${PreStep}
#
# reset_timestep 0
#
# compute pe all pe/atom
# dump 0 all custom ${nevery} ${out_dump_file_0} id type x y z fx fy fz c_pe
# dump_modify 0 format line "%6d %3d %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e" sort id
#
# dump 1 all custom ${nevery} ${out_dump_file} id type x y z fx fy fz c_pe
# dump_modify 1 sort id
#
# run             ${StepNPT}
# '''

# New npt template for Annealing (based on HfNbTiVAl_template.in)
npt = '''variable dt equal 0.001
variable StepNPT equal 8000000
variable StepEquil equal 40000
# T_start comes from external ${T}
variable T_start equal ${T}
variable T_end equal 300
variable P_target equal 0.0
variable nevery equal 4000
variable out_dump_file_0 string force.0.dump
# variable out_dump_file string
# variable mlip_ini string
# variable random equal 12345

units metal
boundary p p p
atom_style atomic


box tilt large
read_data data.in


# 【Core Fix 1】: Manually add masses for HfNbTiVAl (Solves ASE missing masses issue)
mass 1 26.9815 # Al
mass 2 47.867  # Ti
mass 3 50.9415 # V
mass 4 92.9064 # Nb
mass 5 178.49  # Hf



pair_style mlip ${mlip_ini}
pair_coeff * *

neighbor 2.0 bin
neigh_modify delay 10 check yes

timestep ${dt}
variable Tdamp equal "v_dt * 100"
variable Pdamp equal "v_dt * 1000"

# 【Core Fix 4】: Use rot no to avoid angular momentum issues with p p p boundary
velocity all create ${T_start} ${random} mom yes rot no dist gaussian

thermo ${nevery}
thermo_style custom step temp pe press vol lx ly lz

# ----------------- Equilibration Phase (Isothermal-Isobaric at T_start) -----------------
fix 100 all npt temp ${T_start} ${T_start} ${Tdamp} iso ${P_target} ${P_target} ${Pdamp}
run ${StepEquil}

reset_timestep 0

compute pe all pe/atom

# Dump 0 (forces reference)
dump 0 all custom ${nevery} ${out_dump_file_0} id type x y z fx fy fz c_pe
dump_modify 0 format line "%6d %3d %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e" sort id

# Dump 1 (trajectory)
dump 1 all custom ${nevery} ${out_dump_file} id type x y z fx fy fz c_pe
# 【Core Fix 5】: Element mapping for OVITO
dump_modify 1 element Al Ti V Nb Hf format float "%20.15e" sort id

# ----------------- Production Run: Annealing (T_start -> T_end) -----------------
unfix 100
fix 100 all npt temp ${T_start} ${T_end} ${Tdamp} iso ${P_target} ${P_target} ${Pdamp}
run ${StepNPT}

# Output final structure at 300K
write_data final_300K_structure.data
'''

nvt = '''variable dt equal 0.003
variable StepNVT equal 100000
variable PreStep equal 0
variable nevery equal 100
variable out_dump_file_0 string force.0.dump
# variable out_dump_file string
# variable mlip_ini string

units metal
boundary p p p
atom_style atomic

box tilt large
read_data data.in


pair_style  mlip ${mlip_ini}
pair_coeff * *

neighbor 2.0 bin
neigh_modify delay 10 check yes

timestep    ${dt}
variable Tdamp equal "v_dt * 100"
variable Pdamp equal "v_dt * 1000"
velocity        all create ${T} ${random} mom yes rot yes dist gaussian
thermo 100

fix 100 all nvt temp ${T} ${T} ${Tdamp}
velocity all scale ${T}
velocity all zero linear

run             ${PreStep}

reset_timestep 0

compute pe all pe/atom
dump 0 all custom ${nevery} ${out_dump_file_0} id type x y z fx fy fz c_pe
dump_modify 0 format line "%6d %3d %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e %20.15e" sort id

dump 1 all custom ${nevery} ${out_dump_file} id type x y z fx fy fz c_pe
dump_modify 1 sort id

run             ${StepNVT}
'''
